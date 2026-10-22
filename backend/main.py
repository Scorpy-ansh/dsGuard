
import csv
import json
import os
import re
import shutil
import subprocess
import threading
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Annotated
import csv
import io

from dotenv import load_dotenv
from fastapi import (
    BackgroundTasks,
    FastAPI,
    File,
    Form,
    HTTPException,
    UploadFile,
)
from fastapi.middleware.cors import CORSMiddleware
from pymongo import MongoClient
from pymongo.errors import PyMongoError
from fastapi.responses import FileResponse


# ---------------------------------------------------------
# Paths and configuration
# ---------------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent
ENGINE_DIR = ROOT / "Random"
DSGUARD_EXE = ROOT / ".venv" / "Scripts" / "dsguard.exe"

UPLOAD_ROOT = ROOT / "backend" / "uploads"
OUTPUT_ROOT = ROOT / "outputs" / "jobs"

load_dotenv(ROOT / ".env")

MONGODB_URI = os.getenv("MONGODB_URI")
MONGODB_DATABASE = os.getenv("MONGODB_DATABASE", "dsguard")

if not MONGODB_URI:
    raise RuntimeError(
        "MONGODB_URI is missing. Check the project's .env file."
    )

mongo_client = MongoClient(
    MONGODB_URI,
    serverSelectionTimeoutMS=10000,
)

db = mongo_client[MONGODB_DATABASE]
jobs = db["jobs"]

UPLOAD_LIMIT_BYTES = 2 * 1024 * 1024 * 1024

app = FastAPI(
    title="dsGuard API",
    description="Local API for dsGuard sequence screening",
    version="0.3.0",
    openapi_version="3.0.2",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:5173",
        "http://127.0.0.1:5173",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

job_lock = threading.Lock()


# ---------------------------------------------------------
# Utility functions
# ---------------------------------------------------------

def now_utc():
    return datetime.now(timezone.utc).isoformat()


def safe_panel_name(value: str) -> str:
    return re.sub(
        r"[^A-Za-z0-9_-]",
        "_",
        value.strip(),
    )[:40]


def save_upload(upload: UploadFile, destination: Path) -> int:
    """Save and perform basic structural validation on a FASTA upload."""

    if not upload.filename:
        raise ValueError("Uploaded file has no filename.")

    allowed_extensions = (
        ".fa",
        ".fasta",
        ".fna",
        ".fas",
        ".txt",
    )

    if not upload.filename.lower().endswith(allowed_extensions):
        raise ValueError(
            f"{upload.filename}: use a FASTA extension "
            "(.fa, .fasta, .fna, .fas, or .txt)."
        )

    destination.parent.mkdir(parents=True, exist_ok=True)

    total_bytes = 0

    try:
        with destination.open("wb") as output:
            while True:
                chunk = upload.file.read(1024 * 1024)

                if not chunk:
                    break

                total_bytes += len(chunk)

                if total_bytes > UPLOAD_LIMIT_BYTES:
                    raise ValueError(
                        "Each uploaded file must be 2 GiB or smaller."
                    )

                output.write(chunk)

        if total_bytes == 0:
            raise ValueError(f"{upload.filename} is empty.")

        # Basic FASTA check: first non-empty line must be a header.
        first_content = b""

        with destination.open("rb") as source:
            for line in source:
                if line.strip():
                    first_content = line.strip()
                    break

        if not first_content.startswith(b">"):
            raise ValueError(
                f"{upload.filename} is not recognized as FASTA: "
                "the first non-empty line must start with '>'."
            )

        return total_bytes

    except Exception:
        destination.unlink(missing_ok=True)
        raise


# ---------------------------------------------------------
# Screening worker
# ---------------------------------------------------------

def run_screening(
    job_id: str,
    target_path: Path,
    panel_paths: list[tuple[str, Path]],
    params: dict,
):
    """Run the existing dsGuard CLI and persist its results."""

    job_output = OUTPUT_ROOT / job_id

    jobs.update_one(
        {"_id": job_id},
        {
            "$set": {
                "status": "running",
                "started_at": now_utc(),
            }
        },
    )

    try:
        if not DSGUARD_EXE.is_file():
            raise RuntimeError(
                f"dsGuard executable not found: {DSGUARD_EXE}"
            )

        job_output.mkdir(parents=True, exist_ok=True)

        command = [
            str(DSGUARD_EXE),
            "--target",
            str(target_path),
            "--window",
            str(params["window"]),
            "--k",
            str(params["k"]),
            "--l50",
            str(params["l50"]),
            "--kappa",
            str(params["kappa"]),
            "--max-hazard",
            str(params["max_hazard"]),
            "--out",
            str(job_output),
        ]

        for panel_name, panel_path in panel_paths:
            command.extend(
                [
                    "--panel",
                    f"{panel_name}={panel_path}",
                ]
            )

        result = subprocess.run(
            command,
            cwd=str(ENGINE_DIR),
            capture_output=True,
            text=True,
            timeout=24 * 60 * 60,
            check=False,
        )

        (job_output / "engine.log").write_text(
            "COMMAND:\n"
            + subprocess.list2cmdline(command)
            + "\n\nSTDOUT:\n"
            + result.stdout
            + "\n\nSTDERR:\n"
            + result.stderr,
            encoding="utf-8",
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"dsGuard exited with code {result.returncode}. "
                "Inspect the job's engine.log."
            )

        summary_path = job_output / "summary.csv"
        manifest_path = job_output / "manifest.json"

        if not summary_path.is_file():
            raise RuntimeError("dsGuard did not produce summary.csv.")

        if not manifest_path.is_file():
            raise RuntimeError("dsGuard did not produce manifest.json.")

        with summary_path.open(
            "r",
            newline="",
            encoding="utf-8",
        ) as source:
            summary_rows = list(csv.DictReader(source))

        manifest = json.loads(
            manifest_path.read_text(encoding="utf-8")
        )

        jobs.update_one(
            {"_id": job_id},
            {
                "$set": {
                    "status": "completed",
                    "completed_at": now_utc(),
                    "summary": summary_rows,
                    "manifest": manifest,
                    "result_directory": str(job_output),
                    "engine_exit_code": result.returncode,
                }
            },
        )

    except Exception as exc:
        jobs.update_one(
            {"_id": job_id},
            {
                "$set": {
                    "status": "failed",
                    "completed_at": now_utc(),
                    "error": str(exc)[:2000],
                    "result_directory": str(job_output),
                }
            },
        )


# ---------------------------------------------------------
# Basic endpoints
# ---------------------------------------------------------

@app.get("/")
def home():
    return {
        "project": "dsGuard",
        "message": "dsGuard backend is running",
        "docs": "/docs",
    }


@app.get("/api/health")
def health():
    mongo_status = "connected"

    try:
        mongo_client.admin.command("ping")
    except PyMongoError:
        mongo_status = "disconnected"

    return {
        "status": (
            "ok" if mongo_status == "connected" else "degraded"
        ),
        "mongodb": mongo_status,
        "engine_exists": DSGUARD_EXE.is_file(),
        "timestamp": now_utc(),
    }


# ---------------------------------------------------------
# Create screening job: one target + one non-target panel
# ---------------------------------------------------------

@app.post("/api/screen")
async def create_screening_job(
    background_tasks: BackgroundTasks,
    target_file: Annotated[UploadFile, File()],
    panel_file: Annotated[UploadFile, File()],
    panel_name: Annotated[str, Form()],
    window: Annotated[int, Form()] = 400,
    k: Annotated[int, Form()] = 11,
    l50: Annotated[float, Form()] = 21.0,
    kappa: Annotated[float, Form()] = 2.0,
    max_hazard: Annotated[float, Form()] = 0.05,
):
    if not DSGUARD_EXE.is_file():
        raise HTTPException(
            status_code=500,
            detail="dsGuard executable was not found.",
        )

    name = safe_panel_name(panel_name)

    if not name:
        raise HTTPException(
            status_code=400,
            detail="Panel name is required.",
        )

    if window < 1:
        raise HTTPException(
            status_code=400,
            detail="Window must be positive.",
        )

    if not 2 <= k <= 31:
        raise HTTPException(
            status_code=400,
            detail="k must be between 2 and 31.",
        )

    if l50 <= 0 or kappa <= 0:
        raise HTTPException(
            status_code=400,
            detail="l50 and kappa must be positive.",
        )

    if not 0 <= max_hazard <= 1:
        raise HTTPException(
            status_code=400,
            detail="max_hazard must be between 0 and 1.",
        )

    job_id = str(uuid.uuid4())

    job_upload_dir = UPLOAD_ROOT / job_id
    job_output_dir = OUTPUT_ROOT / job_id

    target_path = job_upload_dir / "target.fasta"
    panel_path = job_upload_dir / "panel_1.fasta"

    try:
        save_upload(target_file, target_path)
        save_upload(panel_file, panel_path)

    except ValueError as exc:
        shutil.rmtree(job_upload_dir, ignore_errors=True)

        raise HTTPException(
            status_code=400,
            detail=str(exc),
        ) from exc

    finally:
        await target_file.close()
        await panel_file.close()

    params = {
        "window": window,
        "k": k,
        "l50": l50,
        "kappa": kappa,
        "max_hazard": max_hazard,
    }

    document = {
        "_id": job_id,
        "status": "queued",
        "created_at": now_utc(),
        "parameters": params,
        "panel_names": [name],
        "target_filename": Path(
            target_file.filename or "target.fasta"
        ).name,
        "panel_filenames": [
            Path(panel_file.filename or "panel.fasta").name
        ],
        "upload_directory": str(job_upload_dir),
        "result_directory": str(job_output_dir),
    }

    try:
        jobs.insert_one(document)

    except PyMongoError as exc:
        shutil.rmtree(job_upload_dir, ignore_errors=True)

        raise HTTPException(
            status_code=503,
            detail="Could not save job metadata to MongoDB.",
        ) from exc

    background_tasks.add_task(
        run_screening,
        job_id,
        target_path,
        [(name, panel_path)],
        params,
    )

    return {
        "job_id": job_id,
        "status": "queued",
        "message": "Screening job accepted.",
        "status_url": f"/api/jobs/{job_id}",
    }


# ---------------------------------------------------------
# Job history and results
# ---------------------------------------------------------


@app.get("/api/jobs")
def list_jobs(limit: int = 20):
    limit = max(1, min(limit, 100))

    documents = jobs.find(
        {},
        {
            "_id": 1,
            "status": 1,
            "created_at": 1,
            "started_at": 1,
            "completed_at": 1,
            "target_filename": 1,
            "panel_names": 1,
            "error": 1,
            "summary": 1,
            "manifest": 1,
            "parameters": 1,
        },
    ).sort("created_at", -1).limit(limit)

    return [
        {
            "job_id": document["_id"],
            **{
                key: value
                for key, value in document.items()
                if key != "_id"
            },
        }
        for document in documents
    ]



@app.get("/api/jobs/{job_id}")
def get_job(job_id: str):
    document = jobs.find_one(
        {"_id": job_id},
        {"_id": 0},
    )

    if not document:
        raise HTTPException(
            status_code=404,
            detail="Job not found.",
        )

    return document



@app.delete("/api/jobs/{job_id}")
def delete_job(job_id: str):
    """Delete one completed or failed screening and its stored files."""

    with job_lock:
        document = jobs.find_one({"_id": job_id})

        if not document:
            raise HTTPException(
                status_code=404,
                detail="Job not found.",
            )

        if document.get("status") in {"queued", "running"}:
            raise HTTPException(
                status_code=409,
                detail="Cannot delete a screening while it is queued or running.",
            )

        upload_dir = (UPLOAD_ROOT / job_id).resolve()
        output_dir = (OUTPUT_ROOT / job_id).resolve()

        # Ensure both paths are direct children of their configured roots.
        if (
            upload_dir.parent != UPLOAD_ROOT.resolve()
            or output_dir.parent != OUTPUT_ROOT.resolve()
        ):
            raise HTTPException(
                status_code=400,
                detail="Invalid job path.",
            )

        # Remove files first. Retain the database record if cleanup fails.
        try:
            for directory in (upload_dir, output_dir):
                if directory.exists():
                    shutil.rmtree(directory)
        except OSError as exc:
            raise HTTPException(
                status_code=500,
                detail="Could not remove screening files. The screening record was retained.",
            ) from exc

        try:
            result = jobs.delete_one({"_id": job_id})

            if result.deleted_count != 1:
                raise HTTPException(
                    status_code=404,
                    detail="Job not found.",
                )

        except PyMongoError as exc:
            raise HTTPException(
                status_code=500,
                detail="Files were removed, but the screening record could not be deleted.",
            ) from exc

    return {
        "job_id": job_id,
        "status": "deleted",
        "message": "Screening history entry and associated files deleted.",
    }



@app.get("/api/jobs/{job_id}/results")
def get_job_results(job_id: str):
    document = jobs.find_one(
        {"_id": job_id},
        {
            "_id": 0,
            "status": 1,
            "summary": 1,
            "manifest": 1,
            "error": 1,
        },
    )

    if not document:
        raise HTTPException(
            status_code=404,
            detail="Job not found.",
        )

    return document


@app.get("/api/jobs/{job_id}/reports")
def list_gene_reports(job_id: str):
    document = jobs.find_one(
        {"_id": job_id},
        {"_id": 0, "status": 1},
    )

    if not document:
        raise HTTPException(
            status_code=404,
            detail="Job not found.",
        )

    if document.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Screening job is not completed.",
        )

    job_dir = (OUTPUT_ROOT / job_id).resolve()
    output_root = OUTPUT_ROOT.resolve()

    if not job_dir.is_relative_to(output_root):
        raise HTTPException(
            status_code=400,
            detail="Invalid output path.",
        )

    if not job_dir.is_dir():
        raise HTTPException(
            status_code=404,
            detail="Job output directory not found.",
        )

    reports = [
        {
            "gene": folder.name,
            "display_name": folder.name.split("_accession_")[0],
        }
        for folder in sorted(job_dir.iterdir())
        if folder.is_dir()
        and (folder / "report.md").is_file()
    ]

    return {
        "job_id": job_id,
        "total_reports": len(reports),
        "reports": reports,
    }


@app.get("/api/jobs/{job_id}/reports/{gene_name}")
def get_gene_report(job_id: str, gene_name: str):
    if not re.fullmatch(r"[A-Za-z0-9._-]{1,100}", gene_name):
        raise HTTPException(
            status_code=400,
            detail="Invalid gene name.",
        )

    document = jobs.find_one(
        {"_id": job_id},
        {"_id": 0, "status": 1},
    )

    if not document:
        raise HTTPException(
            status_code=404,
            detail="Job not found.",
        )

    if document.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Screening job is not completed.",
        )

    report_path = (
        OUTPUT_ROOT / job_id / gene_name / "report.md"
    )

    if not report_path.is_file():
        raise HTTPException(
            status_code=404,
            detail="Gene report not found.",
        )

    return {
        "gene": gene_name,
        "report": report_path.read_text(encoding="utf-8"),
    }


from fastapi.responses import StreamingResponse


@app.get("/api/jobs/{job_id}/hits")
def get_job_hits(job_id: str):
    document = jobs.find_one(
        {"_id": job_id},
        {"_id": 0, "status": 1},
    )

    if not document:
        raise HTTPException(
            status_code=404,
            detail="Job not found.",
        )

    if document.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Screening job is not completed.",
        )

    job_dir = (OUTPUT_ROOT / job_id).resolve()
    output_root = OUTPUT_ROOT.resolve()

    if not job_dir.is_relative_to(output_root):
        raise HTTPException(
            status_code=400,
            detail="Invalid output path.",
        )

    if not job_dir.is_dir():
        raise HTTPException(
            status_code=404,
            detail="Job output directory not found.",
        )

    hits = []

    for csv_path in sorted(job_dir.glob("*/hits.csv")):
        gene_name = csv_path.parent.name.split("_accession_")[0]

        with csv_path.open(
            "r",
            newline="",
            encoding="utf-8-sig",
        ) as handle:
            reader = csv.DictReader(handle)

            for row in reader:
                if not row or not any(
                    value and value.strip()
                    for value in row.values()
                ):
                    continue

                hits.append({
                    "gene": gene_name,
                    **row,
                })

    return {
        "job_id": job_id,
        "total_hits": len(hits),
        "hits": hits,
    }

@app.get("/api/jobs/{job_id}/hits/download")
def download_job_hits(job_id: str):
    document = jobs.find_one(
        {"_id": job_id},
        {"_id": 0, "status": 1},
    )

    if not document:
        raise HTTPException(
            status_code=404,
            detail="Job not found.",
        )

    if document.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Screening job is not completed.",
        )

    job_dir = (OUTPUT_ROOT / job_id).resolve()
    output_root = OUTPUT_ROOT.resolve()

    if not job_dir.is_relative_to(output_root):
        raise HTTPException(
            status_code=400,
            detail="Invalid output path.",
        )

    if not job_dir.is_dir():
        raise HTTPException(
            status_code=404,
            detail="Job output directory not found.",
        )

    csv_files = sorted(job_dir.glob("*/hits.csv"))

    if not csv_files:
        raise HTTPException(
            status_code=404,
            detail="No gene-level hit files found.",
        )

    output = io.StringIO(newline="")
    fieldnames = [
        "gene",
        "organism",
        "transcript",
        "transcript_pos_1based",
        "gene_start",
        "gene_end",
        "length",
        "strand",
    ]

    writer = csv.DictWriter(
        output,
        fieldnames=fieldnames,
        extrasaction="ignore",
    )
    writer.writeheader()

    total_hits = 0

    for csv_path in csv_files:
        gene_name = csv_path.parent.name.split("_accession_")[0]

        with csv_path.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)

            for row in reader:
                if not row or not any(
                    value and value.strip()
                    for value in row.values()
                ):
                    continue

                writer.writerow({
                    "gene": gene_name,
                    **row,
                })
                total_hits += 1

    if total_hits == 0:
        raise HTTPException(
            status_code=404,
            detail="No detailed hits were found for this job.",
        )

    csv_content = output.getvalue()
    output.close()

    return StreamingResponse(
        iter([csv_content]),
        media_type="text/csv",
        headers={
            "Content-Disposition":
                f'attachment; filename="dsguard-hits-{job_id}.csv"'
        },
    )


    if not document:
        raise HTTPException(
            status_code=404,
            detail="Job not found.",
        )

    if document.get("status") != "completed":
        raise HTTPException(
            status_code=409,
            detail="Screening job is not completed.",
        )

    hits_path = (OUTPUT_ROOT / job_id / "hits.csv").resolve()
    output_root = OUTPUT_ROOT.resolve()

    if not hits_path.is_relative_to(output_root):
        raise HTTPException(
            status_code=400,
            detail="Invalid output path.",
        )

    if not hits_path.is_file():
        raise HTTPException(
            status_code=404,
            detail="Detailed hits file was not generated for this job.",
        )

    return FileResponse(
        path=str(hits_path),
        media_type="text/csv",
        filename=f"dsguard-hits-{job_id}.csv",
    )
