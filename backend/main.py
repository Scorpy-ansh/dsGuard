

import csv
import hashlib

import io

import json

import logging

import os

import re

import shutil

import subprocess

import sys

import threading

import uuid

import urllib.parse

import urllib.error

import urllib.request

from datetime import datetime, timezone

from pathlib import Path

from typing import Annotated



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

from fastapi.responses import StreamingResponse

from pymongo import MongoClient

from pymongo.errors import PyMongoError


logger = logging.getLogger("dsguard.backend")


try:

    from backend.reporting import (

        build_candidate_fasta,

        build_offtarget_csv,

        build_pdf_report,

        build_results_json,

    )

except ImportError:

    from reporting import (

        build_candidate_fasta,

        build_offtarget_csv,

        build_pdf_report,

        build_results_json,

    )













ROOT = Path(__file__).resolve().parent.parent

ENGINE_DIR = ROOT / "Random"

load_dotenv(ROOT / ".env")

def engine_command() -> list[str] | None:
    interpreter = Path(sys.executable)
    if interpreter.is_file():
        return [str(interpreter), "-m", "dsguard"]
    return None





UPLOAD_ROOT = ROOT / "backend" / "uploads"

OUTPUT_ROOT = ROOT / "outputs" / "jobs"

TARGET_CATALOG_PATH = ROOT / "Random" / "data" / "targets" / "target_catalog.json"



PANEL_DIR = Path(os.getenv("DSGUARD_PANEL_DIR", str(ROOT / "reference_panels"))).resolve()
PANEL_CACHE_DIR = Path(
    os.getenv("DSGUARD_PANEL_CACHE_DIR", str(ROOT / ".cache" / "reference_panels"))
).resolve()
PANEL_CATALOG = {

    "rice": {
        "label": "Rice (Oryza sativa)",
        "organism": "Oryza sativa",
        "filename": "rice.fasta",
    },

    "cyrtorhinus_lividipennis": {
        "label": "Cyrtorhinus lividipennis",
        "organism": "Cyrtorhinus lividipennis",
        "filename": "cyrtorhinus_lividipennis.fasta",
    },

    "honeybee": {
        "label": "Honeybee (Apis mellifera)",
        "organism": "Apis mellifera",
        "filename": "honeybee.fasta",
    },

    "spider": {
        "label": "Spider (reference species panel)",
        "organism": "reference spider species",
        "filename": "spider.fasta",
    },

}

_PANEL_STATE_LOCK = threading.Lock()
_PANEL_STATES: dict[str, dict] = {}


def _panel_state_key(name: str) -> tuple[str, str, str, str, str]:
    return (
        name,
        str(PANEL_DIR),
        str(PANEL_CACHE_DIR),
        os.getenv("RICE_PANEL_URL", ""),
        os.getenv("RICE_PANEL_SHA256", ""),
    )


def _validate_fasta(path: Path) -> None:
    record_count = 0
    sequence_length = 0
    header_seen = False

    with path.open("rb") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(b">"):
                if len(line) == 1:
                    raise ValueError(f"FASTA header is missing a name at line {line_number}.")
                if header_seen and sequence_length == 0:
                    raise ValueError(f"FASTA record has no sequence before line {line_number}.")
                header_seen = True
                record_count += 1
                sequence_length = 0
                continue
            if not header_seen:
                raise ValueError(
                    f"FASTA sequence data appears before the first header at line {line_number}."
                )
            sequence_length += len(line.replace(b" ", b"").replace(b"\t", b""))

    if not header_seen or sequence_length == 0:
        raise ValueError("FASTA contains no non-empty sequence records.")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _url_label(url: str) -> str:
    parsed = urllib.parse.urlsplit(url)
    filename = Path(parsed.path).name or "/"
    return f"{parsed.hostname or '<invalid-host>'}/{filename}"


def _sanitised_response_body(body: bytes) -> str:
    text = body[:500].decode("utf-8", "replace")
    return "".join(character if character.isprintable() else " " for character in text)


def _download_rice_panel() -> Path | None:
    url = os.getenv("RICE_PANEL_URL", "").strip()
    expected_sha256 = os.getenv("RICE_PANEL_SHA256", "").strip().lower()
    if not url or not expected_sha256:
        logger.warning(
            "Rice panel configuration is incomplete: RICE_PANEL_URL and "
            "RICE_PANEL_SHA256 are required."
        )
        return None
    if len(expected_sha256) != 64 or any(
        character not in "0123456789abcdef" for character in expected_sha256
    ):
        raise RuntimeError("RICE_PANEL_SHA256 must be a 64-character hexadecimal SHA-256 checksum.")
    if not url.startswith(("https://", "http://")):
        raise RuntimeError("RICE_PANEL_URL must use HTTP or HTTPS.")

    request = urllib.request.Request(
        url,
        headers={
            "Accept": "*/*",
            "User-Agent": "dsGuard/0.4",
        },
        method="GET",
    )
    request_label = _url_label(url)
    logger.info(
        "Rice panel download started: host_file=%s method=%s",
        request_label,
        request.method,
    )
    PANEL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    destination = PANEL_CACHE_DIR / PANEL_CATALOG["rice"]["filename"]
    if destination.is_file():
        actual_sha256 = _sha256(destination)
        if actual_sha256 == expected_sha256:
            try:
                _validate_fasta(destination)
            except (OSError, ValueError):
                logger.warning("Rice panel cached FASTA validation failed: path=%s", destination)
                destination.unlink(missing_ok=True)
            else:
                logger.info("Rice panel checksum verification passed: sha256=%s", actual_sha256)
                logger.info("Rice panel cache reused: path=%s bytes=%d", destination, destination.stat().st_size)
                return destination
        else:
            logger.error(
                "Rice panel checksum verification failed: expected=%s actual=%s",
                expected_sha256,
                actual_sha256,
            )

    temporary = PANEL_CACHE_DIR / f".{destination.name}.{os.getpid()}.tmp"
    try:
        try:
            response_context = urllib.request.urlopen(request, timeout=60)
            with response_context as response, temporary.open("wb") as output:
                status = getattr(response, "status", None)
                final_url = response.url
                content_type = response.headers.get("Content-Type", "")
                redirected = final_url != url
                logger.info(
                    "Rice panel HTTP response: status=%s final_host_file=%s "
                    "content_type=%s method=%s redirected=%s",
                    status,
                    _url_label(final_url),
                    content_type,
                    request.method,
                    redirected,
                )
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    output.write(chunk)
        except urllib.error.HTTPError as exc:
            response_url = exc.url
            body = _sanitised_response_body(exc.read())
            logger.error(
                "Rice panel HTTP response: status=%s final_host_file=%s "
                "content_type=%s method=%s redirected=%s body=%r",
                exc.code,
                _url_label(response_url),
                exc.headers.get("Content-Type", ""),
                request.method,
                response_url != url,
                body,
            )
            raise RuntimeError(
                f"Rice panel server returned HTTP {exc.code} for {request_label}."
            ) from exc
        byte_size = temporary.stat().st_size
        logger.info("Rice panel download completed: bytes=%d", byte_size)
        actual_sha256 = _sha256(temporary)
        if actual_sha256 != expected_sha256:
            logger.error(
                "Rice panel checksum verification failed: expected=%s actual=%s",
                expected_sha256,
                actual_sha256,
            )
            raise RuntimeError("Downloaded Rice panel failed the RICE_PANEL_SHA256 check.")
        logger.info("Rice panel checksum verification passed: sha256=%s", actual_sha256)
        _validate_fasta(temporary)
        logger.info("Rice panel FASTA validation passed: path=%s", temporary)
        temporary.replace(destination)
        logger.info("Rice panel installed: path=%s bytes=%d", destination, byte_size)
        return destination
    except Exception as exc:
        temporary.unlink(missing_ok=True)
        logger.error("Rice panel reference-loading failure: %s", exc)
        if isinstance(exc, RuntimeError):
            raise
        raise RuntimeError(f"Could not download the Rice reference panel: {exc}") from exc


def panel_path(name: str) -> Path | None:
    config = PANEL_CATALOG[name]
    local_path = PANEL_DIR / config["filename"]
    if local_path.is_file() and local_path.stat().st_size > 0:
        return local_path
    if name == "rice":
        return _download_rice_panel()
    return None


def _panel_status(name: str, path: Path | None, reason: str | None = None) -> dict:
    config = PANEL_CATALOG[name]
    exists = bool(path and path.is_file())
    size = path.stat().st_size if exists else 0
    status = _PANEL_STATES.get(name, {})
    return {
        "key": name,
        "name": config["label"],
        "label": config["label"],
        "organism": config["organism"],
        "installed": bool(status.get("installed", exists and size > 0)),
        "available": bool(status.get("installed", exists and size > 0)),
        "file_exists": exists,
        "validation_status": status.get("validation_status", "not_checked"),
        "last_validation_result": status.get("last_validation_result", reason),
        "filename": config["filename"],
        "path": str(path) if path else None,
        "size_bytes": size,
        "reason": reason,
    }


def ensure_panel(name: str) -> Path | None:
    with _PANEL_STATE_LOCK:
        existing = _PANEL_STATES.get(name)
        if existing is not None and existing.get("state_key") == _panel_state_key(name):
            return existing.get("path")

        config = PANEL_CATALOG[name]
        candidate = PANEL_DIR / config["filename"]
        if name != "rice":
            try:
                installed = candidate.is_file() and candidate.stat().st_size > 0
                if installed:
                    _validate_fasta(candidate)
                validation_status = "passed" if installed else "not_checked"
                validation_result = "FASTA validation passed" if installed else None
            except (OSError, ValueError) as exc:
                installed = False
                validation_status = "failed"
                validation_result = str(exc)
            _PANEL_STATES[name] = {
                "path": candidate if (installed or candidate.exists()) else None,
                "installed": installed,
                "validation_status": validation_status,
                "last_validation_result": validation_result,
                "state_key": _panel_state_key(name),
            }
            logger.info("Panel availability: name=%s installed=%s path=%s", name, installed, candidate)
            return candidate if installed else None

        try:
            path = candidate if candidate.is_file() and candidate.stat().st_size > 0 else None
            if path is None:
                cached = PANEL_CACHE_DIR / config["filename"]
                path = cached if cached.is_file() and cached.stat().st_size > 0 else None
            if path is None:
                path = _download_rice_panel()
            if path is None:
                reason = "Rice panel configuration is missing."
                _PANEL_STATES[name] = {
                    "path": path if path and path.exists() else None,
                    "installed": False,
                    "validation_status": "not_checked",
                    "last_validation_result": reason,
                    "state_key": _panel_state_key(name),
                }
                logger.warning("Panel availability: name=%s installed=false reason=%s", name, reason)
                return None

            _validate_fasta(path)
            _PANEL_STATES[name] = {
                "path": path,
                "installed": True,
                "validation_status": "passed",
                "last_validation_result": "FASTA validation passed",
                "state_key": _panel_state_key(name),
            }
            logger.info("Panel availability: name=%s installed=true path=%s", name, path)
            return path
        except (OSError, ValueError, RuntimeError) as exc:
            reason = str(exc)
            _PANEL_STATES[name] = {
                "path": path if path and path.exists() else candidate if candidate.exists() else None,
                "installed": False,
                "validation_status": "failed",
                "last_validation_result": reason,
                "state_key": _panel_state_key(name),
            }
            logger.error("Reference-loading failure: name=%s reason=%s", name, reason)
            return None


def load_target_catalog() -> dict:

    try:

        catalog = json.loads(TARGET_CATALOG_PATH.read_text(encoding="utf-8"))

    except (OSError, json.JSONDecodeError) as exc:

        raise RuntimeError(f"Target catalog could not be loaded: {exc}") from exc

    organisms = catalog.get("organisms")

    if not isinstance(organisms, list):

        raise RuntimeError("Target catalog must contain an organisms list.")

    for organism in organisms:

        if not isinstance(organism, dict) or not organism.get("key") or not organism.get("scientific_name"):

            raise RuntimeError("Target catalog contains an invalid organism record.")

        if not isinstance(organism.get("genes"), list):

            raise RuntimeError(f"Target catalog organism {organism['key']} has no genes list.")

        for gene in organism["genes"]:

            if not isinstance(gene, dict) or not gene.get("gene_name"):

                raise RuntimeError(f"Target catalog organism {organism['key']} contains an invalid gene.")

            if gene.get("sequence") is not None:

                raise RuntimeError("Target catalog must not store nucleotide sequences.")

    return catalog





TARGET_CATALOG = load_target_catalog()

EXPECTED_TARGET_ORGANISMS = {

    "sogatella_furcifera",

    "chilo_suppressalis",

    "nilaparvata_lugens",

}

if {str(item.get("key", "")).lower() for item in TARGET_CATALOG["organisms"]} != EXPECTED_TARGET_ORGANISMS:

    raise RuntimeError(

        "Target catalog must contain exactly Sogatella furcifera, "

        "Chilo suppressalis, and Nilaparvata lugens."

    )





















WINDOW_MIN = 21

WINDOW_MAX = 501



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

)


def cors_origins(frontend_origin=None):
    origins = [
        "http://localhost:5173",
        "http://127.0.0.1:5173",
    ]
    configured_origin = (
        os.getenv("FRONTEND_ORIGIN")
        if frontend_origin is None
        else frontend_origin
    )
    if configured_origin:
        normalized_origin = configured_origin.strip().rstrip("/")
        if normalized_origin and normalized_origin not in origins:
            origins.append(normalized_origin)
    return origins


app.add_middleware(
    CORSMiddleware,
    allow_origins=cors_origins(),





    allow_origin_regex=r"^http://(localhost|127\.0\.0\.1)(:\d+)?$",

    allow_credentials=True,

    allow_methods=["*"],

    allow_headers=["*"],

)



job_lock = threading.Lock()













def now_utc():

    return datetime.now(timezone.utc).isoformat()





def safe_panel_name(value: str) -> str:

    return re.sub(

        r"[^A-Za-z0-9_-]",

        "_",

        value.strip(),

    )[:40]





def save_upload(upload: UploadFile, destination: Path) -> int:



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













def run_screening(

    job_id: str,

    target_path: Path,

    panel_paths: list[tuple[str, Path]],

    params: dict,

    input_type: str = "fasta_file",

    input_label: str = "",

    skipped_panels: list[tuple[str, str]] | None = None,

):



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

        if not panel_paths:

            raise RuntimeError(

                "No reference panel was available for this job. Nothing was screened."

            )



        engine = engine_command()

        if engine is None:

            raise RuntimeError(

                "The Python interpreter running the backend is not available."

            )



        job_output.mkdir(parents=True, exist_ok=True)



        command = engine + [

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

            "--input-type",

            str(input_type),

            "--input-label",

            str(input_label or target_path.name),

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



        for panel_name, reason in (skipped_panels or []):

            command.extend(

                [

                    "--skipped-panel",

                    f"{panel_name}={reason}",

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

        results_path = job_output / "results.json"



        if not summary_path.is_file():

            raise RuntimeError("dsGuard did not produce summary.csv.")



        if not manifest_path.is_file():

            raise RuntimeError("dsGuard did not produce manifest.json.")



        if not results_path.is_file():

            raise RuntimeError(

                "dsGuard did not produce the structured results object "

                "(results.json). Downloadable reports cannot be generated for "

                "this job; re-run the screening."

            )



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

                    "results_object": str(results_path) if results_path.is_file() else None,

                    "skipped_panels": [

                        {"name": name, "reason": reason}

                        for name, reason in (skipped_panels or [])

                    ],

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

        "engine_exists": engine_command() is not None,

        "engine_command": engine_command(),

        "timestamp": now_utc(),

    }













@app.get("/api/panels")

def list_reference_panels():

    panels = []

    for key, config in PANEL_CATALOG.items():
        path = ensure_panel(key)
        state_path = _PANEL_STATES.get(key, {}).get("path")
        if path is None:
            reason = _PANEL_STATES.get(key, {}).get("last_validation_result")
            if not reason:
                reason = (
                    f"{PANEL_DIR / config['filename']} was not found or is empty."
                )
        else:
            reason = None
        panels.append(_panel_status(key, path or state_path, reason))

    return panels





@app.get("/api/targets")

def list_target_catalog():

    return {

        "organisms": [

            {

                "key": organism["key"],

                "scientific_name": organism["scientific_name"],

                "taxonomy_id": organism.get("taxonomy_id"),

                "genes": organism["genes"],

            }

            for organism in TARGET_CATALOG["organisms"]

        ]

    }





_NCBI_HEADERS = {"User-Agent": "dsGuard/0.4 contact: dsguard-local-app"}







_NCBI_PREFIX_WEIGHTS = (

    ("NM_", 12.0),

    ("NR_", 10.0),

    ("XM_", 8.0),

    ("XR_", 6.0),

    ("NC_", -12.0),

    ("NW_", -12.0),

    ("NG_", -8.0),

    ("NZ_", -6.0),

)





def ncbi_gene_search_score(query: str, accession: str, title: str, length: int) -> float:

    score = 0.0

    for prefix, weight in _NCBI_PREFIX_WEIGHTS:

        if accession.startswith(prefix):

            score += weight

            break

    if query.lower() in (title or "").lower():

        score += 15.0

    else:





        score -= 25.0

    if 100 <= length <= 20_000:

        score += 2.0

    elif length > 20_000:





        score -= 4.0

    return score





def choose_ncbi_record(query: str, records: list[dict]) -> dict | None:

    best, best_key = None, None

    for index, record in enumerate(records):

        key = (

            ncbi_gene_search_score(

                query,

                record.get("accession", ""),

                record.get("title", ""),

                int(record.get("length") or 0),

            ),

            -index,

        )

        if best_key is None or key > best_key:

            best, best_key = record, key

    return best





def parse_ncbi_summary(payload: dict) -> list[dict]:

    result = payload.get("result", {})

    records = []

    for uid in result.get("uids", []):

        item = result.get(uid) or {}

        records.append({

            "uid": uid,

            "accession": item.get("accessionversion") or item.get("caption", ""),

            "title": item.get("title", ""),

            "length": item.get("slen", 0),

        })

    return records





_NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"





def _ncbi_get(url: str, timeout: int) -> bytes:

    request = urllib.request.Request(url, headers=_NCBI_HEADERS)

    with urllib.request.urlopen(request, timeout=timeout) as response:

        return response.read()





def _ncbi_search(term: str, retmax: int) -> list[str]:

    url = _NCBI_EUTILS + "esearch.fcgi?" + urllib.parse.urlencode({

        "db": "nuccore",

        "term": term,

        "retmode": "json",

        "retmax": retmax,

        "sort": "relevance",

    })

    payload = json.loads(_ncbi_get(url, 20).decode("utf-8"))

    return payload.get("esearchresult", {}).get("idlist", [])













PROJECT_TARGET_FASTA = (

    ROOT / "Random" / "data" / "targets" / "nilaparvata_lugens_13_cds.fasta"

)





def project_target_accessions() -> dict[str, str]:

    mapping: dict[str, str] = {}

    if not PROJECT_TARGET_FASTA.is_file():

        return mapping



    try:

        text = PROJECT_TARGET_FASTA.read_text(encoding="utf-8", errors="replace")

    except OSError:

        return mapping

    for line in text.splitlines():

        if not line.startswith(">"):

            continue

        header = line[1:].strip()

        alias = header.split("|", 1)[0].strip()

        match = re.search(r"accession=([A-Za-z0-9_.]+)", header)

        if alias and match:

            mapping.setdefault(alias.lower(), match.group(1))

    return mapping





def catalog_target_accessions() -> dict[str, str]:

    return {

        gene["gene_name"].lower(): gene["accession"]

        for organism in TARGET_CATALOG["organisms"]

        for gene in organism["genes"]

        if gene.get("accession") and gene.get("accession_status") == "verified"

    }





def catalog_gene(query: str, organism: str = "") -> dict | None:

    lowered = query.strip().lower()

    organism_lower = organism.strip().lower()

    for record in TARGET_CATALOG["organisms"]:

        if organism_lower and record["scientific_name"].lower() != organism_lower:

            continue

        for gene in record["genes"]:

            if gene["gene_name"].lower() == lowered:

                return gene

    return None





def fasta_header_line(fasta: str) -> str:

    for line in (fasta or "").splitlines():

        if line.startswith(">"):

            return line[1:].strip()

    return ""





def organism_is_annotated(header_line: str, organism: str) -> bool:

    return organism.strip().lower() in (header_line or "").lower()





def fetch_ncbi_fasta(query: str, organism: str = "") -> tuple[str, str]:

    query = query.strip()

    if not query:

        raise ValueError("Enter an NCBI accession or gene name.")

    base = _NCBI_EUTILS

    alias_accession = {

        **catalog_target_accessions(),

        **project_target_accessions(),

    }.get(query.lower())

    catalog_record = catalog_gene(query, organism)

    if catalog_record and catalog_record.get("accession_status") != "verified":

        raise ValueError(

            f"Catalog target '{query}' has no verified accession and cannot be "

            "screened until its sequence is supplied or resolved from NCBI."

        )

    is_accession = bool(re.fullmatch(r"(?:[A-Za-z]{1,6}_)?[A-Za-z]{0,8}\d{3,}(?:\.\d+)?", query))

    if alias_accession:

        term = alias_accession

        is_accession = True

    elif is_accession:

        term = query

    else:

        term = f'("{query}"[Gene Name] OR "{query}"[All Fields]) AND "{organism or "Insecta"}"[Organism]'



    try:

        if is_accession:

            ids = _ncbi_search(term, 5)

        else:







            ids = _ncbi_search(f'{term} AND "refseq"[Filter]', 20) or _ncbi_search(term, 20)

        if not ids:

            raise ValueError(

                f"NCBI found no nucleotide records for '{query}'. Try an accession ID "

                "or include a more specific gene name/organism."

            )



        chosen_id = ids[0]

        if not is_accession:

            summary_url = base + "esummary.fcgi?" + urllib.parse.urlencode({

                "db": "nuccore",

                "id": ",".join(ids),

                "retmode": "json",

            })

            records = parse_ncbi_summary(json.loads(_ncbi_get(summary_url, 25).decode("utf-8")))

            chosen = choose_ncbi_record(query, records)

            if chosen:

                chosen_id = chosen["uid"]



        fetch_url = base + "efetch.fcgi?" + urllib.parse.urlencode({

            "db": "nuccore",

            "id": chosen_id,

            "rettype": "fasta",

            "retmode": "text",

        })

        fasta = _ncbi_get(fetch_url, 30).decode("utf-8").strip()

    except ValueError:

        raise

    except Exception as exc:

        raise ValueError(f"Could not retrieve sequence from NCBI: {exc}") from exc

    if not fasta.startswith(">") or not re.search(r"\n[ACGTUNacgtun\s]+$", fasta):

        raise ValueError("NCBI did not return a valid nucleotide FASTA sequence.")



    header = fasta_header_line(fasta)

    record_id = header.split()[0] if header else chosen_id

    if organism and not organism_is_annotated(header, organism):

        raise ValueError(

            f"NCBI record {record_id} is not annotated for '{organism}'. Its "

            f"definition line is '{header}'. Refusing to screen an unrelated "

            "record; check the accession/organism combination."

        )

    return fasta + "\n", record_id













@app.post("/api/screen")

async def create_screening_job(

    background_tasks: BackgroundTasks,

    panel_key: Annotated[str, Form()] = "",

    panel_keys: Annotated[str, Form()] = "",

    target_mode: Annotated[str, Form()] = "file",

    target_file: Annotated[UploadFile | None, File()] = None,

    target_query: Annotated[str, Form()] = "",

    target_sequence: Annotated[str, Form()] = "",

    target_organism: Annotated[str, Form()] = "",

    window: Annotated[int, Form()] = 400,

    k: Annotated[int, Form()] = 11,

    l50: Annotated[float, Form()] = 21.0,

    kappa: Annotated[float, Form()] = 2.0,

    max_hazard: Annotated[float, Form()] = 0.05,

):

    if engine_command() is None:

        raise HTTPException(

            status_code=500,

            detail=(

                "The dsGuard CLI was not found because the Python interpreter "
                "running the backend is not available, so screening cannot run."

            ),

        )







    requested = []

    for raw in (panel_keys.split(",") if panel_keys else []) + ([panel_key] if panel_key else []):

        name = raw.strip()

        if name and name not in requested:

            requested.append(name)

    if not requested:

        raise HTTPException(status_code=400, detail="Choose at least one non-target organism panel.")



    unknown = [name for name in requested if name not in PANEL_CATALOG]

    if unknown:

        raise HTTPException(

            status_code=400,

            detail="Unsupported non-target panel(s): " + ", ".join(unknown) + ".",

        )



    if not WINDOW_MIN <= window <= WINDOW_MAX:

        raise HTTPException(

            status_code=400,

            detail=(

                f"Window must be a whole number between {WINDOW_MIN} and {WINDOW_MAX} bp "

                f"(got {window})."

            ),

        )

    if not 2 <= k <= 31:

        raise HTTPException(status_code=400, detail="k must be between 2 and 31.")

    if l50 <= 0 or kappa <= 0:

        raise HTTPException(status_code=400, detail="l50 and kappa must be positive.")

    if not 0 <= max_hazard <= 1:

        raise HTTPException(status_code=400, detail="max_hazard must be between 0 and 1.")







    panel_paths: list[tuple[str, Path]] = []

    skipped_panels: list[tuple[str, str]] = []

    for name in requested:
        config = PANEL_CATALOG[name]
        reason = None
        try:
            path = ensure_panel(name)
        except RuntimeError as exc:
            path = None
            reason = str(exc)
        if path is None:
            if name == "rice" and reason is None:
                reason = (
                    "Rice panel is unavailable. Set RICE_PANEL_URL and "
                    "RICE_PANEL_SHA256 or install rice.fasta locally."
                )
            skipped_panels.append((
                name,
                reason or (
                    f"{PANEL_DIR / config['filename']} was not found or is empty"
                ),
            ))
        else:
            panel_paths.append((name, path))



    if not panel_paths:

        raise HTTPException(

            status_code=400,

            detail=(

                "None of the requested reference panels are installed, so nothing "

                "could be screened. "

                + " ".join(f"{PANEL_CATALOG[n]['label']}: {reason}." for n, reason in skipped_panels)

                + f" Add the FASTA files to {PANEL_DIR} or set DSGUARD_PANEL_DIR."

            ),

        )



    job_id = str(uuid.uuid4())

    job_upload_dir = UPLOAD_ROOT / job_id

    job_output_dir = OUTPUT_ROOT / job_id

    target_path = job_upload_dir / "target.fasta"

    target_label = "target.fasta"

    input_type = "fasta_file"

    try:

        if target_mode == "file":

            if target_file is None:

                raise ValueError("Choose a target FASTA file.")

            save_upload(target_file, target_path)

            target_label = Path(target_file.filename or "target.fasta").name

            input_type = "fasta_file"

        elif target_mode == "sequence":

            sequence = target_sequence.strip()

            if not sequence:

                raise ValueError("Paste a target FASTA or nucleotide sequence.")

            if not sequence.startswith(">"):

                sequence = ">pasted_target\n" + sequence

            sequence = sequence.replace("\r", "")

            lines = sequence.splitlines()

            if len(lines) < 2 or not lines[0].startswith(">"):

                raise ValueError("Sequence must contain a FASTA header and nucleotide sequence.")

            bases = "".join(lines[1:]).replace(" ", "").replace("\t", "")

            if not bases or re.search(r"[^ACGTUNacgtun]", bases):

                raise ValueError("Target sequence contains invalid characters. Use nucleotide bases A, C, G, T, U, or N.")

            job_upload_dir.mkdir(parents=True, exist_ok=True)

            target_path.write_text(lines[0] + "\n" + bases + "\n", encoding="utf-8")

            target_label = lines[0][1:].strip() or "pasted_target"

            input_type = "pasted_sequence"

        elif target_mode == "ncbi":

            fasta, record_id = fetch_ncbi_fasta(target_query, target_organism)

            job_upload_dir.mkdir(parents=True, exist_ok=True)

            target_path.write_text(fasta, encoding="utf-8")

            target_label = f"NCBI:{record_id} ({target_query.strip()})"

            input_type = "ncbi"

        else:

            raise ValueError("Unsupported target input mode.")

    except ValueError as exc:

        shutil.rmtree(job_upload_dir, ignore_errors=True)

        raise HTTPException(status_code=400, detail=str(exc)) from exc

    finally:

        if target_file is not None:

            await target_file.close()



    params = {"window": window, "k": k, "l50": l50, "kappa": kappa, "max_hazard": max_hazard}

    document = {

        "_id": job_id,

        "status": "queued",

        "created_at": now_utc(),

        "parameters": params,

        "panel_names": [name for name, _ in panel_paths],

        "skipped_panels": [{"name": n, "reason": r} for n, r in skipped_panels],

        "target_filename": target_label,

        "target_input_type": input_type,

        "target_organism": target_organism.strip(),

        "panel_filenames": [PANEL_CATALOG[n]["filename"] for n, _ in panel_paths],

        "upload_directory": str(job_upload_dir),

        "result_directory": str(job_output_dir),

    }

    try:

        jobs.insert_one(document)

    except PyMongoError as exc:

        shutil.rmtree(job_upload_dir, ignore_errors=True)

        raise HTTPException(status_code=503, detail="Could not save job metadata to MongoDB.") from exc



    background_tasks.add_task(

        run_screening, job_id, target_path, panel_paths, params,

        input_type, target_label, skipped_panels,

    )

    return {

        "job_id": job_id,

        "status": "queued",

        "message": "Screening job accepted.",

        "target": target_label,

        "panels": [PANEL_CATALOG[name]["label"] for name, _ in panel_paths],

        "skipped_panels": [{"name": n, "reason": r} for n, r in skipped_panels],

        "status_url": f"/api/jobs/{job_id}",

    }















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





        if (

            upload_dir.parent != UPLOAD_ROOT.resolve()

            or output_dir.parent != OUTPUT_ROOT.resolve()

        ):

            raise HTTPException(

                status_code=400,

                detail="Invalid job path.",

            )





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





@app.delete("/api/jobs")

def clear_job_history():

    with job_lock:

        documents = list(jobs.find({}, {"_id": 1, "status": 1}))

        active = [doc["_id"] for doc in documents

                  if doc.get("status") in {"queued", "running"}]

        if active:

            raise HTTPException(

                status_code=409,

                detail=(

                    "Cannot clear history while screenings are queued or running. "

                    "Wait for them to finish and try again."

                ),

            )



        roots = (UPLOAD_ROOT.resolve(), OUTPUT_ROOT.resolve())

        directories = []

        for document in documents:

            job_id = str(document["_id"])

            upload_dir = (UPLOAD_ROOT / job_id).resolve()

            output_dir = (OUTPUT_ROOT / job_id).resolve()

            if upload_dir.parent != roots[0] or output_dir.parent != roots[1]:

                raise HTTPException(status_code=400, detail="Invalid job path.")

            directories.extend((upload_dir, output_dir))



        try:

            for directory in directories:

                if directory.exists():

                    shutil.rmtree(directory)

        except OSError as exc:

            raise HTTPException(

                status_code=500,

                detail="Could not remove screening files. History records were retained.",

            ) from exc



        try:

            result = jobs.delete_many({})

        except PyMongoError as exc:

            raise HTTPException(

                status_code=500,

                detail="Files were removed, but screening history could not be cleared.",

            ) from exc



    return {

        "status": "deleted",

        "deleted_count": result.deleted_count,

        "message": "All application-scoped screening history and owned artifacts were deleted.",

    }







def load_job_results_object(job_id: str) -> dict | None:

    output_root = OUTPUT_ROOT.resolve()

    job_dir = (OUTPUT_ROOT / job_id).resolve()



    if job_dir.parent != output_root:

        raise HTTPException(status_code=400, detail="Invalid job path.")



    results_path = job_dir / "results.json"

    if not results_path.is_file():

        return None



    try:

        return json.loads(results_path.read_text(encoding="utf-8"))

    except (OSError, json.JSONDecodeError) as exc:

        raise HTTPException(

            status_code=500,

            detail=f"Stored analysis results could not be read: {exc}",

        ) from exc





def require_completed_results(job_id: str) -> dict:

    document = jobs.find_one({"_id": job_id}, {"_id": 0, "status": 1})



    if not document:

        raise HTTPException(status_code=404, detail="Job not found.")



    if document.get("status") != "completed":

        raise HTTPException(

            status_code=409,

            detail="Screening job is not completed.",

        )



    results_object = load_job_results_object(job_id)



    if results_object is None:

        raise HTTPException(

            status_code=404,

            detail=(

                "This job's structured results (results.json) were not found on "

                "disk. Jobs created before structured results were introduced do "

                "not have one, and the missing analysis cannot be fabricated or "

                "substituted from another run. Re-run the screening from its "

                "original target and reference panels to regenerate downloadable "

                "reports."

            ),

        )



    return results_object





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

            "target_organism": 1,

        },

    )



    if not document:

        raise HTTPException(

            status_code=404,

            detail="Job not found.",

        )



    document["analysis"] = load_job_results_object(job_id)

    if document["analysis"] is not None and document.get("target_organism"):

        document["analysis"].setdefault("input", {})["organism"] = document["target_organism"]



    return document













def _safe_download_name(job_id: str, suffix: str) -> str:

    return f"dsguard-{safe_panel_name(job_id)[:36]}-{suffix}"





@app.get("/api/jobs/{job_id}/export/results.csv")

def export_offtarget_csv(job_id: str):

    results = require_completed_results(job_id)



    try:

        content = build_offtarget_csv(results)

    except Exception as exc:

        raise HTTPException(

            status_code=500,

            detail=f"Could not generate the CSV export: {exc}",

        ) from exc



    return StreamingResponse(

        iter([content]),

        media_type="text/csv",

        headers={

            "Content-Disposition":

                f'attachment; filename="{_safe_download_name(job_id, "offtargets.csv")}"'

        },

    )





@app.get("/api/jobs/{job_id}/export/candidates.fasta")

def export_candidate_fasta(job_id: str):

    results = require_completed_results(job_id)



    try:

        content = build_candidate_fasta(results)

    except Exception as exc:

        raise HTTPException(

            status_code=500,

            detail=f"Could not generate the FASTA export: {exc}",

        ) from exc



    if not content:





        raise HTTPException(

            status_code=404,

            detail=(

                "This analysis did not record a candidate dsRNA sequence for any "

                "target gene, so there is no FASTA to download. Each gene's "

                "candidate_error states why no candidate was generated. Re-run "

                "the screening with a shorter window or a longer target."

            ),

        )



    return StreamingResponse(

        iter([content]),

        media_type="text/x-fasta",

        headers={

            "Content-Disposition":

                f'attachment; filename="{_safe_download_name(job_id, "candidates.fasta")}"'

        },

    )





@app.get("/api/jobs/{job_id}/export/results.json")

def export_results_json(job_id: str):

    results = require_completed_results(job_id)



    try:

        content = build_results_json(results)

    except Exception as exc:

        raise HTTPException(

            status_code=500,

            detail=f"Could not generate the JSON export: {exc}",

        ) from exc



    return StreamingResponse(

        iter([content]),

        media_type="application/json",

        headers={

            "Content-Disposition":

                f'attachment; filename="{_safe_download_name(job_id, "results.json")}"'

        },

    )





@app.get("/api/jobs/{job_id}/export/report.pdf")

def export_pdf_report(job_id: str):

    results = require_completed_results(job_id)



    try:

        pdf_bytes = build_pdf_report(results)

    except Exception as exc:

        raise HTTPException(

            status_code=500,

            detail=f"Could not generate the PDF report: {exc}",

        ) from exc



    return StreamingResponse(

        iter([pdf_bytes]),

        media_type="application/pdf",

        headers={

            "Content-Disposition":

                f'attachment; filename="{_safe_download_name(job_id, "report.pdf")}"'

        },

    )







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





@app.get("/api/jobs/{job_id}/profiles/{gene_name}")

def get_job_profiles(job_id: str, gene_name: str):

    results = require_completed_results(job_id)

    gene = next((item for item in results.get("genes", [])

                 if item.get("gene") == gene_name), None)

    if gene is None:

        raise HTTPException(status_code=404, detail="Target gene not found.")



    job_dir = (OUTPUT_ROOT / job_id).resolve()

    output_root = OUTPUT_ROOT.resolve()

    if not job_dir.is_relative_to(output_root):

        raise HTTPException(status_code=400, detail="Invalid output path.")



    gene_dir = next(

        (path for path in job_dir.iterdir()

         if path.is_dir() and path.name == re.sub(r"[^A-Za-z0-9._-]+", "_", gene_name)[:60]),

        None,

    )

    if gene_dir is None:





        matches = [path for path in job_dir.iterdir()

                   if path.is_dir() and path.name.startswith(

                       re.sub(r"[^A-Za-z0-9._-]+", "_", gene_name)[:60])]

        gene_dir = matches[0] if len(matches) == 1 else None

    if gene_dir is None:

        raise HTTPException(status_code=404, detail="Positional data unavailable.")



    windows_path = gene_dir / "windows.csv"

    profiles_path = gene_dir / "profiles.csv"

    if not windows_path.is_file() or not profiles_path.is_file():

        return {

            "available": False,

            "reason": "This run did not record positional profile data.",

            "gene": gene_name,

        }



    def read_rows(path: Path) -> list[dict]:

        with path.open("r", newline="", encoding="utf-8-sig") as handle:

            return list(csv.DictReader(

                line for line in handle if not line.startswith("#")

            ))



    try:

        windows = read_rows(windows_path)

        profiles = read_rows(profiles_path)

    except (OSError, csv.Error) as exc:

        raise HTTPException(

            status_code=500, detail=f"Positional data could not be read: {exc}"

        ) from exc



    return {

        "available": bool(windows and profiles),

        "gene": gene_name,

        "target_length": gene.get("length_bp"),

        "candidate": gene.get("candidate"),

        "windows": windows,

        "profiles": profiles,

        "species": [row.get("name") for row in results.get(

            "reference_panels", {}).get("searched", []) if row.get("name")],

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

        "reference_record",

        "reference_pos_1based",

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

                    "organism": row.get("organism", ""),

                    "reference_record": row.get("transcript", ""),

                    "reference_pos_1based": row.get("transcript_pos_1based", ""),

                    "gene_start": row.get("gene_start", ""),

                    "gene_end": row.get("gene_end", ""),

                    "length": row.get("length", ""),

                    "strand": row.get("strand", ""),

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

                f'attachment; filename="{_safe_download_name(job_id, "hits.csv")}"'

        },

    )
