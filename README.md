# dsGuard
### Sequence-Based Off-Target Screening and Risk Assessment

**dsGuard** is a research-oriented software prototype for screening candidate target sequences against reference transcript panels. It identifies exact sequence matches, ranks candidate windows, generates per-gene reports, and exports structured results for further analysis.

> **Research limitation:** dsGuard is a screening aid, not a validated biological safety predictor. Sequence similarity alone cannot establish biological harm or prove non-target safety.

## Features

- **Sequence screening:** Process target coding sequences in FASTA format.
- **Reference-panel comparison:** Compare target sequences against selected reference transcript panels.
- **Candidate-window analysis:** Evaluate candidate windows using configurable parameters.
- **Exact-match detection:** Identify and report exact sequence matches meeting the configured search criteria.
- **Gene-level reports:** Generate individual Markdown reports for each target gene.
- **Structured exports:** Produce summary CSV files, detailed hit records, manifests, window profiles, Pareto results, and primer-related outputs.
- **Web dashboard:** Submit screening jobs, review screening history, inspect reports, and download results.
- **Persistent job metadata:** Store job information and summaries in MongoDB Atlas through a FastAPI backend.

## Technology Stack

| Component | Technology |
|---|---|
| Frontend | React, Vite, CSS |
| Backend API | FastAPI, Uvicorn |
| Screening engine | Python |
| Sequence analysis | Custom dsGuard modules and NumPy |
| Database | MongoDB Atlas |
| API documentation | OpenAPI / Swagger UI |
| Test framework | pytest |

## Project Structure

```text
dsGuard/
├── backend/
│   └── main.py
├── frontend/
│   ├── src/
│   │   ├── App.jsx
│   │   ├── App.css
│   │   ├── index.css
│   │   └── main.jsx
│   └── package.json
├── outputs/
│   └── jobs/
├── Random/
│   ├── dsguard/
│   │   ├── cli.py
│   │   ├── engine.py
│   │   ├── batch.py
│   │   ├── optimise.py
│   │   ├── primers.py
│   │   ├── risk.py
│   │   ├── seq.py
│   │   └── sirna.py
│   ├── data/
│   │   ├── targets/
│   │   └── reference_panels/
│   ├── test_batch.py
│   ├── test_engine.py
│   ├── test_seq.py
│   ├── test_sirna.py
│   ├── pyproject.toml
│   └── README.md
├── .env
└── .gitignore
```

The tree above is a simplified overview. Generated files, local environments, and installed dependencies are intentionally omitted.

## Requirements

- Python 3.12 or a compatible version
- Node.js and npm
- MongoDB Atlas account or a compatible MongoDB deployment
- Valid target FASTA sequences
- An appropriately prepared reference transcript panel

## Configuration

Create a `.env` file in the project root:

```ini
MONGODB_URI="your-mongodb-connection-string"
MONGODB_DATABASE="dsguard"
```

Replace the placeholder connection string with your own MongoDB connection string. Never commit credentials, database passwords, or other secrets to version control.

## Running the Application

Open three terminals from the project directory.

### 1. Start the backend

```bat
cd /d "C:\Users\gfxsc\Desktop\dsGuard"
.venv\Scripts\python.exe -m uvicorn backend.main:app --reload
```

Backend API: `http://127.0.0.1:8000`

Interactive API documentation: `http://127.0.0.1:8000/docs`

Health check: `http://127.0.0.1:8000/api/health`

### 2. Start the frontend

```bat
cd /d "C:\Users\gfxsc\Desktop\dsGuard\frontend"
npm run dev
```

Open the local URL printed by Vite, normally `http://localhost:5173`.

### 3. Run the screening engine from the command line

```bat
cd /d "C:\Users\gfxsc\Desktop\dsGuard\Random"
dsguard --target "data\targets\nilaparvata_lugens_13_cds.fasta" --panel rice="data\reference_panels\rice\rice_transcripts.fna" --window 400 --k 11 --out "rice_screen_test_out"
```

This command uses the available *Nilaparvata lugens* target CDS dataset and the rice reference transcript panel. Ensure both files exist before running it.

## Screening Workflow

1. Upload a target FASTA file and select a reference panel.
2. Configure the screening parameters.
3. Submit a screening job through the dashboard.
4. Monitor the job status.
5. Inspect the gene-level summary and detailed exact-match records.
6. Open individual gene reports.
7. Export summary and detailed-hit CSV files for further analysis.

## Output Files

| File | Purpose |
|---|---|
| `summary.csv` | Per-gene screening summary |
| `manifest.json` | Run metadata and configuration |
| `report.md` | Individual gene report |
| `hits.csv` | Detailed sequence-match records |
| `windows.csv` | Candidate-window results |
| `pareto.csv` | Pareto-analysis output |
| `profiles.csv` | Window/profile data |
| `primers.json` | Primer-related output, where generated |
| `engine.log` | Engine execution log, where generated |

Output availability depends on the job and engine configuration.

## API Endpoints

The current backend includes endpoints for:

- `GET /api/health` — Check backend and engine health.
- `POST /api/screen` — Submit a screening job.
- `GET /api/jobs` — Retrieve job history.
- `GET /api/jobs/{job_id}` — Retrieve job details.
- `GET /api/jobs/{job_id}/results` — Retrieve summary results and manifest.
- `GET /api/jobs/{job_id}/hits` — Retrieve detailed hit records.
- `GET /api/jobs/{job_id}/hits/download` — Download detailed hits as CSV.
- `GET /api/jobs/{job_id}/reports` — List gene reports.
- `GET /api/jobs/{job_id}/reports/{gene_name}` — Retrieve an individual gene report.

Use `/docs` to inspect the live API schemas and endpoint behavior.

## Testing

Run the Python test suite:

```bat
cd /d "C:\Users\gfxsc\Desktop\dsGuard"
.venv\Scripts\python.exe -m pytest -q "Random"
```

The current development baseline previously reported 29 passing tests. Rerun the suite after code changes and record the actual result for each release.

## Scientific Limitations

- Exact sequence matches do not independently establish biological activity or toxicity.
- A lack of detected exact matches does not prove non-target safety.
- Results depend on reference-panel coverage, sequence quality, and screening parameters.
- The current prototype does not comprehensively model uptake, nuclease stability, RISC loading, dose, or miRNA-like partial-complementarity effects.
- The default risk-model parameters are placeholders unless independently fitted and validated.
- Potency methods derived from mammalian siRNA data have not been validated for insects in this prototype.

Scientific conclusions require appropriate reference datasets, independent validation, and experimental evidence.

## Security Notes

- Keep `.env` and credentials out of Git.
- Validate uploaded files and enforce appropriate upload-size limits.
- Restrict access to screening jobs and downloaded results in public deployments.
- Configure production CORS, persistent file storage, and database permissions before deployment.

## Project Status

dsGuard is an actively developed research prototype. Its current implementation includes a command-line screening engine, a web dashboard, a FastAPI backend, MongoDB integration, job history, gene-level reports, and CSV exports.

Further work includes comprehensive edge-case testing, scientific model validation, deployment hardening, and preparation of reproducible evaluation materials.
