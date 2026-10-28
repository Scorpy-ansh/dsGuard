# dsGuard



Non-target-aware dsRNA design. Given a pest target gene, dsGuard ranks candidate

double-stranded RNA windows by predicted potency and screens each window against

curated reference panels of non-target organisms, reporting the longest perfect

sequence match found in any panel, a heuristic hazard score, Pareto-optimal

candidates, preventive (risk-mitigation) measures, and downloadable PDF/CSV/JSON

reports.



dsGuard measures **sequence similarity only**. A low score is necessary, not

sufficient, for non-target safety.



## Repository layout



```

Random/dsguard/          Analysis engine (pure Python + numpy); also the CLI

Random/data/             Demo targets and reference-panel provenance notes

backend/main.py          FastAPI service (jobs, panels, NCBI, exports)

backend/reporting.py     PDF / CSV / JSON report builders

backend/test_*.py        Backend tests (exports, jobs, NCBI, input modes)

Random/test_*.py         Engine tests (engine, batch, ranking, scoring, reproducibility)

frontend/                React + Vite dashboard (single-page)

reference_panels/        Curated non-target FASTA panels (see below)

outputs/jobs/<job_id>/   Per-run results and generated reports

```



## Architecture



- **Engine (`Random/dsguard/`)** — pure-Python/numpy command-line tool. It reads a

  target FASTA and one or more reference panels, performs a strided seed-and-extend

  exact match (`k = 11`, provably complete for matches ≥ `2k-1 = 21 nt`) on both

  strands, computes per-window potency and hazard, extracts the Pareto front, picks

  a recommendation, and writes `summary.csv`, `manifest.json`, `results.json` and a

  per-gene folder (`report.md`, `windows.csv`, `pareto.csv`, `hits.csv`,

  `offtarget_ranked.csv`, `profiles.csv`, `primers.json`).

- **Backend (`backend/`)** — FastAPI. Accepts a screening job, validates inputs,

  resolves installed panels, writes the target, invokes the engine as a subprocess,

  stores job metadata in MongoDB, and generates all three downloads from the run's

  `results.json`.

- **Frontend (`frontend/`)** — React + Vite dashboard. Provides the four input

  methods, panel selection, the ranked off-target table, Pareto plot, preventive

  measures, per-gene reports, and the three download buttons.



## Requirements



- Windows (documented commands use PowerShell). Linux/macOS work with the

  equivalent `python -m venv` / `bin/activate` commands.

- Python 3.10+ and Node.js 18+ (Node 22 was used during development).

- A MongoDB connection string (MongoDB Atlas or a local `mongod`).

- Outbound internet access only for NCBI input mode.



## Installation



```powershell

# 1. Python environment (from the project root)

python -m venv .venv

.\.venv\Scripts\python.exe -m pip install -U pip

.\.venv\Scripts\python.exe -m pip install -e ".[dev]"

.\.venv\Scripts\python.exe -m pip install fastapi "uvicorn[standard]" pymongo python-multipart python-dotenv reportlab



# 2. Frontend dependencies

cd frontend

npm install

cd ..

```



`pip install -e ".[dev]"` installs the `dsguard` console script and pytest from the

root `pyproject.toml` (the package sources live in `Random/dsguard`).



## Configuration



Create `.env` at the project root (it is git-ignored):



```

MONGODB_URI=mongodb+srv://<user>:<password>@<cluster>/...

MONGODB_DATABASE=dsguard

```



Optional: `DSGUARD_PANEL_DIR` overrides the directory that holds reference panels

(default: `reference_panels/`).



## Reference panels (required data files)



Put curated non-target FASTA files in `reference_panels/` with these exact names:



| Panel | File | Notes |

|---|---|---|

| Rice | `rice.fasta` | installed in this checkout (15 records, ~383 Mb) |

| *Cyrtorhinus lividipennis* | `cyrtorhinus_lividipennis.fasta` | not installed |

| Honeybee | `honeybee.fasta` | not installed |

| Spider | `spider.fasta` | not installed; document the species in the FASTA headers |



Panels without a file are disabled in the UI and rejected by the API. When only

some requested panels exist, the missing ones are recorded as `skipped` with a

reason and surfaced as an unresolved risk — a missing panel is **never** treated as

zero risk. `GET /api/panels` reports exactly which files were found and why one is

unavailable.



## Running the application (Windows)



Two terminals, both from the project root.



```powershell

# Terminal 1 — backend

.\.venv\Scripts\python.exe -m uvicorn backend.main:app --host 127.0.0.1 --port 8000

# Health check: http://127.0.0.1:8000/api/health   Docs: http://127.0.0.1:8000/docs



# Terminal 2 — frontend

cd frontend

npm run dev

# Dashboard: http://localhost:5173

```



The backend refuses to start with a clear message if `MONGODB_URI` is missing.



### Supported input methods



1. **Upload FASTA** — a `.fa/.fasta/.fna/.fas/.txt` file whose first non-empty line

   starts with `>`.

2. **Gene / NCBI search** — an NCBI accession (fetched directly) or a gene name

   (searched in NCBI Nucleotide; the candidate records are ranked so a curated

   RefSeq transcript naming the gene is preferred).

3. **Paste sequence** — a FASTA record or a raw nucleotide sequence (A/C/G/T/U/N);

   a header is synthesised when one is omitted.



Each valid input produces a new analysis from the current run; results are never

fabricated or substituted from another job.



### Downloads



Every completed job offers three exports, all generated from that job's

`results.json`:



- `results.csv` — ranked off-target matches (job-wide rank, target ID, window,

  coordinates, matched organism, reference accession/record ID, match length,

  identity/mismatch/gap, hazard and category, threshold status, evidence).

- `results.json` — complete structured export with parameters, provenance,

  scoring definitions and preventive measures.

- `report.pdf` — formatted human-readable report.



### CLI usage (engine only)



```powershell

.\.venv\Scripts\python.exe -m dsguard `

  --target Random\data\targets\nilaparvata_lugens_13_cds.fasta `

  --panel rice=reference_panels\rice.fasta `

  --out outputs\cli_run

```



## Testing



```powershell

# Python (engine + backend)

.\.venv\Scripts\python.exe -m pytest -q



# Frontend

cd frontend

npm run lint

npm test

npm run build

```



## Troubleshooting



- **`MONGODB_URI is missing`** — create `.env` at the project root (see

  Configuration). The backend reads it relative to the project root, not the shell

  working directory.

- **`The dsGuard CLI was not found`** — run

  `.\.venv\Scripts\python.exe -m pip install -e ".[dev]"` so

  `.venv\Scripts\dsguard.exe` exists (or so `python -m dsguard` resolves).

- **`None of the requested reference panels are installed`** — add the panel FASTA

  to `reference_panels/` or point `DSGUARD_PANEL_DIR` at the folder that has it.

  `GET /api/panels` lists every panel's expected absolute path and why it is

  unavailable.

- **CORS errors in the browser** — the backend accepts any loopback dev origin, so

  a Vite fallback port (5174, 5175, …) is fine; a non-loopback origin is not.

- **Port already in use (8000 / 5173)** — a previous dev server is still listening.

  Identify it and confirm it is a stale dsGuard process before stopping anything:



  ```powershell

  # Who owns the default ports? (shows LocalAddress, LocalPort, OwningProcess)

  Get-NetTCPConnection -LocalPort 8000,5173 -State Listen |

    Select-Object LocalAddress, LocalPort, OwningProcess



  # Show the process behind a PID (name; CommandLine may be blank if elevated)

  Get-Process -Id <PID> | Select-Object Id, ProcessName

  Get-CimInstance Win32_Process -Filter "ProcessId=<PID>" | Select-Object Name, CommandLine



  # Stop ONLY the stale dsGuard PIDs, from an ELEVATED PowerShell (Run as Administrator):

  Stop-Process -Id 25132,34956 -Force          # replace with the PIDs found above

  ```



  `Stop-Process` / `taskkill` fail with **"Access is denied"** from a normal shell

  when the old servers were started elevated; use an Administrator shell. Vite binds

  to `::1` (IPv6 loopback) by default, so the frontend answers on

  `http://localhost:5173` but not `http://127.0.0.1:5173` — that is expected.

  To avoid the conflict entirely, use another port: `--port 8001` for uvicorn or

  `npm run dev -- --port 5174` for Vite.

- **NCBI retrieval fails** — the API returns a 400 with the transport error; check

  connectivity or use an accession, which is fetched directly.

- **A job reports missing `results.json`** — it predates structured results; the

  API returns a 404 explaining that the analysis must be re-run from its original

  inputs. It is not fabricated.

- **19–21 nt matches and 17-vs-9 counts** — `HITS` counts exact matches ≥ 21 nt

  (the provably complete floor). The ranked table additionally lists the longest

  match inside the recommended window even when shorter, because that match sets

  the window's hazard.



## Scoring model: what is measured and what is assumed



The engine reports three distinct things. Keep them separate when reading any result:



- **Observed exact match (nt).** A perfect, gapless match between a candidate window

  and a reference record on either strand. Enumeration is provably complete for

  matches of at least `2k-1 = 21 nt` (default `k=11`); shorter matches are best-effort

  lower bounds. This is a measured sequence fact.

- **Hazard score (0–1).** `hazard(m) = 1 / (1 + exp(-(m - L50) / kappa))`. A monotone

  heuristic that rises with the longest perfect match. It is **not** a probability of

  biological harm.

- **Potency score (0–1).** The mean over the sense-strand 21-mers of

  `0.5 * (reynolds/7 + ui_tei/4)`, combining published mammalian siRNA design rules

  (Reynolds 2004, 7 of 8 criteria; Ui-Tei 2004, 4 rules). Used only as a weak ordering

  prior.



Neither score is a validated prediction of biological effect, and no calibrated

validation data is bundled with this build.



### Candidate selection and screening status



The configured `max_hazard` value is retained as an informational comparison and

for compatibility with older API clients and result files; it is not a mandatory

candidate eligibility or pass/fail gate. For each valid candidate window, dsGuard

uses the existing Pareto front and selects the highest-potency member, breaking

ties by lower hazard and then deterministic candidate coordinate/order. New

results record `selection_policy = pareto_highest_potency_then_lowest_hazard`.



New per-gene results use neutral `screening_status` values:



- `screened_candidate_available` — screening completed and a real candidate

  sequence was recorded.

- `screened_no_candidate` — screening completed but no valid full-length window

  could be produced.

- `screening_incomplete` — screening did not produce enough information to

  determine a result.



Legacy fields including `max_hazard`, `threshold_status`,

`recommendation_status`, `is_fallback`, and

`hazard_threshold_equivalent_match_bp` remain available for compatibility. They

must not be interpreted as validated biological safety or efficacy.



### Placeholder parameters (L50 and kappa)



`L50` and `kappa` are **placeholder constants, not fitted values**:



- `L50 = 21.0 nt` (the midpoint where the score is 0.5)

- `kappa = 2.0 nt` (the steepness)



They were chosen so the score behaves sensibly between the guaranteed detection floor

(21 nt) and the window length; they were **not** fitted to published non-target outcome

data, and the run records `l50_kappa_fitted: false`. The default hazard threshold

`max_hazard = 0.05` corresponds to an equivalent perfect-match length of roughly

`L50 + kappa * ln(0.05/0.95) ≈ 15.1 nt` under these placeholders. **Change the

threshold, not the model, if you re-fit:** re-fitting L50/kappa changes what a given

score means and would invalidate comparison with prior runs. Every report labels these

values "(placeholder)".



### Coordinate conventions



Internally matches are 0-based, half-open; every file and report states coordinates as

**1-based inclusive** (`target_start = q_start + 1`, `target_end = q_end`). Ambiguous

bases (anything other than A/C/G/T/U) never match and cannot be spanned by a match.



## Reproducibility and provenance



Each run writes `manifest.json` and `results.json` recording, for the exact run: the

dsGuard version, Python/numpy versions and platform; all analysis parameters (`window`,

`k`, `L50`, `kappa`, `max_hazard`, `min_report_match`, and the derived equivalent-match

length with `l50_kappa_fitted`); the target file path and its SHA-256; and, per searched

panel, the file path, SHA-256, size in Mb, record count, a sample of record accessions,

the number of records lacking a header annotation, and the first header annotation as

auditable provenance.



Runs are deterministic. Identical target, panel and parameters produce identical ranked

tables, recommendations and provenance (only the run id and wall-clock timings differ);

this is covered by an automated test. Panels that are configured but not installed are

reported as `skipped` with a reason and surfaced as an unresolved risk — a missing panel

is never treated as zero risk.



## Limitations



- Hazard is a heuristic over perfect sequence similarity. It does not model uptake,

  nuclease stability, RISC loading or dose, and cannot detect miRNA-like

  partial-complementarity effects.

- Potency rules were derived from synthetic siRNA in mammalian cells and are not

  validated for insect Dicer-2 processing of long dsRNA.

- Only exact, gapless matches are detected; mismatched or gapped off-targets are not.

- Reference coverage equals the installed panels; absence of a match in a panel is not

  evidence of absence in the organism.
