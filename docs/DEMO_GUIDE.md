# dsGuard — Live Demonstration Guide



This guide is for **running a live demonstration** of dsGuard in front of an audience.

It lists the exact commands that were executed against this repository, the URLs to

open, three concrete inputs to screen, what each screen means, and the scientific

caveats you must state out loud.



It deliberately does **not** repeat the material that already exists elsewhere:



| Read this | For |

|---|---|

| [../README.md](../README.md) | installation, `.env` configuration, panel layout, scoring model |

| [OVERVIEW.md](OVERVIEW.md) | architecture, methodology, module list, full test-area breakdown |

| [DEMO_READINESS_REPORT.md](DEMO_READINESS_REPORT.md) | the verification history and known limitations |



Everything below was re-checked against this checkout on Windows (PowerShell) while

writing it. Where a number is quoted, it came from a command that was actually run.



---



## 1. Pre-flight checklist (do this before the audience arrives)



### 1.1 Prerequisites



Installation and `.env` creation are covered in [../README.md](../README.md) §Installation

and §Configuration. Make sure both are already done: the `.venv` exists with the backend

dependencies, and `.env` exists at the project root with `MONGODB_URI` and

`MONGODB_DATABASE`.



### 1.2 MongoDB must be reachable — this is the #1 demo risk



The dashboard stores each screening job in MongoDB. If the database is unreachable the

app still starts, but `/api/health` reports `"status":"degraded"` and `"mongodb":"disconnected"`,

and creating a job fails with HTTP 500 — the demo appears to "break" for a reason that has

nothing to do with the software.



> **Verified on this machine:** the Atlas cluster in `.env` is **not reachable** —

> `pymongo.errors.ServerSelectionTimeoutError: SSL handshake failed ...

> [SSL: TLSV1_ALERT_INTERNAL_ERROR]`. This is a network/TLS block, not a code bug.



For a live demo, the reliable option is a **local MongoDB**:



```powershell

# 1. Start a local MongoDB (install MongoDB Community Server first if needed).

#    Use any existing empty folder as the data directory.

mongod --dbpath "$env:USERPROFILE\dsguard-db" --port 27017

```



```ini

# 2. Point .env at it (project root, git-ignored). Keep the database name.

MONGODB_URI=mongodb://127.0.0.1:27017

MONGODB_DATABASE=dsguard

```



If you must use Atlas instead, confirm it works **from the demo network** (a campus

network or a phone hotspot may behave differently) before you present — see the health

check in §4.



### 1.3 Free the default ports



```powershell

Get-NetTCPConnection -LocalPort 8000,5173 -State Listen -ErrorAction SilentlyContinue |

  Select-Object LocalPort, OwningProcess

```



No output = both ports are free. If something is listed, see §10.



### 1.4 Run the test suite and the production build (once, before the demo)



These prove the *software* runs. They do **not** prove any biological claim.



```powershell

# From the project root

.\.venv\Scripts\python.exe -m pytest -q        # verified: 148 passed

```



```powershell

# From frontend/

cd frontend

npm test                                        # verified: 28 passed, 0 failed

npm run lint                                    # verified: clean (no ESLint errors)

npm run build                                   # verified: success, 1932 modules, ~268 kB JS

```



---



## 2. Start the application



Two terminals, both started from the project root.



```powershell

# Terminal 1 — backend

.\.venv\Scripts\python.exe -m uvicorn backend.main:app --host 127.0.0.1 --port 8000

```



```powershell

# Terminal 2 — frontend

cd frontend

npm run dev

```



| URL | What it is |

|---|---|

| http://localhost:5173 | The dashboard (open this) |

| http://127.0.0.1:8000/api/health | Backend health |

| http://127.0.0.1:8000/api/panels | Which reference panels are installed |

| http://127.0.0.1:8000/docs | Interactive API documentation |



> The dashboard talks to `http://127.0.0.1:8000` by default. If you start the backend on

> another port (`--port 8001`), also set `VITE_API_BASE=http://127.0.0.1:8001` before

> `npm run dev`. Vite serves the dashboard on `http://localhost:5173` (IPv6 loopback), so

> use `localhost`, not `127.0.0.1`, in the browser. The backend accepts any loopback dev

> origin, so a Vite fallback port is fine.



---



## 3. Thirty-second smoke test



```powershell

# Should be status "ok" and "mongodb":"connected" once the database is reachable.

Invoke-WebRequest -UseBasicParsing http://127.0.0.1:8000/api/health | Select-Object -Expand Content

```



Open http://127.0.0.1:8000/api/panels — you should see `rice` as `available: true` and the

other panels (`cyrtorhinus_lividipennis`, `honeybee`, `spider`) as `available: false` with a

`reason`. Then open the dashboard and confirm the "System health" card matches.



If health says `degraded` / `disconnected`, fix MongoDB (§1.2) before continuing.



---



## 4. The demonstration (~10 minutes, three runs)



Use the project's own dataset for all three runs:

`Random/data/targets/nilaparvata_lugens_13_cds.fasta` (13 *Nilaparvata lugens* CDS records,

all verified to parse).



### Run 1 — Upload the 13-gene FASTA file (≈2 minutes)



1. In **Configure a screening job**, click **Upload FASTA**.

2. Choose `Random\data\targets\nilaparvata_lugens_13_cds.fasta`.

3. Leave the reference panel on **Rice (Oryza sativa)** → **Start screening**.

4. Wait for the job to finish (the rice panel is ~383 Mb, so the scan takes roughly

   1–2 minutes).



**Expected result (verified):** 13 genes screened; **all 13** report

`FALLBACK — THRESHOLD EXCEEDED`; an amber banner reads **"Configured hazard threshold

exceeded"**. In the verified run there were **8 exact matches ≥ 21 nt** (5 in

`NlCHSA_NlCHS1`, 1 each in `NlSRP54`, `NlEcR`, `NlRan`), and the recommended window's

longest match was **17–19 nt** for every gene.



Point out the honesty of this: the tool did **not** declare success. It found candidate

windows, but none of them met the configured risk limit, and it says so.



### Run 2 — One gene alias in the search box (≈1 minute)



1. Click **Gene / NCBI search**.

2. Type the alias **`NlRan`** → **Start screening**.



**Expected result (verified):** the alias resolves to the intended *N. lugens* record

`NlRan|accession=KT313028.1|cds_length=645`. The table shows it under the clean name

`NlRan`, and the report keeps the full identifier and accession. If NCBI is offline, use

the accession instead (Run 2b).



### Run 3 — One pasted sequence (≈1 minute)



1. Click **Paste sequence**.

2. Paste the **complete `NlRan` record** (header + sequence, 645 nt) shown below — or copy

   the first record straight out of the FASTA file:



```text

>NlRan|accession=KT313028.1|cds_length=645

ATGGAAGAACAGGATATGCCAACATTCAAGTGTGTGCTGGTTGGAGATGGTGGTACCGGCAAAACTACAT

TTGTCAAACGACATCTTACCGGAGAATTTGAAAAGAAGTATGTTGCCACCCTTGGAGTTGAAGTTCACCC

CCTTGTATTTCACACAAACAGAGGTGTGATTAGGTTCAATGTGTGGGACACAGCTGGCCAGGAAAAGTTC

GGTGGACTTCGTGATGGATATTACATTCAGGGACAATGCGCCATCATTATGTTTGACGTAACGTCAAGAG

TCACCTACAAGAACGTTCCCAACTGGCACAGAGATTTAGTGAGGGTTTGCGAAAACATTCCCATTGTACT

ATGCGGCAACAAAGTAGACATCAAGGACAGGAAAGTCAAGGCCAAGAGCATAGTCTTCCATAGGAAGAAG

AACCTTCAGTACTACGACATCAGTGCGAAAAGCAACTACAACTTCGAGAAGCCGTTCCTGTGGTTGGCAA

AGAAGCTGATCGGTGACCCCAACCTGGAGTTCGTCGCCATGCCCGCCCTCCTCCCTCCCGAGGTCACAAT

GGACCCCCAATGGCAGTCGCAGATCGAAAAAGATCTCAAGGAAGCTCAGGCCACCGCCCTTCCCGAGGAT

GATGAAGATCTGTAG

```



3. Leave the panel on **Rice** → **Start screening**.



**Expected result (verified):** one gene, displayed as `NlRan`, with the full header and

accession preserved in the report and exports. To demonstrate the *raw* variant instead,

delete the first line and paste only the bases — the app then synthesises a record named

`pasted_target`.



Avoid pasting a short fragment: a window is 400 nt, so a sequence shorter than that

produces a "no window" result rather than a ranked report.



### What to show in the results



- **Ranked off-target matches** — sorted by hazard (worst first) with the `THRESHOLD`

  column showing `meets_threshold` / `exceeds_threshold`.

- **Pareto front (potency vs hazard)** — the non-dominated windows for one gene

  (pick a gene from the dropdown).

- **Preventive measures & risk mitigation** — findings, explanations and actions derived

  from *this* run; heuristic scores and observed matches are kept in separate categories.

- **Download report** — the three exports, all generated from the same run:

  `results.csv` (ranked matches), `results.json` (full structured export) and `report.pdf`.



---



## 5. What each dashboard section means



| Section | Plain-language meaning |

|---|---|

| **Screen with confidence** (top) | Summary strip: the most recent completed job, its panel and its headline numbers. |

| **Recent screening jobs** | The last few runs, newest first. `TARGET` shows the job's target (a file name, an accession, or a pasted record's alias). |

| **System health** | Whether the backend is reachable and the database is connected. `degraded` here usually means MongoDB is down. |

| **Configure a screening job** | The form. Three input modes — **Upload FASTA**, **Gene / NCBI search**, **Paste sequence** — plus the reference panel and **Start screening**. |

| **Screening history** | Every job with its ID, target, panel, status and creation time. Click a row to load that job's results. Statuses: `queued` → `running` → `completed` / `failed`. |

| **Screening results** | Why: everything below describes the single job you selected. |

| **Run parameters and reference coverage** | The exact settings used (`window`, `k`, `min_report_match`, `max_hazard`, L50/κ) and, per searched panel, its size, record count and SHA-256. Panel checksums are how you prove which data was used. |

| **Summary table** | One row per gene: `GENE`, `SCREENING STATUS`, `RECOMMENDATION`, `HAZARD SCORE`, `MAX HAZARD`, `LENGTH`, `GC %`, `MAX MATCH` and `HITS (≥21 nt)`. |

| **`HITS (≥21 nt)`** | Count of **confirmed** exact matches at or above the guaranteed-detection floor (21 nt). It is deliberately narrower than the ranked list below. |

| **Ranked off-target matches** | The full candidate list with target window, matched organism, reference ID, match length, identity, strand, coordinates, hazard and `THRESHOLD` status. |

| **Pareto front (potency vs hazard)** | A scatter of the non-dominated candidates for the selected gene: up-and-left (high potency, low hazard) is better. |

| **Preventive measures & risk mitigation** | Concrete follow-up actions for this run, grouped as findings / explanations / actions / evidence. |

| **Download report** | The three exports for the selected job (`results.csv`, `results.json`, `report.pdf`). |

| **Detailed screening hits** | Collapsible raw list of every confirmed match ≥ the detection floor, exportable as the detailed CSV. |

| **Gene report viewer** | The engine's per-gene `report.md` for the selected gene. |



---



## 6. The hazard limit, and why all 13 genes fail it



**The configured limit is `max_hazard = 0.05`** (recorded in every run's parameters). Under

the current placeholder model that corresponds to an equivalent perfect-match length of

**≈15.1 nt** — the run records this as `hazard_threshold_equivalent_match_bp: 15.11`.



Why the 13 *N. lugens* genes fail:



- The rice reference panel is a **whole genome** (~382.8 Mb, 15 records). In a genome that

  size, a 400 nt window very often contains a **17–19 nt** perfect match to *something*.

- A 17–19 nt match sits **above** the ~15.1 nt equivalent of the 0.05 limit, so the

  window's hazard exceeds the limit even though it is well below the 21 nt

  guaranteed-detection floor used by `HITS`.

- The result: every gene's best window is a **fallback** — the lowest-risk window the

  engine could find — and it still fails the configured criterion.



**Why the fallback must stay labelled as a failure.** The engine never silently promotes a

fallback to a pass. A fallback means *"no window met the criterion you configured"*, so the

dashboard shows **FALLBACK — THRESHOLD EXCEEDED** and the amber **"Configured hazard

threshold exceeded"** banner. Never describe these results as "safe", "passing" or

"acceptable"; the correct phrasing is *"no window satisfied the configured hazard limit"*.

If the criterion itself is wrong for your organism, change the **threshold** — but then

state clearly that you relaxed it, and note that this changes the meaning of the result.



---



## 7. Limitations you must state out loud



1. **The hazard score is a heuristic, not a probability.**

   `hazard(m) = 1 / (1 + exp(-(m - L50) / kappa))` is a monotone rescaling of the longest

   perfect match `m`. A score of 0.4 does **not** mean a 40 % chance of harm. Each run

   records `calibrated: false`.

2. **`L50 = 21.0` and `kappa = 2.0` are placeholders, not fitted values.**

   They were chosen to behave sensibly; they were **not** fitted to outcome data, and the

   run records `l50_kappa_fitted: false`. The **potency** score is borrowed from mammalian

   siRNA design rules (Reynolds 2004 / Ui-Tei 2004) and is not validated for insect

   Dicer-2 processing of long dsRNA.

3. **Reference coverage equals the installed panels.** Only `rice.fasta` is installed in

   this checkout; `cyrtorhinus_lividipennis`, `honeybee` and `spider` are **not installed

   and are rejected by the API**. A missing panel is recorded as `skipped` with a reason and

   surfaced as an unresolved risk — it is **never** treated as zero risk. Absence of a match

   in the panels you did search is not evidence of absence in the organism.

4. **Exact, gapless matches only.** Mismatched or gapped off-targets, and miRNA-like

   partial-complementarity effects, are not detected.



---



## 8. Software success is not biological validation



Keep these two statements separate, every time:



- **What a green run proves:** the pipeline executed correctly — inputs were parsed, the

  engine scanned the installed panels deterministically, the numbers are reproducible, and

  the reports agree with each other.

- **What it does not prove:** that the dsRNA is safe or effective in any organism. dsGuard

  measures **sequence similarity**. A low hazard score is *necessary, not sufficient*, for

  non-target safety, and no wet-lab result is claimed or implied by this tool.



Do **not** say, at any point: "all candidates passed", "these genes are safe", "biological

safety is proven", or "the tool validated the design". Say instead: "the software ran

correctly and reported that no candidate met the configured similarity limit."



---



## 9. Quick reference — verified commands



| Purpose | Command (from project root) | Verified result |

|---|---|---|

| Backend tests + engine tests | `.\.venv\Scripts\python.exe -m pytest -q` | **148 passed** |

| Frontend tests | `cd frontend; npm test` | **28 passed** |

| Frontend lint | `cd frontend; npm run lint` | clean |

| Production build | `cd frontend; npm run build` | success, 1932 modules, ~268 kB JS |

| Engine only (CLI) | `.\.venv\Scripts\python.exe -m dsguard --target Random\data\targets\nilaparvata_lugens_13_cds.fasta --panel rice=reference_panels\rice.fasta --out outputs\cli_run` | 13 genes ranked; per-gene files written |



---



## 10. Troubleshooting (demo-specific)



| Symptom | Cause and fix |

|---|---|

| `/api/health` shows `"mongodb":"disconnected"`, jobs return HTTP 500 | MongoDB is unreachable. Use a local `mongod` and set `MONGODB_URI=mongodb://127.0.0.1:27017` in `.env` (§1.2). Atlas TLS failures (`TLSV1_ALERT_INTERNAL_ERROR`) are a network block, not a bug. |

| Dashboard says "Cannot reach the backend. Check that FastAPI is running at 127.0.0.1:8000." | The backend is not running, or it runs on a different port. Start it on 8000, or set `VITE_API_BASE` to the port you used, then restart `npm run dev`. |

| Backend exits with `MONGODB_URI is missing` | `.env` is missing or not at the project root. |

| Backend exits with `The dsGuard CLI was not found` | Run `.\.venv\Scripts\python.exe -m pip install -e ".[dev]"` so `.venv\Scripts\dsguard.exe` exists. |

| `"None of the requested reference panels are installed"` (HTTP 400) | Add the panel FASTA to `reference_panels/`, or point `DSGUARD_PANEL_DIR` at the folder that has it. `GET /api/panels` lists each expected path. |

| Port 8000 / 5173 already in use | A previous dev server is still listening. Identify it — `Get-NetTCPConnection -LocalPort 8000,5173 -State Listen \| Select-Object LocalPort,OwningProcess`, then `Get-Process -Id <PID>`. Stop **only** the stale dsGuard PIDs from an **elevated** PowerShell (`Stop-Process -Id <PID> -Force`); stopping a process started by an elevated shell fails from a normal shell with "Access is denied". Or just switch ports (see §2). |

| NCBI search returns HTTP 400 | No internet, or the record does not exist. Use the **accession** (fetched directly) or the **Upload FASTA** / **Paste sequence** modes, which need no network. |

| A job's results show a 404 asking to re-run | That job predates structured results. Re-run it from its original inputs; the app never fabricates results. |

| The audience asks why `HITS` is smaller than the ranked list | `HITS (≥21 nt)` counts only confirmed matches at or above the guaranteed-detection floor; the ranked table also lists shorter matches inside the recommended window because that shorter match is what sets the window's hazard. |



---



*All commands in this guide were executed against this checkout on Windows (PowerShell)

while writing it. Where a result is stated, it is the result that command produced.*
