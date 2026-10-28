# dsGuard — Technical Overview



This document describes what dsGuard is, how it works, what has been verified, and

what its limitations are. Every number and behaviour below comes from the code in

this repository and from test runs against it; no research findings, performance

benchmarks or novelty claims are asserted.



---



## 1. Problem statement and objectives



RNA interference (RNAi) can be triggered in an insect by a double-stranded RNA

(dsRNA) that matches a pest gene, which makes ingested dsRNA a candidate pest-control

agent. The same dsRNA can also match a gene in a non-target organism and cause an

unintended effect. The practical design question is therefore: *among the possible

dsRNA windows of a pest gene, which ones are least likely to match any curated

non-target sequence?*



dsGuard addresses the **sequence-similarity part** of that question. Its objectives,

as implemented, are:



1. Detect, deterministically and completely, every perfect (gapless) match of at

   least 21 nt between any candidate window of a target gene and any installed

   non-target reference panel, on either strand.

2. Rank candidate windows by a composite potency prior and by the hazardousness of

   their longest panel match.

3. Select the highest-potency window whose hazard is at or below a configured

   threshold, and expose a Pareto front of potency-vs-risk candidates for review.

4. Derive preventive (risk-mitigation) measures from the findings of the specific

   run, keeping an observed match, a heuristic score and experimental validation

   strictly apart.

5. Record enough provenance (versions, parameters, file checksums, panel metadata)

   that a run can be audited and reproduced.

6. Deliver the results as a dashboard plus PDF, CSV and JSON downloads.



dsGuard does **not** model biological efficacy or safety; it ranks sequence

similarity.



---



## 2. System architecture and methodology



### 2.1 Components



| Layer | Location | Responsibility |

|---|---|---|

| Analysis engine | `Random/dsguard/` | Exact-match scan, scoring, Pareto selection, recommendation, per-gene and job-level output files |

| Backend | `backend/main.py`, `backend/reporting.py` | HTTP API, input handling, panel resolution, job execution and storage, report generation |

| Frontend | `frontend/src/` | React dashboard: input, results, downloads |

| Data | `reference_panels/`, `Random/data/` | Curated non-target panels and demo/target data |



The backend does not re-implement any analysis: it validates input, resolves panels,

invokes the engine as a subprocess, and generates reports from the engine's

`results.json`. This keeps a single source of truth for every number that reaches the

dashboard or a download.



### 2.2 Matching methodology



- **Strided seed-and-extend exact match.** A k-mer seed is placed on every k-th

  reference position; both the query and the reference are encoded into 2-bit arrays

  so a seed comparison is a packed integer comparison. Candidate hits are extended in

  both directions while bases are identical.

- **Guaranteed detection threshold.** With seed length `k`, every exact match of

  length `≥ 2k-1` is found. The default `k = 11` therefore guarantees complete

  enumeration of matches `≥ 21 nt`. Shorter matches are best-effort lower bounds and

  are labelled as such.

- **Both strands.** The reference is scanned on the forward and reverse-complement

  strands.

- **Ambiguity.** Bases other than A/C/G/T/U never match and cannot be spanned by a

  match, so an ambiguous base splits a potential match instead of being scored.

- **Batching.** Multiple target genes are concatenated with a separator so each panel

  is scanned once for all genes; matches are mapped back per gene.



### 2.3 Pipeline



1. Read and encode the target FASTA; count ambiguous bases.

2. For each requested panel, load the reference and scan all genes on both strands.

3. Slide a window (default 400 bp) over each gene; compute per-window potency and the

   longest panel match (per organism and overall).

4. Convert the longest match into a hazard score; flag windows whose match reaches the

   detection threshold.

5. Extract the Pareto front (maximise potency, minimise hazard) over non-flagged,

   non-ambiguous windows.

6. Recommend a window; design T7 primers for it.

7. Emit per-gene files, the ranked off-target table, the preventive measures, and the

   job-level `summary.csv` / `manifest.json` / `results.json`.



---



## 3. Implementation and scoring approach



### 3.1 Engine modules



| Module | Role |

|---|---|

| `seq.py` | 2-bit encoding, FASTA parsing, `Reference` (records + annotations), SHA-256 |

| `engine.py` | Seed-and-extend matching, both-strand scan, per-window longest match, coverage profile |

| `batch.py` | Multi-gene concatenation and match splitting |

| `sirna.py` | Reynolds (2004) and Ui-Tei (2004) design rules; per-window potency |

| `risk.py` | Logistic hazard (`L50`, `kappa`), graded hazard |

| `optimise.py` | Pareto front and recommendation selection |

| `primers.py` | T7 primer design for the recommended window |

| `cli.py` | Orchestration, ranking, provenance, all output files |



### 3.2 Hazard score



```

hazard(m) = 1 / (1 + exp(-(m - L50) / kappa))

```



where `m` is the longest perfect match (in nt) between a candidate window and any

reference record on either strand, clipped to the window. The score is in `[0, 1]` and

rises monotonically with `m`. It is a **heuristic**, explicitly recorded as

`calibrated: false`; it is not a probability of biological harm.



`L50 = 21.0` nt and `kappa = 2.0` nt are **placeholder constants**, not fitted values

(`l50_kappa_fitted: false`). The default threshold `max_hazard = 0.05` corresponds to

an equivalent match length of `L50 + kappa * ln(0.05/0.95) ≈ 15.1 nt`.



### 3.3 Potency score



```

potency = mean over sense-strand 21-mers in the window of

          0.5 * (reynolds / 7 + ui_tei / 4)

```



Reynolds (2004) contributes 7 of the 8 published criteria (the no-internal-repeat

criterion is not implemented); Ui-Tei (2004) contributes 4 rules. The score is used

only as a weak ordering prior and is not validated for insect Dicer-2 processing of

long dsRNA.



### 3.4 Pareto selection and recommendation



- A window is **flagged** when its longest match reaches the detection threshold.

- The Pareto front is computed over non-flagged, non-ambiguous windows, maximising

  potency and minimising hazard.

- The recommendation is the highest-potency window whose hazard is `≤ max_hazard`. If

  none qualifies, the lowest-risk window is returned as an explicit **fallback** and

  the report states that the threshold was exceeded.



### 3.5 Ranking and reporting



Ranked off-target rows are ordered by hazard (descending) with a deterministic

tie-break `(reference_panel, matched_organism, transcript_or_reference_id,

target_start, strand)`, so the dashboard, the ranked CSV and the CLI produce an

identical order. Coordinates are 0-based half-open internally and 1-based inclusive in

every output. The PDF, CSV and JSON exports are all generated from the same

`results.json`, so they cannot describe different runs.



### 3.6 Input methods



The API accepts an uploaded FASTA file, a pasted sequence, an NCBI accession (fetched

directly), or an NCBI gene name (searched, then candidate records are ranked so a

curated RefSeq transcript naming the gene is preferred). All four paths converge on a

validated FASTA that the engine consumes.



---



## 4. Test results and demonstrated functionality



Commands run from the project root (PowerShell, Windows 11):



| Command | Result |

|---|---|

| `.\.venv\Scripts\python.exe -m pytest -q` | **148 passed** in ~15 s |

| `npm run lint` (in `frontend/`) | clean (no ESLint errors) |

| `npm test` (in `frontend/`) | **28 passed**, 0 failed |

| `npm run build` (in `frontend/`) | success (`dist/` produced, 1932 modules) |



Test areas covered by the Python suite:



- **Matching correctness** — the strided scan is checked against a brute-force

  implementation, and multi-gene batch scanning against per-gene scanning.

- **Ranking** — hazard-descending order, deterministic tie-break, de-duplication of

  identical rows, and the counted-vs-listed match distinction.

- **Scoring** — potency rule evaluation and hazard behaviour, including that ambiguous

  windows are excluded.

- **Preventive measures** — each finding carries finding/why/action/priority/evidence;

  observed match, heuristic score and validation are kept in separate categories;

  missing panels and unannotated records are reported.

- **Reports and exports** — the ranked CSV exposes the required fields, the JSON

  export carries parameters/provenance/scoring, the PDF is a valid document, and the

  threshold agrees across CSV, JSON and the summary.

- **Job handling** — deletion guards (active jobs, unsafe paths, cleanup failure), and

  actionable 404s when structured results are missing.

- **Input modes (end-to-end)** — uploaded FASTA, pasted sequence, NCBI accession and

  NCBI gene name each create a job and schedule the worker with the correct input type;

  invalid input, missing panels, an unknown panel, missing engine and out-of-range

  parameters each produce a clear 4xx/5xx.

- **NCBI behaviour** — ranking of candidate records, accession-direct lookups,

  no-results vs transport-failure errors, and rejection of non-FASTA payloads, all

  offline with canned payloads.

- **Reproducibility** — two full engine runs on identical inputs produce byte-identical

  `summary.csv` and ranked CSVs, and identical recommendations, Pareto sets, parameters

  and panel provenance.



### Demonstrated functionality (verified live, this checkout)



- The engine was run against the installed rice panel

  (`reference_panels/rice.fasta`, 15 records, ~383 Mb, SHA-256

  `61856c0b…ba1e`) with `Random/data/targets/nilaparvata_lugens_13_cds.fasta`:

  **13 genes** were screened, the panel scan took ~78 s, and per-gene recommendations,

  Pareto fronts and ranked tables were produced. The recorded panel SHA-256 matched the

  file on disk, and the run recorded `l50_kappa_fitted: false`, the parameters, and the

  panel header annotation.

- The three download endpoints were exercised against a completed job and each returned

  a parseable file of the correct type (PDF, CSV, JSON).

- A run that had a configured-but-missing panel reported it as `skipped` with a reason

  and raised it as an unresolved risk, rather than reporting a clean screen.



---



## 5. Scientific limitations and future improvements



### Limitations (as implemented)



1. **Hazard is a heuristic.** `L50`/`kappa` are placeholders and no calibration data

   is bundled. The score orders candidates; it does not predict harm.

2. **Potency is borrowed.** Reynolds/Ui-Tei were derived from synthetic siRNA in

   mammalian cells and are not validated for insect RNAi.

3. **Exact matches only.** Mismatched and gapped off-targets, and miRNA-like

   partial-complementarity effects, are not detected.

4. **Coverage equals the installed panels.** Only `rice.fasta` is installed in this

   checkout; other configured panels are reported as skipped. Absence of a match is not

   evidence of absence in the organism.

5. **No experimental validation.** No wet-lab result is claimed or implied.



### Future improvements



- Fit `L50`/`kappa` to published non-target outcome data and record a fit version, so

  scores become calibrated rather than placeholder.

- Implement the remaining Reynolds criterion (no internal repeats) and validate the

  potency prior against insect RNAi data.

- Add mismatch/gap-tolerant matching (e.g. seed-and-extend with an edit budget) and

  report identity for near-matches.

- Install the remaining non-target panels (honeybee, *Cyrtorhinus lividipennis*,

  spider) and record their species and accessions per panel.

- Expose the equivalent-match length and calibration status more prominently in the

  dashboard.
