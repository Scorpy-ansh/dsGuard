# dsGuard — Demonstration Readiness Report



**Audit date:** 2026-10-24

**Dataset:** `Random/data/targets/nilaparvata_lugens_13_cds.fasta` (13 *Nilaparvata lugens* CDS records)

**Reference panel installed in the dashboard:** `reference_panels/rice.fasta`

**Verdict:** The three input methods (Upload FASTA, Gene/NCBI search, Paste sequence) work correctly and

reproducibly on the project dataset. One real, reproducible defect was found and fixed (gene-alias search

returned an unrelated/absent record). **All 13 genes legitimately fail the configured hazard threshold and are

correctly reported as failed-threshold fallbacks.** No score, threshold, or warning was altered to change this.



---



## 1. What was run and how it was verified



| Activity | Command / artefact | Result |

|---|---|---|

| Full automated test suite | `pytest -q` (testpaths `Random`, `backend`) | **148 passed** |

| Frontend unit tests | `npm test` (`node --test`) | **28 passed** |

| Frontend lint | `npm run lint` (`eslint .`) | clean |

| Production build | `npm run build` (`vite build`) | built, 1932 modules, 268.18 kB JS |

| Live end-to-end (3 input methods + exports) | `e2e_api.py` against a live backend | **9/9 checks PASS, 3/3 exports PASS** |

| Per-job result isolation & display names | `job_isolation.py` against a live backend | **7/7 checks PASS** |

| Cross-format consistency | `consistency_api.py` on the live job | **20/20 checks PASS** |

| Live NCBI alias resolution | `alias_live.py` | **13/13 resolved to the correct accession** |



---



## 2. Data provenance (checksums)



| File | SHA-256 | Size |

|---|---|---|

| `Random/data/targets/nilaparvata_lugens_13_cds.fasta` | `9641e30fa0e04e2125455f5a10a07d1934e3613eb4650062d2b00e1b74c830bf` | — |

| `reference_panels/rice.fasta` (dashboard panel) | `61856c0bbd63757b2f6a40295f42b37e9cbdde8559364cefaac22be1ed82ba1e` | 387,564,294 B (~382.8 Mb, 15 chromosomes) |

| `Random/data/reference_panels/rice/rice_transcripts.fna` (validation panel) | `47e1568777b0e259393588a0e36335fff9f61034cf9cbd70d8ce3712be63bba9` | 115,778,255 B (107.9 Mb, 55,675 transcripts) |



The target checksum is byte-identical to the checksum recorded in `dsguard_test_output/manifest.json`

(2026-10-23), so the current results are directly comparable to the earlier project validation.



**Important:** the dashboard and the earlier validation use **two different rice panels**. This is the sole

reason the reported hit counts differ (§7).



---



## 3. Item 1 — FASTA parsing: PASS (13/13)



Every record was parsed and checked against the authoritative mapping in

`Random/data/targets/target_cds_retrieval_report.json`: alias, accession, sequence length and nucleotide

alphabet all agree; no record is missing, empty, duplicated, or contains invalid characters.



| # | Alias | Accession | Length (nt) | Bases valid | Report agrees |

|---|---|---|---|---|---|

| 1 | NlRan | KT313028.1 | 645 | yes | yes |

| 2 | NlSRP54 | XM_022328877.2 | 1530 | yes | yes |

| 3 | NlEcR | FJ263048.1 | 2067 | yes | yes |

| 4 | NlCHSA_NlCHS1 | JQ040014.1 | 4737 | yes | yes |

| 5 | Nlcar | JF915744.1 | 675 | yes | yes |

| 6 | Nltry | JF915745.1 | 825 | yes | yes |

| 7 | Nlsid-1 | JF915743.1 | 1875 | yes | yes |

| 8 | Nlaub | JF915742.1 | 2085 | yes | yes |

| 9 | NlRab2 | KT984796.1 | 639 | yes | yes |

| 10 | NlRab7 | KT984797.1 | 621 | yes | yes |

| 11 | NlRab6 | KT984812.1 | 627 | yes | yes |

| 12 | NlSar1 | KT984793.1 | 582 | yes | yes |

| 13 | NlSRbeta | KU568374.1 | 762 | yes | yes |



The engine independently reported **13 genes screened** for the uploaded file, confirming the app's own

parser agrees with this count.



**Note on display names.** The project headers are written as `>NlRan|accession=KT313028.1|cds_length=645` with

no space. Per the FASTA convention the identifier is the first whitespace-delimited token, so the whole

pipe-delimited string is the sequence name — correct parsing, but unreadable in a table. The dashboard now

shows the **clean alias** (`NlRan`) on screens, while **all reports and exports keep the full header with its

accession** (see §12a). Nothing about parsing changed; only the on-screen label.



---



## 4. Item 2 — Gene / NCBI search: one real bug found and fixed



### The defect

Searching a **project gene alias** (`NlRan`, `Nlcar`, `Nlsid-1`, `NlCHSA_NlCHS1`, …) returned **0 records**,

because these aliases are not present as gene names in NCBI Nucleotide. A naive esearch returned nothing, the

broad fallback could surface an unrelated record, and the previous code also returned NCBI's internal UID

instead of the accession. A teacher typing `NlRan` would therefore get no result, or the wrong one.



### The fix (`backend/main.py`)

1. `project_target_accessions()` reads the curated alias → accession map directly from the project FASTA headers.

2. `fetch_ncbi_fasta()` resolves a project alias to its curated accession first; genuine accessions are still

   fetched directly; other gene names still go through NCBI esearch.

3. The returned record id is now taken from the FASTA **header** (e.g. `KT313028.1`), not the internal UID.

4. `organism_is_annotated()` + an explicit rejection: if an organism was supplied and the fetched definition

   line does not name it, the request is **refused** with an explanatory error instead of silently screening an

   unrelated record.



### Live verification — 13/13 correct



| Alias typed | Resolved accession | Definition line (truncated) |

|---|---|---|

| NlRan | KT313028.1 | `KT313028.1 Nilaparvata lugens RanGTPase (Ran) mRNA, complete cds` |

| NlSRP54 | XM_022328877.2 | `XM_022328877.2 PREDICTED: Nilaparvata lugens signal recognition particle…` |

| NlEcR | FJ263048.1 | `FJ263048.1 Nilaparvata lugens ecdysteroid receptor A isoform (EcR-A) mRNA` |

| NlCHSA_NlCHS1 | JQ040014.1 | `JQ040014.1 Nilaparvata lugens chitin synthase 1 variant a (CHS1) mRNA` |

| Nlcar | JF915744.1 | `JF915744.1 Nilaparvata lugens carboxypeptidase mRNA, complete cds` |

| Nltry | JF915745.1 | `JF915745.1 Nilaparvata lugens trypsin-like serine protein mRNA` |

| Nlsid-1 | JF915743.1 | `JF915743.1 Nilaparvata lugens SID-1-like protein mRNA, complete cds` |

| Nlaub | JF915742.1 | `JF915742.1 Nilaparvata lugens aubergine (aub) mRNA, complete cds` |

| NlRab2 | KT984796.1 | `KT984796.1 Nilaparvata lugens Ras superfamily small GTPase Rab2 (Rab2)` |

| NlRab7 | KT984797.1 | `KT984797.1 Nilaparvata lugens Ras superfamily small GTPase Rab7 (Rab7)` |

| NlRab6 | KT984812.1 | `KT984812.1 Nilaparvata lugens Ras superfamily small GTPase (Rab6) mRNA` |

| NlSar1 | KT984793.1 | `KT984793.1 Nilaparvata lugens Ras superfamily small GTPase Sar1 (Sar1)` |

| NlSRbeta | KU568374.1 | `KU568374.1 Nilaparvata lugens SRP receptor beta subunit mRNA, complete cds` |



Both `NlRan` (alias) and `KT313028.1` (accession) were also exercised live through the API; both produced a

job titled `NCBI:KT313028.1 (…)` screening the single gene `KT313028.1`. 6 new regression tests cover this

behaviour, including a test that an unrelated organism is rejected.



---



## 5. Item 3 — Pasted sequence: PASS



* **Raw nucleotide bases** (no header) → accepted, labelled `pasted_target`, screened normally (1 gene, 1 match).

* **Pasted FASTA with header** → accepted, header preserved, screened normally.

* Invalid characters are still rejected with a clear message; the checksum-verified file path is unaffected.



---



## 6. Item 4 — Engine run and threshold diagnosis (dashboard rice panel)



**Parameters (defaults):** window 400 nt, k = 11, min_report_match 21 nt, L50 = 21.0, κ = 2.0, **max_hazard = 0.05**.

Hazard is `hazard(m) = 1 / (1 + exp(-(m − 21) / 2))`, where `m` is the longest perfect (gapless) match in the window.



**Result: 13/13 genes → `recommendation_status = fallback_threshold_exceeded`. 69 ranked rows, all `exceeds_threshold`.**



| Gene | Recommended window | Best match `m` (nt) | Hazard | Hazard ≤ 0.05? | Hits ≥ 21 nt | Matches below floor |

|---|---|---|---|---|---|---|

| Nlsid-1 | 505–904 | 17 | 0.119203 | no | 0 | 46,071 |

| NlCHSA_NlCHS1 | 1099–1498 | 17 | 0.119203 | no | 5 | 107,571 |

| NlRab7 | 14–413 | 17 | 0.119203 | no | 0 | 11,645 |

| Nlaub | 866–1265 | 17 | 0.119203 | no | 0 | 37,093 |

| NlRab6 | 2–401 | 18 | 0.182426 | no | 0 | 14,423 |

| NlSRP54 | 285–684 | 18 | 0.182426 | no | 1 | 30,881 |

| NlEcR | 55–454 | 18 | 0.182426 | no | 1 | 43,647 |

| Nltry | 426–825 | 18 | 0.182426 | no | 0 | 16,940 |

| NlRab2 | 1–400 | 18 | 0.182426 | no | 0 | 15,597 |

| Nlcar | 270–669 | 19 | 0.268941 | no | 0 | 17,962 |

| NlSRbeta | 355–754 | 19 | 0.268941 | no | 0 | 21,413 |

| NlRan | 45–444 | 19 | 0.268941 | no | 1 | 14,620 |

| NlSar1 | 1–400 | 19 | 0.268941 | no | 0 | 10,820 |



Ranked-row match-length distribution across the job: 17 nt → 27, 18 nt → 27, 19 nt → 7, 21 nt → 2,

22 nt → 4, 23 nt → 1, 24 nt → 1. Total guaranteed (≥ 21 nt) hits: **8**, in 4 genes

(NlCHSA_NlCHS1 = 5, NlSRP54 = 1, NlEcR = 1, NlRan = 1).



### Exact cause of every failure

`max_hazard = 0.05` is equivalent to a perfect match of only **≈ 15.11 nt**

(inverting the sigmoid: `m = 21 + 2·ln(19)`). Any window containing a **> 15.11 nt** perfect match therefore

exceeds the threshold by construction. Against a **382.8 Mb** whole-genome reference, a 400 nt window almost

inevitably contains a chance 16–19 nt perfect match: the expected number of length-`L` matches is

≈ `2 · (400 − L + 1) · 3.828×10⁸ / 4^L`, which is ≫ 1 for `L = 17` and still ≳ 1 for `L = 19`.

The observed best matches (17–19 nt) are exactly what this expectation predicts.



**This is not a defect.** The optimizer already chose, for every gene, a window whose hazard is *lower* than

every ≥ 21 nt hit window (e.g. `hazard(17) = 0.119` vs `hazard(24) = 0.818`). No window can reach 0.05;

the correct behaviour is the transparent fallback that is reported.



---



## 7. Item 5 — Comparison with the existing project validation run



| | Existing validation (`dsguard_test_output/`, 2026-10-23) | Current run (dashboard) |

|---|---|---|

| Panel | `rice_transcripts.fna` — 55,675 transcripts, 107.9 Mb | `rice.fasta` — 15 chromosomes, ~382.8 Mb |

| Target checksum | `9641e30f…` | `9641e30f…` (identical) |

| Genes screened | 13/13 `ok` | 13/13 `ok` |

| Guaranteed hits (≥ 21 nt) | **3** across 2 genes (NlCHSA_NlCHS1 = 2, NlRan = 1) | **8** across 4 genes (NlCHSA_NlCHS1 = 5, NlSRP54 = 1, NlEcR = 1, NlRan = 1) |

| Recommended-window best match | 16–19 nt | 17–19 nt |

| Recommendation status | all `fallback_threshold_exceeded` | all `fallback_threshold_exceeded` |

| `DATASET_STATUS.md` claim | "3 matches across two target genes" | reproduced exactly (see the transcript-panel run) |



Re-running the engine on the validation transcript panel reproduced **3 matches across two genes**, matching

`Random/data/reference_panels/rice/DATASET_STATUS.md` exactly — the engine is cross-validated. The extra hits in

the dashboard run are fully explained by the panel: the genome adds introns, intergenic and repeat sequence

(≈ 3.5× more bases), so more chance matches are expected. **No ranking, scoring or parsing bug was demonstrated.**



---



## 8. Items 6 & 7 — Integrity of the result (no fabrication)



* The hazard threshold was **not** raised; scoring formula, L50/κ and defaults are untouched.

* No candidate was re-labelled as safe. All 13 genes remain `fallback_threshold_exceeded`, and the dashboard

  labels them as failed-threshold fallbacks.

* No warning was suppressed and no result was invented. `matches_below_report_floor` (10,820–107,571 per gene)

  is reported honestly, and the UI explains that `HITS` counts only exact matches ≥ 21 nt while ranked rows may

  include shorter sub-floor matches.

* Fallback behaviour is preserved and explicit rather than silently reporting success.



---



## 9. Item 8 — Cross-format consistency: 20/20 PASS



Verified on the **actual downloaded bytes** of the live upload job (13 genes):



* CSV header: 18 unique, non-empty columns; first column `gene`; documented comment preamble declaring

  `configured_max_hazard = 0.05` and the exact-match/heuristic caveats.

* CSV rows = dashboard ranked rows = **69**; match length, hazard, threshold status, reference record and strand

  form an **identical multiset**.

* JSON export is tagged `dsguard-results-v1`; its `summary`, `genes`, `parameters`, `preventive_measures` and

  `limitations` are **identical** to the dashboard payload.

* CSV hazard values equal dashboard hazard values; PDF is a valid document (33,901 B).

* Every gene carries an explicit `recommendation_status`; every ranked row carries an explicit

  `threshold_status`; preventive measures and limitations are present in the payload and rendered in the PDF.



---



## 10. Item 9 — Tests and live end-to-end



**Automated:** `pytest` **148 passed** (backend + engine); `npm test` **28 passed**; `eslint` clean;

`vite build` succeeded (1932 modules, 268.18 kB JS / 82.15 kB gzip, 13.73 kB CSS).



**Live end-to-end (all three input methods, against a real backend):**



| # | Case | Result |

|---|---|---|

| 0 | Parse project FASTA | PASS — 13 records, names/accessions/lengths verified |

| 1 | Upload FASTA | PASS — 13 genes, 8 hits ≥ 21 nt, all fallback; rec best match 17/18/19 |

| 2 | Paste raw nucleotide bases | PASS — 1 gene `pasted_target`, 1 match |

| 3 | Gene search by alias `NlRan` | PASS — resolved live to `KT313028.1` |

| 4 | Search by accession `KT313028.1` | PASS — resolved to the same record |

| 5 | Paste FASTA with header | PASS — header preserved, screened |

| 6–8 | Export `results.csv`, `results.json`, `report.pdf` | PASS — correct content and attachment filenames |



---



## 11. Item 10 — Per-gene pass/fail



"Parsing/search" = the record is parsed and (where applicable) resolved to the correct *N. lugens* accession.

"Eligibility" = whether the recommended window satisfies `hazard ≤ 0.05`. **An eligibility FAIL is a correct

result, not a bug** — it means no candidate met the threshold.



| Gene | Parsing/search | Eligibility (hazard ≤ 0.05) | Cause |

|---|---|---|---|

| NlRan | PASS | FAIL | best window match 19 nt → hazard 0.268941 |

| NlSRP54 | PASS | FAIL | best window match 18 nt → hazard 0.182426 |

| NlEcR | PASS | FAIL | best window match 18 nt → hazard 0.182426 |

| NlCHSA_NlCHS1 | PASS | FAIL | best window match 17 nt → hazard 0.119203 |

| Nlcar | PASS | FAIL | best window match 19 nt → hazard 0.268941 |

| Nltry | PASS | FAIL | best window match 18 nt → hazard 0.182426 |

| Nlsid-1 | PASS | FAIL | best window match 17 nt → hazard 0.119203 |

| Nlaub | PASS | FAIL | best window match 17 nt → hazard 0.119203 |

| NlRab2 | PASS | FAIL | best window match 18 nt → hazard 0.182426 |

| NlRab7 | PASS | FAIL | best window match 17 nt → hazard 0.119203 |

| NlRab6 | PASS | FAIL | best window match 18 nt → hazard 0.182426 |

| NlSar1 | PASS | FAIL | best window match 19 nt → hazard 0.268941 |

| NlSRbeta | PASS | FAIL | best window match 19 nt → hazard 0.268941 |



**Summary:** parsing/search 13/13 PASS · eligibility 0/13 (all correctly reported as failed-threshold fallbacks).



---



## 12. Changes made



| File | Change |

|---|---|

| `backend/main.py` | Curated project-alias → accession resolution; organism verification with explicit rejection of unrelated records; return the FASTA accession (not the internal UID) as the record id. |

| `frontend/src/App.jsx` | NCBI input help text now documents that project gene aliases resolve to their curated record. |

| `backend/test_ncbi.py` | Updated 3 existing tests to the header-based record id; added 6 tests (alias map parsing, missing-list case, alias→accession resolution, unrelated-organism rejection, case-insensitivity). |



No changes were made to the hazard formula, L50, κ, the threshold, ranking logic, reference data or exports.



### 12a. Demo-readiness pass (2026-10-24, this revision)



A second, minimal pass targeted only presentation and result-isolation robustness. No scoring, threshold,

ranking, reference data or export content was touched.



| Item | File | Change |

|---|---|---|

| Clean gene names on screens | `frontend/src/analysis.js` | Added `geneDisplayName()` (takes the alias before the first `|`) and `jobDisplayName()` (same rule for job titles). Pure functions. |

| Clean gene names on screens | `frontend/src/App.jsx` | Results summary, ranked table, Pareto options/labels and hits now render `geneDisplayName(...)`; job/history labels render `jobDisplayName(job)`. |

| Reports keep provenance | `frontend/src/App.jsx`, `backend/main.py` | The exported gene column and the backend CSV/JSON/PDF keep the **full** identifier/accession unchanged. |

| Stale results (item 4) | `frontend/src/App.jsx` | Monotonic `resultsRequestRef` / `reportRequestRef`; a response is applied only if it belongs to the latest request. `deleteJob` bumps the counter so a deleted job cannot reappear. |

| UI / error text (item 3) | `frontend/src/App.jsx`, `frontend/index.html` | Backend-unreachable message now derives the host from `API` instead of a hardcoded port; page `<title>` set to `dsGuard - dsRNA off-target screening`. |

| Regression tests | `frontend/src/analysis.test.js` | +6 tests for `geneDisplayName` / `jobDisplayName` (pipe alias, plain id, empty/null, pasted header, file name, NCBI label). |

| Regression tests | `backend/test_screen_inputs.py` | +1 test: a pasted FASTA keeps its own header (previously only pasted raw bases were covered). |

| Backend URL (demo blocker) | `frontend/src/App.jsx` | The dashboard's API base had drifted to `http://127.0.0.1:8001`, which does not match the documented startup port (8000), so a demo started per the README showed "Cannot reach the backend". Restored the default to `http://127.0.0.1:8000` and made it overridable via `VITE_API_BASE`. |



Verified after these changes: `pytest` 148 passed · `npm test` 28 passed · `eslint` clean ·

`vite build` OK · live E2E 9/9 · per-job isolation & display 7/7. The item-4 isolation check interleaves

requests for two different jobs and confirms each response still carries its **own** gene set (upload job =

13 genes; paste job = the single `NlRan` record), with the full header `NlRan|accession=KT313028.1|cds_length=645`

kept in the analysis/report payloads and the clean alias `NlRan` served to the screens.



### 12b. Final verification run (2026-10-24)



Re-verified end to end after the demo-readiness pass. **No reproducible software defect was found; no score,

threshold or reference data was changed.**



| Check | Result |

|---|---|

| 13 genes parsed / identified | PASS — 13 FASTA headers are the 13 `summary` genes verbatim; alias **and** accession sets identical; `input.records = 13` |

| Determinism (two engine runs, same inputs) | PASS — `summary.csv` byte-identical; `summary`, `ranked_offtargets`, `recommendation`, `pareto` identical; distinct `analysis_id` |

| Ranked matches / threshold / summary / recommendations consistent | PASS — 69 ranked rows; per-gene `rank` = 1..n; hazard descending; `recommendation` (by gene) == `summary`; `threshold_status` ∈ {meets, exceeds}; 13/13 `fallback_threshold_exceeded`; rec best match 17–19 nt (0.05 unchanged) |

| Missing panels reported as missing | PASS — live `/api/panels` marks rice available, and `cyrtorhinus_lividipennis` / `honeybee` / `spider` `available=false` with a "was not found" reason; engine records a skipped panel as `not_searched`; selecting only an uninstalled panel → HTTP 400 |

| Invalid inputs → clear error, no crash | PASS — live backend: 6 malformed requests each returned HTTP 400 with a specific message; `/api/health` still 200 afterwards |

| Input modes | PASS — automated: uploaded FASTA, pasted raw bases, pasted FASTA (new test), gene alias, accession, plus invalid-input cases (`backend/test_screen_inputs.py`, `backend/test_ncbi.py`) |

| Exports correspond to the selected job | PASS — `backend/test_exports.py` drives `/results`, CSV, JSON and PDF from the job's own `results.json`; 409 while running, 404 for unknown/missing |

| Tests / lint / build | PASS — `pytest` 148, `npm test` 28, `eslint` clean, `vite build` 268.18 kB |



**Known environmental blocker (not a code defect):** MongoDB Atlas is unreachable from this network

(`SSL: TLSV1_ALERT_INTERNAL_ERROR`), so a job cannot be *created* on the real backend here. The live checks

above therefore cover every path that runs before persistence (health, panels, request validation) plus a real

offline engine run; the DB-backed job flow (create → poll → results → exports) is covered by the automated

backend tests, which call the production handlers with a stubbed job collection. Once Atlas is reachable,

re-run the live `create → results → export` flow to confirm the same behaviour on the wire.



---



## 13. Remaining scientific limitations



1. **Placeholder parameters.** `L50 = 21`, `κ = 2` are not fitted to *N. lugens*; the hazard value is a

   heuristic ranking aid, not a calibrated biological probability.

2. **Exact-match model only.** The engine detects perfect, gapless matches; mismatches and gaps score 0.

   miRNA-like partial-complementarity effects are not modelled.

3. **Potency rules.** Potency is derived from mammalian siRNA rules and is unvalidated for insect Dicer-2.

4. **Short-match detection is best-effort.** Matches of 21 nt (2k−1) or more are guaranteed; the 16–19 nt

   matches driving the recommended-window hazard are best-effort finds, not exhaustive.

5. **Sequence similarity is not biological harm.** A low score is necessary, not sufficient, for non-target

   safety; uptake, nuclease stability, RISC loading and dose are not modelled.

6. **Panel scope.** The installed dashboard panel is *Oryza sativa* only; the other catalogued panels

   (*Cyrtorhinus lividipennis*, honeybee, spider) are reported as unavailable and must not be read as zero risk.

7. **Two rice panels exist** (whole-genome vs transcriptome). They give different hit counts; neither is

   "wrong", but state which one is used when comparing runs.

8. **Above 0.05 does not mean unsafe**, and below 0.05 would not mean safe — it means the window met the

   configured, admittedly placeholder, criterion.



---



## 14. Recommended cleanup before the demonstration



The repository root contains leftover scratch artefacts that are not part of the deliverable, e.g.

`_tmp_cli_log.txt`, `_tmp_cli_out/`, `_tmp_demo_out/`, `_tmp_report.pdf`, `_tmp_verify.py`,

`_tmp_old_cli.py`, `_tmp_old_optimise.py`, `_tmp_old_risk.py`, `_tmp_orig_readme.md`, and a stray

`dsguard_test_output/`. They are harmless but untidy; they were **not** deleted automatically because

`_tmp_orig_readme.md` looks like a deliberate backup. Remove them only if you no longer need them.



Also note there is a **duplicate module set at the repository root** (`engine.py`, `cli.py`, `risk.py`, `seq.py`,

`sirna.py`, `batch.py`, `brute.py`, `optimise.py`, `primers.py`, `test_*.py`) alongside the real package in

`Random/dsguard/`. The test configuration only collects `Random` and `backend`, so the root copies are dead

weight and could confuse a reviewer.
