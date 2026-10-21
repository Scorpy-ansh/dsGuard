# dsGuard Rice Reference Dataset

## Dataset provenance

- **Organism:** *Oryza sativa* Japonica Group
- **Assembly accession:** GCF_034140825.1
- **Assembly name:** AGIS1.0
- **Annotation release:** RS_2024_06
- **Source:** Official NCBI Datasets RNA FASTA
- **Project file:** `rice_transcripts.fna`
- **Transcript records:** 55,675
- **Total sequence bases:** 107,882,961
- **Minimum transcript length:** 52 nt
- **Maximum transcript length:** 23,949 nt
- **Duplicate IDs:** 0
- **Empty sequences:** 0
- **Invalid nucleotide characters:** 0
- **SHA-256:** `47e1568777b0e259393588a0e36335fff9f61034cf9cbd70d8ce3712be63bba9`

## Sequence integrity observations

- 1,255 records have trailing A-runs of at least 10 bases.
- 400 records have trailing A-runs of at least 20 bases.
- Maximum observed trailing A-run: 61 bases.
- No sequences were trimmed.
- The integrity report is stored in `integrity_report.json`.

## dsGuard screening integration

- The rice transcript FASTA was used as a non-target reference panel.
- Screening run parameters: window size 400 nt, k = 11, exact-match threshold = 21 nt.
- Target input: `data/targets/nilaparvata_lugens_13_cds.fasta`.
- Target records processed successfully: 13 of 13.
- Reported rice matches: 3 across two target genes.
- The three reported matches were independently checked against the original target and rice FASTA sequences and their recorded coordinates.
- All 13 recommended-window records passed validation against the stored `windows.csv` results: their reported `max_match` was below 21 and `binary_flag_ge21` was zero.
- The matching engine test suite passed 29 tests, including randomized comparison with a brute-force reference implementation.

## Known limitations

- Sequence integrity checks do not establish biological completeness or transcript annotation correctness.
- The recommended-window validation checks stored screening output; it is not an independent full-panel rescan of every recommended window.
- The exact-match screen does not establish biological non-target safety.
- The screening model does not account for uptake, nuclease stability, RISC loading, dose, or miRNA-like partial-complementarity effects.
- Potency parameters L50 = 21 and kappa = 2 are placeholders and have not been fitted to the target organism.
- A low hazard score is not sufficient evidence of biological safety.

## Intended use

This dataset is a non-target transcript reference panel for dsGuard screening. It is not the brown planthopper target-gene dataset.

## Reproducibility

The screening output is stored in `rice_screen_test_out/`. The run parameters, input paths, input checksums, software versions, and stated limitations are recorded in `rice_screen_test_out/manifest.json`.

The target retrieval and validation reports are stored under `data/targets/`.
