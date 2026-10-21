# dsGuard (7-day MVP)

Designs a dsRNA window against a pest gene and certifies it against a panel of
non-target organisms, using exact continuous-match profiling (strided seeding, k = 11,
provably finds every match >= 21 nt).

## Set up in VS Code
1. Install Python 3.10+ and the VS Code **Python** extension.
2. Open this folder in VS Code. Terminal -> New Terminal, then:
   ```
   python -m venv .venv
   .venv\Scripts\activate        (Windows)   |   source .venv/bin/activate (Mac/Linux)
   pip install -e ".[dev]"
   pytest
   ```
   Pick the `.venv` interpreter (Ctrl+Shift+P -> "Python: Select Interpreter").

## Try it with synthetic data (no downloads)
```
python scripts/make_demo_data.py demo
dsguard --target demo/target.fa --panel crop=demo/crop.fa --panel predator=demo/predator.fa --panel pollinator=demo/pollinator.fa --panel human=demo/human.fa --out demo_out
```

## Many genes at once (batch mode)
Put all the genes your sir gives you in ONE multi-FASTA file (one `>name` record per gene) and run:
```
dsguard --target demo/genes.fa --panel crop=demo/crop.fa --panel predator=demo/predator.fa --panel pollinator=demo/pollinator.fa --panel human=demo/human.fa --out demo_out
```
Each organism is scanned once for all genes together. You get `summary.csv` (one ranked row per gene)
and a folder per gene with `report.md` (the detail), `windows.csv`, `pareto.csv`, `hits.csv`, `profiles.csv`, `primers.json`.
Status values: `ok`, `no_window_free_of_21nt_match` (every window hits a non-target), `too_short`.

## Real run
Download from NCBI (RefSeq RNA / TSA FASTA): the target gene CDS (*N. lugens*) and a
transcriptome per organism (rice, predator *Propylea japonica* or a substitute, parasitoid,
honeybee, *Daphnia magna*, human). Record accession + licence for each. Then:
```
dsguard --target nlug_gene.fa --panel rice=rice_rna.fa.gz --panel human=human_rna.fa.gz --out results
```

## What is real and what is a placeholder
- Exactness (>= 21 nt) is tested against a brute-force implementation (tests/test_engine.py).
- `--l50 / --kappa` are PLACEHOLDERS (21 / 2). Fit them to published non-target data on Day 5-6, or say clearly they are defaults.
- Potency uses Reynolds (7 of 8 criteria) + Ui-Tei; not validated for insects. Say so in the report.
- Matches shorter than 21 nt are best-effort, so `max_match` is exact only when >= 21.
- Coordinates: 0-based half-open in code, 1-based inclusive in output files.

## Limits (print on every report)
dsGuard ranks hazard by perfect sequence similarity only. It does not model uptake, nuclease
stability, RISC loading or dose, and cannot detect miRNA-like partial-complementarity effects.
A low score is necessary, not sufficient, for non-target safety.
