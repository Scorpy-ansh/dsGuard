
import csv
from pathlib import Path

BASE = Path(__file__).resolve().parent
OUT_DIR = BASE / "rice_screen_test_out"
SUMMARY = OUT_DIR / "summary.csv"
THRESHOLD = 21

# Load the summary containing each gene's recommended coordinates.
with SUMMARY.open(newline="", encoding="utf-8-sig") as f:
    summary_rows = list(csv.DictReader(f))

passed = 0
failed = 0

for summary in summary_rows:
    gene = summary["gene"].split("|")[0]
    rec_start = int(summary["rec_start"])
    rec_end = int(summary["rec_end"])

    folder = next(
        (p for p in OUT_DIR.iterdir()
         if p.is_dir() and p.name.split("_accession_")[0] == gene),
        None,
    )

    if folder is None:
        print(f"FAIL {gene}: per-gene output folder not found")
        failed += 1
        continue

    windows_file = folder / "windows.csv"
    if not windows_file.exists():
        print(f"FAIL {gene}: windows.csv not found")
        failed += 1
        continue

    # The first line is a comment, so skip it before parsing CSV columns.
    with windows_file.open(newline="", encoding="utf-8-sig") as f:
        next(f, None)
        window_rows = list(csv.DictReader(f))

    matching = [
        row for row in window_rows
        if int(row["start"]) == rec_start
        and int(row["end"]) == rec_end
    ]

    if len(matching) != 1:
        print(
            f"FAIL {gene}: expected exactly one row for "
            f"{rec_start}-{rec_end}, found {len(matching)}"
        )
        failed += 1
        continue

    row = matching[0]
    max_match = int(row["max_match"])
    flag = int(row["binary_flag_ge21"])

    if max_match < THRESHOLD and flag == 0:
        print(
            f"PASS {gene}: {rec_start}-{rec_end}, "
            f"max_match={max_match}, flag={flag}"
        )
        passed += 1
    else:
        print(
            f"FAIL {gene}: {rec_start}-{rec_end}, "
            f"max_match={max_match}, flag={flag}"
        )
        failed += 1

print(f"\nRESULT: {passed}/{len(summary_rows)} recommended windows passed; "
      f"{failed} failed.")

if failed:
    raise SystemExit(1)
