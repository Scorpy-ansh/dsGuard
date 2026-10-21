
import csv
from pathlib import Path
from Bio import SeqIO

BASE = Path(__file__).resolve().parent
TARGET_FASTA = BASE / "data" / "targets" / "nilaparvata_lugens_13_cds.fasta"
RICE_FASTA = BASE / "data" / "reference_panels" / "rice" / "rice_transcripts.fna"
OUT_DIR = BASE / "rice_screen_test_out"

# Target records are identified by the gene name before the first "|".
targets = {}
for record in SeqIO.parse(str(TARGET_FASTA), "fasta"):
    gene = record.id.split("|")[0]
    targets[gene] = str(record.seq).upper()

# Rice records are identified by their accession in the FASTA header.
rice = {
    record.id: str(record.seq).upper()
    for record in SeqIO.parse(str(RICE_FASTA), "fasta")
}

total = 0
passed = 0
failed = 0

for csv_path in sorted(OUT_DIR.glob("*/hits.csv")):
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        continue

    gene = csv_path.parent.name.split("_accession_")[0]

    for row in rows:
        total += 1
        transcript = row["transcript"]
        start = int(row["gene_start"])
        end = int(row["gene_end"])
        rice_start = int(row["transcript_pos_1based"])
        length = int(row["length"])
        strand = row["strand"]

        errors = []

        if gene not in targets:
            errors.append(f"target gene not found: {gene}")
        if transcript not in rice:
            errors.append(f"rice transcript not found: {transcript}")

        if not errors:
            target_seq = targets[gene][start - 1:end]
            rice_seq = rice[transcript][rice_start - 1:rice_start - 1 + length]

            if end - start + 1 != length:
                errors.append("target coordinates do not match reported length")
            if len(target_seq) != length:
                errors.append("target slice has unexpected length")
            if len(rice_seq) != length:
                errors.append("rice slice has unexpected length")
            if strand == "+" and target_seq != rice_seq:
                errors.append("target and rice sequences differ on reported + strand")
            elif strand == "-" and target_seq != str(
                __import__("Bio.Seq", fromlist=["Seq"]).Seq(rice_seq).reverse_complement()
            ):
                errors.append("target and rice sequences differ on reported - strand")
            elif strand not in ("+", "-"):
                errors.append(f"unrecognized strand: {strand}")

        if errors:
            failed += 1
            print(f"FAIL: {gene} vs {transcript}: {'; '.join(errors)}")
        else:
            passed += 1
            print(
                f"PASS: {gene} vs {transcript} | "
                f"target {start}-{end} | rice {rice_start}-"
                f"{rice_start + length - 1} | {length} nt | {strand}"
            )
            print(f"      Sequence: {target_seq}")

print(f"\nRESULT: {passed}/{total} hits verified; {failed} failed.")

if failed:
    raise SystemExit(1)
