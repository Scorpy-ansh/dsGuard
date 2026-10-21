
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parent
RAW_DIR = ROOT / "data" / "targets" / "ncbi_records"
FASTA = ROOT / "data" / "targets" / "nilaparvata_lugens_13_cds.fasta"
REPORT = ROOT / "data" / "targets" / "target_cds_validation.json"

VALID_BASES = set("ACGTNRYKMSWBDHV")


def find_genbank_file(accession):
    """Find the saved GenBank record with or without accession version."""
    base_accession = accession.split(".")[0]

    matches = [
        path
        for path in RAW_DIR.glob("*.gb")
        if path.stem == accession or path.stem == base_accession
    ]

    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one saved GenBank file for {accession}; "
            f"found {len(matches)}"
        )

    return matches[0]


def validate_record(fasta_record):
    parts = fasta_record.id.split("|")
    gene = parts[0]

    accession_part = next(
        (part for part in parts if part.startswith("accession=")),
        None,
    )

    item = {
        "gene": gene,
        "fasta_id": fasta_record.id,
        "length": len(fasta_record.seq),
    }

    if accession_part is None:
        item["status"] = "VALIDATION_ERROR"
        item["error"] = "FASTA header has no accession field"
        return item

    accession = accession_part.split("=", 1)[1]
    sequence = str(fasta_record.seq).upper()

    item["accession"] = accession
    item["valid_alphabet"] = (
        bool(sequence) and not (set(sequence) - VALID_BASES)
    )
    item["length_divisible_by_3"] = len(sequence) % 3 == 0

    try:
        gb_path = find_genbank_file(accession)
        record = SeqIO.read(gb_path, "genbank")

        item["source_file"] = gb_path.name
        item["source_record_id"] = record.id
        item["organism"] = record.annotations.get(
            "organism", "unknown"
        )

        cds_features = [
            feature
            for feature in record.features
            if feature.type == "CDS"
        ]
        item["cds_feature_count"] = len(cds_features)

        if len(cds_features) != 1:
            item["status"] = "REVIEW_CDS_FEATURE_COUNT"
            item["message"] = (
                "Expected exactly one CDS feature in the source record."
            )
            return item

        feature = cds_features[0]
        extracted = str(feature.extract(record.seq)).upper()
        item["matches_source_cds"] = extracted == sequence

        codon_start = int(
            feature.qualifiers.get("codon_start", ["1"])[0]
        )
        table = int(
            feature.qualifiers.get("transl_table", ["1"])[0]
        )

        item["codon_start"] = codon_start
        item["translation_table"] = table

        coding_sequence = extracted[codon_start - 1:]
        translated = str(
            Seq(coding_sequence).translate(table=table)
        )

        annotated_translation = "".join(
            feature.qualifiers.get("translation", [""])[0].split()
        ).upper()

        item["has_translation_annotation"] = bool(annotated_translation)
        item["translation_matches_annotation"] = (
            translated.rstrip("*") == annotated_translation.rstrip("*")
            if annotated_translation
            else None
        )
        item["internal_stop_codons"] = translated.rstrip("*").count("*")

        location_text = str(feature.location)
        item["partial_location"] = (
            "<" in location_text or ">" in location_text
        )

        issues = []

        if item["organism"] != "Nilaparvata lugens":
            issues.append("ORGANISM_MISMATCH")

        if not item["valid_alphabet"]:
            issues.append("INVALID_BASES_OR_EMPTY_SEQUENCE")

        if not item["matches_source_cds"]:
            issues.append("SOURCE_SEQUENCE_MISMATCH")

        if not item["length_divisible_by_3"]:
            issues.append("LENGTH_NOT_DIVISIBLE_BY_3")

        if item["partial_location"]:
            issues.append("PARTIAL_CDS_LOCATION")

        if not annotated_translation:
            issues.append("MISSING_TRANSLATION_ANNOTATION")
        elif not item["translation_matches_annotation"]:
            issues.append("TRANSLATION_MISMATCH")

        if item["internal_stop_codons"] > 0:
            issues.append("INTERNAL_STOP_CODON")

        item["issues"] = issues
        item["status"] = "PASS_BASIC_CHECKS" if not issues else "REVIEW"

    except Exception as exc:
        item["status"] = "VALIDATION_ERROR"
        item["error"] = f"{type(exc).__name__}: {exc}"

    return item


def main():
    if not FASTA.is_file():
        raise FileNotFoundError(f"Target FASTA not found: {FASTA}")

    if not RAW_DIR.is_dir():
        raise FileNotFoundError(
            f"NCBI record folder not found: {RAW_DIR}"
        )

    fasta_records = list(SeqIO.parse(FASTA, "fasta"))
    results = []
    seen_genes = set()

    for fasta_record in fasta_records:
        gene = fasta_record.id.split("|")[0]

        if gene in seen_genes:
            results.append({
                "gene": gene,
                "status": "REVIEW_DUPLICATE_GENE",
                "error": "Duplicate gene identifier in target FASTA",
            })
            continue

        seen_genes.add(gene)
        results.append(validate_record(fasta_record))

    report = {
        "target_species": "Nilaparvata lugens",
        "records_checked": len(fasta_records),
        "unique_gene_identifiers": len(seen_genes),
        "pass_basic_checks": sum(
            result.get("status") == "PASS_BASIC_CHECKS"
            for result in results
        ),
        "review_or_error_count": sum(
            result.get("status") != "PASS_BASIC_CHECKS"
            for result in results
        ),
        "results": results,
        "limitations": (
            "Basic sequence checks only. Passing does not prove that an "
            "accession is the intended RNAi target region or establish "
            "biological efficacy or safety."
        ),
    }

    REPORT.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )

    print(f"Records checked: {len(fasta_records)}")
    print(f"Passed basic checks: {report['pass_basic_checks']}")
    print(f"Review/errors: {report['review_or_error_count']}")

    for result in results:
        print(
            f"{result.get('gene', 'UNKNOWN'):16} "
            f"length={result.get('length', 0):5} "
            f"translation_match="
            f"{result.get('translation_matches_annotation', 'N/A')} "
            f"status={result['status']}"
        )

        if result.get("issues"):
            print("  Issues:", ", ".join(result["issues"]))

        if result.get("error"):
            print("  Error:", result["error"])

    print(f"\nReport saved to: {REPORT}")


if __name__ == "__main__":
    main()
