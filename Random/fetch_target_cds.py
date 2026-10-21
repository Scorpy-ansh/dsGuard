
import json
import time
from io import StringIO
from pathlib import Path

from Bio import Entrez, SeqIO

PROJECT = Path(__file__).resolve().parent
OUT_DIR = PROJECT / "data" / "targets"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_FASTA = OUT_DIR / "nilaparvata_lugens_13_cds.fasta"
REPORT = OUT_DIR / "target_cds_retrieval_report.json"
RAW_DIR = OUT_DIR / "ncbi_records"
RAW_DIR.mkdir(parents=True, exist_ok=True)

ACCESSIONS = [
    ("NlRan", "KT313028"),
    ("NlSRP54", "XM_022328877.2"),
    ("NlEcR", "FJ263048.1"),
    ("NlCHSA_NlCHS1", "JQ040014.1"),
    ("Nlcar", "JF915744"),
    ("Nltry", "JF915745"),
    ("Nlsid-1", "JF915743"),
    ("Nlaub", "JF915742"),
    ("NlRab2", "KT984796"),
    ("NlRab7", "KT984797"),
    ("NlRab6", "KT984812"),
    ("NlSar1", "KT984793"),
    ("NlSRbeta", "KU568374"),
]


def main():
    email = input("Enter your email for NCBI requests: ").strip()
    if "@" not in email:
        raise ValueError("Enter a valid email address.")

    Entrez.email = email

    output_records = []
    results = []

    for gene, accession in ACCESSIONS:
        entry = {
            "gene": gene,
            "requested_accession": accession,
            "status": "not_processed",
        }

        try:
            print(f"Fetching {gene}: {accession} ...")

            with Entrez.efetch(
                db="nuccore",
                id=accession,
                rettype="gb",
                retmode="text",
            ) as handle:
                text = handle.read()

            if not text.strip() or text.lstrip().startswith("Error"):
                raise ValueError("NCBI returned no GenBank record.")

            # Preserve the source record for audit and later review.
            raw_path = RAW_DIR / f"{accession}.gb"
            raw_path.write_text(text, encoding="utf-8")
            entry["source_record_file"] = str(raw_path)

            records = list(SeqIO.parse(StringIO(text), "genbank"))
            if len(records) != 1:
                raise ValueError(
                    f"Expected one GenBank record; got {len(records)}."
                )

            record = records[0]
            entry["returned_accession"] = record.id
            entry["record_length"] = len(record.seq)
            entry["organism"] = record.annotations.get(
                "organism", "unknown"
            )

            if "Nilaparvata lugens" not in entry["organism"]:
                entry["status"] = "REVIEW_ORGANISM_MISMATCH"
                entry["message"] = (
                    "Record organism does not match the target species."
                )
                results.append(entry)
                time.sleep(0.4)
                continue

            cds_features = [
                feature
                for feature in record.features
                if feature.type == "CDS"
            ]

            if not cds_features:
                entry["status"] = "NO_ANNOTATED_CDS"
                entry["message"] = (
                    "No CDS feature found. Full sequence was not "
                    "substituted for a CDS."
                )
                results.append(entry)
                time.sleep(0.4)
                continue

            if len(cds_features) != 1:
                entry["status"] = "REVIEW_MULTIPLE_CDS"
                entry["cds_feature_count"] = len(cds_features)
                entry["message"] = (
                    "Multiple CDS features found; selection requires review."
                )
                results.append(entry)
                time.sleep(0.4)
                continue

            feature = cds_features[0]
            sequence = str(feature.extract(record.seq)).upper()

            # Keep partial CDS visible for review; do not silently
            # label it a complete coding sequence.
            qualifiers = feature.qualifiers
            codon_start = qualifiers.get("codon_start", ["1"])[0]
            translation = qualifiers.get("translation", [None])[0]

            partial = (
                feature.location.start.__class__.__name__ == "BeforePosition"
                or feature.location.end.__class__.__name__ == "AfterPosition"
                or "<" in str(feature.location)
                or ">" in str(feature.location)
            )

            if not sequence or set(sequence) - set("ACGTNRYKMSWBDHV"):
                entry["status"] = "INVALID_CDS_SEQUENCE"
                entry["message"] = (
                    "Empty sequence or unsupported nucleotide characters."
                )
                results.append(entry)
                time.sleep(0.4)
                continue

            record_id = (
                f"{gene}|accession={record.id}|"
                f"cds_length={len(sequence)}"
            )
            output_records.append((record_id, sequence))

            entry["cds_length"] = len(sequence)
            entry["codon_start"] = codon_start
            entry["has_translation_annotation"] = translation is not None
            entry["partial_location"] = partial
            entry["status"] = (
                "REVIEW_PARTIAL_CDS" if partial else "CDS_EXTRACTED_REVIEW"
            )
            entry["message"] = (
                "Extracted annotated CDS; verify completeness and target "
                "identity before biological interpretation."
            )

        except Exception as exc:
            entry["status"] = "FETCH_OR_PARSE_ERROR"
            entry["message"] = str(exc)

        results.append(entry)
        time.sleep(0.4)

    with OUTPUT_FASTA.open("w", encoding="utf-8", newline="\n") as handle:
        for record_id, sequence in output_records:
            handle.write(f">{record_id}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start:start + 70] + "\n")

    report = {
        "target_species": "Nilaparvata lugens",
        "requested_gene_count": len(ACCESSIONS),
        "cds_records_written": len(output_records),
        "fasta_path": str(OUTPUT_FASTA),
        "results": results,
        "warning": (
            "CDS extraction is not proof of target suitability or "
            "biological efficacy. Review every status before screening."
        ),
    }

    REPORT.write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )

    print("\nRetrieval summary")
    print(f"Requested genes: {len(ACCESSIONS)}")
    print(f"CDS records written: {len(output_records)}")
    print(f"FASTA: {OUTPUT_FASTA}")
    print(f"Report: {REPORT}")
    for item in results:
        print(
            f"{item['gene']:16} {item['requested_accession']:16} "
            f"{item['status']}"
        )


if __name__ == "__main__":
    main()
