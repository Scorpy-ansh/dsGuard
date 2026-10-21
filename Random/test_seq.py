
import pytest
from dsguard.seq import read_fasta, encode, Reference


def test_read_valid_fasta(tmp_path):
    fasta = tmp_path / "valid.fa"
    fasta.write_text(">gene1 description\nACGT\nTGCA\n>gene2\nAAAA\n")

    records = list(read_fasta(fasta))

    assert records == [
        ("gene1", b"ACGTTGCA"),
        ("gene2", b"AAAA"),
    ]


def test_empty_fasta_rejected(tmp_path):
    fasta = tmp_path / "empty.fa"
    fasta.write_text("")

    with pytest.raises(ValueError):
        Reference.from_fasta(fasta)


def test_sequence_before_header_rejected(tmp_path):
    fasta = tmp_path / "invalid.fa"
    fasta.write_text("ACGT\n>gene1\nACGT\n")

    with pytest.raises(ValueError):
        list(read_fasta(fasta))


def test_header_without_name_rejected(tmp_path):
    fasta = tmp_path / "invalid.fa"
    fasta.write_text(">\nACGT\n")

    with pytest.raises(ValueError):
        list(read_fasta(fasta))


def test_ambiguous_bases_are_encoded(tmp_path):
    encoded = encode("ACGTNacgt")
    assert encoded.tolist() == [0, 1, 2, 3, 4, 0, 1, 2, 3]
