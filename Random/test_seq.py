

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





def test_reference_records_carry_their_header_annotation(tmp_path):

    fasta = tmp_path / "annotated.fa"

    fasta.write_text(

        ">NC_008394.4 Oryza sativa Japonica Group chromosome 1\nACGT\n"

        ">NC_008395.2\nTTTT\n"

    )



    reference = Reference.from_fasta(fasta)



    assert reference.names == ["NC_008394.4", "NC_008395.2"]



    assert reference.descriptions[0] == "Oryza sativa Japonica Group chromosome 1"

    assert reference.descriptions[1] == ""

    assert reference.records_without_annotation == 1





def test_fully_annotated_reference_reports_no_missing_annotation(tmp_path):

    fasta = tmp_path / "annotated.fa"

    fasta.write_text(">gene1 Oryza sativa\nACGT\n>gene2 Zea mays\nGGGG\n")



    reference = Reference.from_fasta(fasta)



    assert reference.records_without_annotation == 0
