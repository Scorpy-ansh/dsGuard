import json

import urllib.parse



import pytest



import backend.main as main





def _summary_payload(records):

    result = {"uids": [r["uid"] for r in records]}

    for r in records:

        result[r["uid"]] = {

            "accessionversion": r["accession"],

            "title": r["title"],

            "slen": r["length"],

        }

    return json.dumps({"result": result}).encode("utf-8")





REFSEQ_TRANSCRIPT = {

    "uid": "1",

    "accession": "NM_002046.7",

    "title": "Homo sapiens glyceraldehyde-3-phosphate dehydrogenase (GAPDH), mRNA",

    "length": 1285,

}

WRONG_REFSEQ = {

    "uid": "2",

    "accession": "NM_002740.6",

    "title": "Homo sapiens protein kinase C iota (PRKCI), mRNA",

    "length": 4887,

}

GENBANK_EXON = {

    "uid": "3",

    "accession": "MH011443.1",

    "title": "Homo sapiens TP53 (TP53) gene, exon 5 and partial cds",

    "length": 123,

}

WHOLE_GENE = {

    "uid": "4",

    "accession": "NG_017013.2",

    "title": "Homo sapiens tumor protein p53 (TP53), RefSeqGene (LRG_321)",

    "length": 32772,

}





def test_score_prefers_curated_transcript_over_whole_gene():

    transcript = main.ncbi_gene_search_score("TP53", "NM_000546.6", "tumor protein p53 (TP53), transcript variant 1", 2512)

    whole_gene = main.ncbi_gene_search_score("TP53", "NG_017013.2", "tumor protein p53 (TP53), RefSeqGene", 32772)

    chromosome = main.ncbi_gene_search_score("TP53", "NC_000017.11", "Homo sapiens chromosome 17", 83_000_000)

    assert transcript > whole_gene > chromosome





def test_score_penalises_records_that_only_match_free_text():

    matching = main.ncbi_gene_search_score("GAPDH", "NM_002046.7", "glyceraldehyde-3-phosphate dehydrogenase (GAPDH), mRNA", 1285)

    free_text = main.ncbi_gene_search_score("GAPDH", "NM_002740.6", "protein kinase C iota (PRKCI), mRNA", 4887)

    assert matching > free_text





def test_choose_record_rejects_the_wrong_gene():

    chosen = main.choose_ncbi_record("GAPDH", [WRONG_REFSEQ, REFSEQ_TRANSCRIPT])

    assert chosen["accession"] == "NM_002046.7"





def test_choose_record_prefers_transcript_over_exon_fragment():

    chosen = main.choose_ncbi_record("TP53", [GENBANK_EXON, WHOLE_GENE, {

        "uid": "5",

        "accession": "NM_000546.6",

        "title": "Homo sapiens tumor protein p53 (TP53), transcript variant 1, mRNA",

        "length": 2512,

    }])

    assert chosen["accession"] == "NM_000546.6"





def test_choose_record_breaks_ties_with_ncbi_relevance_order():

    first = dict(REFSEQ_TRANSCRIPT, uid="10", accession="NM_001407264.1")

    second = dict(REFSEQ_TRANSCRIPT, uid="11", accession="NM_001407263.1")

    assert main.choose_ncbi_record("GAPDH", [first, second])["uid"] == "10"





def test_choose_record_returns_none_without_candidates():

    assert main.choose_ncbi_record("TP53", []) is None





def test_parse_summary_skips_missing_entries_and_falls_back_to_caption():

    payload = {

        "result": {

            "uids": ["7", "8"],

            "7": {"accessionversion": "NM_1.1", "title": "t", "slen": 10},

            "8": {"caption": "XM_9", "title": "u", "slen": 20},

        }

    }

    records = main.parse_ncbi_summary(payload)

    assert records == [

        {"uid": "7", "accession": "NM_1.1", "title": "t", "length": 10},

        {"uid": "8", "accession": "XM_9", "title": "u", "length": 20},

    ]





def test_empty_query_is_rejected_without_network(monkeypatch):

    monkeypatch.setattr(main, "_ncbi_get", lambda *a, **k: pytest.fail("network call"))

    with pytest.raises(ValueError, match="Enter an NCBI accession or gene name"):

        main.fetch_ncbi_fasta("   ")





def test_no_results_is_distinguished_from_retrieval_failure(monkeypatch):

    monkeypatch.setattr(main, "_ncbi_search", lambda term, retmax: [])

    with pytest.raises(ValueError, match="found no nucleotide records"):

        main.fetch_ncbi_fasta("NOTAGENE12345")





def test_retrieval_failure_reports_a_transport_error(monkeypatch):

    def boom(url, timeout):

        raise OSError("connection reset")



    monkeypatch.setattr(main, "_ncbi_get", boom)

    with pytest.raises(ValueError, match="Could not retrieve sequence from NCBI"):

        main.fetch_ncbi_fasta("NM_000546.6")





def test_gene_search_prefers_refseq_candidates(monkeypatch):

    seen = []



    def fake_get(url, timeout):

        seen.append(url)

        if "esearch" in url:

            term = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)["term"][0]

            if "refseq" in term:

                return json.dumps({"esearchresult": {"idlist": ["1"]}}).encode()

            return json.dumps({"esearchresult": {"idlist": ["3"]}}).encode()

        if "esummary" in url:

            return _summary_payload([REFSEQ_TRANSCRIPT])

        return b">NM_002046.7 Homo sapiens glyceraldehyde-3-phosphate dehydrogenase\nACGTACGT\n"



    monkeypatch.setattr(main, "_ncbi_get", fake_get)

    fasta, record_id = main.fetch_ncbi_fasta("GAPDH", "Homo sapiens")



    assert record_id == "NM_002046.7"

    assert fasta.startswith(">NM_002046.7")

    assert any("refseq" in urllib.parse.unquote(u) for u in seen)





def test_gene_search_falls_back_to_broad_search(monkeypatch):

    def fake_get(url, timeout):

        if "esearch" in url:

            term = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)["term"][0]

            if "refseq" in term:

                return json.dumps({"esearchresult": {"idlist": []}}).encode()

            return json.dumps({"esearchresult": {"idlist": ["3"]}}).encode()

        if "esummary" in url:

            return _summary_payload([GENBANK_EXON])

        return b">MH011443.1 Nilaparvata lugens chitin synthase, partial cds\nACGTACGT\n"



    monkeypatch.setattr(main, "_ncbi_get", fake_get)

    fasta, record_id = main.fetch_ncbi_fasta("chitin synthase", "Nilaparvata lugens")

    assert record_id == "MH011443.1"

    assert fasta.startswith(">MH011443.1")





def test_accession_mode_fetches_directly_without_ranking(monkeypatch):

    def fake_get(url, timeout):

        if "esearch" in url:

            return json.dumps({"esearchresult": {"idlist": ["99"]}}).encode()

        if "esummary" in url:

            pytest.fail("accession lookups must not need a summary call")

        return b">NM_000546.6 Homo sapiens tumor protein p53\nACGTACGT\n"



    monkeypatch.setattr(main, "_ncbi_get", fake_get)

    fasta, record_id = main.fetch_ncbi_fasta("NM_000546.6")

    assert record_id == "NM_000546.6"

    assert fasta.startswith(">NM_000546.6")





def test_non_fasta_payload_is_rejected(monkeypatch):

    def fake_get(url, timeout):

        if "esearch" in url:

            return json.dumps({"esearchresult": {"idlist": ["99"]}}).encode()

        return b"<html>error page</html>"



    monkeypatch.setattr(main, "_ncbi_get", fake_get)

    with pytest.raises(ValueError, match="did not return a valid nucleotide FASTA"):

        main.fetch_ncbi_fasta("NM_000546.6")













def test_project_target_accessions_parses_curated_headers(tmp_path, monkeypatch):

    fasta = tmp_path / "targets.fasta"

    fasta.write_text(

        ">NlRan|accession=KT313028.1|cds_length=645\nACGT\n"

        ">NlSRP54|accession=XM_022328877.2|cds_length=1530\nACGT\n"

        ">alias_without_accession\nACGT\n",

        encoding="utf-8",

    )

    monkeypatch.setattr(main, "PROJECT_TARGET_FASTA", fasta)

    assert main.project_target_accessions() == {

        "nlran": "KT313028.1",

        "nlsrp54": "XM_022328877.2",

    }





def test_project_target_accessions_is_empty_when_list_is_missing(tmp_path, monkeypatch):

    monkeypatch.setattr(main, "PROJECT_TARGET_FASTA", tmp_path / "absent.fasta")

    assert main.project_target_accessions() == {}





def test_project_gene_alias_resolves_to_its_curated_accession(monkeypatch, tmp_path):

    fasta_file = tmp_path / "targets.fasta"

    fasta_file.write_text(

        ">NlRan|accession=KT313028.1|cds_length=645\nACGT\n", encoding="utf-8"

    )

    monkeypatch.setattr(main, "PROJECT_TARGET_FASTA", fasta_file)



    def fake_get(url, timeout):

        if "esearch" in url:

            term = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)["term"][0]

            assert term == "KT313028.1"

            return json.dumps({"esearchresult": {"idlist": ["42"]}}).encode()

        if "esummary" in url:

            pytest.fail("alias lookups must not need a summary call")

        return (

            b">KT313028.1 Nilaparvata lugens RanGTPase (Ran) mRNA, complete cds\n"

            b"ACGTACGT\n"

        )



    monkeypatch.setattr(main, "_ncbi_get", fake_get)

    fasta, record_id = main.fetch_ncbi_fasta("NlRan", "Nilaparvata lugens")

    assert record_id == "KT313028.1"

    assert fasta.startswith(">KT313028.1")





def test_unrelated_record_is_rejected_when_organism_does_not_match(monkeypatch):

    def fake_get(url, timeout):

        if "esearch" in url:

            return json.dumps({"esearchresult": {"idlist": ["7"]}}).encode()

        return b">NM_000546.6 Homo sapiens tumor protein p53\nACGTACGT\n"



    monkeypatch.setattr(main, "_ncbi_get", fake_get)

    with pytest.raises(ValueError, match="not annotated for 'Nilaparvata lugens'"):

        main.fetch_ncbi_fasta("NM_000546.6", "Nilaparvata lugens")





def test_organism_check_is_case_insensitive():

    assert main.organism_is_annotated("XM_1.1 PREDICTED: Nilaparvata lugens Ran", "nilaparvata lugens")

    assert not main.organism_is_annotated("NM_1.1 Homo sapiens TP53", "Nilaparvata lugens")





def test_target_catalog_contains_all_target_organisms_and_no_sequences():

    response = main.list_target_catalog()

    organisms = {item["key"]: item for item in response["organisms"]}



    assert set(organisms) == {

        "nilaparvata_lugens",

        "chilo_suppressalis",

        "sogatella_furcifera",

    }

    assert len(organisms["nilaparvata_lugens"]["genes"]) == 13

    assert len(organisms["chilo_suppressalis"]["genes"]) == 33

    assert len(organisms["sogatella_furcifera"]["genes"]) == 33

    assert all(

        gene["sequence"] is None

        for item in organisms.values()

        for gene in item["genes"]

    )

    assert "chilo_suppressalis" not in main.PANEL_CATALOG

    assert "sogatella_furcifera" not in main.PANEL_CATALOG





def test_catalog_preserves_shared_mitochondrial_coordinates_and_unresolved_accessions():

    response = main.list_target_catalog()

    chilo = next(

        item for item in response["organisms"] if item["key"] == "chilo_suppressalis"

    )

    mitochondrial = [

        gene for gene in chilo["genes"]

        if gene["source_workbook"] == "Chilo_suppressalis_13_genes.xlsx"

    ]

    assert len(mitochondrial) == 13

    assert {gene["accession"] for gene in mitochondrial} == {"NC_015612.1"}

    assert len({

        (gene["gene_name"], gene["start_bp"], gene["end_bp"])

        for gene in mitochondrial

    }) == 13

    unresolved = next(gene for gene in chilo["genes"] if gene["gene_name"] == "CsCHS1")

    assert unresolved["accession"] is None

    assert unresolved["accession_status"] == "unresolved"





def test_unresolved_catalog_gene_is_not_sent_to_ncbi():

    with pytest.raises(ValueError, match="no verified accession"):

        main.fetch_ncbi_fasta("CsCHS1", "Chilo suppressalis")
