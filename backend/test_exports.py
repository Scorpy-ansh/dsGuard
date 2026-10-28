

import asyncio

import contextlib

import csv

import io

import json

from unittest.mock import patch



import pytest

from fastapi import HTTPException

from reportlab import rl_config



from backend import main, reporting





@contextlib.contextmanager

def _uncompressed_pdf():

    previous = rl_config.pageCompression

    rl_config.pageCompression = 0



    try:

        yield

    finally:

        rl_config.pageCompression = previous





def _ranked_row(**overrides):

    row = {

        "rank": 1,

        "target_id": "testgene",

        "target_window": "1-30",

        "matched_organism": "rice",

        "transcript_or_reference_id": "rice.1",

        "reference_panel": "rice",

        "transcript_pos_1based": 5,

        "target_start": 1,

        "target_end": 30,

        "match_length_bp": 30,

        "identity_pct": 100.0,

        "mismatches": 0,

        "gaps": 0,

        "strand": "+",

        "potency": 0.5,

        "hazard": 0.989013,

        "risk_category": "Critical",

        "threshold_status": "exceeds_threshold",

        "hazard_is_exact": True,

        "evidence": "confirmed exact match",

        "in_recommended_window": 1,

        "governs_recommendation": True,

    }

    row.update(overrides)

    return row





@pytest.fixture

def sample_results():

    return {

        "project": "dsGuard",

        "analysis_id": "analysis-123",

        "created_utc": "2026-01-01T00:00:00+00:00",

        "input": {

            "label": "testgene",

            "type": "fasta_file",

            "sha256": "deadbeef",

            "genes": 1,

            "ambiguous_bases": 0,

            "records": [{"name": "testgene", "length_bp": 40}],

        },

        "parameters": {

            "window": 30,

            "k": 11,

            "threshold": 21,

            "l50": 21.0,

            "kappa": 2.0,

            "max_hazard": 0.05,

            "min_report_match": 21,

            "hazard_threshold_equivalent_match_bp": 15.11,

            "l50_kappa_fitted": False,

        },

        "reporting": {

            "min_report_match_bp": 21,

            "guaranteed_detection_threshold_bp": 21,

            "note": "Only perfect gapless matches at or above the floor are listed.",

        },

        "scoring": {

            "hazard": {

                "formula": "hazard(m) = 1 / (1 + exp(-(m - L50) / kappa))",

                "input": "longest perfect match",

                "range": "[0, 1]",

                "meaning": "heuristic",

                "calibrated": False,

                "caveat": "placeholder constants",

            },

            "potency": {"formula": "p", "meaning": "m", "caveat": "c"},

            "matching_method": {"description": "d", "completeness": "c"},

        },

        "reference_panels": {

            "searched": [

                {"name": "rice", "path": "p", "sha256": "s", "Mb": 1.0,

                 "records": 1, "seconds": 0.1}

            ],

            "skipped": [

                {"name": "honeybee", "status": "not_searched",

                 "detail": "FASTA not installed"}

            ],

        },

        "summary": [

            {

                "gene": "testgene",

                "status": "ok",

                "recommendation_status": "fallback_threshold_exceeded",

                "rec_start": 11,

                "rec_end": 40,

                "rec_potency": 0.5,

                "rec_max_match": 20,

                "rec_hazard": 0.377541,

                "total_hits": 1,

                "max_hazard": 0.05,

            }

        ],

        "genes": [

            {

                "gene": "testgene",

                "row": {"gene": "testgene"},

                "length_bp": 40,

                "window": 30,

                "total_matches_found": 2,

                "matches_below_report_floor": 0,

                "pareto": [

                    {"start": 11, "end": 40, "potency": 0.5, "max_match": 20,

                     "hazard": 0.377541, "threshold_status": "exceeds_threshold",

                     "hazard_is_exact": True}

                ],

                "recommendation": {

                    "start": 11, "end": 40, "hazard": 0.377541,

                    "recommendation_status": "fallback_threshold_exceeded",

                    "is_fallback": True,

                },

                "ranked_offtargets": [

                    _ranked_row(),

                    _ranked_row(

                        rank=2,

                        target_window="11-40",

                        transcript_pos_1based=25,

                        target_start=11,

                        target_end=40,

                        match_length_bp=20,

                        hazard=0.377541,

                        risk_category="High",

                        in_recommended_window=1,

                        governs_recommendation=False,

                    ),

                ],





                "candidate": {

                    "candidate_id": "testgene:12-41",

                    "gene": "testgene",

                    "requested_length_bp": 30,

                    "actual_length_bp": 30,

                    "sequence": "ACGTACGTACACGTACGTACACGTACGTAC",

                    "start_1based": 12,

                    "end_1based": 41,

                    "strand": "sense",

                    "ambiguous_bases": 0,

                },

                "candidate_error": None,

            },

            {





                "gene": "shortgene",

                "row": {"gene": "shortgene"},

                "length_bp": 15,

                "window": 30,

                "total_matches_found": 0,

                "matches_below_report_floor": 0,

                "pareto": [],

                "recommendation": {

                    "start": None, "end": None, "hazard": None,

                    "recommendation_status": "too_short",

                    "is_fallback": False,

                },

                "ranked_offtargets": [],

                "candidate": None,

                "candidate_error": (

                    "target is 15 bp, shorter than the requested 30 bp window"

                ),

            },

        ],

        "preventive_measures": {

            "findings": [

                {"finding": "f", "why": "w", "action": "a", "priority": "High",

                 "evidence": "rec_hazard=0.377541 > max_hazard=0.05"}

            ],

            "categories": {"before_accepting_a_candidate": ["check it"]},

        },

        "limitations": ["heuristic only"],

    }













def test_offtarget_csv_is_ranked_by_hazard_descending(sample_results):

    text = reporting.build_offtarget_csv(sample_results)







    parsed = [

        row for row in csv.reader(io.StringIO(text))

        if row and not row[0].startswith("#") and row[0] != "gene"

    ]



    assert len(parsed) == 2



    hazard_column = reporting.RANKED_HEADERS.index("Risk score (hazard)") + 1

    assert float(parsed[0][hazard_column]) >= float(parsed[1][hazard_column])

    assert parsed[0][1] == "1"

    assert parsed[1][1] == "2"



    assert "configured_max_hazard,0.05" in text



    assert parsed[0][0] == "testgene"

    assert parsed[0][reporting.RANKED_HEADERS.index("Threshold status") + 1] == "exceeds_threshold"





def test_results_json_includes_provenance_and_scoring(sample_results):

    payload = json.loads(reporting.build_results_json(sample_results))



    assert payload["analysis_id"] == "analysis-123"

    assert payload["parameters"]["max_hazard"] == 0.05

    assert payload["scoring"]["hazard"]["calibrated"] is False

    assert payload["reference_panels"]["skipped"][0]["name"] == "honeybee"

    assert payload["export"]["format"] == "dsguard-results-v1"





def test_pdf_report_is_a_valid_pdf_document(sample_results):

    pdf = reporting.build_pdf_report(sample_results)



    assert isinstance(pdf, bytes)

    assert pdf.startswith(b"%PDF")

    assert len(pdf) > 1000













def _collect(response):

    async def gather():

        chunks = []

        async for chunk in response.body_iterator:

            chunks.append(chunk if isinstance(chunk, bytes) else chunk.encode("utf-8"))

        return b"".join(chunks)



    return asyncio.run(gather())





@pytest.fixture

def job_setup(tmp_path, sample_results):

    output_root = tmp_path / "outputs"

    job_dir = output_root / "job-1"

    job_dir.mkdir(parents=True)

    (job_dir / "results.json").write_text(

        json.dumps(sample_results), encoding="utf-8"

    )



    with (

        patch.object(main, "OUTPUT_ROOT", output_root),

        patch.object(main, "jobs") as jobs_mock,

    ):

        yield jobs_mock





def test_results_endpoint_includes_current_analysis_object(job_setup):

    job_setup.find_one.return_value = {"status": "completed"}



    payload = main.get_job_results("job-1")



    assert payload["status"] == "completed"

    assert payload["analysis"]["analysis_id"] == "analysis-123"

    assert payload["analysis"]["genes"][0]["recommendation"]["is_fallback"] is True





def test_results_endpoint_analysis_is_none_without_results_file(job_setup):

    job_setup.find_one.return_value = {"status": "completed"}



    payload = main.get_job_results("job-without-results")



    assert payload["analysis"] is None











CANDIDATE_FIELDS = (

    "candidate_id",

    "sequence",

    "requested_length_bp",

    "actual_length_bp",

    "start_1based",

    "end_1based",

    "strand",

    "ambiguous_bases",

)





def test_results_endpoint_returns_the_candidate_object_unchanged(job_setup, sample_results):

    job_setup.find_one.return_value = {"status": "completed"}



    payload = main.get_job_results("job-1")



    recorded = sample_results["genes"][0]["candidate"]

    returned = payload["analysis"]["genes"][0]["candidate"]







    assert returned == recorded

    for field in CANDIDATE_FIELDS:

        assert field in returned, f"the API dropped {field}"





def test_results_endpoint_candidate_lengths_are_self_consistent(job_setup):

    job_setup.find_one.return_value = {"status": "completed"}



    candidate = main.get_job_results("job-1")["analysis"]["genes"][0]["candidate"]







    assert candidate["actual_length_bp"] == len(candidate["sequence"])

    assert candidate["actual_length_bp"] == candidate["requested_length_bp"]





def test_results_endpoint_preserves_null_candidate_and_its_error(job_setup):

    job_setup.find_one.return_value = {"status": "completed"}



    gene = main.get_job_results("job-1")["analysis"]["genes"][1]



    assert gene["candidate"] is None

    assert "shorter than the requested" in gene["candidate_error"]





def test_results_json_export_keeps_the_candidate_object(job_setup, sample_results):

    job_setup.find_one.return_value = {"status": "completed"}



    payload = json.loads(_collect(main.export_results_json("job-1")).decode("utf-8"))



    assert payload["genes"][0]["candidate"] == sample_results["genes"][0]["candidate"]

    assert payload["genes"][1]["candidate"] is None













def _fasta_records(text):

    records = []

    header = None

    chunks = []



    for line in text.splitlines():

        if line.startswith(">"):

            if header is not None:

                records.append((header, "".join(chunks)))

            header, chunks = line, []

        elif header is not None:

            chunks.append(line)



    if header is not None:

        records.append((header, "".join(chunks)))



    return records





def test_candidate_fasta_has_one_sanitised_header_with_the_exported_length(sample_results):

    text = reporting.build_candidate_fasta(sample_results)

    records = _fasta_records(text)







    assert len(records) == 1

    assert text.count(">") == 1



    header, sequence = records[0]

    candidate = sample_results["genes"][0]["candidate"]



    assert header.startswith(">testgene|candidate_id=testgene:12-41|")

    assert f"|length_bp={candidate['actual_length_bp']}|" in header

    assert "|target_position=12-41|" in header

    assert "|strand=sense" in header

    assert "\n" not in header





def test_candidate_fasta_sequence_matches_the_candidate_object_exactly(sample_results):

    _, sequence = _fasta_records(reporting.build_candidate_fasta(sample_results))[0]

    candidate = sample_results["genes"][0]["candidate"]







    assert sequence == candidate["sequence"]

    assert len(sequence) == candidate["actual_length_bp"]

    assert len(sequence) != sample_results["genes"][0]["length_bp"]





def test_candidate_fasta_wraps_long_sequences_to_standard_lines(sample_results):

    sequence = "ACGT" * 33

    candidate = sample_results["genes"][0]["candidate"]

    candidate["sequence"] = sequence

    candidate["requested_length_bp"] = len(sequence)

    candidate["actual_length_bp"] = len(sequence)



    text = reporting.build_candidate_fasta(sample_results)

    header, exported = _fasta_records(text)[0]



    assert [len(line) for line in text.splitlines()[1:]] == [60, 60, 12]

    assert exported == sequence

    assert f"|length_bp={len(sequence)}|" in header





def test_candidate_fasta_header_metadata_cannot_inject_an_extra_record(sample_results):

    gene = sample_results["genes"][0]

    gene["gene"] = "evil\n>fake_record|length_bp=4"

    gene["candidate"]["candidate_id"] = "id\n>ACGT\nGTAC"

    gene["candidate"]["strand"] = "sense>"



    text = reporting.build_candidate_fasta(sample_results)

    records = _fasta_records(text)



    assert len(records) == 1, "header metadata must not create a second record"

    header, sequence = records[0]

    assert header.count(">") == 1

    assert "\n" not in header

    assert sequence == "ACGTACGTACACGTACGTACACGTACGTAC"





def test_candidate_fasta_refuses_a_non_nucleotide_sequence(sample_results):

    sample_results["genes"][0]["candidate"]["sequence"] = "ACGT>fake_record"





    with pytest.raises(ValueError):

        reporting.build_candidate_fasta(sample_results)





def test_candidate_fasta_includes_organism_only_when_the_run_recorded_one(sample_results):





    assert "organism=" not in reporting.build_candidate_fasta(sample_results)



    sample_results["input"]["organism"] = "Nilaparvata lugens"

    header = _fasta_records(reporting.build_candidate_fasta(sample_results))[0][0]



    assert "|organism=Nilaparvata lugens|" in header





def test_candidate_fasta_is_empty_when_no_gene_has_a_candidate(sample_results):

    for gene in sample_results["genes"]:

        gene["candidate"] = None



    assert reporting.build_candidate_fasta(sample_results) == ""





def test_candidate_fasta_handles_results_that_predate_candidate_output(sample_results):

    for gene in sample_results["genes"]:

        gene.pop("candidate", None)

        gene.pop("candidate_error", None)



    assert reporting.build_candidate_fasta(sample_results) == ""

    assert [(candidate, sequence) for _, candidate, sequence

            in reporting.candidate_gene_rows(sample_results)] == [

        (None, ""),

        (None, ""),

    ]





def test_export_fasta_endpoint_returns_the_recorded_candidate(job_setup, sample_results):

    job_setup.find_one.return_value = {"status": "completed"}



    response = main.export_candidate_fasta("job-1")

    body = _collect(response).decode("utf-8")



    assert response.media_type == "text/x-fasta"

    assert "candidates.fasta" in response.headers["content-disposition"]

    assert body.startswith(">")



    header, sequence = _fasta_records(body)[0]

    assert sequence == sample_results["genes"][0]["candidate"]["sequence"]

    assert "testgene" in header





def test_export_fasta_endpoint_404_when_no_candidate_was_recorded(tmp_path, sample_results):

    output_root = tmp_path / "outputs"

    job_dir = output_root / "job-no-candidate"

    job_dir.mkdir(parents=True)

    for gene in sample_results["genes"]:

        gene["candidate"] = None

    (job_dir / "results.json").write_text(

        json.dumps(sample_results), encoding="utf-8"

    )



    with (

        patch.object(main, "OUTPUT_ROOT", output_root),

        patch.object(main, "jobs") as jobs_mock,

    ):

        jobs_mock.find_one.return_value = {"status": "completed"}



        with pytest.raises(HTTPException) as exc:

            main.export_candidate_fasta("job-no-candidate")





    assert exc.value.status_code == 404

    assert "candidate" in exc.value.detail

    assert "candidate_error" in exc.value.detail





def test_export_fasta_endpoint_500_for_a_corrupt_stored_sequence(tmp_path, sample_results):

    output_root = tmp_path / "outputs"

    job_dir = output_root / "job-corrupt"

    job_dir.mkdir(parents=True)

    sample_results["genes"][0]["candidate"]["sequence"] = "ACGT<bad>"

    (job_dir / "results.json").write_text(

        json.dumps(sample_results), encoding="utf-8"

    )



    with (

        patch.object(main, "OUTPUT_ROOT", output_root),

        patch.object(main, "jobs") as jobs_mock,

    ):

        jobs_mock.find_one.return_value = {"status": "completed"}



        with pytest.raises(HTTPException) as exc:

            main.export_candidate_fasta("job-corrupt")



    assert exc.value.status_code == 500

    assert "FASTA" in exc.value.detail



















def test_pdf_report_renders_the_candidate_sequence_with_its_metadata(sample_results):

    with _uncompressed_pdf():

        pdf = reporting.build_pdf_report(sample_results)



    candidate = sample_results["genes"][0]["candidate"]



    assert pdf.startswith(b"%PDF")



    assert b"dsRNA candidate sequences" in pdf

    assert candidate["sequence"].encode("ascii") in pdf

    assert b"testgene:12-41" in pdf

    assert b"Actual length" in pdf



    assert b"not the 40 bp target gene" in pdf



    assert str(sample_results["genes"][0]["length_bp"]).encode() in pdf





def test_pdf_report_states_why_a_gene_has_no_candidate_sequence(sample_results):

    sample_results["genes"][1]["candidate_error"] = "TooShortForWindow"



    with _uncompressed_pdf():

        pdf = reporting.build_pdf_report(sample_results)



    assert pdf.startswith(b"%PDF")

    assert b"shortgene" in pdf

    assert b"no candidate sequence was generated" in pdf



    assert b"TooShortForWindow" in pdf





def test_pdf_report_handles_a_run_that_predates_candidate_output(sample_results):

    for gene in sample_results["genes"]:

        gene.pop("candidate", None)

        gene.pop("candidate_error", None)



    with _uncompressed_pdf():

        pdf = reporting.build_pdf_report(sample_results)



    assert pdf.startswith(b"%PDF")



    assert b"no candidate sequence was generated" in pdf

    assert b"shortgene" in pdf

    assert b"predates candidate" in pdf





def test_export_csv_endpoint_returns_ranked_rows(job_setup):

    job_setup.find_one.return_value = {"status": "completed"}



    response = main.export_offtarget_csv("job-1")

    body = _collect(response).decode("utf-8")



    assert response.media_type == "text/csv"

    assert "offtargets.csv" in response.headers["content-disposition"]

    assert "rice.1" in body





def test_export_json_endpoint_returns_reproducible_results(job_setup):

    job_setup.find_one.return_value = {"status": "completed"}



    response = main.export_results_json("job-1")

    payload = json.loads(_collect(response).decode("utf-8"))



    assert response.media_type == "application/json"

    assert payload["analysis_id"] == "analysis-123"





def test_export_pdf_endpoint_returns_pdf(job_setup):

    job_setup.find_one.return_value = {"status": "completed"}



    response = main.export_pdf_report("job-1")

    body = _collect(response)



    assert response.media_type == "application/pdf"

    assert body.startswith(b"%PDF")





def test_export_endpoints_reject_incomplete_job(job_setup):

    job_setup.find_one.return_value = {"status": "running"}



    for endpoint in (

        main.export_candidate_fasta,

        main.export_offtarget_csv,

        main.export_results_json,

        main.export_pdf_report,

    ):

        with pytest.raises(HTTPException) as exc:

            endpoint("job-1")



        assert exc.value.status_code == 409





def test_export_endpoints_404_for_unknown_job(job_setup):

    job_setup.find_one.return_value = None



    with pytest.raises(HTTPException) as exc:

        main.export_results_json("missing")



    assert exc.value.status_code == 404





def test_export_endpoint_404_when_results_json_missing(job_setup):

    job_setup.find_one.return_value = {"status": "completed"}



    with pytest.raises(HTTPException) as exc:

        main.export_pdf_report("job-without-results")



    assert exc.value.status_code == 404

    assert "results.json" in exc.value.detail





def test_load_job_results_object_rejects_path_traversal(tmp_path):

    with patch.object(main, "OUTPUT_ROOT", tmp_path / "outputs"):

        (tmp_path / "outputs").mkdir()

        (tmp_path / "outside").mkdir()



        with pytest.raises(HTTPException) as exc:

            main.load_job_results_object("../outside")



    assert exc.value.status_code == 400













REQUIRED_RANKED_FIELDS = (

    "Rank",

    "Target ID",

    "Target window",

    "Matched organism",

    "Transcript / reference ID",

    "Match (bp)",

    "Identity (%)",

    "Mismatch / gap",

    "Risk score (hazard)",

    "Risk category",

    "Threshold status",

)





def test_ranked_csv_exposes_the_required_fields(sample_results):

    text = reporting.build_offtarget_csv(sample_results)

    header = next(

        row for row in csv.reader(io.StringIO(text)) if row and row[0] == "gene"

    )



    for field in REQUIRED_RANKED_FIELDS:

        assert field in header, f"{field} missing from the ranked CSV header"





def test_ranked_headers_use_reference_record_wording(sample_results):



    assert "Transcript / reference ID" in reporting.RANKED_HEADERS

    assert not any(

        header.lower() == "reference transcript" for header in reporting.RANKED_HEADERS

    )





def test_ranked_csv_rank_column_is_a_job_wide_position(sample_results):



    second_gene = json.loads(json.dumps(sample_results["genes"][0]))

    second_gene["gene"] = "othergene"

    for row in second_gene["ranked_offtargets"]:

        row["target_id"] = "othergene"

        row["hazard"] = 0.9139

    sample_results["genes"].append(second_gene)



    text = reporting.build_offtarget_csv(sample_results)

    parsed = [

        row for row in csv.reader(io.StringIO(text))

        if row and not row[0].startswith("#") and row[0] != "gene"

    ]



    assert [row[1] for row in parsed] == ["1", "2", "3", "4"]

    hazards = [float(row[reporting.RANKED_HEADERS.index("Risk score (hazard)") + 1])

               for row in parsed]

    assert hazards == sorted(hazards, reverse=True)





def test_flatten_is_deterministic_for_equal_hazard():

    def results(rows):

        return {"genes": [{"gene": "G1", "ranked_offtargets": rows}]}



    a = _ranked_row(rank=None, hazard=0.5, transcript_or_reference_id="b.1")

    b = _ranked_row(rank=None, hazard=0.5, transcript_or_reference_id="a.1")



    first = reporting.flatten_ranked_offtargets(results([a, b]))

    second = reporting.flatten_ranked_offtargets(results([b, a]))



    expected = ["a.1", "b.1"]

    assert [r["transcript_or_reference_id"] for r in first] == expected

    assert [r["transcript_or_reference_id"] for r in second] == expected

    assert [r["rank"] for r in first] == [1, 2]













def test_threshold_agrees_across_csv_json_and_summary(sample_results):

    csv_text = reporting.build_offtarget_csv(sample_results)

    payload = json.loads(reporting.build_results_json(sample_results))



    assert "configured_max_hazard,0.05" in csv_text

    assert payload["parameters"]["max_hazard"] == 0.05

    assert payload["summary"][0]["max_hazard"] == 0.05





def test_json_export_preserves_findings_evidence_and_priority(sample_results):

    payload = json.loads(reporting.build_results_json(sample_results))



    finding = payload["preventive_measures"]["findings"][0]

    assert finding["priority"] == "High"

    assert finding["evidence"] == "rec_hazard=0.377541 > max_hazard=0.05"





def test_pdf_renders_findings_with_and_without_evidence(sample_results):

    sample_results["preventive_measures"]["findings"].append(

        {"finding": "no evidence yet", "why": "w", "action": "a"}

    )



    pdf = reporting.build_pdf_report(sample_results)



    assert pdf.startswith(b"%PDF")













def test_missing_results_detail_is_actionable(job_setup):

    job_setup.find_one.return_value = {"status": "completed"}



    with pytest.raises(HTTPException) as exc:

        main.require_completed_results("job-without-results")



    detail = exc.value.detail

    assert exc.value.status_code == 404

    assert "results.json" in detail

    assert "cannot be fabricated" in detail

    assert "Re-run the screening" in detail













def test_ranked_csv_header_is_unique_and_metadata_is_marked(sample_results):

    text = reporting.build_offtarget_csv(sample_results)

    rows = list(csv.reader(io.StringIO(text)))



    header_index = next(

        index for index, row in enumerate(rows) if row and row[0] == "gene"

    )

    header = rows[header_index]



    assert len(header) == len(set(header)), "ranked CSV header must be unique"





    for row in rows[:header_index]:

        assert row and row[0].startswith("#")





    data = [row for row in rows[header_index + 1:] if row]

    assert data

    for row in data:

        assert len(row) == len(header)





@pytest.fixture

def hits_setup(tmp_path):

    output_root = tmp_path / "outputs"

    gene_dir = output_root / "job-hits" / "geneA_accession_X"

    gene_dir.mkdir(parents=True)

    (gene_dir / "hits.csv").write_text(

        "organism,transcript,transcript_pos_1based,gene_start,gene_end,"

        "length,strand\r\n"

        "rice,NC_008394.4,25463079,571,593,23,+\r\n"

        "rice,XM_001.2,100,1,21,21,-\r\n",

        encoding="utf-8",

    )



    with (

        patch.object(main, "OUTPUT_ROOT", output_root),

        patch.object(main, "jobs") as jobs_mock,

    ):

        jobs_mock.find_one.return_value = {"status": "completed"}

        yield jobs_mock





def test_detailed_hits_csv_has_no_preamble_and_is_well_formed(hits_setup):

    response = main.download_job_hits("job-hits")

    body = _collect(response).decode("utf-8")



    assert response.media_type == "text/csv"



    rows = [row for row in csv.reader(io.StringIO(body)) if row]

    header = rows[0]





    assert not header[0].startswith("#"), "detailed-hits CSV must not start with prose"

    assert header == [

        "gene",

        "organism",

        "reference_record",

        "reference_pos_1based",

        "gene_start",

        "gene_end",

        "length",

        "strand",

    ]

    assert len(header) == len(set(header)), "detailed-hits CSV header must be unique"



    data = rows[1:]

    assert len(data) == 2

    for row in data:

        assert len(row) == len(header), "every hits row must match the header width"



    assert data[0][:3] == ["geneA", "rice", "NC_008394.4"]
