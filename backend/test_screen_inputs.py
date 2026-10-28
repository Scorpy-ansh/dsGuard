

import asyncio

import inspect

import json

from io import BytesIO

from unittest.mock import patch



import pytest

from fastapi import HTTPException

from starlette.background import BackgroundTasks

from starlette.datastructures import UploadFile



from backend import main



VALID_FASTA = b">target_gene test\nACGTACGTACGTACGTACGTACGT\n"

RICE_FASTA = ">NC_TEST.1 Oryza sativa Japonica Group test chromosome\n" + "ACGT" * 50 + "\n"





@pytest.fixture

def screen_env(tmp_path):

    panel_dir = tmp_path / "reference_panels"

    panel_dir.mkdir()

    (panel_dir / "rice.fasta").write_text(RICE_FASTA, encoding="utf-8")



    upload_root = tmp_path / "uploads"

    output_root = tmp_path / "outputs"

    upload_root.mkdir()

    output_root.mkdir()



    with (

        patch.object(main, "PANEL_DIR", panel_dir),

        patch.object(main, "UPLOAD_ROOT", upload_root),

        patch.object(main, "OUTPUT_ROOT", output_root),

        patch.object(main, "jobs") as jobs_mock,

        patch.object(main, "engine_command", return_value=["python", "-m", "dsguard"]),

        patch.object(main, "run_screening") as worker,

    ):

        yield jobs_mock, worker, upload_root, panel_dir





def call_screen(background=None, **kwargs):

    background = background if background is not None else BackgroundTasks()

    result = asyncio.run(main.create_screening_job(background, **kwargs))

    return result, background





def make_upload(name: str, content: bytes) -> UploadFile:

    return UploadFile(filename=name, file=BytesIO(content))





def scheduled_args(background, worker):

    assert len(background.tasks) == 1

    task = background.tasks[0]

    assert task.func is worker

    return task.args













def test_uploaded_fasta_is_accepted_and_scheduled(screen_env):

    jobs_mock, worker, upload_root, _ = screen_env



    result, background = call_screen(

        panel_keys="rice",

        target_mode="file",

        target_file=make_upload("target.fa", VALID_FASTA),

    )



    assert result["status"] == "queued"

    assert result["panels"] == [main.PANEL_CATALOG["rice"]["label"]]

    jobs_mock.insert_one.assert_called_once()

    assert len(background.tasks) == 1



    job_id, target_path, panel_paths, params, input_type, label, skipped = \
        scheduled_args(background, worker)

    assert job_id == result["job_id"]

    assert input_type == "fasta_file"

    assert label == "target.fa"

    assert skipped == []

    assert panel_paths[0][0] == "rice"

    assert target_path.is_file()

    assert target_path.read_bytes().startswith(b">")

    assert (upload_root / job_id).is_dir()





def test_file_mode_without_a_file_is_rejected(screen_env):

    with pytest.raises(HTTPException) as exc:

        call_screen(panel_keys="rice", target_mode="file")



    assert exc.value.status_code == 400

    assert "Choose a target FASTA file" in exc.value.detail





def test_upload_with_wrong_extension_is_rejected(screen_env):

    with pytest.raises(HTTPException) as exc:

        call_screen(

            panel_keys="rice",

            target_mode="file",

            target_file=make_upload("target.csv", VALID_FASTA),

        )



    assert exc.value.status_code == 400

    assert "FASTA extension" in exc.value.detail





def test_upload_that_is_not_fasta_is_rejected(screen_env):

    with pytest.raises(HTTPException) as exc:

        call_screen(

            panel_keys="rice",

            target_mode="file",

            target_file=make_upload("target.fa", b"ACGT\nACGT\n"),

        )



    assert exc.value.status_code == 400

    assert "not recognized as FASTA" in exc.value.detail













def test_pasted_sequence_is_accepted(screen_env):

    _, worker, _, _ = screen_env



    result, background = call_screen(

        panel_keys="rice",

        target_mode="sequence",

        target_sequence="ACGTACGTACGTACGTACGTACGT",

    )



    assert result["status"] == "queued"

    _, target_path, _, _, input_type, label, _ = scheduled_args(background, worker)

    assert input_type == "pasted_sequence"

    assert label == "pasted_target"



    written = target_path.read_text(encoding="utf-8")

    assert written.startswith(">pasted_target")

    assert "ACGTACGT" in written





def test_pasted_fasta_keeps_its_own_header(screen_env):

    _, worker, _, _ = screen_env



    pasted = ">NlRan|accession=KT313028.1|cds_length=645\nACGTACGTACGTACGTACGTACGT\n"



    result, background = call_screen(

        panel_keys="rice",

        target_mode="sequence",

        target_sequence=pasted,

    )



    assert result["status"] == "queued"

    _, target_path, _, _, input_type, label, _ = scheduled_args(background, worker)

    assert input_type == "pasted_sequence"

    assert label == "NlRan|accession=KT313028.1|cds_length=645"

    written = target_path.read_text(encoding="utf-8")

    assert written.startswith(">NlRan|accession=KT313028.1|cds_length=645")

    assert "ACGTACGT" in written





def test_pasted_sequence_with_invalid_characters_is_rejected(screen_env):

    with pytest.raises(HTTPException) as exc:

        call_screen(panel_keys="rice", target_mode="sequence", target_sequence="ACGTZZZ")



    assert exc.value.status_code == 400

    assert "invalid characters" in exc.value.detail





def test_empty_pasted_sequence_is_rejected(screen_env):

    with pytest.raises(HTTPException) as exc:

        call_screen(panel_keys="rice", target_mode="sequence", target_sequence="   ")



    assert exc.value.status_code == 400

    assert "Paste a target" in exc.value.detail













def test_ncbi_accid_is_fetched_and_recorded(screen_env):

    _, worker, _, _ = screen_env



    with patch.object(

        main, "fetch_ncbi_fasta",

        return_value=(">NM_000546.6 Homo sapiens tumor protein p53\nACGTACGT\n", "12345"),

    ):

        result, background = call_screen(

            panel_keys="rice",

            target_mode="ncbi",

            target_query="NM_000546.6",

        )



    assert result["status"] == "queued"

    _, target_path, _, _, input_type, label, _ = scheduled_args(background, worker)

    assert input_type == "ncbi"

    assert label.startswith("NCBI:12345")

    assert target_path.read_text(encoding="utf-8").startswith(">NM_000546.6")





def test_ncbi_gene_name_uses_the_same_mode(screen_env):

    _, worker, _, _ = screen_env



    with patch.object(

        main, "fetch_ncbi_fasta",

        return_value=(">XM_022328877.2 Nilaparvata lugens\nACGTACGT\n", "999"),

    ):

        result, background = call_screen(

            panel_keys="rice",

            target_mode="ncbi",

            target_query="NlRan",

            target_organism="Nilaparvata lugens",

        )



    assert result["status"] == "queued"

    _, _, _, _, input_type, label, _ = scheduled_args(background, worker)

    assert input_type == "ncbi"

    assert "NlRan" in label





def test_ncbi_network_failure_is_surfaced_as_400(screen_env):

    with patch.object(

        main, "fetch_ncbi_fasta",

        side_effect=ValueError("Could not retrieve sequence from NCBI: timed out"),

    ):

        with pytest.raises(HTTPException) as exc:

            call_screen(panel_keys="rice", target_mode="ncbi", target_query="NM_000546.6")



    assert exc.value.status_code == 400

    assert "Could not retrieve sequence from NCBI" in exc.value.detail





def test_ncbi_empty_query_is_rejected_without_a_job(screen_env):

    jobs_mock, worker, _, _ = screen_env

    background = BackgroundTasks()



    with patch.object(

        main, "fetch_ncbi_fasta",

        side_effect=ValueError("Enter an NCBI accession or gene name."),

    ):

        with pytest.raises(HTTPException) as exc:

            call_screen(

                background,

                panel_keys="rice",

                target_mode="ncbi",

                target_query="",

            )



    assert exc.value.status_code == 400

    jobs_mock.insert_one.assert_not_called()

    assert background.tasks == []

    worker.assert_not_called()













def test_missing_panel_is_reported_not_silently_skipped(screen_env):

    _, worker, _, _ = screen_env



    result, background = call_screen(

        panel_keys="rice,honeybee",

        target_mode="sequence",

        target_sequence="ACGTACGTACGTACGTACGTACGT",

    )





    assert [p["name"] for p in result["skipped_panels"]] == ["honeybee"]

    assert "not found" in result["skipped_panels"][0]["reason"]

    *_, skipped = scheduled_args(background, worker)

    assert [name for name, _ in skipped] == ["honeybee"]





def test_no_installed_panel_rejects_the_job(screen_env):

    with pytest.raises(HTTPException) as exc:

        call_screen(

            panel_keys="honeybee",

            target_mode="sequence",

            target_sequence="ACGTACGTACGTACGTACGTACGT",

        )



    assert exc.value.status_code == 400

    assert "None of the requested reference panels are installed" in exc.value.detail





def test_no_panel_selected_is_rejected(screen_env):

    with pytest.raises(HTTPException) as exc:

        call_screen(target_mode="sequence", target_sequence="ACGTACGTACGT")



    assert exc.value.status_code == 400

    assert "Choose at least one" in exc.value.detail





def test_unknown_panel_is_rejected(screen_env):

    with pytest.raises(HTTPException) as exc:

        call_screen(

            panel_keys="dragon",

            target_mode="sequence",

            target_sequence="ACGTACGTACGT",

        )



    assert exc.value.status_code == 400

    assert "Unsupported non-target panel" in exc.value.detail













@pytest.mark.parametrize(

    "params, fragment",

    [

        ({"window": 0}, "Window must be a whole number between"),

        ({"window": 20}, "Window must be a whole number between"),

        ({"window": -3}, "Window must be a whole number between"),

        ({"window": 502}, "Window must be a whole number between"),

        ({"k": 1}, "k must be between"),

        ({"l50": 0}, "l50 and kappa must be positive"),

        ({"kappa": -1}, "l50 and kappa must be positive"),

        ({"max_hazard": 2}, "max_hazard must be between"),

    ],

)

def test_invalid_parameters_are_rejected(screen_env, params, fragment):

    with pytest.raises(HTTPException) as exc:

        call_screen(

            panel_keys="rice",

            target_mode="sequence",

            target_sequence="ACGTACGTACGT",

            **params,

        )



    assert exc.value.status_code == 400

    assert fragment in exc.value.detail













WINDOW_MIN = main.WINDOW_MIN

WINDOW_MAX = main.WINDOW_MAX









SUPPORTED_WINDOWS = [

    WINDOW_MIN,

    50, 100, 101, 200, 201, 300, 301, 400, 401, WINDOW_MAX,

]





@pytest.mark.parametrize("window", SUPPORTED_WINDOWS)

def test_supported_window_lengths_are_accepted_and_passed_through(screen_env, window):

    _, worker, _, _ = screen_env



    result, background = call_screen(

        panel_keys="rice",

        target_mode="sequence",

        target_sequence="ACGTACGTACGTACGTACGTACGT",

        window=window,

    )



    assert result["status"] == "queued"

    _, _, _, params, *_ = scheduled_args(background, worker)

    assert params["window"] == window





@pytest.mark.parametrize("window", [0, 1, WINDOW_MIN - 1, WINDOW_MAX + 1, 5000])

def test_out_of_range_window_is_rejected_not_replaced(screen_env, window):

    jobs_mock, worker, _, _ = screen_env

    background = BackgroundTasks()



    with pytest.raises(HTTPException) as exc:

        call_screen(

            background,

            panel_keys="rice",

            target_mode="sequence",

            target_sequence="ACGTACGTACGTACGTACGTACGT",

            window=window,

        )



    assert exc.value.status_code == 400

    assert f"between {WINDOW_MIN} and {WINDOW_MAX} bp" in exc.value.detail



    jobs_mock.insert_one.assert_not_called()

    worker.assert_not_called()

    assert background.tasks == []





def test_window_default_stays_backward_compatible(screen_env):

    _, worker, _, _ = screen_env



    result, background = call_screen(

        panel_keys="rice",

        target_mode="sequence",

        target_sequence="ACGTACGTACGTACGTACGTACGT",

    )



    assert result["status"] == "queued"

    _, _, _, params, *_ = scheduled_args(background, worker)

    assert params["window"] == 400





def test_window_form_declaration_is_the_backward_compatible_default():

    default = inspect.signature(main.create_screening_job).parameters["window"].default

    assert default == 400





def test_missing_engine_returns_a_clear_500(screen_env):

    with patch.object(main, "engine_command", return_value=None):

        with pytest.raises(HTTPException) as exc:

            call_screen(

                panel_keys="rice",

                target_mode="sequence",

                target_sequence="ACGTACGTACGT",

            )



    assert exc.value.status_code == 500

    assert "dsGuard CLI was not found" in exc.value.detail













def test_job_document_records_input_provenance(screen_env):

    jobs_mock, _, _, _ = screen_env



    call_screen(

        panel_keys="rice",

        target_mode="sequence",

        target_sequence="ACGTACGTACGTACGTACGTACGT",

    )



    document = jobs_mock.insert_one.call_args.args[0]

    assert document["status"] == "queued"

    assert document["target_input_type"] == "pasted_sequence"

    assert document["panel_names"] == ["rice"]

    assert document["panel_filenames"] == ["rice.fasta"]

    assert document["parameters"]["max_hazard"] == 0.05



    json.dumps({k: v for k, v in document.items() if k != "_id"})
