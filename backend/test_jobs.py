

from unittest.mock import Mock, patch



import pytest

from fastapi import HTTPException



from backend import main


def test_fastapi_app_uses_supported_constructor_options():
    schema = main.app.openapi()

    assert schema["info"]["title"] == "dsGuard API"
    assert schema["info"]["version"] == "0.3.0"
    assert schema["openapi"]


def test_cors_origins_include_configured_frontend_origin():
    assert main.cors_origins("https://dsguard-ui.onrender.com/") == [
        "http://localhost:5173",
        "http://127.0.0.1:5173",
        "https://dsguard-ui.onrender.com",
    ]


@pytest.fixture
def deletion_setup(tmp_path):

    upload_root = tmp_path / "uploads"

    output_root = tmp_path / "outputs"

    upload_root.mkdir()

    output_root.mkdir()



    upload_dir = upload_root / "job-test"

    output_dir = output_root / "job-test"

    upload_dir.mkdir()

    output_dir.mkdir()



    with (

        patch.object(main, "UPLOAD_ROOT", upload_root),

        patch.object(main, "OUTPUT_ROOT", output_root),

        patch.object(main, "jobs") as jobs_mock,

    ):

        yield jobs_mock, upload_dir, output_dir





def test_delete_job_removes_record_and_files(deletion_setup):

    jobs_mock, upload_dir, output_dir = deletion_setup

    jobs_mock.find_one.return_value = {

        "_id": "job-test",

        "status": "completed",

    }

    jobs_mock.delete_one.return_value = Mock(deleted_count=1)



    result = main.delete_job("job-test")



    assert result["status"] == "deleted"

    jobs_mock.delete_one.assert_called_once_with({"_id": "job-test"})

    assert not upload_dir.exists()

    assert not output_dir.exists()





def test_delete_job_returns_404_for_missing_job(deletion_setup):

    jobs_mock, _, _ = deletion_setup

    jobs_mock.find_one.return_value = None



    with pytest.raises(HTTPException) as exc:

        main.delete_job("missing-job")



    assert exc.value.status_code == 404

    jobs_mock.delete_one.assert_not_called()





@pytest.mark.parametrize("status", ["queued", "running"])

def test_delete_job_rejects_active_job(deletion_setup, status):

    jobs_mock, upload_dir, output_dir = deletion_setup

    jobs_mock.find_one.return_value = {

        "_id": "job-test",

        "status": status,

    }



    with pytest.raises(HTTPException) as exc:

        main.delete_job("job-test")



    assert exc.value.status_code == 409

    jobs_mock.delete_one.assert_not_called()

    assert upload_dir.exists()

    assert output_dir.exists()





def test_delete_job_rejects_unsafe_path(deletion_setup):

    jobs_mock, _, _ = deletion_setup

    jobs_mock.find_one.return_value = {

        "_id": "..",

        "status": "completed",

    }



    with pytest.raises(HTTPException) as exc:

        main.delete_job("..")



    assert exc.value.status_code == 400

    jobs_mock.delete_one.assert_not_called()



def test_delete_job_retains_record_when_file_cleanup_fails(deletion_setup):

    jobs_mock, upload_dir, output_dir = deletion_setup

    jobs_mock.find_one.return_value = {

        "_id": "job-test",

        "status": "completed",

    }



    with patch.object(main.shutil, "rmtree", side_effect=OSError("Permission denied")):

        with pytest.raises(HTTPException) as exc:

            main.delete_job("job-test")



    assert exc.value.status_code == 500

    assert "record was retained" in exc.value.detail

    jobs_mock.delete_one.assert_not_called()

    assert upload_dir.exists()

    assert output_dir.exists()





def test_clear_job_history_removes_owned_artifacts_only(tmp_path):

    upload_root = tmp_path / "uploads"

    output_root = tmp_path / "outputs"

    panel_root = tmp_path / "reference_panels"

    catalog = tmp_path / "target_catalog.json"

    upload_root.mkdir()

    output_root.mkdir()

    panel_root.mkdir()

    catalog.write_text("catalog", encoding="utf-8")

    (panel_root / "rice.fasta").write_text(">rice\nACGT\n", encoding="utf-8")

    (upload_root / "job-a").mkdir()

    (output_root / "job-a").mkdir()



    with (

        patch.object(main, "UPLOAD_ROOT", upload_root),

        patch.object(main, "OUTPUT_ROOT", output_root),

        patch.object(main, "jobs") as jobs_mock,

    ):

        jobs_mock.find.return_value = iter([{"_id": "job-a", "status": "completed"}])

        jobs_mock.delete_many.return_value = Mock(deleted_count=1)



        result = main.clear_job_history()



    assert result["deleted_count"] == 1

    jobs_mock.delete_many.assert_called_once_with({})

    assert not (upload_root / "job-a").exists()

    assert not (output_root / "job-a").exists()

    assert (panel_root / "rice.fasta").exists()

    assert catalog.exists()





def test_clear_job_history_rejects_active_jobs(tmp_path):

    upload_root = tmp_path / "uploads"

    output_root = tmp_path / "outputs"

    upload_root.mkdir()

    output_root.mkdir()



    with (

        patch.object(main, "UPLOAD_ROOT", upload_root),

        patch.object(main, "OUTPUT_ROOT", output_root),

        patch.object(main, "jobs") as jobs_mock,

    ):

        jobs_mock.find.return_value = iter([{"_id": "job-a", "status": "running"}])

        with pytest.raises(HTTPException) as exc:

            main.clear_job_history()



    assert exc.value.status_code == 409

    jobs_mock.delete_many.assert_not_called()
