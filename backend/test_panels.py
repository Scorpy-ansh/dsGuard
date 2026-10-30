import hashlib

import io

import urllib.error

from email.message import Message

from backend import main


class FakeResponse:
    def __init__(self, content, url="https://example.invalid/rice.fasta", status=200):
        self.content = content
        self.status = status
        self.url = url
        self.headers = Message()
        self.headers["Content-Type"] = "application/octet-stream"

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return None

    def read(self, _size):
        content, self.content = self.content, b""
        return content

    def geturl(self):
        return self.url


class InterruptedResponse(FakeResponse):
    def read(self, _size):
        raise OSError("connection interrupted")

def test_rice_panel_download_is_verified_and_cached(tmp_path, monkeypatch):
    fasta = b">rice\nACGTACGT\n"
    checksum = hashlib.sha256(fasta).hexdigest()
    monkeypatch.setattr(main, "PANEL_CACHE_DIR", tmp_path)
    monkeypatch.setenv("RICE_PANEL_URL", "https://example.invalid/rice.fasta")
    monkeypatch.setenv("RICE_PANEL_SHA256", checksum)
    monkeypatch.setattr(main.urllib.request, "urlopen", lambda *_args, **_kwargs: FakeResponse(fasta))

    path = main._download_rice_panel()

    assert path == tmp_path / "rice.fasta"
    assert path.read_bytes() == fasta

    def fail_download(*_args, **_kwargs):
        raise AssertionError("cached panel should not be downloaded")

    monkeypatch.setattr(main.urllib.request, "urlopen", fail_download)
    assert main._download_rice_panel() == path


def test_rice_panel_download_rejects_bad_checksum(tmp_path, monkeypatch):
    monkeypatch.setattr(main, "PANEL_CACHE_DIR", tmp_path)
    monkeypatch.setenv("RICE_PANEL_URL", "https://example.invalid/rice.fasta")
    monkeypatch.setenv("RICE_PANEL_SHA256", "0" * 64)
    monkeypatch.setattr(
        main.urllib.request,
        "urlopen",
        lambda *_args, **_kwargs: FakeResponse(b">rice\nACGT\n"),
    )

    try:
        main._download_rice_panel()
    except RuntimeError as exc:
        assert "SHA256" in str(exc)
    else:
        raise AssertionError("checksum failure was not reported")
    assert not list(tmp_path.glob("*.fasta"))


def test_rice_panel_download_rejects_invalid_fasta(tmp_path, monkeypatch):
    fasta = b"not-a-fasta\n"
    monkeypatch.setattr(main, "PANEL_CACHE_DIR", tmp_path)
    monkeypatch.setenv("RICE_PANEL_URL", "https://example.invalid/rice.fasta")
    monkeypatch.setenv("RICE_PANEL_SHA256", hashlib.sha256(fasta).hexdigest())
    monkeypatch.setattr(
        main.urllib.request,
        "urlopen",
        lambda *_args, **_kwargs: FakeResponse(fasta),
    )

    try:
        main._download_rice_panel()
    except RuntimeError as exc:
        assert "FASTA" in str(exc)
    else:
        raise AssertionError("invalid FASTA was not reported")
    assert not list(tmp_path.glob("*.fasta"))


def test_rice_panel_download_requires_configuration(tmp_path, monkeypatch):
    monkeypatch.setattr(main, "PANEL_CACHE_DIR", tmp_path)
    monkeypatch.delenv("RICE_PANEL_URL", raising=False)
    monkeypatch.delenv("RICE_PANEL_SHA256", raising=False)

    assert main._download_rice_panel() is None


def test_rice_panel_download_reports_network_failure(tmp_path, monkeypatch):
    monkeypatch.setattr(main, "PANEL_CACHE_DIR", tmp_path)
    monkeypatch.setenv("RICE_PANEL_URL", "https://example.invalid/rice.fasta")
    monkeypatch.setenv("RICE_PANEL_SHA256", "0" * 64)

    def fail_download(*_args, **_kwargs):
        raise OSError("network unavailable")

    monkeypatch.setattr(main.urllib.request, "urlopen", fail_download)

    try:
        main._download_rice_panel()
    except RuntimeError as exc:
        assert "Could not download" in str(exc)
    else:
        raise AssertionError("download failure was not reported")
    assert not list(tmp_path.iterdir())


def test_rice_panel_download_reports_http_403_and_response_body(
    tmp_path, monkeypatch, caplog
):
    monkeypatch.setattr(main, "PANEL_CACHE_DIR", tmp_path)
    monkeypatch.setenv("RICE_PANEL_URL", "https://example.invalid/rice.fasta")
    monkeypatch.setenv("RICE_PANEL_SHA256", "0" * 64)
    error = urllib.error.HTTPError(
        "https://example.invalid/rice.fasta",
        403,
        "Forbidden",
        {"Content-Type": "text/plain"},
        io.BytesIO(b"blocked by edge policy"),
    )
    monkeypatch.setattr(main.urllib.request, "urlopen", lambda *_args, **_kwargs: (_ for _ in ()).throw(error))

    with caplog.at_level("ERROR"):
        try:
            main._download_rice_panel()
        except RuntimeError as exc:
            assert "HTTP 403" in str(exc)
        else:
            raise AssertionError("HTTP 403 was not reported")

    assert "blocked by edge policy" in caplog.text
    assert "example.invalid/rice.fasta" in caplog.text
    assert not list(tmp_path.iterdir())


def test_rice_panel_download_logs_status_headers_and_redirect(
    tmp_path, monkeypatch, caplog
):
    fasta = b">rice\nACGTACGT\n"
    monkeypatch.setattr(main, "PANEL_CACHE_DIR", tmp_path)
    monkeypatch.setenv("RICE_PANEL_URL", "https://example.invalid/rice.fasta")
    monkeypatch.setenv("RICE_PANEL_SHA256", hashlib.sha256(fasta).hexdigest())
    response = FakeResponse(
        fasta,
        url="https://cdn.example.invalid/rice.fasta",
    )
    requests = []

    def open_request(request, **_kwargs):
        requests.append(request)
        return response

    monkeypatch.setattr(main.urllib.request, "urlopen", open_request)

    with caplog.at_level("INFO"):
        path = main._download_rice_panel()

    assert path == tmp_path / "rice.fasta"
    assert "status=200" in caplog.text
    assert "content_type=application/octet-stream" in caplog.text
    assert "method=GET" in caplog.text
    assert "redirected=True" in caplog.text
    assert "cdn.example.invalid/rice.fasta" in caplog.text
    assert requests[0].method == "GET"
    assert requests[0].get_header("Accept") == "*/*"
    assert requests[0].get_header("User-agent") == "dsGuard/0.4"


def test_rice_panel_download_removes_partial_file_on_interruption(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(main, "PANEL_CACHE_DIR", tmp_path)
    monkeypatch.setenv("RICE_PANEL_URL", "https://example.invalid/rice.fasta")
    monkeypatch.setenv("RICE_PANEL_SHA256", "0" * 64)
    monkeypatch.setattr(
        main.urllib.request,
        "urlopen",
        lambda *_args, **_kwargs: InterruptedResponse(b"partial"),
    )

    try:
        main._download_rice_panel()
    except RuntimeError as exc:
        assert "Could not download" in str(exc)
    else:
        raise AssertionError("interrupted download was not reported")
    assert not list(tmp_path.iterdir())


def test_missing_rice_panel_configuration_is_explained(tmp_path, monkeypatch):
    monkeypatch.setattr(main, "PANEL_DIR", tmp_path)
    monkeypatch.delenv("RICE_PANEL_URL", raising=False)
    monkeypatch.delenv("RICE_PANEL_SHA256", raising=False)

    assert main.panel_path("rice") is None


def test_panel_path_prefers_existing_local_rice_panel(tmp_path, monkeypatch):
    panel_dir = tmp_path / "local"
    panel_dir.mkdir()
    local_path = panel_dir / "rice.fasta"
    local_path.write_text(">local\nACGT\n", encoding="ascii")
    monkeypatch.setattr(main, "PANEL_DIR", panel_dir)
    monkeypatch.setattr(main, "_download_rice_panel", lambda: (_ for _ in ()).throw(
        AssertionError("local panel should be preferred")
    ))

    assert main.panel_path("rice") == local_path


def test_panel_status_endpoint_initializes_rice_once(tmp_path, monkeypatch):
    fasta = b">rice Oryza sativa\nACGTACGT\n"
    checksum = hashlib.sha256(fasta).hexdigest()
    monkeypatch.setattr(main, "PANEL_DIR", tmp_path / "local")
    monkeypatch.setattr(main, "PANEL_CACHE_DIR", tmp_path / "cache")
    monkeypatch.setenv("RICE_PANEL_URL", "https://example.invalid/rice.fasta")
    monkeypatch.setenv("RICE_PANEL_SHA256", checksum)
    calls = []

    def download(*_args, **_kwargs):
        calls.append(True)
        return FakeResponse(fasta)

    monkeypatch.setattr(main.urllib.request, "urlopen", download)
    main._PANEL_STATES.clear()

    first = {panel["key"]: panel for panel in main.list_reference_panels()}
    second = {panel["key"]: panel for panel in main.list_reference_panels()}

    assert len(calls) == 1
    assert first["rice"]["installed"] is True
    assert first["rice"]["file_exists"] is True
    assert first["rice"]["validation_status"] == "passed"
    assert first["rice"]["organism"] == "Oryza sativa"
    assert first["rice"]["size_bytes"] == len(fasta)
    assert first["rice"]["path"] == second["rice"]["path"]


def test_panel_status_reports_invalid_local_fasta(tmp_path, monkeypatch):
    panel_dir = tmp_path / "local"
    panel_dir.mkdir()
    (panel_dir / "rice.fasta").write_bytes(b"not fasta")
    monkeypatch.setattr(main, "PANEL_DIR", panel_dir)
    monkeypatch.delenv("RICE_PANEL_URL", raising=False)
    monkeypatch.delenv("RICE_PANEL_SHA256", raising=False)
    main._PANEL_STATES.clear()

    panels = {panel["key"]: panel for panel in main.list_reference_panels()}

    assert panels["rice"]["installed"] is False
    assert panels["rice"]["file_exists"] is True
    assert panels["rice"]["validation_status"] == "failed"
    assert "FASTA" in panels["rice"]["last_validation_result"]
