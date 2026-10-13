"""scripts/download_tcia.py against a local fake NBIA server: corrupt, truncated and hostile
downloads must never be mistaken for complete ones."""

import importlib.util
import io
import json
import threading
import zipfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parent.parent / "scripts" / "download_tcia.py"
spec = importlib.util.spec_from_file_location("download_tcia", SCRIPT)
dl = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dl)

SERIES = [
    {
        "PatientID": "P1",
        "Modality": "CT",
        "SeriesInstanceUID": "1.2.3.4.5678",
        "StudyInstanceUID": "1.2.3",
        "ImageCount": 3,
        "FileSize": "3000",
        "LicenseName": "CC BY 4.0",
        "LicenseURI": "u",
        "Collection": "Fake",
    },
    {
        "PatientID": "P1",
        "Modality": "RTDOSE",
        "SeriesInstanceUID": "1.2.3.4.9999",
        "StudyInstanceUID": "1.2.3",
        "ImageCount": 1,
        "FileSize": "1000",
        "LicenseName": "CC BY 4.0",
        "LicenseURI": "u",
        "Collection": "Fake",
    },
]


def make_zip(n_dcm: int, extra: dict[str, bytes] | None = None) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for i in range(n_dcm):
            z.writestr(f"1-{i:02d}.dcm", b"x" * 100)
        z.writestr("LICENSE", "CC BY 4.0")
        for name, data in (extra or {}).items():
            z.writestr(name, data)
    return buf.getvalue()


@pytest.fixture
def server():
    state = {"plan": []}  # queue of behaviours for successive getImage calls; default "ok"

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def _send(self, body: bytes, ctype="application/json"):
            self.send_response(200)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            if self.path.startswith("/getCollectionValues"):
                return self._send(json.dumps([{"Collection": "Fake"}]).encode())
            if self.path.startswith("/getSeries"):
                return self._send(json.dumps(SERIES).encode())
            n = 1 if SERIES[1]["SeriesInstanceUID"] in self.path else 3
            mode = state["plan"].pop(0) if state["plan"] else "ok"
            if mode == "html":
                return self._send(b"<html>Service unavailable</html>", "text/html")
            if mode == "truncated":
                return self._send(make_zip(n)[:150], "application/zip")
            if mode == "short":
                return self._send(make_zip(max(n - 1, 0)), "application/zip")
            if mode == "traversal":
                return self._send(make_zip(n, {"../evil.txt": b"pwned"}), "application/zip")
            return self._send(make_zip(n), "application/zip")

    srv = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    state["url"] = f"http://127.0.0.1:{srv.server_address[1]}"
    yield state
    srv.shutdown()


def run(server, out, *extra):
    return dl.main(
        ["Fake", str(out), "P1", "--base-url", server["url"], "--min-free-gb", "0", *extra]
    )


def test_dry_run_creates_nothing(server, tmp_path):
    out = tmp_path / "data" / "x"
    assert run(server, out) == 0
    assert not out.exists() and not (tmp_path / "data").exists()


def test_happy_path_is_complete_and_rerun_skips(server, tmp_path, capsys):
    out = tmp_path / "out"
    assert run(server, out, "--yes") == 0
    dirs = sorted(d for d in (out / "P1").iterdir())
    assert len(dirs) == 2 and all((d / dl.COMPLETE).exists() for d in dirs)
    assert not list(out.rglob("*.partial")) and not list(out.rglob("*.zip"))
    assert json.loads((out / "MANIFEST_P1.json").read_text())[0]["LicenseName"] == "CC BY 4.0"
    capsys.readouterr()
    server["plan"] = ["html"] * 10  # any request now would fail: a rerun must not make one
    assert run(server, out, "--yes") == 0
    assert capsys.readouterr().out.count("skip (complete)") == 2


@pytest.mark.parametrize("bad", ["html", "truncated", "short"])
def test_bad_download_is_not_kept_and_a_rerun_recovers(server, tmp_path, bad):
    out = tmp_path / "out"
    server["plan"] = [bad]
    with pytest.raises(SystemExit, match="re-run"):
        run(server, out, "--yes", "--retries", "1")
    assert not list(out.rglob(dl.COMPLETE)) and not list(out.rglob("*.partial"))
    assert not list(out.rglob("*.zip")) and not (out / "MANIFEST_P1.json").exists()
    assert run(server, out, "--yes", "--retries", "1") == 0  # same command resumes cleanly
    assert len(list(out.rglob(dl.COMPLETE))) == 2


def test_retry_recovers_from_a_transient_failure(server, tmp_path, monkeypatch):
    monkeypatch.setattr(dl.time, "sleep", lambda s: None)
    server["plan"] = ["truncated", "html"]
    assert run(server, tmp_path / "out", "--yes", "--retries", "3") == 0


def test_path_traversal_in_archive_is_rejected(server, tmp_path):
    out = tmp_path / "nested" / "out"
    server["plan"] = ["traversal"]
    with pytest.raises(SystemExit):
        run(server, out, "--yes", "--retries", "1")
    assert not (tmp_path / "nested" / "evil.txt").exists() and not (out / "evil.txt").exists()
    assert not list(tmp_path.rglob("evil.txt"))


def test_unknown_collection_and_patient_exit_with_a_message(server, tmp_path):
    with pytest.raises(SystemExit, match="public listing"):
        dl.main(["Nope", str(tmp_path), "P1", "--base-url", server["url"]])
    with pytest.raises(SystemExit, match="no series"):
        dl.main(["Fake", str(tmp_path), "P404", "--base-url", server["url"]])


def test_legacy_directory_without_marker_counts_as_complete_only_if_file_count_matches(tmp_path):
    d = tmp_path / "CT_0_x"
    d.mkdir()
    for i in range(3):
        (d / f"{i}.dcm").write_bytes(b"x")
    (d / "LICENSE").write_text("cc")
    assert dl.is_complete(d, 3) and (d / dl.COMPLETE).exists()
    d2 = tmp_path / "CT_1_y"
    d2.mkdir()
    (d2 / "0.dcm").write_bytes(b"x")
    assert not dl.is_complete(d2, 3)
