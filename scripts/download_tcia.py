#!/usr/bin/env python3
"""Download chosen patients of a *public* TCIA collection through the NBIA API (stdlib only).

Dry run by default: prints what would be fetched (series, size, licence) and the free disk space,
and downloads nothing (and creates nothing) until you add --yes.

    python scripts/download_tcia.py Vestibular-Schwannoma-SEG data/raw/vs-seg VS-SEG-154 VS-SEG-135
    python scripts/download_tcia.py Pancreatic-CT-CBCT-SEG data/raw/pancreas-ct-cbct-seg \
        Pancreas-CT-CB_014 --yes

Each series is extracted into ``<out>/<patient>/<MODALITY>_<n>_<uid tail>/``. Extraction happens in
a ``*.partial`` directory that is renamed only after the zip verified, the file count matched the
series' ``ImageCount`` and a ``.complete`` marker was written, so an interrupted or corrupt download
is never mistaken for a finished one and simply restarts. A ``MANIFEST_<patient>.json`` (series
UIDs, sizes, licence) is written only after every series of that patient succeeded.

Collections that need a login or a signed data-use agreement are not reachable this way; accept
those terms yourself on the TCIA site. Downloaded data is licensed by its collection (the licence
text is inside every series); cite the dataset, see data/README.md.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath

BASE = "https://services.cancerimagingarchive.net/nbia-api/services/v1/"
COMPLETE = ".complete"
MANIFEST_KEYS = (
    "PatientID",
    "Modality",
    "SeriesDescription",
    "SeriesInstanceUID",
    "StudyInstanceUID",
    "ImageCount",
    "FileSize",
    "LicenseName",
    "LicenseURI",
    "Collection",
)


def api(base: str, endpoint: str, **params):
    url = base + endpoint + "?" + urllib.parse.urlencode({**params, "format": "json"})
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)


def _dcm_count(directory: Path) -> int:
    return sum(1 for f in directory.rglob("*") if f.is_file() and f.suffix.lower() == ".dcm")


def is_complete(directory: Path, expected_images: int) -> bool:
    """Finished series: marker present, or (older downloads) the DICOM count matches exactly."""
    if not directory.is_dir():
        return False
    if (directory / COMPLETE).exists():
        return True
    if _dcm_count(directory) == expected_images and not list(directory.glob("*.zip")):
        (directory / COMPLETE).write_text("verified by file count\n")
        return True
    return False


def _safe_extract(zpath: Path, dest: Path, max_bytes: float) -> None:
    with zipfile.ZipFile(zpath) as z:
        bad = z.testzip()
        if bad is not None:
            raise zipfile.BadZipFile(f"corrupt member {bad}")
        total = 0
        for info in z.infolist():
            p = PurePosixPath(info.filename)
            if p.is_absolute() or ".." in p.parts:
                raise zipfile.BadZipFile(f"unsafe path in archive: {info.filename}")
            total += info.file_size
        if total > max_bytes:
            raise zipfile.BadZipFile(f"archive unpacks to {total / 1e6:.0f} MB, far above expected")
        z.extractall(dest)


def fetch_series(base: str, series: dict, dest: Path, retries: int = 3) -> None:
    """Download one series atomically into ``dest`` (raises after the last failed attempt)."""
    expected = int(series["ImageCount"])
    max_bytes = 4 * float(series["FileSize"]) + 50e6
    url = (
        base
        + "getImage?"
        + urllib.parse.urlencode({"SeriesInstanceUID": series["SeriesInstanceUID"]})
    )
    partial = dest.with_name(dest.name + ".partial")
    last: Exception | None = None
    for attempt in range(1, retries + 1):
        shutil.rmtree(partial, ignore_errors=True)
        partial.mkdir(parents=True)
        try:
            zpath = partial / "_series.zip"
            with urllib.request.urlopen(url, timeout=300) as r, open(zpath, "wb") as f:
                shutil.copyfileobj(r, f, 1 << 20)
            _safe_extract(zpath, partial, max_bytes)
            zpath.unlink()
            n = _dcm_count(partial)
            if n != expected:
                raise OSError(f"got {n} DICOM files, expected {expected}")
            (partial / COMPLETE).write_text("ok\n")
            shutil.rmtree(dest, ignore_errors=True)
            partial.rename(dest)
            return
        except (OSError, zipfile.BadZipFile, urllib.error.URLError) as exc:
            last = exc
            shutil.rmtree(partial, ignore_errors=True)
            if attempt < retries:
                print(f"  attempt {attempt} failed ({exc}); retrying", file=sys.stderr)
                time.sleep(min(2**attempt, 10))
    raise RuntimeError(f"download failed for {dest.name}: {last}") from last


def _existing_parent(path: Path) -> Path:
    path = path.resolve()
    while not path.exists():
        path = path.parent
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("collection")
    ap.add_argument("outdir", type=Path)
    ap.add_argument("patients", nargs="+")
    ap.add_argument("--yes", action="store_true", help="actually download (default: dry run)")
    ap.add_argument(
        "--min-free-gb",
        type=float,
        default=2.0,
        help="refuse to go below this much free disk (default 2)",
    )
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--base-url", default=BASE, help=argparse.SUPPRESS)  # for tests
    args = ap.parse_args(argv)
    base = args.base_url if args.base_url.endswith("/") else args.base_url + "/"

    if args.collection not in {c["Collection"] for c in api(base, "getCollectionValues")}:
        sys.exit(f"{args.collection!r} is not in the public listing (login/agreement required?)")
    series = api(base, "getSeries", Collection=args.collection)
    chosen = {
        p: sorted(
            (s for s in series if s["PatientID"] == p),
            key=lambda s: (s["Modality"], s["SeriesInstanceUID"]),
        )
        for p in args.patients
    }
    missing = [p for p, rows in chosen.items() if not rows]
    if missing:
        sys.exit(f"no series for: {missing}")

    def dest_of(p: str, i: int, s: dict) -> Path:
        return args.outdir / p / f"{s['Modality']}_{i}_{s['SeriesInstanceUID'][-8:]}"

    todo_bytes = sum(
        float(s["FileSize"])
        for p, rows in chosen.items()
        for i, s in enumerate(rows)
        if not is_complete(dest_of(p, i, s), int(s["ImageCount"]))
    )
    licences = {s.get("LicenseName") for rows in chosen.values() for s in rows}
    free = shutil.disk_usage(_existing_parent(args.outdir)).free
    for p, rows in chosen.items():
        mods: dict[str, int] = {}
        for s in rows:
            mods[s["Modality"]] = mods.get(s["Modality"], 0) + 1
        size = sum(float(s["FileSize"]) for s in rows) / 1e6
        print(f"{p}: {len(rows)} series {mods}, {size:.0f} MB")
    print(
        f"to download {todo_bytes / 1e6:.0f} MB | free {free / 1e9:.1f} GB | "
        f"licence: {', '.join(map(str, licences))}"
    )
    if free - todo_bytes < args.min_free_gb * 1e9:
        sys.exit(f"would leave less than {args.min_free_gb} GB free; aborting")
    if not args.yes:
        print("dry run: nothing downloaded (add --yes to proceed)")
        return 0

    try:
        for p, rows in chosen.items():
            manifest = []
            for i, s in enumerate(rows):
                d = dest_of(p, i, s)
                manifest.append({k: s.get(k) for k in MANIFEST_KEYS} | {"dir": str(d)})
                if is_complete(d, int(s["ImageCount"])):
                    print("skip (complete)", d.name)
                    continue
                d.parent.mkdir(parents=True, exist_ok=True)
                t0 = time.time()
                fetch_series(base, s, d, args.retries)
                print(
                    f"{p} {s['Modality']:9s} {s['ImageCount']:>4} images "
                    f"{float(s['FileSize']) / 1e6:6.1f} MB {time.time() - t0:4.1f}s"
                )
            (args.outdir / f"MANIFEST_{p}.json").write_text(json.dumps(manifest, indent=1))
    except RuntimeError as exc:
        sys.exit(f"error: {exc}\nre-run the same command to resume; finished series are kept")
    return 0


if __name__ == "__main__":
    sys.exit(main())
