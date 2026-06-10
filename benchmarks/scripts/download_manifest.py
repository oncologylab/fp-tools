#!/usr/bin/env python3
"""Download files listed in an fp-tools public-data manifest.

The script supports dry-runs and writes a TSV report. It is intentionally
conservative: it skips existing files unless --force is supplied and verifies
MD5 checksums when the manifest provides them.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib.request import urlretrieve


@dataclass
class DownloadResult:
    file_accession: str
    url: str
    local_path: str
    status: str
    message: str
    bytes: int
    checksum_ok: str


def read_manifest(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="	"))


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_one(row: dict[str, str], *, dry_run: bool = False, force: bool = False, downloader: str = "auto") -> DownloadResult:
    url = row.get("url", "")
    local_path = row.get("local_path", "")
    accession = row.get("file_accession", "")
    expected_md5 = row.get("checksum", "")
    if not url or not local_path:
        return DownloadResult(accession, url, local_path, "skipped", "missing url or local_path", 0, "not_checked")

    path = Path(local_path)
    if dry_run:
        return DownloadResult(accession, url, str(path), "dry_run", "planned", 0, "not_checked")
    if path.exists() and not force:
        checksum_ok = "not_checked"
        if expected_md5:
            checksum_ok = str(md5sum(path) == expected_md5).lower()
        return DownloadResult(accession, url, str(path), "skipped", "exists", path.stat().st_size, checksum_ok)

    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        if downloader == "aria2c" or (downloader == "auto" and shutil.which("aria2c")):
            subprocess.run(["aria2c", "-c", "-d", str(path.parent), "-o", path.name, url], check=True)
        elif downloader == "wget" or (downloader == "auto" and shutil.which("wget")):
            subprocess.run(["wget", "-c", "-O", str(path), url], check=True)
        else:
            urlretrieve(url, path)
    except Exception as exc:
        return DownloadResult(accession, url, str(path), "failed", str(exc), path.stat().st_size if path.exists() else 0, "not_checked")

    checksum_ok = "not_checked"
    if expected_md5:
        checksum_ok = str(md5sum(path) == expected_md5).lower()
        if checksum_ok == "false":
            return DownloadResult(accession, url, str(path), "failed", "md5 mismatch", path.stat().st_size, checksum_ok)
    return DownloadResult(accession, url, str(path), "downloaded", "ok", path.stat().st_size, checksum_ok)


def write_report(results: list[DownloadResult], path: str | Path) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(DownloadResult.__dataclass_fields__.keys()), delimiter="	")
        writer.writeheader()
        for result in results:
            writer.writerow(result.__dict__)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Input manifest TSV.")
    parser.add_argument("--report", default="benchmarks/download_reports/download_report.tsv", help="Output report TSV.")
    parser.add_argument("--dry-run", action="store_true", help="Write planned downloads without downloading files.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing files.")
    parser.add_argument("--downloader", choices=["auto", "aria2c", "wget", "urllib"], default="auto")
    args = parser.parse_args()

    rows = read_manifest(args.manifest)
    results = [download_one(row, dry_run=args.dry_run, force=args.force, downloader=args.downloader) for row in rows]
    write_report(results, args.report)
    failed = sum(result.status == "failed" for result in results)
    print(f"wrote {len(results)} download records to {args.report}; failed={failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
