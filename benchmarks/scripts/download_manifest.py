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


def download_one(
    row: dict[str, str],
    *,
    dry_run: bool = False,
    force: bool = False,
    downloader: str = "auto",
) -> DownloadResult:
    url = row.get("url", "")
    local_path = row.get("local_path", "")
    accession = row.get("file_accession", "")
    expected_md5 = row.get("checksum", "")
    expected_bytes = int(row.get("expected_bytes") or 0)
    if not url or not local_path:
        return DownloadResult(
            accession,
            url,
            local_path,
            "skipped",
            "missing url or local_path",
            0,
            "not_checked",
        )

    path = Path(local_path)
    if dry_run:
        return DownloadResult(
            accession, url, str(path), "dry_run", "planned", 0, "not_checked"
        )
    if path.exists() and not force:
        size_ok = not expected_bytes or path.stat().st_size == expected_bytes
        checksum_ok = not expected_md5 or md5sum(path) == expected_md5
        if size_ok and checksum_ok:
            return DownloadResult(
                accession,
                url,
                str(path),
                "skipped",
                "verified existing file",
                path.stat().st_size,
                str(checksum_ok).lower() if expected_md5 else "not_checked",
            )
        path.unlink()

    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        if downloader == "aria2c" or (downloader == "auto" and shutil.which("aria2c")):
            subprocess.run(
                [
                    "aria2c",
                    "-c",
                    "--max-tries=5",
                    "--retry-wait=5",
                    "--file-allocation=none",
                    "-x",
                    "4",
                    "-s",
                    "4",
                    "-d",
                    str(path.parent),
                    "-o",
                    path.name,
                    url,
                ],
                check=True,
            )
        elif downloader == "wget" or (downloader == "auto" and shutil.which("wget")):
            subprocess.run(["wget", "-c", "-O", str(path), url], check=True)
        else:
            urlretrieve(url, path)
    except Exception as exc:
        return DownloadResult(
            accession,
            url,
            str(path),
            "failed",
            str(exc),
            path.stat().st_size if path.exists() else 0,
            "not_checked",
        )

    checksum_ok = "not_checked"
    if expected_bytes and path.stat().st_size != expected_bytes:
        observed = path.stat().st_size
        path.unlink(missing_ok=True)
        return DownloadResult(
            accession,
            url,
            str(path),
            "failed",
            f"size mismatch: expected {expected_bytes}, found {observed}",
            observed,
            checksum_ok,
        )
    if expected_md5:
        checksum_ok = str(md5sum(path) == expected_md5).lower()
        if checksum_ok == "false":
            observed = path.stat().st_size
            path.unlink(missing_ok=True)
            return DownloadResult(
                accession,
                url,
                str(path),
                "failed",
                "md5 mismatch",
                observed,
                checksum_ok,
            )
    return DownloadResult(
        accession, url, str(path), "downloaded", "ok", path.stat().st_size, checksum_ok
    )


def write_report(results: list[DownloadResult], path: str | Path) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(DownloadResult.__dataclass_fields__.keys()),
            delimiter="	",
        )
        writer.writeheader()
        for result in results:
            writer.writerow(result.__dict__)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, help="Input manifest TSV.")
    parser.add_argument(
        "--report",
        default="benchmarks/download_reports/download_report.tsv",
        help="Output report TSV.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write planned downloads without downloading files.",
    )
    parser.add_argument(
        "--force", action="store_true", help="Overwrite existing files."
    )
    parser.add_argument(
        "--include",
        action="append",
        default=[],
        help="Download only matching file or run accessions; repeat as needed.",
    )
    parser.add_argument(
        "--downloader", choices=["auto", "aria2c", "wget", "urllib"], default="auto"
    )
    args = parser.parse_args()

    rows = read_manifest(args.manifest)
    if args.include:
        wanted = set(args.include)
        rows = [
            row
            for row in rows
            if row.get("file_accession") in wanted or row.get("run_accession") in wanted
        ]
        found = {
            value
            for row in rows
            for value in (row.get("file_accession"), row.get("run_accession"))
            if value in wanted
        }
        missing = sorted(wanted - found)
        if missing:
            parser.error(f"unknown requested accession(s): {', '.join(missing)}")
    results = [
        download_one(
            row, dry_run=args.dry_run, force=args.force, downloader=args.downloader
        )
        for row in rows
    ]
    write_report(results, args.report)
    failed = sum(result.status == "failed" for result in results)
    print(f"wrote {len(results)} download records to {args.report}; failed={failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
