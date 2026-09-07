"""Managed reference genomes shared by analysis and preprocessing workflows."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import tempfile
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REFERENCE_MANIFEST: dict[str, dict[str, str]] = {
    "hg38": {
        "fasta_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        "fasta_md5": "1c9dcaddfa41027f17cd8f7a82c7293b",
        "blacklist_url": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/61a04d2c5e49341d76735d485c61f0d1177d08a8/lists/hg38-blacklist.v2.bed.gz",
        "blacklist_md5": "83fe6bf8187a64dee8079b80f75ba289",
        "tss_gtf_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz",
        "tss_gtf_md5": "8486a6bdcd27a8a7a08232d01cc13b77",
        "macs_genome_size": "hs",
    },
    "mm10": {
        "fasta_url": "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz",
        "fasta_md5": "db005b65828db31735f384e4c5787be5",
        "blacklist_url": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/61a04d2c5e49341d76735d485c61f0d1177d08a8/lists/mm10-blacklist.v2.bed.gz",
        "blacklist_md5": "4ae47e40309533c2a71de55494cda9bc",
        "tss_gtf_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz",
        "tss_gtf_md5": "c5125258a0a2c5250ddb4c192abbf4e8",
        "macs_genome_size": "mm",
    },
}


@dataclass(frozen=True)
class AnalysisReference:
    """Resolved reference inputs for footprint-analysis commands."""

    assembly: str | None
    fasta: Path
    blacklist: Path | None
    managed: bool


def default_reference_dir() -> Path:
    """Return the default managed-reference cache without creating it."""

    return Path.home() / ".cache" / "fp-tools" / "references"


def is_managed_assembly(value: str | os.PathLike[str]) -> bool:
    return str(value).strip().lower() in REFERENCE_MANIFEST


def _digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _metadata_path(path: Path) -> Path:
    return path.with_name(path.name + ".fp-tools.json")


def _read_metadata(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return {}
    return value if isinstance(value, dict) else {}


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


class _DirectoryLock:
    """Small cross-platform exclusive lock for one reference directory."""

    def __init__(self, path: Path, timeout: float = 300.0):
        self.path = path
        self.timeout = timeout
        self.descriptor: int | None = None

    def __enter__(self):
        deadline = time.monotonic() + self.timeout
        while self.descriptor is None:
            try:
                self.descriptor = os.open(
                    self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600
                )
                os.write(self.descriptor, f"{os.getpid()}\n".encode("ascii"))
            except FileExistsError:
                try:
                    stale = time.time() - self.path.stat().st_mtime > self.timeout
                except FileNotFoundError:
                    continue
                if stale:
                    try:
                        self.path.unlink()
                    except FileNotFoundError:
                        pass
                    continue
                if time.monotonic() >= deadline:
                    raise TimeoutError(
                        f"Timed out waiting for managed-reference lock: {self.path}"
                    )
                time.sleep(0.05)
        return self

    def __exit__(self, exc_type, exc, traceback):
        if self.descriptor is not None:
            os.close(self.descriptor)
            self.descriptor = None
        self.path.unlink(missing_ok=True)
        return False


def _installed_asset_is_valid(
    output: Path, *, source_url: str, compressed_md5: str
) -> bool:
    if not output.is_file() or output.stat().st_size == 0:
        return False
    metadata = _read_metadata(_metadata_path(output))
    stat = output.stat()
    if (
        metadata.get("source_url") != source_url
        or metadata.get("compressed_md5") != compressed_md5
        or not metadata.get("sha256")
    ):
        return False
    if metadata.get("size") == stat.st_size and metadata.get("mtime_ns") == stat.st_mtime_ns:
        return True
    return _digest(output, "sha256") == metadata["sha256"]


def _install_gzip_asset(
    source_url: str,
    compressed_md5: str,
    output: Path,
    *,
    timeout: int = 120,
) -> bool:
    """Install a compressed remote asset atomically; return whether it changed."""

    if _installed_asset_is_valid(
        output, source_url=source_url, compressed_md5=compressed_md5
    ):
        return False

    download_descriptor, download_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".download.gz", dir=output.parent
    )
    os.close(download_descriptor)
    unpack_descriptor, unpack_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".unpacked", dir=output.parent
    )
    os.close(unpack_descriptor)
    download = Path(download_name)
    unpacked = Path(unpack_name)
    try:
        with (
            urllib.request.urlopen(source_url, timeout=timeout) as response,
            download.open("wb") as target,
        ):
            shutil.copyfileobj(response, target)
        if _digest(download, "md5") != compressed_md5:
            raise RuntimeError(
                f"Checksum mismatch while downloading managed reference asset: {source_url}"
            )
        with gzip.open(download, "rb") as source, unpacked.open("wb") as target:
            shutil.copyfileobj(source, target)
        if unpacked.stat().st_size == 0:
            raise RuntimeError(
                f"Managed reference download produced an empty asset: {source_url}"
            )
        installed_sha256 = _digest(unpacked, "sha256")
        os.replace(unpacked, output)
        stat = output.stat()
        _write_json_atomic(
            _metadata_path(output),
            {
                "compressed_md5": compressed_md5,
                "mtime_ns": stat.st_mtime_ns,
                "sha256": installed_sha256,
                "size": stat.st_size,
                "source_url": source_url,
            },
        )
        return True
    finally:
        download.unlink(missing_ok=True)
        unpacked.unlink(missing_ok=True)


def _write_fasta_index(fasta: Path, output: Path) -> None:
    """Write a samtools-compatible FAI without requiring an external program."""

    records: list[tuple[str, int, int, int, int]] = []
    name: str | None = None
    length = 0
    offset = 0
    line_bases = 0
    line_width = 0
    with fasta.open("rb") as handle:
        while True:
            position = handle.tell()
            raw = handle.readline()
            if not raw:
                break
            if raw.startswith(b">"):
                if name is not None:
                    records.append((name, length, offset, line_bases, line_width))
                header = raw[1:].strip().split(maxsplit=1)
                if not header:
                    raise ValueError(f"FASTA contains an empty sequence name: {fasta}")
                name = header[0].decode("utf-8")
                length = 0
                offset = handle.tell()
                line_bases = 0
                line_width = 0
                continue
            if name is None:
                if raw.strip():
                    raise ValueError(f"FASTA sequence appears before its header: {fasta}")
                continue
            sequence = raw.rstrip(b"\r\n")
            if not sequence:
                continue
            if line_bases == 0:
                offset = position
                line_bases = len(sequence)
                line_width = len(raw)
            length += len(sequence)
    if name is not None:
        records.append((name, length, offset, line_bases, line_width))
    if not records or any(record[1] <= 0 or record[3] <= 0 for record in records):
        raise ValueError(f"FASTA contains no nonempty sequences: {fasta}")
    with output.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write("\t".join(str(value) for value in record) + "\n")


def _fasta_index_is_valid(fasta: Path, index: Path) -> bool:
    if not index.is_file() or index.stat().st_size == 0:
        return False
    metadata = _read_metadata(_metadata_path(index))
    fasta_stat = fasta.stat()
    index_stat = index.stat()
    if (
        metadata.get("fasta_size") != fasta_stat.st_size
        or metadata.get("fasta_mtime_ns") != fasta_stat.st_mtime_ns
    ):
        return False
    if (
        metadata.get("size") != index_stat.st_size
        or metadata.get("mtime_ns") != index_stat.st_mtime_ns
        or metadata.get("sha256") != _digest(index, "sha256")
    ):
        return False
    try:
        rows = [line.rstrip().split("\t") for line in index.read_text().splitlines()]
        return bool(rows) and all(
            len(row) == 5 and all(int(value) >= 0 for value in row[1:])
            for row in rows
        )
    except (OSError, UnicodeError, ValueError):
        return False


def _ensure_fasta_index(fasta: Path) -> Path:
    index = Path(str(fasta) + ".fai")
    if _fasta_index_is_valid(fasta, index):
        return index
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{index.name}.", suffix=".tmp", dir=index.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        _write_fasta_index(fasta, temporary)
        index_sha256 = _digest(temporary, "sha256")
        os.replace(temporary, index)
        fasta_stat = fasta.stat()
        index_stat = index.stat()
        _write_json_atomic(
            _metadata_path(index),
            {
                "fasta_sha256": _digest(fasta, "sha256"),
                "fasta_mtime_ns": fasta_stat.st_mtime_ns,
                "fasta_size": fasta_stat.st_size,
                "mtime_ns": index_stat.st_mtime_ns,
                "sha256": index_sha256,
                "size": index_stat.st_size,
            },
        )
    finally:
        temporary.unlink(missing_ok=True)
    return index


def resolve_analysis_reference(
    genome: str | os.PathLike[str],
    *,
    reference_dir: str | os.PathLike[str] | None = None,
    blacklist: str | os.PathLike[str] | None = None,
    no_blacklist: bool = False,
    dry_run: bool = False,
) -> AnalysisReference:
    """Resolve a managed assembly name or a user-owned FASTA path."""

    if blacklist and no_blacklist:
        raise ValueError("--blacklist and --no-blacklist are mutually exclusive")
    genome_text = str(genome).strip()
    assembly = genome_text.lower()
    custom_blacklist = Path(blacklist).expanduser().resolve() if blacklist else None
    if custom_blacklist is not None and not custom_blacklist.is_file():
        raise ValueError(f"--blacklist file does not exist: {blacklist}")

    if assembly not in REFERENCE_MANIFEST:
        fasta = Path(genome_text).expanduser().resolve()
        if not fasta.is_file():
            raise ValueError(
                f"--genome must be hg38, mm10, or an existing FASTA file: {genome}"
            )
        return AnalysisReference(
            assembly=None,
            fasta=fasta,
            blacklist=None if no_blacklist else custom_blacklist,
            managed=False,
        )

    cache_root = Path(reference_dir).expanduser().resolve() if reference_dir else default_reference_dir()
    root = cache_root / assembly
    fasta = root / f"{assembly}.fa"
    managed_blacklist = root / f"{assembly}.blacklist.bed"
    selected_blacklist = None if no_blacklist else custom_blacklist or managed_blacklist
    if dry_run:
        return AnalysisReference(
            assembly=assembly,
            fasta=fasta,
            blacklist=selected_blacklist,
            managed=True,
        )

    root.mkdir(parents=True, exist_ok=True)
    manifest = REFERENCE_MANIFEST[assembly]
    with _DirectoryLock(root / ".reference.lock"):
        _install_gzip_asset(
            manifest["fasta_url"], manifest["fasta_md5"], fasta
        )
        _ensure_fasta_index(fasta)
        if selected_blacklist == managed_blacklist:
            _install_gzip_asset(
                manifest["blacklist_url"],
                manifest["blacklist_md5"],
                managed_blacklist,
            )
    return AnalysisReference(
        assembly=assembly,
        fasta=fasta,
        blacklist=selected_blacklist,
        managed=True,
    )
