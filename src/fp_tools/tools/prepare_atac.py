"""Prepare raw ATAC-seq reads as stable inputs for fp-tools.

The module intentionally remains an orchestration layer.  It validates sample
metadata, resolves public FASTQs, invokes established command-line tools, and
writes a portable project layout; scientific footprint analysis stays in the
existing fp-tools commands.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed

from fp_tools.utils import bigwig as pyBigWig
try:
    import pysam
except ImportError:  # Raw-read preparation is documented through WSL on Windows.
    pysam = None
import yaml


DEFAULTS: dict[str, Any] = {
    "profile": "modern",
    "download": {"provider": "auto", "retries": 5, "keep_fastq": True},
    "trim": {"enabled": True, "min_length": 20, "extra_args": []},
    "align": {
        "mapq": 30,
        "max_insert": 2000,
        "extra_args": [
            "--very-sensitive",
            "--no-mixed",
            "--no-discordant",
            "--dovetail",
        ],
    },
    "filter": {"remove_duplicates": True, "remove_mito": True, "extra_args": []},
    "peaks": {
        "qvalue": 0.01,
        "single_end_shift": -75,
        "single_end_extsize": 150,
        "extra_args": [],
    },
    "tracks": {"normalization_reads": 10_000_000, "extra_args": []},
    "qc": {
        "fastqc": True,
        "multiqc": True,
        "warn_frip_below": 0.2,
        "warn_tss_below": {"hg38": 5.0, "mm10": 10.0},
        "fail_on_warning": False,
    },
    "cleanup": {"keep_intermediates": False},
    "resources": {
        "cores": max(1, os.cpu_count() or 1),
        "max_parallel_samples": 1,
        "memory_gb": None,
    },
}

PUBLIC_PROFILES = ("modern", "homer-atac")
PROFILE_ALIASES = {"legacy-atac": "homer-atac"}

PROFILE_DEFAULTS: dict[str, dict[str, Any]] = {
    "modern": {},
    "homer-atac": {
        "align": {
            "mapq": 0,
            "max_insert": 1000,
            "extra_args": [
                "--very-sensitive-local",
                "--no-mixed",
                "--dovetail",
                "--phred33",
            ],
        },
        "filter": {"remove_duplicates": True, "remove_mito": False, "extra_args": []},
        "peaks": {
            "qvalue": 0.00001,
            "single_end_shift": -75,
            "single_end_extsize": 150,
            "extra_args": [],
            "homer_style": "factor",
            "homer_local_fold": 15,
            "homer_local_size": 150000,
        },
        "resources": {
            "cores": max(1, os.cpu_count() or 1),
            "max_parallel_samples": 1,
            "memory_gb": 24,
            "sample_memory_gb": 16,
        },
    },
}

REFERENCE_MANIFEST = {
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

MITO_CHROMS = {"chrM", "chrMT", "M", "MT", "Mito"}
PROFILE_EXECUTABLES = {
    "modern": ("fastp", "bowtie2", "bowtie2-build", "samtools", "bedtools", "macs3"),
    "homer-atac": (
        "trim_galore",
        "fastqc",
        "bowtie2",
        "bowtie2-build",
        "samtools",
        "bedtools",
        "picard",
        "makeTagDirectory",
        "makeUCSCfile",
        "findPeaks",
        "pos2bed.pl",
    ),
}


def normalize_profile(profile: str) -> str:
    """Return the canonical public name for an ATAC preprocessing profile."""
    canonical = PROFILE_ALIASES.get(str(profile), str(profile))
    if canonical not in PROFILE_DEFAULTS:
        raise ValueError(
            f"Unknown ATAC preprocessing profile: {profile}; "
            f"choose one of {', '.join(PUBLIC_PROFILES)}"
        )
    return canonical


def _profile_argument(value: str) -> str:
    try:
        return normalize_profile(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


@dataclass(frozen=True)
class RunInput:
    accession: str
    sample: str
    condition: str
    replicate: str = "1"
    fastq_1: str = ""
    fastq_2: str = ""
    md5_1: str = ""
    md5_2: str = ""


@dataclass
class SampleInput:
    sample: str
    condition: str
    replicate: str
    runs: list[RunInput] = field(default_factory=list)


@dataclass(frozen=True)
class ReferenceBundle:
    assembly: str
    fasta: Path
    index_prefix: Path
    blacklist: Path | None
    chrom_sizes: Path
    tss: Path | None
    macs_genome_size: str


def _deep_merge(base: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    result = json.loads(json.dumps(base))
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_settings(
    path: str | Path | None = None, overrides: dict[str, Any] | None = None
) -> dict[str, Any]:
    loaded: dict[str, Any] = {}
    if path:
        with Path(path).expanduser().open(encoding="utf-8") as handle:
            loaded = yaml.safe_load(handle) or {}
        if not isinstance(loaded, dict):
            raise ValueError("ATAC preprocessing config must be a YAML mapping")
        unknown = set(loaded) - set(DEFAULTS)
        if unknown:
            raise ValueError(
                f"Unknown preprocessing config section(s): {', '.join(sorted(unknown))}"
            )
    profile = normalize_profile(
        str(
            (overrides or {}).get("profile")
            or loaded.get("profile")
            or DEFAULTS["profile"]
        )
    )
    settings = _deep_merge(DEFAULTS, PROFILE_DEFAULTS[profile])
    settings = _deep_merge(settings, loaded)
    if overrides:
        settings = _deep_merge(settings, overrides)
    settings["profile"] = profile
    resources = settings["resources"]
    if "legacy_sample_memory_gb" in resources:
        resources["sample_memory_gb"] = resources["legacy_sample_memory_gb"]
        resources.pop("legacy_sample_memory_gb", None)
    provider = str(settings["download"].get("provider", "auto"))
    if provider not in {"auto", "ena", "sra"}:
        raise ValueError("download.provider must be auto, ena, or sra")
    if int(settings["resources"].get("cores") or 0) < 1:
        raise ValueError("resources.cores must be at least 1")
    if int(settings["resources"].get("max_parallel_samples") or 0) < 1:
        raise ValueError("resources.max_parallel_samples must be at least 1")
    for section in ("trim", "align", "filter", "peaks", "tracks"):
        if not isinstance(settings[section].get("extra_args", []), list):
            raise ValueError(f"{section}.extra_args must be a YAML list")
    return settings


def write_default_config(path: str | Path, profile: str = "modern") -> Path:
    output = Path(path).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    profile = normalize_profile(profile)
    settings = _deep_merge(DEFAULTS, PROFILE_DEFAULTS[profile])
    settings["profile"] = profile
    output.write_text(yaml.safe_dump(settings, sort_keys=False), encoding="utf-8")
    return output


def _first(row: dict[str, str], names: Iterable[str], default: str = "") -> str:
    lower = {str(key).strip().lower(): key for key in row}
    for name in names:
        key = lower.get(name.lower())
        if key is not None and str(row.get(key) or "").strip():
            return str(row[key]).strip()
    return default


def _metadata_reader(path: Path) -> csv.DictReader:
    handle = path.open(encoding="utf-8", newline="")
    sample = handle.read(8192)
    handle.seek(0)
    try:
        delimiter = csv.Sniffer().sniff(sample, delimiters="\t,").delimiter
    except csv.Error:
        delimiter = "\t" if "\t" in sample.splitlines()[0] else ","
    reader = csv.DictReader(handle, delimiter=delimiter)
    reader._fp_tools_handle = handle  # type: ignore[attr-defined]
    return reader


def read_preprocess_metadata(
    path: str | Path,
    id_column: str | None = None,
    sample_column: str | None = None,
    condition_column: str | None = None,
) -> list[SampleInput]:
    metadata = Path(path).expanduser()
    reader = _metadata_reader(metadata)
    try:
        if not reader.fieldnames:
            raise ValueError(f"{metadata} is missing a header row")
        runs: list[RunInput] = []
        for line_no, raw in enumerate(reader, start=2):
            row = {str(key): str(value or "").strip() for key, value in raw.items()}
            accession = _first(
                row,
                (id_column,)
                if id_column
                else ("id", "run_accession", "accession", "run"),
            )
            fastq_1 = _first(row, ("fastq_1", "fastq1", "r1", "read1", "fq1"))
            fastq_2 = _first(row, ("fastq_2", "fastq2", "r2", "read2", "fq2"))
            if not accession and not fastq_1:
                raise ValueError(
                    f"{metadata}:{line_no} requires an accession/ID or fastq_1"
                )
            if accession and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", accession):
                raise ValueError(
                    f"{metadata}:{line_no} accession is not filesystem-safe: {accession!r}"
                )
            sample = _first(
                row,
                (sample_column,)
                if sample_column
                else ("sample", "sample_name", "sample_id", "name"),
                accession
                or Path(urllib.parse.urlparse(fastq_1).path).name.split(".")[0],
            )
            condition = _first(
                row,
                (condition_column,)
                if condition_column
                else ("condition", "condition_label", "group"),
                sample,
            )
            replicate = _first(row, ("replicate", "rep", "biological_replicate"), "1")
            if not sample:
                raise ValueError(f"{metadata}:{line_no} has an empty sample name")
            if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
                raise ValueError(
                    f"{metadata}:{line_no} sample name is not filesystem-safe: {sample!r}"
                )
            if any(value in condition for value in ("\t", "\n", "\r")):
                raise ValueError(
                    f"{metadata}:{line_no} condition contains a tab or newline"
                )
            runs.append(
                RunInput(
                    accession=accession,
                    sample=sample,
                    condition=condition,
                    replicate=replicate,
                    fastq_1=fastq_1,
                    fastq_2=fastq_2,
                    md5_1=_first(row, ("fastq_1_md5", "md5_1", "fastq_md5_1")),
                    md5_2=_first(row, ("fastq_2_md5", "md5_2", "fastq_md5_2")),
                )
            )
    finally:
        reader._fp_tools_handle.close()  # type: ignore[attr-defined]
    if not runs:
        raise ValueError(f"{metadata} does not contain any sample rows")
    grouped: dict[str, SampleInput] = {}
    seen_accessions: set[str] = set()
    for run in runs:
        if run.accession and run.accession in seen_accessions:
            raise ValueError(f"Duplicate run accession in {metadata}: {run.accession}")
        seen_accessions.add(run.accession)
        current = grouped.setdefault(
            run.sample, SampleInput(run.sample, run.condition, run.replicate)
        )
        if current.condition != run.condition or current.replicate != run.replicate:
            raise ValueError(
                f"Technical runs for sample {run.sample} disagree on condition or replicate"
            )
        current.runs.append(run)
    return list(grouped.values())


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_ena_fastqs(
    accession: str, timeout: int = 60
) -> list[tuple[str, str, int | None]]:
    fields = "run_accession,fastq_ftp,fastq_md5,fastq_bytes,library_layout"
    query = urllib.parse.urlencode(
        {"accession": accession, "result": "read_run", "fields": fields}
    )
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?{query}"
    with urllib.request.urlopen(url, timeout=timeout) as response:
        text = response.read().decode("utf-8")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    if len(rows) != 1 or not rows[0].get("fastq_ftp"):
        raise RuntimeError(f"ENA did not return downloadable FASTQs for {accession}")
    urls = rows[0]["fastq_ftp"].split(";")
    md5s = (rows[0].get("fastq_md5") or "").split(";")
    sizes = (rows[0].get("fastq_bytes") or "").split(";")
    if len(urls) not in (1, 2):
        raise RuntimeError(
            f"{accession} resolved to {len(urls)} FASTQs; expected one or two"
        )
    result = []
    for idx, remote in enumerate(urls):
        remote = remote if "://" in remote else f"https://{remote}"
        size = int(sizes[idx]) if idx < len(sizes) and sizes[idx].isdigit() else None
        result.append((remote, md5s[idx] if idx < len(md5s) else "", size))
    return result


def _run(
    command: list[str], log: Path | None = None, cwd: Path | None = None, stdout=None
) -> None:
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        with log.open("a", encoding="utf-8") as handle:
            handle.write("$ " + " ".join(command) + "\n")
            result = subprocess.run(
                command,
                cwd=cwd,
                stdout=stdout or handle,
                stderr=handle,
                text=stdout is None,
            )
    else:
        result = subprocess.run(command, cwd=cwd, stdout=stdout)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)


def _run_pipeline(first: list[str], second: list[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as handle:
        handle.write("$ " + " ".join(first) + " | " + " ".join(second) + "\n")
        left = subprocess.Popen(
            first, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=handle
        )
        right = subprocess.Popen(
            second, stdin=left.stdout, stdout=handle, stderr=handle
        )
        if left.stdout:
            left.stdout.close()
        right_rc = right.wait()
        left_rc = left.wait()
    if left_rc:
        raise subprocess.CalledProcessError(left_rc, first)
    if right_rc:
        raise subprocess.CalledProcessError(right_rc, second)


def download_file(
    url: str,
    output: Path,
    expected_md5: str = "",
    retries: int = 5,
    *,
    expected_size: int = 0,
) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    expected_size = max(0, int(expected_size or 0))

    def is_valid(path: Path) -> bool:
        if not path.is_file() or path.stat().st_size == 0:
            return False
        if expected_size and path.stat().st_size != expected_size:
            return False
        if expected_md5 and _md5(path) != expected_md5:
            return False
        return bool(expected_size or expected_md5)

    if is_valid(output):
        return output

    partial = output.with_suffix(output.suffix + ".partial")
    if output.exists():
        if partial.exists():
            partial.unlink()
        output.replace(partial)

    if shutil.which("aria2c"):
        _run(
            [
                "aria2c",
                "--continue=true",
                f"--max-tries={retries}",
                "--file-allocation=none",
                "--dir",
                str(partial.parent),
                "--out",
                partial.name,
                url,
            ]
        )
    elif shutil.which("curl"):
        _run(
            [
                "curl",
                "--fail",
                "--location",
                "--continue-at",
                "-",
                "--retry",
                str(retries),
                "--output",
                str(partial),
                url,
            ]
        )
    else:
        with (
            urllib.request.urlopen(url, timeout=120) as response,
            partial.open("wb") as handle,
        ):
            shutil.copyfileobj(response, handle)

    if expected_size and partial.stat().st_size != expected_size:
        partial.unlink(missing_ok=True)
        raise RuntimeError(
            f"Size mismatch for {output}: expected {expected_size} bytes"
        )
    if expected_md5 and _md5(partial) != expected_md5:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"MD5 mismatch for {output}")
    if not partial.is_file() or partial.stat().st_size == 0:
        partial.unlink(missing_ok=True)
        raise RuntimeError(f"Download produced an empty file for {output}")
    partial.replace(output)
    return output


def _gzip_file(path: Path) -> Path:
    output = path.with_suffix(path.suffix + ".gz")
    if shutil.which("pigz"):
        _run(["pigz", "-f", str(path)])
    else:
        with path.open("rb") as source, gzip.open(output, "wb") as target:
            shutil.copyfileobj(source, target)
        path.unlink()
    return output


def materialize_run_fastqs(
    run: RunInput, fastq_dir: Path, settings: dict[str, Any]
) -> tuple[Path, Path | None]:
    run_id = run.accession or run.sample
    target = fastq_dir / run_id
    target.mkdir(parents=True, exist_ok=True)
    if run.fastq_1:
        values = []
        for idx, (source, checksum) in enumerate(
            ((run.fastq_1, run.md5_1), (run.fastq_2, run.md5_2)), start=1
        ):
            if not source:
                continue
            parsed = urllib.parse.urlparse(source)
            if parsed.scheme in {"http", "https", "ftp"}:
                name = Path(parsed.path).name or f"{run_id}_{idx}.fastq.gz"
                values.append(
                    download_file(
                        source,
                        target / name,
                        checksum,
                        int(settings["download"]["retries"]),
                    )
                )
            else:
                local = Path(source).expanduser().resolve()
                if not local.exists():
                    raise FileNotFoundError(local)
                if checksum and _md5(local) != checksum:
                    raise RuntimeError(f"MD5 mismatch for {local}")
                values.append(local)
        return values[0], values[1] if len(values) > 1 else None
    provider = str(settings["download"].get("provider", "auto"))
    if provider in {"auto", "ena"}:
        try:
            resolved = resolve_ena_fastqs(run.accession)
            files = [
                download_file(
                    url,
                    target / Path(urllib.parse.urlparse(url).path).name,
                    md5,
                    int(settings["download"]["retries"]),
                    expected_size=size,
                )
                for url, md5, size in resolved
            ]
            return files[0], files[1] if len(files) > 1 else None
        except Exception:
            if provider == "ena":
                raise
    for executable in ("prefetch", "vdb-validate", "fasterq-dump"):
        if not shutil.which(executable):
            raise RuntimeError(
                f"ENA resolution failed and NCBI fallback requires {executable}"
            )
    sra_dir = target / run.accession
    _run(["prefetch", run.accession, "--max-size", "u", "-O", str(sra_dir)])
    _run(["vdb-validate", str(sra_dir)])
    _run(
        [
            "fasterq-dump",
            "--split-files",
            "--threads",
            str(settings["resources"]["cores"]),
            "--outdir",
            str(target),
            str(sra_dir),
        ]
    )
    r1 = target / f"{run.accession}_1.fastq"
    r2 = target / f"{run.accession}_2.fastq"
    single = target / f"{run.accession}.fastq"
    if r1.exists():
        return _gzip_file(r1), _gzip_file(r2) if r2.exists() else None
    if single.exists():
        return _gzip_file(single), None
    raise RuntimeError(f"fasterq-dump did not create FASTQs for {run.accession}")


def _gunzip_to(source: Path, target: Path) -> None:
    opener = gzip.open if source.suffix == ".gz" else open
    with opener(source, "rb") as input_handle, target.open("wb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle)


def _download_reference_asset(url: str, output: Path, expected_md5: str = "") -> Path:
    compressed = output.with_suffix(output.suffix + ".gz")
    download_file(url, compressed, expected_md5)
    if not output.exists() or output.stat().st_mtime < compressed.stat().st_mtime:
        tmp = output.with_suffix(output.suffix + ".tmp")
        _gunzip_to(compressed, tmp)
        tmp.replace(output)
    return output


def _download_tss_bed(url: str, output: Path, expected_md5: str = "") -> Path:
    gtf_gz = output.with_suffix(".gtf.gz")
    download_file(url, gtf_gz, expected_md5)
    if output.exists() and output.stat().st_mtime >= gtf_gz.stat().st_mtime:
        return output
    tmp = output.with_suffix(output.suffix + ".tmp")
    seen: set[tuple[str, int, str]] = set()
    with (
        gzip.open(gtf_gz, "rt", encoding="utf-8") as source,
        tmp.open("w", encoding="utf-8") as target,
    ):
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            chrom, start, end, strand = (
                fields[0],
                int(fields[3]) - 1,
                int(fields[4]),
                fields[6],
            )
            tss = start if strand == "+" else end - 1
            key = (chrom, tss, strand)
            if key in seen:
                continue
            seen.add(key)
            target.write(f"{chrom}\t{tss}\t{tss + 1}\tTSS_{len(seen)}\t0\t{strand}\n")
    tmp.replace(output)
    return output


def _index_complete(prefix: Path) -> bool:
    normal = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
    large = [value + "l" for value in normal]
    return all(Path(str(prefix) + suffix).exists() for suffix in normal) or all(
        Path(str(prefix) + suffix).exists() for suffix in large
    )


def prepare_reference(
    genome: str,
    reference_dir: str | Path,
    fasta: str | Path | None = None,
    bowtie2_index: str | Path | None = None,
    blacklist: str | Path | None = None,
    tss: str | Path | None = None,
    macs_genome_size: str | None = None,
    dry_run: bool = False,
    cores: int | None = None,
) -> ReferenceBundle:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", genome):
        raise ValueError(f"Genome label is not filesystem-safe: {genome!r}")
    assembly = genome
    root = Path(reference_dir).expanduser() / genome
    if not dry_run:
        root.mkdir(parents=True, exist_ok=True)
    manifest = REFERENCE_MANIFEST.get(genome, {})
    fasta_path = Path(fasta).expanduser().resolve() if fasta else root / f"{genome}.fa"
    if not fasta and not fasta_path.exists():
        if not manifest:
            raise ValueError("Custom genomes require --fasta")
        if not dry_run:
            _download_reference_asset(
                manifest["fasta_url"], fasta_path, manifest["fasta_md5"]
            )
    blacklist_path = (
        Path(blacklist).expanduser().resolve()
        if blacklist
        else (root / f"{genome}.blacklist.bed" if manifest else None)
    )
    if blacklist_path and not blacklist and not blacklist_path.exists() and not dry_run:
        _download_reference_asset(
            manifest["blacklist_url"], blacklist_path, manifest["blacklist_md5"]
        )
    if (
        not dry_run
        and fasta_path.exists()
        and not Path(str(fasta_path) + ".fai").exists()
    ):
        _run(["samtools", "faidx", str(fasta_path)])
    chrom_sizes = root / f"{genome}.chrom.sizes"
    if not dry_run and not chrom_sizes.exists():
        with (
            Path(str(fasta_path) + ".fai").open(encoding="utf-8") as source,
            chrom_sizes.open("w", encoding="utf-8") as target,
        ):
            for line in source:
                fields = line.split("\t")
                target.write(f"{fields[0]}\t{fields[1]}\n")
    index_prefix = (
        Path(bowtie2_index).expanduser().resolve()
        if bowtie2_index
        else root / "bowtie2" / genome
    )
    if not dry_run and not _index_complete(index_prefix):
        index_prefix.parent.mkdir(parents=True, exist_ok=True)
        _run(
            [
                "bowtie2-build",
                "--threads",
                str(max(1, cores or os.cpu_count() or 1)),
                str(fasta_path),
                str(index_prefix),
            ]
        )
    resolved_macs_size = macs_genome_size or str(manifest.get("macs_genome_size") or "")
    if not resolved_macs_size:
        raise ValueError("Custom genomes require --macs-genome-size")
    tss_path = (
        Path(tss).expanduser().resolve()
        if tss
        else (root / f"{genome}.tss.bed" if manifest else None)
    )
    if tss_path and not tss and not tss_path.exists() and not dry_run:
        _download_tss_bed(manifest["tss_gtf_url"], tss_path, manifest["tss_gtf_md5"])
    return ReferenceBundle(
        assembly=assembly,
        fasta=fasta_path,
        index_prefix=index_prefix,
        blacklist=blacklist_path,
        chrom_sizes=chrom_sizes,
        tss=tss_path,
        macs_genome_size=resolved_macs_size,
    )


def dependency_report(profile: str = "modern") -> list[dict[str, str]]:
    profile = normalize_profile(profile)
    required = PROFILE_EXECUTABLES[profile]
    rows = []
    optional = (
        "multiqc",
        "aria2c",
        "curl",
        "prefetch",
        "vdb-validate",
        "fasterq-dump",
        "pigz",
    )
    for name in (*required, *optional):
        rows.append(
            {
                "tool": name,
                "path": shutil.which(name) or "",
                "required": "yes" if name in required else "optional",
            }
        )
    return rows


def _tool_version(path: str) -> str:
    try:
        result = subprocess.run(
            [path, "--version"],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=10,
        )
        text = (result.stdout or result.stderr).strip().splitlines()
        return text[0] if text else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


@lru_cache(maxsize=2)
def _software_identity(profile: str = "modern") -> list[dict[str, str | int]]:
    rows: list[dict[str, str | int]] = []
    for row in dependency_report(profile):
        path = row["path"]
        if path:
            stat = Path(path).stat()
            rows.append(
                {
                    "tool": row["tool"],
                    "path": path,
                    "mtime_ns": stat.st_mtime_ns,
                    "version": _tool_version(path),
                }
            )
    return rows


def doctor(profile: str = "modern") -> int:
    profile = normalize_profile(profile)
    rows = dependency_report(profile)
    print(f"profile\t{profile}")
    print("tool\trequired\tpath")
    for row in rows:
        print(f"{row['tool']}\t{row['required']}\t{row['path'] or 'MISSING'}")
    return 1 if any(row["required"] == "yes" and not row["path"] for row in rows) else 0


def _fingerprint(
    sample: SampleInput,
    reference: ReferenceBundle,
    settings: dict[str, Any],
    fastqs: list[tuple[Path, Path | None]],
) -> str:
    reference_paths = [
        reference.fasta,
        reference.chrom_sizes,
        reference.blacklist,
        reference.tss,
    ]
    reference_paths.extend(
        sorted(
            reference.index_prefix.parent.glob(reference.index_prefix.name + "*.bt2*")
        )
    )
    payload = {
        "sample": asdict(sample),
        "reference": asdict(reference),
        "settings": settings,
        "fastqs": [
            [
                str(r1),
                r1.stat().st_size,
                r1.stat().st_mtime_ns,
                str(r2 or ""),
                r2.stat().st_size if r2 else 0,
                r2.stat().st_mtime_ns if r2 else 0,
            ]
            for r1, r2 in fastqs
        ],
        "software": _software_identity(str(settings.get("profile", "modern"))),
        "pipeline_source": [
            [str(path), path.stat().st_mtime_ns, path.stat().st_size]
            for path in (
                Path(__file__),
                Path(__file__).with_name("prepare_atac_legacy.py"),
            )
            if path.exists()
        ],
        "reference_files": [
            [str(path), path.stat().st_size, path.stat().st_mtime_ns]
            for path in reference_paths
            if path is not None and path.exists()
        ],
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, default=str).encode()
    ).hexdigest()


def _read_chrom_sizes(path: Path) -> list[tuple[str, int]]:
    with path.open(encoding="utf-8") as handle:
        return [
            (line.split()[0], int(line.split()[1])) for line in handle if line.strip()
        ]


def _bedgraph_to_bigwig(bedgraph: Path, chrom_sizes: Path, output: Path) -> None:
    header = _read_chrom_sizes(chrom_sizes)
    order = {chrom: idx for idx, (chrom, _) in enumerate(header)}
    bw = pyBigWig.open(str(output), "w")
    try:
        bw.addHeader(header)
        entries: list[tuple[str, int, int, float]] = []
        last_key = (-1, -1)
        with bedgraph.open(encoding="utf-8") as handle:
            for line in handle:
                chrom, start, end, value = line.rstrip().split("\t")[:4]
                if chrom not in order:
                    continue
                entry = (chrom, int(start), int(end), float(value))
                key = (order[chrom], entry[1])
                if key < last_key:
                    raise ValueError(
                        f"bedGraph is not sorted by reference order: {bedgraph}"
                    )
                last_key = key
                entries.append(entry)
                if len(entries) >= 100_000:
                    bw.addEntries(
                        [x[0] for x in entries],
                        [x[1] for x in entries],
                        ends=[x[2] for x in entries],
                        values=[x[3] for x in entries],
                    )
                    entries.clear()
        if entries:
            bw.addEntries(
                [x[0] for x in entries],
                [x[1] for x in entries],
                ends=[x[2] for x in entries],
                values=[x[3] for x in entries],
            )
    finally:
        bw.close()


def _bam_count(path: Path) -> int:
    with pysam.AlignmentFile(path, "rb") as bam:
        return sum(1 for _ in bam.fetch(until_eof=True))


def _samtools_sort_memory(settings: dict[str, Any]) -> str:
    """Divide half of the run memory budget across samtools sort threads."""
    cores = max(1, int(settings["resources"].get("cores") or 1))
    configured = settings["resources"].get("memory_gb")
    if configured is None:
        configured = max(2.0, min(24.0, _available_memory_gb() - 8.0))
    per_thread_mb = int(float(configured) * 1024 * 0.5 / cores)
    return f"{max(64, min(768, per_thread_mb))}M"


def _write_chromosome_subset(
    source: Path, output: Path, chromosomes: list[str], cores: str, log: Path
) -> None:
    """Index a BAM before asking samtools to extract named chromosomes."""
    _run(["samtools", "index", "-@", cores, str(source)], log)
    _run(
        [
            "samtools",
            "view",
            "-@",
            cores,
            "-b",
            "-o",
            str(output),
            str(source),
            *chromosomes,
        ],
        log,
    )


def _fragment_metrics(
    bam_path: Path, output: Path, max_length: int = 1000
) -> dict[str, float | int | None]:
    counts = [0] * (max_length + 1)
    total = 0
    nfr = 0
    mono = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_duplicate
            ):
                continue
            if read.is_paired:
                if not read.is_read1 or read.template_length <= 0:
                    continue
                length = read.template_length
            else:
                length = read.query_alignment_length or 0
            if length <= 0:
                continue
            total += 1
            nfr += int(length < 100)
            mono += int(180 <= length <= 247)
            counts[min(length, max_length)] += 1
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        handle.write("fragment_length\tcount\n")
        for length, count in enumerate(counts):
            if count:
                handle.write(f"{length}\t{count}\n")
    median = None
    if total:
        midpoint = (total + 1) // 2
        cumulative = 0
        for length, count in enumerate(counts):
            cumulative += count
            if cumulative >= midpoint:
                median = length
                break
    return {
        "fragment_count": total,
        "median_fragment_length": median,
        "nucleosome_free_fraction": round(nfr / total, 6) if total else None,
        "mononucleosome_fraction": round(mono / total, 6) if total else None,
    }


def _frip(bam: Path, peaks: Path, log: Path) -> float:
    total = _bam_count(bam)
    if total == 0:
        return 0.0
    first = subprocess.Popen(
        ["bedtools", "intersect", "-u", "-abam", str(bam), "-b", str(peaks)],
        stdout=subprocess.PIPE,
        stderr=log.open("a"),
    )
    second = subprocess.run(
        ["samtools", "view", "-c", "-"],
        stdin=first.stdout,
        capture_output=True,
        text=True,
    )
    if first.stdout:
        first.stdout.close()
    rc = first.wait()
    if rc or second.returncode:
        raise RuntimeError("Failed to calculate FRiP")
    return int(second.stdout.strip() or 0) / total


def _tss_enrichment(
    bam_path: Path, tss_path: Path | None, flank: int = 2000
) -> float | None:
    if tss_path is None or not tss_path.exists():
        return None
    profile = [0.0] * (flank * 2 + 1)
    used = 0
    with (
        pysam.AlignmentFile(bam_path, "rb") as bam,
        tss_path.open(encoding="utf-8") as handle,
    ):
        references = set(bam.references)
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3 or fields[0] not in references:
                continue
            chrom = fields[0]
            start, end = int(fields[1]), int(fields[2])
            strand = fields[5] if len(fields) > 5 else "+"
            center = start if strand != "-" else max(start, end - 1)
            lo, hi = max(0, center - flank), center + flank + 1
            for read in bam.fetch(chrom, lo, hi):
                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read.is_duplicate
                ):
                    continue
                cut = (
                    (read.reference_end - 1 - 5)
                    if read.is_reverse
                    else (read.reference_start + 4)
                )
                offset = cut - center
                if strand == "-":
                    offset = -offset
                if -flank <= offset <= flank:
                    profile[offset + flank] += 1
            used += 1
    if used == 0:
        return None
    edge = profile[:100] + profile[-100:]
    background = sum(edge) / max(1, len(edge))
    if background <= 0:
        return None
    center_values = profile[flank - 50 : flank + 51]
    return round(max(center_values) / background, 4)


def _relative_link(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.unlink(missing_ok=True)
    try:
        destination.symlink_to(os.path.relpath(source, destination.parent))
    except OSError:
        if source.stat().st_dev == destination.parent.stat().st_dev:
            os.link(source, destination)
        else:
            shutil.copy2(source, destination)


def process_sample(
    sample: SampleInput,
    root: Path,
    reference: ReferenceBundle,
    settings: dict[str, Any],
    resume: bool = True,
) -> dict[str, str]:
    if settings.get("profile") == "homer-atac":
        from fp_tools.tools.prepare_atac_legacy import process_legacy_sample

        return process_legacy_sample(sample, root, reference, settings, resume)
    fastq_dir = root / "fastq"
    sample_root = root / "samples" / sample.sample
    alignment = sample_root / "alignment"
    peaks_dir = sample_root / "peaks"
    tracks = sample_root / "tracks"
    qc = sample_root / "qc"
    work = sample_root / ".work"
    for directory in (alignment, peaks_dir, tracks, qc, work):
        directory.mkdir(parents=True, exist_ok=True)
    fastqs = [materialize_run_fastqs(run, fastq_dir, settings) for run in sample.runs]
    paired = all(r2 is not None for _, r2 in fastqs)
    if not paired and any(r2 is not None for _, r2 in fastqs):
        raise ValueError(f"Sample {sample.sample} mixes paired- and single-end runs")
    final_bam = alignment / f"{sample.sample}.filtered.bam"
    final_bai = Path(str(final_bam) + ".bai")
    peak_bed = peaks_dir / f"{sample.sample}.narrowPeak"
    track_bw = tracks / f"{sample.sample}.rp10m.bw"
    state_path = sample_root / "state.json"
    fingerprint = _fingerprint(sample, reference, settings, fastqs)
    if (
        resume
        and state_path.exists()
        and all(path.exists() for path in (final_bam, final_bai, peak_bed, track_bw))
    ):
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if (
            state.get("fingerprint") == fingerprint
            and state.get("status") == "complete"
        ):
            return {
                "sample": sample.sample,
                "condition": sample.condition,
                "bam": str(final_bam),
                "peaks": str(peak_bed),
                "bigwig": str(track_bw),
                "status": "cached",
            }
    cores = str(max(1, int(settings["resources"]["cores"])))
    sort_memory = _samtools_sort_memory(settings)
    log = qc / "commands.log"
    trimmed: list[tuple[Path, Path | None]] = []
    for idx, (r1, r2) in enumerate(fastqs, start=1):
        prefix = f"{sample.sample}.{idx}"
        if settings["qc"].get("fastqc") and shutil.which("fastqc"):
            _run(
                [
                    "fastqc",
                    "--threads",
                    cores,
                    "--outdir",
                    str(qc),
                    str(r1),
                    *([str(r2)] if r2 else []),
                ],
                log,
            )
        if settings["trim"].get("enabled"):
            out1 = work / f"{prefix}.R1.trim.fastq.gz"
            out2 = work / f"{prefix}.R2.trim.fastq.gz" if r2 else None
            command = [
                "fastp",
                "--in1",
                str(r1),
                "--out1",
                str(out1),
                "--thread",
                cores,
                "--length_required",
                str(settings["trim"]["min_length"]),
                "--html",
                str(qc / f"{prefix}.fastp.html"),
                "--json",
                str(qc / f"{prefix}.fastp.json"),
            ]
            if r2:
                command.extend(
                    ["--in2", str(r2), "--out2", str(out2), "--detect_adapter_for_pe"]
                )
            command.extend(str(x) for x in settings["trim"].get("extra_args", []))
            _run(command, log)
            trimmed.append((out1, out2))
        else:
            trimmed.append((r1, r2))
    command = [
        "bowtie2",
        "-x",
        str(reference.index_prefix),
        "-p",
        cores,
        "-X",
        str(settings["align"]["max_insert"]),
        *[str(x) for x in settings["align"].get("extra_args", [])],
    ]
    if paired:
        command.extend(
            [
                "-1",
                ",".join(str(x[0]) for x in trimmed),
                "-2",
                ",".join(str(x[1]) for x in trimmed),
            ]
        )
    else:
        command.extend(["-U", ",".join(str(x[0]) for x in trimmed)])
    initial = work / "initial.bam"
    flags = "2828" if paired else "2820"
    view = ["samtools", "view", "-@", cores, "-b"]
    if paired:
        view.extend(["-f", "2"])
    view.extend(
        [
            "-F",
            flags,
            "-q",
            str(settings["align"]["mapq"]),
            *[str(x) for x in settings["filter"].get("extra_args", [])],
            "-o",
            str(initial),
            "-",
        ]
    )
    _run_pipeline(command, view, log)
    name_sorted = work / "name_sorted.bam"
    fixmate = work / "fixmate.bam"
    coord = work / "coord.bam"
    dedup = work / "dedup.bam"
    _run(
        [
            "samtools",
            "sort",
            "-@",
            cores,
            "-m",
            sort_memory,
            "-n",
            "-o",
            str(name_sorted),
            str(initial),
        ],
        log,
    )
    _run(
        ["samtools", "fixmate", "-@", cores, "-m", str(name_sorted), str(fixmate)], log
    )
    _run(
        [
            "samtools",
            "sort",
            "-@",
            cores,
            "-m",
            sort_memory,
            "-o",
            str(coord),
            str(fixmate),
        ],
        log,
    )
    markdup = ["samtools", "markdup", "-@", cores]
    if settings["filter"].get("remove_duplicates"):
        markdup.append("-r")
    markdup.append("-s")
    markdup.extend([str(coord), str(dedup)])
    _run(markdup, log)
    keep_chroms = [
        chrom
        for chrom, _ in _read_chrom_sizes(reference.chrom_sizes)
        if not (settings["filter"].get("remove_mito") and chrom in MITO_CHROMS)
    ]
    nomito = work / "nomito.bam"
    _write_chromosome_subset(dedup, nomito, keep_chroms, cores, log)
    dedup_reads = _bam_count(dedup)
    nomito_reads = _bam_count(nomito)
    if reference.blacklist and reference.blacklist.exists():
        unfiltered = work / "blacklist_filtered.unsorted.bam"
        with unfiltered.open("wb") as output:
            _run(
                [
                    "bedtools",
                    "intersect",
                    "-v",
                    "-abam",
                    str(nomito),
                    "-b",
                    str(reference.blacklist),
                ],
                log,
                stdout=output,
            )
        _run(
            [
                "samtools",
                "sort",
                "-@",
                cores,
                "-m",
                sort_memory,
                "-o",
                str(final_bam),
                str(unfiltered),
            ],
            log,
        )
    else:
        shutil.copy2(nomito, final_bam)
    _run(["samtools", "index", "-@", cores, str(final_bam)], log)
    with (qc / "flagstat.tsv").open("w", encoding="utf-8") as output:
        _run(
            ["samtools", "flagstat", "-@", cores, "-O", "tsv", str(final_bam)],
            log,
            stdout=output,
        )
    macs_name = sample.sample
    command = [
        "macs3",
        "callpeak",
        "-t",
        str(final_bam),
        "-n",
        macs_name,
        "--outdir",
        str(peaks_dir),
        "-g",
        reference.macs_genome_size,
        "--keep-dup",
        "all",
        "-q",
        str(settings["peaks"]["qvalue"]),
    ]
    if paired:
        command.extend(["-f", "BAMPE"])
    else:
        command.extend(
            [
                "-f",
                "BAM",
                "--nomodel",
                "--shift",
                str(settings["peaks"]["single_end_shift"]),
                "--extsize",
                str(settings["peaks"]["single_end_extsize"]),
            ]
        )
    command.extend(str(x) for x in settings["peaks"].get("extra_args", []))
    _run(command, log)
    macs_peak = peaks_dir / f"{macs_name}_peaks.narrowPeak"
    if not macs_peak.exists():
        raise RuntimeError(f"MACS3 did not create {macs_peak}")
    shutil.copy2(macs_peak, peak_bed)
    filtered_peak = peaks_dir / f"{sample.sample}.narrowPeak.filtered.bed"
    with (
        peak_bed.open(encoding="utf-8") as source,
        filtered_peak.open("w", encoding="utf-8") as target,
    ):
        for line in source:
            if line.split("\t", 1)[0] not in MITO_CHROMS:
                target.write(line)
    usable_reads = max(1, _bam_count(final_bam) // (2 if paired else 1))
    scale = int(settings["tracks"]["normalization_reads"]) / usable_reads
    bedgraph = work / "coverage.rp10m.bedgraph"
    command = [
        "bedtools",
        "genomecov",
        "-ibam",
        str(final_bam),
        "-bg",
        "-scale",
        str(scale),
    ]
    if paired:
        command.append("-pc")
    command.extend(str(x) for x in settings["tracks"].get("extra_args", []))
    with bedgraph.open("w", encoding="utf-8") as output:
        _run(command, log, stdout=output)
    _bedgraph_to_bigwig(bedgraph, reference.chrom_sizes, track_bw)
    metrics = {
        "sample": sample.sample,
        "condition": sample.condition,
        "paired_end": paired,
        "usable_reads": usable_reads,
        "peaks": sum(1 for line in peak_bed.open(encoding="utf-8") if line.strip()),
        "frip": round(_frip(final_bam, peak_bed, log), 6),
        "tss_enrichment": _tss_enrichment(final_bam, reference.tss),
        "mitochondrial_fraction": round((dedup_reads - nomito_reads) / dedup_reads, 6)
        if dedup_reads
        else None,
    }
    metrics.update(_fragment_metrics(final_bam, qc / "fragment_lengths.tsv"))
    warnings = []
    if metrics["frip"] < float(settings["qc"]["warn_frip_below"]):
        warnings.append(
            f"FRiP {metrics['frip']:.3f} is below {settings['qc']['warn_frip_below']}"
        )
    tss_threshold = settings["qc"].get("warn_tss_below", {}).get(reference.assembly)
    if (
        metrics["tss_enrichment"] is not None
        and tss_threshold is not None
        and metrics["tss_enrichment"] < float(tss_threshold)
    ):
        warnings.append(
            f"TSS enrichment {metrics['tss_enrichment']:.3f} is below {tss_threshold}"
        )
    metrics["warnings"] = warnings
    (qc / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    if warnings and settings["qc"].get("fail_on_warning"):
        raise RuntimeError(
            f"{sample.sample} failed configured QC thresholds: {'; '.join(warnings)}"
        )
    legacy = root / sample.sample
    suffix = reference.assembly
    _relative_link(final_bam, legacy / f"{sample.sample}.{suffix}.filtered.bam")
    _relative_link(final_bai, legacy / f"{sample.sample}.{suffix}.filtered.bam.bai")
    _relative_link(track_bw, legacy / f"{sample.sample}.{suffix}.rp10m.bw")
    _relative_link(peak_bed, legacy / f"{sample.sample}.{suffix}.rp10m.narrowpeaks.bed")
    _relative_link(
        filtered_peak,
        legacy / f"{sample.sample}.{suffix}.rp10m.narrowpeaks.filtered.bed",
    )
    state_path.write_text(
        json.dumps(
            {"status": "complete", "fingerprint": fingerprint, "metrics": metrics},
            indent=2,
        ),
        encoding="utf-8",
    )
    if not settings.get("cleanup", {}).get("keep_intermediates"):
        shutil.rmtree(work, ignore_errors=True)
    if not settings["download"].get("keep_fastq"):
        for run in sample.runs:
            shutil.rmtree(fastq_dir / (run.accession or run.sample), ignore_errors=True)
    return {
        "sample": sample.sample,
        "condition": sample.condition,
        "bam": str(final_bam),
        "peaks": str(peak_bed),
        "bigwig": str(track_bw),
        "status": "complete",
    }


def _merge_peaks(results: list[dict[str, str]], root: Path) -> Path:
    peak_dir = root / "peaks"
    peak_dir.mkdir(parents=True, exist_ok=True)
    combined = peak_dir / "all_sample_peaks.bed"
    with combined.open("w", encoding="utf-8") as target:
        for result in results:
            with Path(result["peaks"]).open(encoding="utf-8") as source:
                for line in source:
                    fields = line.rstrip().split("\t")
                    if (
                        len(fields) >= 3
                        and fields[1].isdigit()
                        and fields[2].isdigit()
                        and fields[0] not in MITO_CHROMS
                    ):
                        target.write("\t".join(fields[:3]) + "\n")
    sorted_bed = peak_dir / "all_sample_peaks.sorted.bed"
    merged = peak_dir / "merged_peaks.bed"
    with sorted_bed.open("w", encoding="utf-8") as output:
        _run(["bedtools", "sort", "-i", str(combined)], stdout=output)
    with merged.open("w", encoding="utf-8") as output:
        _run(["bedtools", "merge", "-i", str(sorted_bed)], stdout=output)
    shutil.copy2(merged, peak_dir / "merged_peaks_filtered.bed")
    return merged


def _write_project_metadata(
    results: list[dict[str, str]],
    samples_path: Path,
    root: Path,
    reference: ReferenceBundle,
    settings: dict[str, Any],
    samples: list[SampleInput] | None = None,
) -> None:
    metadata = root / "metadata"
    metadata.mkdir(parents=True, exist_ok=True)
    (root / "reports").mkdir(parents=True, exist_ok=True)
    shutil.copy2(samples_path, metadata / f"input_samples{samples_path.suffix}")
    if samples is not None:
        with (metadata / "resolved_runs.tsv").open("w", encoding="utf-8") as handle:
            handle.write("accession\tsample\tcondition\treplicate\tfastq_1\tfastq_2\n")
            for sample in samples:
                for run in sample.runs:
                    handle.write(
                        f"{run.accession}\t{run.sample}\t{run.condition}\t{run.replicate}\t{run.fastq_1}\t{run.fastq_2}\n"
                    )
    with (metadata / "samples.tsv").open("w", encoding="utf-8") as handle:
        handle.write("sample\tcondition\tbam\tpeaks\tbigwig\n")
        for row in results:
            handle.write(
                f"{row['sample']}\t{row['condition']}\t{row['bam']}\t{row['peaks']}\t{row['bigwig']}\n"
            )
    (metadata / "resolved_settings.yml").write_text(
        yaml.safe_dump(settings, sort_keys=False), encoding="utf-8"
    )
    (metadata / "reference.json").write_text(
        json.dumps(asdict(reference), indent=2, default=str), encoding="utf-8"
    )
    with (metadata / "software_versions.tsv").open("w", encoding="utf-8") as handle:
        handle.write("tool\tversion\tpath\n")
        for row in _software_identity(str(settings.get("profile", "modern"))):
            handle.write(f"{row['tool']}\t{row['version']}\t{row['path']}\n")
    with (root / "reports" / "qc_summary.tsv").open("w", encoding="utf-8") as handle:
        handle.write(
            "sample\tcondition\tusable_reads\tpeaks\tfrip\ttss_enrichment\tmitochondrial_fraction\tmedian_fragment_length\tnucleosome_free_fraction\tmononucleosome_fraction\twarnings\n"
        )
        for row in results:
            metrics_path = root / "samples" / row["sample"] / "qc" / "metrics.json"
            metrics = (
                json.loads(metrics_path.read_text(encoding="utf-8"))
                if metrics_path.exists()
                else {}
            )
            handle.write(
                f"{row['sample']}\t{row['condition']}\t{metrics.get('usable_reads', '')}\t{metrics.get('peaks', '')}\t"
                f"{metrics.get('frip', '')}\t{metrics.get('tss_enrichment') or ''}\t{metrics.get('mitochondrial_fraction') or ''}\t"
                f"{metrics.get('median_fragment_length') or ''}\t{metrics.get('nucleosome_free_fraction') or ''}\t"
                f"{metrics.get('mononucleosome_fraction') or ''}\t{'; '.join(metrics.get('warnings', []))}\n"
            )
    downstream = {
        "version": 1,
        "run_mode": "single",
        "samples": [
            {
                "sample_id": "atac_correct",
                "tool": "atac-correct",
                "sample_table": str(metadata / "samples.tsv"),
                "genome": str(reference.fasta),
                "peaks": str(root / "peaks" / "merged_peaks.bed"),
                "blacklist": str(reference.blacklist or ""),
                "outdir": str(root),
                "layout": "project",
            }
        ],
        "comparisons": [],
    }
    (metadata / "atac_correct.yml").write_text(
        yaml.safe_dump(downstream, sort_keys=False), encoding="utf-8"
    )


def _available_memory_gb() -> float:
    """Return currently available physical memory without another dependency."""
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / (1024**2)
    except (OSError, ValueError, IndexError):
        pass
    try:
        return os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / (1024**3)
    except (AttributeError, OSError, ValueError):
        return 8.0


def run_preprocessing(args: argparse.Namespace) -> int:
    if os.name == "nt":
        raise SystemExit(
            "prepare-atac requires Unix bioinformatics executables. On Windows, run it in WSL; "
            "the downstream fp-tools commands and GUI run natively."
        )
    overrides: dict[str, Any] = {"profile": args.profile} if args.profile else {}
    resource_overrides = {
        key: value
        for key, value in {
            "cores": args.cores,
            "max_parallel_samples": args.max_parallel_samples,
            "memory_gb": args.memory_gb,
        }.items()
        if value is not None
    }
    if resource_overrides:
        overrides["resources"] = resource_overrides
    if args.keep_intermediates:
        overrides["cleanup"] = {"keep_intermediates": True}
    settings = load_settings(args.config, overrides or None)
    samples = read_preprocess_metadata(
        args.samples, args.id_column, args.sample_column, args.condition_column
    )
    if args.include:
        wanted = set(args.include)
        samples = [
            sample
            for sample in samples
            if sample.sample in wanted
            or any(run.accession in wanted for run in sample.runs)
        ]
        missing = (
            wanted
            - {sample.sample for sample in samples}
            - {run.accession for sample in samples for run in sample.runs}
        )
        if missing:
            raise ValueError(
                f"Unknown requested sample/accession(s): {', '.join(sorted(missing))}"
            )
    root = Path(args.outdir).expanduser().resolve()
    reference_dir = (
        args.reference_dir or Path.home() / ".cache" / "fp-tools" / "references"
    )
    if args.dry_run:
        prepare_reference(
            args.genome,
            reference_dir,
            args.fasta,
            args.bowtie2_index,
            args.blacklist,
            args.tss,
            args.macs_genome_size,
            True,
            int(settings["resources"]["cores"]),
        )
        print(
            f"prepare-atac: {len(samples)} sample(s), genome={args.genome}, outdir={root}"
        )
        for sample in samples:
            print(
                f"{sample.sample}\t{sample.condition}\t{','.join(run.accession or run.fastq_1 for run in sample.runs)}"
            )
        return 0
    profile = str(settings.get("profile", "modern"))
    missing = [
        row["tool"]
        for row in dependency_report(profile)
        if row["required"] == "yes" and not row["path"]
    ]
    if missing:
        raise RuntimeError(
            f"Missing required preprocessing tools: {', '.join(missing)}; run prepare-atac --doctor"
        )
    root.mkdir(parents=True, exist_ok=True)
    reference = prepare_reference(
        args.genome,
        reference_dir,
        args.fasta,
        args.bowtie2_index,
        args.blacklist,
        args.tss,
        args.macs_genome_size,
        False,
        int(settings["resources"]["cores"]),
    )
    results = []
    failures = []
    parallel = min(
        len(samples), int(settings["resources"].get("max_parallel_samples") or 1)
    )
    memory_gb = settings["resources"].get("memory_gb")
    if memory_gb:
        available_gb = _available_memory_gb()
        required_gb = float(memory_gb) + 8.0
        if available_gb < required_gb:
            raise RuntimeError(
                f"Only {available_gb:.1f} GiB RAM is available; this run requires {required_gb:.1f} GiB including the host reserve"
            )
        if profile == "homer-atac":
            per_sample = float(
                settings["resources"].get("sample_memory_gb") or 16
            )
            if float(memory_gb) < per_sample:
                raise RuntimeError(
                    f"homer-atac reserves {per_sample:g} GiB per active sample; increase --memory-gb"
                )
            parallel = min(parallel, max(1, int(float(memory_gb) // per_sample)))
    total_cores = int(settings["resources"]["cores"])
    worker_settings = _deep_merge(
        settings, {"resources": {"cores": max(1, total_cores // max(1, parallel))}}
    )
    if parallel <= 1 or args.fail_fast:
        for sample in samples:
            try:
                results.append(
                    process_sample(
                        sample,
                        root,
                        reference,
                        worker_settings,
                        resume=not args.no_resume,
                    )
                )
            except Exception as exc:
                failures.append((sample.sample, str(exc)))
                if args.fail_fast:
                    raise
    else:
        with ThreadPoolExecutor(max_workers=parallel) as executor:
            future_map = {
                executor.submit(
                    process_sample,
                    sample,
                    root,
                    reference,
                    worker_settings,
                    not args.no_resume,
                ): sample
                for sample in samples
            }
            for future in as_completed(future_map):
                sample = future_map[future]
                try:
                    results.append(future.result())
                except Exception as exc:
                    failures.append((sample.sample, str(exc)))
        order = {sample.sample: idx for idx, sample in enumerate(samples)}
        results.sort(key=lambda row: order[row["sample"]])
    if results:
        _merge_peaks(results, root)
        (root / "reports").mkdir(parents=True, exist_ok=True)
        _write_project_metadata(
            results, Path(args.samples).expanduser(), root, reference, settings, samples
        )
    if settings["qc"].get("multiqc") and shutil.which("multiqc") and results:
        _run(
            [
                "multiqc",
                "--force",
                "--outdir",
                str(root / "reports"),
                str(root / "samples"),
            ]
        )
    if failures:
        print("Failed samples:", file=sys.stderr)
        for sample, message in failures:
            print(f"  {sample}: {message}", file=sys.stderr)
        return 1
    print(f"Prepared {len(results)} ATAC sample(s) in {root}")
    print(f"Downstream sample table: {root / 'metadata' / 'samples.tsv'}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="prepare-atac",
        description="Download single- or paired-end ATAC-seq reads and prepare filtered BAM, peak BED, sequencing-depth-normalized alignment coverage bigWig, and QC files.",
    )
    parser.add_argument(
        "--samples",
        help="TSV/CSV sample metadata with an ID/run_accession or fastq_1 column.",
    )
    parser.add_argument(
        "--genome", help="Named genome (hg38/mm10) or custom genome label."
    )
    parser.add_argument("--outdir", help="Output project directory.")
    parser.add_argument(
        "--config",
        help="Optional preprocessing YAML; CLI values override packaged defaults.",
    )
    parser.add_argument(
        "--profile",
        choices=PUBLIC_PROFILES,
        type=_profile_argument,
        help="Processing method: modern uses fastp, samtools, and MACS3; homer-atac uses Trim Galore, Picard, and HOMER (default: modern).",
    )
    parser.add_argument("--id-column", help="Explicit accession column name.")
    parser.add_argument("--sample-column", help="Explicit sample-name column.")
    parser.add_argument("--condition-column", help="Explicit condition column.")
    parser.add_argument(
        "--include", nargs="+", help="Only process these sample names or accessions."
    )
    parser.add_argument(
        "--reference-dir",
        help="Reference cache root (default: ~/.cache/fp-tools/references).",
    )
    parser.add_argument("--fasta", help="Custom reference FASTA.")
    parser.add_argument("--bowtie2-index", help="Existing Bowtie2 index prefix.")
    parser.add_argument("--blacklist", help="Custom blacklist BED.")
    parser.add_argument("--tss", help="Optional TSS BED for enrichment QC.")
    parser.add_argument(
        "--macs-genome-size",
        help="MACS3 genome size or hs/mm shorthand for custom genomes.",
    )
    parser.add_argument("--cores", type=int, help="Total core budget.")
    parser.add_argument(
        "--max-parallel-samples",
        type=int,
        help="Maximum samples processed concurrently (default: config value).",
    )
    parser.add_argument(
        "--memory-gb",
        type=float,
        help="Enforced total memory budget in GiB; reserves 8 GiB for the host.",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep trimmed FASTQs and intermediate alignment files under each sample .work directory.",
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Recompute completed samples even when fingerprints match.",
    )
    parser.add_argument(
        "--fail-fast", action="store_true", help="Stop after the first failed sample."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and list the planned samples without downloading or processing.",
    )
    parser.add_argument(
        "--doctor",
        action="store_true",
        help="Report external preprocessing dependencies and exit.",
    )
    parser.add_argument(
        "--write-default-config",
        metavar="PATH",
        help="Write the fully documented default YAML and exit.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.doctor:
        raise SystemExit(doctor(args.profile or "modern"))
    if args.write_default_config:
        print(write_default_config(args.write_default_config, args.profile or "modern"))
        return
    missing = [
        flag for flag in ("samples", "genome", "outdir") if not getattr(args, flag)
    ]
    if missing:
        parser.error(
            "the following arguments are required for a run: "
            + ", ".join("--" + value for value in missing)
        )
    try:
        raise SystemExit(run_preprocessing(args))
    except (
        ValueError,
        RuntimeError,
        FileNotFoundError,
        subprocess.CalledProcessError,
    ) as exc:
        parser.exit(2, f"prepare-atac: error: {exc}\n")


if __name__ == "__main__":
    main()
