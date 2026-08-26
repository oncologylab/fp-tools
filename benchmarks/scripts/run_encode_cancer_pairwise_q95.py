#!/usr/bin/env python3
"""Rebuild the seven-line ENCODE browser with the reference Q95 workflow.

Each unordered comparison is an independent, resumable analysis.  Its peak
universe is the union of every released GRCh38 ``IDR thresholded peaks`` file
from the two locked ENCODE experiments.  Corrected cut-site tracks are scaled
together by their peak-bin q95 values before footprint scoring.  This is the
same design used by the preserved K562 versus HepG2 demonstration.

Downloaded BAMs are retained until a separate, checksum-verified Box archive
has completed.  Pair-local generated bigWigs are removed only after the result
table and compact browser payload have both passed validation.
"""

from __future__ import annotations

import argparse
import base64
import concurrent.futures
from datetime import datetime, timezone
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import urllib.request

import pandas as pd

try:
    import fcntl
except ImportError:  # Windows can inspect fixtures but cannot run concurrent jobs.
    fcntl = None


def _lock_file(handle) -> None:
    if fcntl is not None:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814.tsv"
PEAK_MANIFEST = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814_peaks.tsv"
COMPARISONS = ROOT / "benchmarks/manifests/encode_cancer_7line_20260814_comparisons.tsv"
PROJECT = ROOT / "data/public/processed/encode_cancer_pairwise_q95_20260814"
RAW = ROOT / "data/public/raw/encode_cancer_pairwise_q95_20260814"
GENOME = ROOT / "data/public/raw/genome/hg38.fa"
REFERENCE_PROJECT = ROOT / "data/public/processed/encode_k562_hepg2_atac_replicates"
REFERENCE_REPORT = ROOT / "docs/demos/reports/diff_footprints_K562_HepG2.html"
REFERENCE_RESULTS = ROOT / "benchmarks/fixtures/encode_k562_hepg2_q95_diff_footprints_results.txt.gz"
BLACKLIST = REFERENCE_PROJECT / "reference/hg38-blacklist.v2.bed"
CHROM_SIZES = REFERENCE_PROJECT / "reference/hg38.chrom.sizes"
BEDTOOLS = Path("/home/exouser/miniforge3/envs/fp-tools-atac/bin/bedtools")
SAMTOOLS = Path("/home/exouser/miniforge3/envs/fp-tools-atac/bin/samtools")
EXPECTED_CONDITIONS = ("A549", "HCT116", "HepG2", "K562", "MCF-7", "PC-3", "Panc1")
EXPECTED_REPLICATES = {"A549": 3, "HCT116": 2, "HepG2": 3, "K562": 3, "MCF-7": 2, "PC-3": 2, "Panc1": 2}
EXPECTED_MOTIFS = 1019
EXPECTED_PAIRS = 21
EXPECTED_AGGREGATE_SITE_SET = "all"
REFERENCE_JSON_SHA256 = "761181a913f6f538aa47c3af07d005fa34f30f38f986ded0152f2316fc40ad6e"
REPORT_LABEL = "Normalization: Q95-scale; aggregate sites: all motif matches; FDR < 0.001; |Δ FP score| > 0.1; Bound Sites > 500"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.{os.getpid()}.part")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def md5(path: Path) -> str:
    value = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def run(command: list[str], log: Path, *, stdout=None) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("a", encoding="utf-8") as handle:
        handle.write(f"\n[{utc_now()}] $ {' '.join(command)}\n")
        handle.flush()
        result = subprocess.run(
            command,
            cwd=ROOT,
            stdout=stdout if stdout is not None else handle,
            stderr=handle,
            check=False,
        )
    if result.returncode:
        raise RuntimeError(f"Command failed with exit status {result.returncode}; see {log}")


def load_design() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    samples = pd.read_csv(MANIFEST, sep="\t", dtype=str, keep_default_na=False)
    peaks = pd.read_csv(PEAK_MANIFEST, sep="\t", dtype=str, keep_default_na=False)
    comparisons = pd.read_csv(COMPARISONS, sep="\t", dtype=str, keep_default_na=False)
    if len(samples) != 17 or samples["sample"].duplicated().any() or samples["bam_accession"].duplicated().any():
        raise ValueError("The locked sample manifest must contain exactly 17 unique samples and BAMs")
    observed = samples.groupby("condition")["sample"].count().to_dict()
    if observed != EXPECTED_REPLICATES:
        raise ValueError(f"Unexpected replicate design: {observed}")
    for condition, count in EXPECTED_REPLICATES.items():
        subset = samples.loc[samples["condition"].eq(condition)]
        expected_numbers = {str(value) for value in range(1, count + 1)}
        if set(subset["display_order"]) != expected_numbers:
            raise ValueError(f"Unexpected display order for {condition}")
        if set(subset["biological_replicate"]) != expected_numbers:
            raise ValueError(f"Unexpected ENCODE biological replicates for {condition}")
    if tuple(sorted(samples["condition"].unique())) != tuple(sorted(EXPECTED_CONDITIONS)):
        raise ValueError("The sample manifest does not contain the exact seven cancer cell lines")
    if len(peaks) != 30 or peaks["peak_accession"].duplicated().any():
        raise ValueError("The locked peak manifest must contain exactly 30 unique files")
    if set(peaks["output_type"]) != {"IDR thresholded peaks"}:
        raise ValueError("Only ENCODE IDR thresholded peak files are eligible")
    sample_experiments = samples.groupby("condition")["experiment"].first().to_dict()
    if peaks.groupby("condition")["experiment"].first().to_dict() != sample_experiments:
        raise ValueError("Peak and BAM experiments differ")
    expected_pairs = {frozenset(pair) for pair in itertools.combinations(EXPECTED_CONDITIONS, 2)}
    observed_pairs = {frozenset((row.cond1, row.cond2)) for row in comparisons.itertuples(index=False)}
    if len(comparisons) != EXPECTED_PAIRS or observed_pairs != expected_pairs:
        raise ValueError("The comparison manifest must contain all 21 unordered pairs exactly once")
    return samples, peaks, comparisons


def verify_file(path: Path, size: int, checksum: str) -> bool:
    if not path.is_file() or path.stat().st_size != size:
        return False
    stat = path.stat()
    cache_key = hashlib.sha256(str(path.resolve()).encode("utf-8")).hexdigest()
    cache_dir = PROJECT / "verification"
    cache_dir.mkdir(parents=True, exist_ok=True)
    marker = cache_dir / f"{cache_key}.json"
    with (cache_dir / f"{cache_key}.lock").open("a", encoding="utf-8") as lock:
        _lock_file(lock)
        expected = {
            "path": str(path.resolve()),
            "size": size,
            "mtime_ns": stat.st_mtime_ns,
            "md5": checksum,
        }
        if marker.is_file():
            try:
                if json.loads(marker.read_text(encoding="utf-8")) == expected:
                    return True
            except (OSError, ValueError):
                pass
        if md5(path) != checksum:
            return False
        temporary = marker.with_name(f"{marker.name}.{os.getpid()}.part")
        temporary.write_text(json.dumps(expected, sort_keys=True) + "\n", encoding="utf-8")
        temporary.replace(marker)
        return True


def download(accession: str, suffix: str, size: int, checksum: str, destination: Path, log: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    lock_path = destination.with_name(destination.name + ".lock")
    with lock_path.open("a", encoding="utf-8") as lock:
        _lock_file(lock)
        if verify_file(destination, size, checksum):
            return destination
        partial = destination.with_name(destination.name + ".part")
        url = f"https://www.encodeproject.org/files/{accession}/@@download/{accession}.{suffix}"
        command = [
            "curl", "--fail", "--location", "--retry", "12", "--retry-all-errors",
            "--continue-at", "-", "--output", str(partial), url,
        ]
        run(command, log)
        if partial.stat().st_size != size or md5(partial) != checksum:
            raise ValueError(f"Downloaded file failed ENCODE size/MD5 verification: {accession}")
        partial.replace(destination)
    return destination


def existing_bam(accession: str) -> Path | None:
    roots = (
        ROOT / "data/public/raw/encode_k562_hepg2_atac/bam",
        ROOT / "data/public/raw/encode_hct116_reviewer_revision",
        RAW / "bam",
    )
    for root in roots:
        if not root.is_dir():
            continue
        matches = sorted(root.glob(f"*{accession}*.bam"))
        if matches:
            return matches[0]
    return None


def ensure_bam(row, *, allow_download: bool) -> Path:
    path = existing_bam(row.bam_accession)
    if path is None:
        if not allow_download:
            raise ValueError(f"Missing BAM {row.bam_accession}; rerun with --download")
        path = download(
            row.bam_accession,
            "bam",
            int(row.bam_size),
            row.bam_md5,
            RAW / "bam" / f"{row.bam_accession}.bam",
            PROJECT / "logs/download_bams.log",
        )
    if not verify_file(path, int(row.bam_size), row.bam_md5):
        raise ValueError(f"Local BAM does not match the locked ENCODE record: {path}")
    index = Path(str(path) + ".bai")
    if not index.is_file() or not index.stat().st_size:
        run([str(SAMTOOLS), "index", "-@", "4", str(path)], PROJECT / "logs/index_bams.log")
    return path


def ensure_peak(row, *, allow_download: bool) -> Path:
    path = RAW / "peaks" / f"{row.peak_accession}.bed.gz"
    if verify_file(path, int(row.peak_size), row.peak_md5):
        return path
    historical = ROOT / "data/public/raw/encode_k562_hepg2_atac/peaks"
    matches = sorted(historical.glob(f"*{row.peak_accession}*.bed.gz")) if historical.is_dir() else []
    if matches and verify_file(matches[0], int(row.peak_size), row.peak_md5):
        return matches[0]
    if not allow_download:
        raise ValueError(f"Missing peak file {row.peak_accession}; rerun with --download")
    return download(
        row.peak_accession,
        "bed.gz",
        int(row.peak_size),
        row.peak_md5,
        path,
        PROJECT / "logs/download_peaks.log",
    )


def extract_payload(report: Path) -> tuple[dict, str]:
    match = re.search(r'const reportPayloadB64="([^"]+)"', report.read_text(encoding="utf-8"))
    if not match:
        raise ValueError(f"No compressed report payload found in {report}")
    uncompressed = gzip.decompress(base64.b64decode(match.group(1)))
    return json.loads(uncompressed), hashlib.sha256(uncompressed).hexdigest()


def compact_payload(payload: dict, output: Path) -> dict:
    compact = dict(payload)
    compact["logos"] = {}
    raw = json.dumps(compact, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".part")
    with temporary.open("wb") as handle:
        handle.write(gzip.compress(raw, compresslevel=9, mtime=0))
    decoded = json.loads(gzip.decompress(temporary.read_bytes()))
    validate_payload(decoded)
    temporary.replace(output)
    return {
        "payload_sha256": sha256(output),
        "scientific_json_sha256": scientific_digest(decoded),
        "payload_bytes": output.stat().st_size,
    }


def clean_pair_outputs(results: Path, work: Path) -> None:
    """Retain only compact scientific results after a pair passes validation."""
    keep = {"report_payload.json.gz", "diff_footprints_results.txt"}
    if work.is_dir():
        shutil.rmtree(work)
    if not results.is_dir():
        return
    for path in results.iterdir():
        if path.name in keep:
            continue
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def validate_result_table(path: Path) -> None:
    if not path.is_file() or not path.stat().st_size:
        raise ValueError(f"Differential result table is missing: {path}")
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if len(frame) != EXPECTED_MOTIFS or "output_prefix" not in frame:
        raise ValueError(f"Differential result table does not contain {EXPECTED_MOTIFS} motifs: {path}")
    if frame["output_prefix"].nunique() != EXPECTED_MOTIFS:
        raise ValueError(f"Differential result table motif identifiers are not unique: {path}")


def scientific_digest(payload: dict) -> str:
    scientific = {
        key: payload[key]
        for key in ("title", "report_label", "conditions", "groups", "colors", "points", "aggregate", "change_label")
    }
    raw = json.dumps(scientific, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def validate_payload(payload: dict, expected_samples: list[str] | None = None) -> None:
    if len(payload.get("conditions", [])) != 2 or len(payload.get("points", [])) != EXPECTED_MOTIFS:
        raise ValueError("Report payload does not contain two conditions and 1,019 motifs")
    points = payload["points"]
    if len({point["prefix"] for point in points}) != EXPECTED_MOTIFS:
        raise ValueError("Report payload motif prefixes are not unique")
    motifs = payload.get("aggregate", {}).get("motifs", [])
    if not motifs or len({motif["prefix"] for motif in motifs}) != len(motifs):
        raise ValueError("Report payload lacks valid aggregate profiles")
    if payload.get("aggregate", {}).get("site_set") != EXPECTED_AGGREGATE_SITE_SET:
        raise ValueError("Aggregate profiles were not calculated from all motif matches")
    for motif in motifs:
        site_counts = {
            int(condition.get("n_sites", -1))
            for condition in motif.get("conditions", [])
        }
        if len(site_counts) != 1 or next(iter(site_counts), 0) <= 0:
            raise ValueError(
                f"All-site aggregate count differs between conditions: {motif.get('prefix')}"
            )
    observed_samples = {
        row["name"]
        for motif in motifs
        for condition in motif.get("conditions", [])
        for row in condition.get("samples", [])
    }
    if len(observed_samples) != sum(EXPECTED_REPLICATES[name] for name in payload["conditions"]):
        raise ValueError(f"Aggregate sample count differs from the locked design: {sorted(observed_samples)}")
    if expected_samples is not None and observed_samples != set(expected_samples):
        raise ValueError(
            f"Aggregate samples differ from the pair manifest: observed={sorted(observed_samples)}, "
            f"expected={sorted(expected_samples)}"
        )


def report_outputs_valid(report: Path, result_table: Path, expected_samples: list[str]) -> bool:
    """Return whether a differential run reached both validated final outputs."""
    try:
        payload, _source_digest = extract_payload(report)
        validate_payload(payload, expected_samples)
        validate_result_table(result_table)
    except (KeyError, OSError, TypeError, ValueError):
        return False
    return True


def seed_reference(project: Path = PROJECT) -> Path:
    pair = project / "pairs/HepG2_vs_K562"
    pair.mkdir(parents=True, exist_ok=True)
    with (pair / "analysis.lock").open("a", encoding="utf-8") as seed_lock:
        _lock_file(seed_lock)
        return seed_reference_locked(pair)


def seed_reference_locked(pair: Path) -> Path:
    completed = pair / "complete.json"
    payload_path = pair / "results/report_payload.json.gz"
    if completed.is_file() and payload_path.is_file():
        marker = json.loads(completed.read_text(encoding="utf-8"))
        validate_payload(json.loads(gzip.decompress(payload_path.read_bytes())), marker["samples"])
        result_table = pair / "results/diff_footprints_results.txt"
        validate_result_table(result_table)
        marker["results_sha256"] = sha256(result_table)
        write_json_atomic(completed, marker)
        return payload_path
    payload, raw_digest = extract_payload(REFERENCE_REPORT)
    if raw_digest != REFERENCE_JSON_SHA256:
        raise ValueError("The preserved K562/HepG2 report payload has changed")
    validate_payload(
        payload,
        ["K562_rep1", "K562_rep2", "K562_rep3", "HepG2_rep1", "HepG2_rep2", "HepG2_rep3"],
    )
    if payload["conditions"] != ["K562", "HepG2"]:
        raise ValueError("Unexpected condition order in the preserved reference")
    results = pair / "results"
    results.mkdir(parents=True, exist_ok=True)
    metrics = compact_payload(payload, payload_path)
    result_table = results / "diff_footprints_results.txt"
    if not REFERENCE_RESULTS.is_file():
        raise ValueError(f"The compressed reference result fixture is missing: {REFERENCE_RESULTS}")
    with gzip.open(REFERENCE_RESULTS, "rb") as source, result_table.open("wb") as destination:
        shutil.copyfileobj(source, destination)
    validate_result_table(result_table)
    marker = {
        "schema": "fp-tools.encode-cancer-pair.v1",
        "comparison": "HepG2_vs_K562",
        "conditions": ["K562", "HepG2"],
        "source": str(REFERENCE_REPORT.relative_to(ROOT)),
        "source_uncompressed_json_sha256": raw_digest,
        "points": len(payload["points"]),
        "aggregate_motifs": len(payload["aggregate"]["motifs"]),
        "aggregate_site_set": EXPECTED_AGGREGATE_SITE_SET,
        "samples": ["K562_rep1", "K562_rep2", "K562_rep3", "HepG2_rep1", "HepG2_rep2", "HepG2_rep3"],
        "results_sha256": sha256(result_table),
        "completed_at": utc_now(),
        **metrics,
    }
    write_json_atomic(completed, marker)
    return payload_path


def merge_pair_peaks(pair_dir: Path, conditions: tuple[str, str], peaks: pd.DataFrame, *, allow_download: bool) -> tuple[Path, Path]:
    peak_dir = pair_dir / "peaks"
    peak_dir.mkdir(parents=True, exist_ok=True)
    merged = peak_dir / "merged_peaks.bed"
    bins = peak_dir / "merged_peaks.50bp_bins.bed"
    if merged.is_file() and bins.is_file() and merged.stat().st_size and bins.stat().st_size:
        return merged, bins
    selected = peaks[peaks["condition"].isin(conditions)].sort_values(["condition", "peak_accession"])
    sources = [ensure_peak(row, allow_download=allow_download) for row in selected.itertuples(index=False)]
    unsorted = peak_dir / "pair_peaks.unsorted.bed"
    with unsorted.open("w", encoding="utf-8") as handle:
        for source in sources:
            with gzip.open(source, "rt", encoding="utf-8") as incoming:
                for line in incoming:
                    fields = line.split("\t", 3)
                    if len(fields) < 3 or fields[0] in {"chrM", "chrMT", "M", "MT"} or "_" in fields[0]:
                        continue
                    handle.write("\t".join(fields[:3]) + "\n")
    sorted_bed = peak_dir / "pair_peaks.sorted.bed"
    with sorted_bed.open("wb") as output:
        run([str(BEDTOOLS), "sort", "-i", str(unsorted)], pair_dir / "logs/peaks.log", stdout=output)
    temporary = merged.with_name(merged.name + ".part")
    with temporary.open("wb") as output:
        run([str(BEDTOOLS), "merge", "-i", str(sorted_bed)], pair_dir / "logs/peaks.log", stdout=output)
    temporary.replace(merged)
    temporary_bins = bins.with_name(bins.name + ".part")
    with temporary_bins.open("wb") as output:
        run([str(BEDTOOLS), "makewindows", "-b", str(merged), "-w", "50"], pair_dir / "logs/peaks.log", stdout=output)
    temporary_bins.replace(bins)
    unsorted.unlink()
    sorted_bed.unlink()
    return merged, bins


def locate_track(root: Path, sample: str, suffix: str) -> Path:
    matches = sorted(path for path in root.rglob(f"*{suffix}") if sample in path.name or sample in str(path.parent))
    if len(matches) != 1 or not matches[0].stat().st_size:
        raise ValueError(f"Expected one {suffix} track for {sample}, found {matches}")
    return matches[0]


def run_pair(comparison: str, *, cores: int, allow_download: bool, keep_work: bool) -> Path:
    samples, peaks, comparisons = load_design()
    selected_pair = comparisons.loc[comparisons["comparison"].eq(comparison)]
    if len(selected_pair) != 1:
        raise ValueError(f"Unknown comparison: {comparison}")
    row = selected_pair.iloc[0]
    conditions = (str(row.cond1), str(row.cond2))
    if frozenset(conditions) == {"K562", "HepG2"}:
        return seed_reference()
    pair_dir = PROJECT / "pairs" / comparison
    pair_dir.mkdir(parents=True, exist_ok=True)
    pair_lock = (pair_dir / "analysis.lock").open("a", encoding="utf-8")
    _lock_file(pair_lock)
    completed = pair_dir / "complete.json"
    payload_path = pair_dir / "results/report_payload.json.gz"
    if completed.is_file() and payload_path.is_file():
        marker = json.loads(completed.read_text(encoding="utf-8"))
        validate_payload(json.loads(gzip.decompress(payload_path.read_bytes())), marker["samples"])
        result_table = pair_dir / "results/diff_footprints_results.txt"
        validate_result_table(result_table)
        result_digest = sha256(result_table)
        if marker.get("results_sha256") not in (None, result_digest):
            raise ValueError(f"Result-table checksum mismatch: {comparison}")
        marker["results_sha256"] = result_digest
        write_json_atomic(completed, marker)
        if not keep_work:
            clean_pair_outputs(pair_dir / "results", pair_dir / "work")
        return payload_path
    pair_samples = samples[samples["condition"].isin(conditions)].sort_values(["condition", "display_order"])
    bam_paths = [ensure_bam(sample, allow_download=allow_download) for sample in pair_samples.itertuples(index=False)]
    merged, bins = merge_pair_peaks(pair_dir, conditions, peaks, allow_download=allow_download)
    work = pair_dir / "work"
    corrected_dir = work / "atac_correct"
    corrected_dir.mkdir(parents=True, exist_ok=True)
    sample_names = pair_samples["sample"].tolist()
    corrected = []
    if not all(any(sample in path.name or sample in str(path.parent) for path in corrected_dir.rglob("*_corrected.bw")) for sample in sample_names):
        run([
            str(ROOT / ".venv/bin/atac-correct"), "--bams", *map(str, bam_paths),
            "--sample-names", *sample_names, "--genome", str(GENOME), "--peaks", str(merged),
            "--blacklist", str(BLACKLIST), "--outdir", str(corrected_dir), "--cores", str(cores),
        ], pair_dir / "logs/atac_correct.log")
    corrected = [locate_track(corrected_dir, sample, "_corrected.bw") for sample in sample_names]
    normalized_dir = work / "normalized_corrected_q95"
    normalized_dir.mkdir(parents=True, exist_ok=True)
    if len(list(normalized_dir.glob("*.background_scale_q95.bw"))) != len(sample_names):
        run([
            str(ROOT / ".venv/bin/normalize-bigwig"), "--bigwigs", *map(str, corrected),
            "--background", str(bins), "--method", "background-scale", "--stat", "q95",
            "--target", "median", "--chrom-sizes", str(CHROM_SIZES), "--outdir", str(normalized_dir),
        ], pair_dir / "logs/normalize_q95.log")
    normalized = [locate_track(normalized_dir, sample, ".background_scale_q95.bw") for sample in sample_names]
    footprint_dir = work / "footprints_corrected_q95"
    footprint_dir.mkdir(parents=True, exist_ok=True)
    footprints = [footprint_dir / f"{sample}.footprints.bw" for sample in sample_names]
    if not all(output.is_file() and output.stat().st_size for output in footprints):
        run([
            str(ROOT / ".venv/bin/call-footprints"), "--signals", *map(str, normalized),
            "--regions", str(merged), "--outputs", *map(str, footprints),
            "--score", "footprint", "--cores", str(cores),
        ], pair_dir / "logs/call_footprints.log")
    results = pair_dir / "results"
    results.mkdir(parents=True, exist_ok=True)
    report = results / f"diff_footprints_{conditions[0]}_{conditions[1]}.html"
    result_table = results / "diff_footprints_results.txt"
    if not report_outputs_valid(report, result_table, sample_names):
        for stale in results.iterdir():
            if stale.is_dir():
                shutil.rmtree(stale)
            else:
                stale.unlink()
        run([
            str(ROOT / ".venv/bin/diff-footprints"), "--motif-db", "jaspar2026_vertebrates",
            "--signals", *map(str, footprints), "--sample-names", *sample_names,
            "--cond-names", *pair_samples["condition"].tolist(), "--genome", str(GENOME),
            "--peaks", str(merged), "--outdir", str(results), "--prefix", "diff_footprints",
            "--normalization", "none", "--replicate-report", "auto",
            "--aggregate-signals", *map(str, normalized), "--aggregate-normalization", "none",
            "--aggregate-site-set", "all", "--plot-aggregate", "all", "--aggregate-flank", "100",
            "--motif-outputs", "summary", "--report-label", REPORT_LABEL, "--skip-excel", "--cores", str(cores),
        ], pair_dir / "logs/diff_footprints.log")
    payload, source_digest = extract_payload(report)
    validate_payload(payload, sample_names)
    validate_result_table(result_table)
    metrics = compact_payload(payload, payload_path)
    marker = {
        "schema": "fp-tools.encode-cancer-pair.v1",
        "comparison": comparison,
        "conditions": payload["conditions"],
        "samples": sample_names,
        "bam_accessions": pair_samples["bam_accession"].tolist(),
        "peak_accessions": peaks.loc[peaks["condition"].isin(conditions), "peak_accession"].tolist(),
        "merged_peaks_sha256": sha256(merged),
        "source_uncompressed_json_sha256": source_digest,
        "points": len(payload["points"]),
        "aggregate_motifs": len(payload["aggregate"]["motifs"]),
        "aggregate_site_set": EXPECTED_AGGREGATE_SITE_SET,
        "results_sha256": sha256(result_table),
        "completed_at": utc_now(),
        **metrics,
    }
    write_json_atomic(completed, marker)
    if not keep_work:
        clean_pair_outputs(results, work)
    return payload_path


def verify_project() -> list[dict]:
    _samples, _peaks, comparisons = load_design()
    records = []
    for comparison in comparisons["comparison"]:
        pair = PROJECT / "pairs" / comparison
        payload_path = pair / "results/report_payload.json.gz"
        marker_path = pair / "complete.json"
        if not payload_path.is_file() or not marker_path.is_file():
            raise ValueError(f"Incomplete pair: {comparison}")
        marker = json.loads(marker_path.read_text(encoding="utf-8"))
        payload = json.loads(gzip.decompress(payload_path.read_bytes()))
        validate_payload(payload, marker["samples"])
        result_table = pair / "results/diff_footprints_results.txt"
        validate_result_table(result_table)
        if marker.get("results_sha256") != sha256(result_table):
            raise ValueError(f"Result-table checksum mismatch: {comparison}")
        if marker["payload_sha256"] != sha256(payload_path):
            raise ValueError(f"Payload checksum mismatch: {comparison}")
        records.append(marker)
    return records


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("preflight")
    sub.add_parser("seed-reference")
    download_inputs = sub.add_parser("download-inputs")
    download_inputs.add_argument("--workers", type=int, default=3)
    run_one = sub.add_parser("run-pair")
    run_one.add_argument("comparison")
    run_all = sub.add_parser("run")
    for item in (run_one, run_all):
        item.add_argument("--cores", type=int, default=max(1, min(16, os.cpu_count() or 1)))
        item.add_argument("--download", action="store_true")
        item.add_argument("--keep-work", action="store_true")
    run_all.add_argument("--shard-index", type=int, default=0)
    run_all.add_argument("--shard-count", type=int, default=1)
    sub.add_parser("verify")
    return parser.parse_args(argv)


def preflight() -> None:
    samples, peaks, comparisons = load_design()
    required = [GENOME, Path(str(GENOME) + ".fai"), BLACKLIST, CHROM_SIZES, BEDTOOLS, SAMTOOLS]
    required += [ROOT / f".venv/bin/{name}" for name in ("atac-correct", "normalize-bigwig", "call-footprints", "diff-footprints")]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise ValueError(f"Missing required inputs/tools: {missing}")
    print(json.dumps({
        "samples": len(samples), "replicates": EXPECTED_REPLICATES, "peak_files": len(peaks),
        "comparisons": len(comparisons), "new_bams_to_download": [
            row.bam_accession for row in samples.itertuples(index=False) if existing_bam(row.bam_accession) is None
        ],
    }, indent=2))


def download_inputs(workers: int) -> None:
    samples, peaks, _comparisons = load_design()
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
        bam_futures = [
            executor.submit(ensure_bam, row, allow_download=True)
            for row in samples.itertuples(index=False)
        ]
        peak_futures = [
            executor.submit(ensure_peak, row, allow_download=True)
            for row in peaks.itertuples(index=False)
        ]
        for future in concurrent.futures.as_completed(bam_futures + peak_futures):
            print(future.result(), flush=True)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if fcntl is None:
        raise RuntimeError("The concurrent ENCODE benchmark runner requires POSIX file locking")
    if args.command == "preflight":
        preflight()
    elif args.command == "seed-reference":
        print(seed_reference())
    elif args.command == "download-inputs":
        download_inputs(args.workers)
    elif args.command == "run-pair":
        print(run_pair(args.comparison, cores=args.cores, allow_download=args.download, keep_work=args.keep_work))
    elif args.command == "run":
        _samples, _peaks, comparisons = load_design()
        if args.shard_count < 1 or not 0 <= args.shard_index < args.shard_count:
            raise ValueError("--shard-index must be in [0, --shard-count)")
        seed_reference()
        selected = comparisons.iloc[args.shard_index :: args.shard_count]
        for comparison in selected["comparison"]:
            print(f"[{utc_now()}] {comparison}", flush=True)
            run_pair(comparison, cores=args.cores, allow_download=args.download, keep_work=args.keep_work)
    else:
        records = verify_project()
        print(json.dumps({"comparisons": len(records), "verified_at": utc_now()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
