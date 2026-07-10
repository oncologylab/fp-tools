#!/usr/bin/env python3
"""Memory-bounded comparison of two processed ATAC sample directories."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import tempfile
from pathlib import Path

import pyBigWig
import pysam

from fp_tools.tools.prepare_atac import _bedgraph_to_bigwig


def bam_metrics(path: Path) -> dict[str, float | int]:
    metrics = {
        "records": 0,
        "mapped": 0,
        "proper_pair": 0,
        "read1": 0,
        "mitochondrial": 0,
        "mapq_sum": 0,
    }
    lengths = [0] * 1001
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            metrics["records"] += 1
            if not read.is_unmapped:
                metrics["mapped"] += 1
                metrics["mapq_sum"] += read.mapping_quality
                if read.reference_name in {"chrM", "chrMT", "M", "MT"}:
                    metrics["mitochondrial"] += 1
            metrics["proper_pair"] += int(read.is_proper_pair)
            metrics["read1"] += int(read.is_read1)
            if (
                read.is_read1
                and read.is_proper_pair
                and 0 < abs(read.template_length) <= 1000
            ):
                lengths[abs(read.template_length)] += 1
    mapped = max(1, metrics["mapped"])
    total = max(1, metrics["records"])
    metrics["mean_mapq"] = metrics.pop("mapq_sum") / mapped
    metrics["proper_pair_fraction"] = metrics["proper_pair"] / total
    metrics["mitochondrial_fraction"] = metrics["mitochondrial"] / mapped
    fragments = sum(lengths)
    midpoint = (fragments + 1) // 2
    cumulative = 0
    median = 0
    for length, count in enumerate(lengths):
        cumulative += count
        if cumulative >= midpoint:
            median = length
            break
    metrics["fragment_count_le_1000"] = fragments
    metrics["median_fragment_length"] = median
    metrics["nucleosome_free_fraction"] = sum(lengths[1:101]) / max(1, fragments)
    metrics["mononucleosome_fraction"] = sum(lengths[180:248]) / max(1, fragments)
    return metrics


def frip(bam: Path, peaks: Path) -> float:
    total = int(
        subprocess.check_output(["samtools", "view", "-c", str(bam)], text=True)
    )
    left = subprocess.Popen(
        ["bedtools", "intersect", "-u", "-abam", str(bam), "-b", str(peaks)],
        stdout=subprocess.PIPE,
    )
    right = subprocess.run(
        ["samtools", "view", "-c", "-"],
        stdin=left.stdout,
        capture_output=True,
        text=True,
    )
    if left.stdout:
        left.stdout.close()
    if left.wait() or right.returncode:
        raise RuntimeError("Failed to calculate FRiP")
    return int(right.stdout.strip() or 0) / max(1, total)


def write_name_hashes(bam_path: Path, output: Path) -> None:
    with (
        pysam.AlignmentFile(bam_path, "rb") as bam,
        output.open("w", encoding="ascii") as handle,
    ):
        for read in bam.fetch(until_eof=True):
            digest = hashlib.blake2b(
                read.query_name.encode(), digest_size=8
            ).hexdigest()
            handle.write(digest + "\n")


def sorted_hashes(bam_path: Path, output: Path, cores: int) -> Path:
    raw = output.with_suffix(".raw")
    write_name_hashes(bam_path, raw)
    with output.open("w", encoding="ascii") as handle:
        subprocess.run(
            ["sort", "-u", "-S", "2G", f"--parallel={cores}", str(raw)],
            check=True,
            stdout=handle,
        )
    raw.unlink()
    return output


def hash_overlap(left: Path, right: Path) -> dict[str, float | int]:
    intersection = union = 0
    with left.open() as a, right.open() as b:
        va, vb = a.readline(), b.readline()
        while va or vb:
            if not vb or (va and va < vb):
                union += 1
                va = a.readline()
            elif not va or vb < va:
                union += 1
                vb = b.readline()
            else:
                union += 1
                intersection += 1
                va, vb = a.readline(), b.readline()
    return {
        "read_name_intersection": intersection,
        "read_name_union": union,
        "read_name_jaccard": intersection / max(1, union),
    }


def peak_intervals(path: Path) -> list[tuple[str, int, int]]:
    rows = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) >= 3 and fields[1].isdigit() and fields[2].isdigit():
                rows.append((fields[0], int(fields[1]), int(fields[2])))
    return rows


def merged_bp(
    rows: list[tuple[str, int, int]],
) -> tuple[int, dict[str, list[tuple[int, int]]]]:
    grouped: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in rows:
        grouped.setdefault(chrom, []).append((start, end))
    total = 0
    for chrom, intervals in grouped.items():
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
            else:
                merged.append((start, end))
        grouped[chrom] = merged
        total += sum(end - start for start, end in merged)
    return total, grouped


def peak_metrics(left: Path, right: Path) -> dict[str, float | int]:
    a, b = peak_intervals(left), peak_intervals(right)
    bp_a, ma = merged_bp(a)
    bp_b, mb = merged_bp(b)
    overlap = 0
    for chrom in set(ma) & set(mb):
        i = j = 0
        while i < len(ma[chrom]) and j < len(mb[chrom]):
            x, y = ma[chrom][i], mb[chrom][j]
            overlap += max(0, min(x[1], y[1]) - max(x[0], y[0]))
            if x[1] < y[1]:
                i += 1
            else:
                j += 1
    widths_a = sorted(end - start for _, start, end in a)
    widths_b = sorted(end - start for _, start, end in b)
    return {
        "baseline_peak_count": len(a),
        "legacy_peak_count": len(b),
        "baseline_peak_bp": bp_a,
        "legacy_peak_bp": bp_b,
        "baseline_median_width": widths_a[len(widths_a) // 2] if widths_a else 0,
        "legacy_median_width": widths_b[len(widths_b) // 2] if widths_b else 0,
        "peak_overlap_bp": overlap,
        "peak_bp_jaccard": overlap / max(1, bp_a + bp_b - overlap),
        "baseline_bp_recovered": overlap / max(1, bp_a),
        "legacy_bp_recovered": overlap / max(1, bp_b),
    }


def make_rp10m(bam: Path, chrom_sizes: Path, output: Path) -> Path:
    records = int(
        subprocess.check_output(["samtools", "view", "-c", str(bam)], text=True)
    )
    scale = 10_000_000 / max(1, records // 2)
    bedgraph = output.with_suffix(".bedgraph")
    with bedgraph.open("w") as handle:
        subprocess.run(
            [
                "bedtools",
                "genomecov",
                "-ibam",
                str(bam),
                "-bg",
                "-pc",
                "-scale",
                str(scale),
            ],
            check=True,
            stdout=handle,
        )
    _bedgraph_to_bigwig(bedgraph, chrom_sizes, output)
    bedgraph.unlink()
    return output


def bigwig_metrics(
    left: Path, right: Path, bin_size: int = 10_000
) -> dict[str, float | int]:
    n = 0
    sx = sy = sxx = syy = sxy = 0.0
    lsx = lsy = lsxx = lsyy = lsxy = 0.0
    with pyBigWig.open(str(left)) as a, pyBigWig.open(str(right)) as b:
        for chrom, size in a.chroms().items():
            if chrom not in b.chroms():
                continue
            size = min(size, b.chroms(chrom))
            for start in range(0, size, bin_size):
                end = min(size, start + bin_size)
                x = a.stats(chrom, start, end, type="mean")[0] or 0.0
                y = b.stats(chrom, start, end, type="mean")[0] or 0.0
                lx, ly = math.log1p(max(0.0, x)), math.log1p(max(0.0, y))
                n += 1
                sx += x
                sy += y
                sxx += x * x
                syy += y * y
                sxy += x * y
                lsx += lx
                lsy += ly
                lsxx += lx * lx
                lsyy += ly * ly
                lsxy += lx * ly

    def corr(x, y, xx, yy, xy):
        return (n * xy - x * y) / math.sqrt(
            max(1e-30, (n * xx - x * x) * (n * yy - y * y))
        )

    return {
        "bigwig_bins": n,
        "bigwig_pearson": corr(sx, sy, sxx, syy, sxy),
        "bigwig_log1p_pearson": corr(lsx, lsy, lsxx, lsyy, lsxy),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-bam", type=Path, required=True)
    parser.add_argument("--baseline-peaks", type=Path, required=True)
    parser.add_argument("--legacy-bam", type=Path, required=True)
    parser.add_argument("--legacy-peaks", type=Path, required=True)
    parser.add_argument("--legacy-bigwig", type=Path, required=True)
    parser.add_argument("--baseline-bigwig", type=Path)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cores", type=int, default=1)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    baseline_bw = args.baseline_bigwig or make_rp10m(
        args.baseline_bam, args.chrom_sizes, args.output_dir / "baseline.rp10m.bw"
    )
    report = {
        "baseline_bam": bam_metrics(args.baseline_bam),
        "legacy_bam": bam_metrics(args.legacy_bam),
    }
    report["baseline_frip"] = frip(args.baseline_bam, args.baseline_peaks)
    report["legacy_frip"] = frip(args.legacy_bam, args.legacy_peaks)
    report.update(peak_metrics(args.baseline_peaks, args.legacy_peaks))
    with tempfile.TemporaryDirectory(dir=args.output_dir) as temp:
        temp = Path(temp)
        left = sorted_hashes(args.baseline_bam, temp / "baseline.hashes", args.cores)
        right = sorted_hashes(args.legacy_bam, temp / "legacy.hashes", args.cores)
        report.update(hash_overlap(left, right))
    report.update(bigwig_metrics(baseline_bw, args.legacy_bigwig))
    (args.output_dir / "comparison.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    with (args.output_dir / "comparison.tsv").open("w", encoding="utf-8") as handle:
        handle.write("metric\tvalue\n")
        for key, value in report.items():
            if isinstance(value, dict):
                for nested, nested_value in value.items():
                    handle.write(f"{key}.{nested}\t{nested_value}\n")
            else:
                handle.write(f"{key}\t{value}\n")
    lines = ["# ATAC profile comparison", "", "| Metric | Value |", "|---|---:|"]
    for key, value in report.items():
        if not isinstance(value, dict):
            lines.append(f"| {key} | {value} |")
    (args.output_dir / "comparison.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
