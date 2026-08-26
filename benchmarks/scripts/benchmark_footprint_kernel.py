#!/usr/bin/env python
"""Benchmark fast versus legacy call-footprints kernels and compare outputs."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from fp_tools.utils import bigwig as pyBigWig

try:
    import resource
except ImportError:  # Windows does not provide the POSIX resource module.
    resource = None


def run_timed(command: list[str], log_path: Path) -> dict[str, object]:
    start = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    elapsed = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss if resource else 0
    return {
        "command": " ".join(shlex.quote(str(part)) for part in command),
        "exit_code": completed.returncode,
        "wall_seconds": round(elapsed, 3),
        "peak_rss_kb": int(after),
        "log": str(log_path),
    }


def compare_bigwigs(a_path: Path, b_path: Path, chunk_size: int = 1_000_000) -> dict[str, object]:
    a = pyBigWig.open(str(a_path))
    b = pyBigWig.open(str(b_path))
    try:
        if a.chroms() != b.chroms():
            raise ValueError(f"bigWig chromosome headers differ: {a_path} vs {b_path}")
        max_abs = 0.0
        max_rel = 0.0
        sum_abs = 0.0
        count = 0
        for chrom, length in a.chroms().items():
            for start in range(0, length, chunk_size):
                end = min(start + chunk_size, length)
                av = np.asarray(a.values(chrom, start, end, numpy=True), dtype="float64")
                bv = np.asarray(b.values(chrom, start, end, numpy=True), dtype="float64")
                mask = np.isfinite(av) | np.isfinite(bv)
                if not mask.any():
                    continue
                aa = np.nan_to_num(av[mask], nan=0.0, posinf=0.0, neginf=0.0)
                bb = np.nan_to_num(bv[mask], nan=0.0, posinf=0.0, neginf=0.0)
                diff = np.abs(aa - bb)
                rel = diff / np.maximum(np.abs(aa), 1e-12)
                max_abs = max(max_abs, float(diff.max(initial=0.0)))
                max_rel = max(max_rel, float(rel.max(initial=0.0)))
                sum_abs += float(diff.sum())
                count += int(diff.size)
        return {
            "bigwig_count": count,
            "bigwig_max_abs_diff": max_abs,
            "bigwig_mean_abs_diff": sum_abs / max(count, 1),
            "bigwig_max_rel_diff": max_rel,
        }
    finally:
        a.close()
        b.close()


def _read_candidate_bed(path: Path) -> tuple[dict[tuple[str, ...], float], int]:
    rows: dict[tuple[str, ...], float] = {}
    total = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            total += 1
            # Ignore generated footprint_N names and rounded BED score. Keep
            # genomic interval, strand, source region, and footprint center.
            key = (parts[0], parts[1], parts[2], parts[5], parts[6], parts[7])
            rows[key] = float(parts[8])
    return rows, total


def compare_candidate_beds(a_path: Path, b_path: Path) -> dict[str, object]:
    a, a_total = _read_candidate_bed(a_path)
    b, b_total = _read_candidate_bed(b_path)
    a_keys = set(a)
    b_keys = set(b)
    common = a_keys & b_keys
    diffs = [abs(a[key] - b[key]) for key in common]
    union = a_keys | b_keys
    return {
        "bed_legacy_rows": a_total,
        "bed_fast_rows": b_total,
        "bed_common_coords": len(common),
        "bed_legacy_only_coords": len(a_keys - b_keys),
        "bed_fast_only_coords": len(b_keys - a_keys),
        "bed_coordinate_jaccard": len(common) / max(len(union), 1),
        "bed_max_shared_raw_score_diff": max(diffs) if diffs else 0.0,
        "bed_mean_shared_raw_score_diff": sum(diffs) / max(len(diffs), 1),
    }


def compare_tables(a_path: Path, b_path: Path, id_col: str = "output_prefix") -> dict[str, object]:
    a = pd.read_csv(a_path, sep="\t")
    b = pd.read_csv(b_path, sep="\t")
    common_cols = [col for col in a.columns if col in b.columns]
    if id_col not in common_cols:
        id_col = common_cols[0]
    a = a[common_cols].sort_values(id_col).reset_index(drop=True)
    b = b[common_cols].sort_values(id_col).reset_index(drop=True)
    if len(a) != len(b):
        return {"table_rows_a": len(a), "table_rows_b": len(b), "table_same_rows": False}
    numeric_cols = [col for col in common_cols if pd.api.types.is_numeric_dtype(a[col]) and pd.api.types.is_numeric_dtype(b[col])]
    max_abs = 0.0
    for col in numeric_cols:
        diff = (a[col].astype(float) - b[col].astype(float)).abs().replace([np.inf, -np.inf], np.nan).dropna()
        if not diff.empty:
            max_abs = max(max_abs, float(diff.max()))
    non_numeric_cols = [col for col in common_cols if col not in numeric_cols]
    same_labels = a[non_numeric_cols].fillna("").astype(str).equals(b[non_numeric_cols].fillna("").astype(str))
    return {
        "table_rows_a": len(a),
        "table_rows_b": len(b),
        "table_same_rows": True,
        "table_same_non_numeric": bool(same_labels),
        "table_numeric_max_abs_diff": max_abs,
    }


def _run_call_footprints(args: argparse.Namespace, signal: Path, kernel: str, prefix: str, output_bed: bool = True) -> tuple[Path, Path | None, dict[str, object]]:
    out_bw = args.outdir / f"{prefix}_{kernel}_footprints.bw"
    out_bed = args.outdir / f"{prefix}_{kernel}_candidates.bed"
    command = [
        args.call_footprints,
        "--signal",
        str(signal),
        "--regions",
        str(args.regions),
        "--output",
        str(out_bw),
        "--score",
        "footprint",
        "--footprint-kernel",
        kernel,
        "--cores",
        str(args.cores),
        "--verbosity",
        str(args.verbosity),
    ]
    if output_bed:
        command.extend(["--output-bed", str(out_bed)])
    timing = run_timed(command, args.outdir / f"{prefix}_{kernel}.log")
    if timing["exit_code"] != 0:
        raise RuntimeError(f"{prefix} {kernel} call-footprints failed; see {timing['log']}")
    return out_bw, out_bed if output_bed else None, timing


def _motif_argument(args: argparse.Namespace) -> list[str]:
    if args.motif_db:
        return ["--motif-db", args.motif_db]
    if args.motifs:
        return ["--motifs", str(args.motifs)]
    raise ValueError("end-to-end workflow check requires --motifs or --motif-db")


def run_workflow_check(args: argparse.Namespace, first_legacy_bw: Path, first_fast_bw: Path) -> dict[str, object]:
    if not args.workflow_second_signal:
        return {}
    if not args.genome:
        raise ValueError("end-to-end workflow check requires --genome")

    second_legacy_bw, _, second_legacy = _run_call_footprints(args, args.workflow_second_signal, "legacy", args.workflow_second_name, output_bed=False)
    second_fast_bw, _, second_fast = _run_call_footprints(args, args.workflow_second_signal, "fast", args.workflow_second_name, output_bed=False)

    results: dict[str, object] = {
        "workflow_second_legacy_seconds": second_legacy["wall_seconds"],
        "workflow_second_fast_seconds": second_fast["wall_seconds"],
    }
    for kernel, signals in {
        "legacy": [first_legacy_bw, second_legacy_bw],
        "fast": [first_fast_bw, second_fast_bw],
    }.items():
        match_dir = args.outdir / f"{kernel}_match_motifs"
        diff_dir = args.outdir / f"{kernel}_diff_footprints"
        match_cmd = [
            args.match_motifs,
            "--signals",
            *[str(path) for path in signals],
            "--sample-names",
            args.workflow_first_name,
            args.workflow_second_name,
            "--genome",
            str(args.genome),
            "--peaks",
            str(args.regions),
            "--outdir",
            str(match_dir),
            "--cores",
            str(args.cores),
            "--skip-excel",
            "--verbosity",
            str(args.verbosity),
            *_motif_argument(args),
        ]
        match = run_timed(match_cmd, args.outdir / f"{kernel}_match_motifs.log")
        if match["exit_code"] != 0:
            raise RuntimeError(f"{kernel} match-motifs failed; see {match['log']}")
        diff_cmd = [
            args.diff_footprints,
            "--signals",
            *[str(path) for path in signals],
            "--sample-names",
            args.workflow_first_name,
            args.workflow_second_name,
            "--cond-names",
            args.workflow_first_condition,
            args.workflow_second_condition,
            "--genome",
            str(args.genome),
            "--peaks",
            str(args.regions),
            "--outdir",
            str(diff_dir),
            "--cores",
            str(args.cores),
            "--skip-excel",
            "--plot-aggregate",
            "off",
            "--motif-outputs",
            "summary",
            "--verbosity",
            str(args.verbosity),
            *_motif_argument(args),
        ]
        diff = run_timed(diff_cmd, args.outdir / f"{kernel}_diff_footprints.log")
        if diff["exit_code"] != 0:
            raise RuntimeError(f"{kernel} diff-footprints failed; see {diff['log']}")
        results[f"workflow_{kernel}_match_seconds"] = match["wall_seconds"]
        results[f"workflow_{kernel}_diff_seconds"] = diff["wall_seconds"]
        results[f"workflow_{kernel}_diff_results"] = str(diff_dir / "diff_footprints_results.txt")

    results.update(
        {
            f"workflow_{key}": value
            for key, value in compare_tables(
                args.outdir / "legacy_diff_footprints" / "diff_footprints_results.txt",
                args.outdir / "fast_diff_footprints" / "diff_footprints_results.txt",
            ).items()
        }
    )
    return results


def run_kernel_benchmark(args: argparse.Namespace) -> dict[str, object]:
    args.outdir.mkdir(parents=True, exist_ok=True)
    legacy_bw, legacy_bed, legacy = _run_call_footprints(args, args.signal, "legacy", args.workflow_first_name)
    fast_bw, fast_bed, fast = _run_call_footprints(args, args.signal, "fast", args.workflow_first_name)

    summary: dict[str, object] = {
        "legacy_seconds": legacy["wall_seconds"],
        "fast_seconds": fast["wall_seconds"],
        "speedup": round(float(legacy["wall_seconds"]) / max(float(fast["wall_seconds"]), 1e-9), 3),
        "legacy_peak_rss_kb": legacy["peak_rss_kb"],
        "fast_peak_rss_kb": fast["peak_rss_kb"],
        "legacy_bigwig": str(legacy_bw),
        "fast_bigwig": str(fast_bw),
        "legacy_bed": str(legacy_bed),
        "fast_bed": str(fast_bed),
    }
    summary.update(compare_bigwigs(legacy_bw, fast_bw, chunk_size=args.chunk_size))
    assert legacy_bed is not None and fast_bed is not None
    summary.update(compare_candidate_beds(legacy_bed, fast_bed))
    summary.update(run_workflow_check(args, legacy_bw, fast_bw))
    return summary


def write_summary(summary: dict[str, object], outdir: Path) -> None:
    json_path = outdir / "kernel_benchmark_summary.json"
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tsv_path = outdir / "kernel_benchmark_summary.tsv"
    with tsv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary), delimiter="\t")
        writer.writeheader()
        writer.writerow(summary)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signal", type=Path, required=True, help="Corrected cut-site bigWig input.")
    parser.add_argument("--regions", type=Path, required=True, help="BED regions for call-footprints.")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory for benchmark files.")
    parser.add_argument("--call-footprints", default=str(Path(sys.executable).with_name("call-footprints")), help="call-footprints executable.")
    parser.add_argument("--cores", type=int, default=1, help="Cores passed to each call-footprints run.")
    parser.add_argument("--chunk-size", type=int, default=1_000_000, help="bigWig comparison chunk size.")
    parser.add_argument("--verbosity", type=int, default=1, help="fp-tools command verbosity.")
    parser.add_argument("--workflow-second-signal", type=Path, help="Optional second corrected bigWig for end-to-end match/diff consistency.")
    parser.add_argument("--workflow-first-name", default="sample1", help="Sample name for --signal.")
    parser.add_argument("--workflow-second-name", default="sample2", help="Sample name for --workflow-second-signal.")
    parser.add_argument("--workflow-first-condition", default="condition1", help="Condition label for --signal in the workflow check.")
    parser.add_argument("--workflow-second-condition", default="condition2", help="Condition label for --workflow-second-signal in the workflow check.")
    parser.add_argument("--genome", type=Path, help="Genome FASTA for optional workflow check.")
    parser.add_argument("--motifs", type=Path, help="Motif file for optional workflow check.")
    parser.add_argument("--motif-db", help="Built-in motif database name for optional workflow check.")
    parser.add_argument("--match-motifs", default=str(Path(sys.executable).with_name("match-motifs")), help="match-motifs executable.")
    parser.add_argument("--diff-footprints", default=str(Path(sys.executable).with_name("diff-footprints")), help="diff-footprints executable.")
    args = parser.parse_args(argv)

    summary = run_kernel_benchmark(args)
    write_summary(summary, args.outdir)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
