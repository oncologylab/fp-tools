#!/usr/bin/env python
"""Plot coverage of bound JASPAR sites by all de novo motif hits."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

REPO_ROOT = SCRIPT_DIR.parents[1]

from figure_style import apply_style, bold_all_text  # noqa: E402


DEFAULT_VALIDATION_DIR = (
    REPO_ROOT
    / "data/public/processed/encode_k562_hepg2_atac_replicates/fp_tools/denovo_motif_validation_maxcover_n250"
)
DEFAULT_OUT_PREFIX = REPO_ROOT / "manuscript/figures/denovo_all_vs_bound_jaspar_coverage"
WINDOW_HALF_WIDTH = 10


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def merge_interval_dict(data: dict[str, list[tuple[int, int]]]) -> dict[str, list[tuple[int, int]]]:
    return {chrom: merge_intervals(intervals) for chrom, intervals in data.items() if intervals}


def interval_count(data: dict[str, list[tuple[int, int]]]) -> int:
    return sum(len(intervals) for intervals in data.values())


def interval_bp(data: dict[str, list[tuple[int, int]]]) -> int:
    return sum(end - start for intervals in data.values() for start, end in intervals)


def read_center_windows(paths: list[Path], half_width: int = WINDOW_HALF_WIDTH) -> tuple[dict[str, list[tuple[int, int]]], int]:
    windows: dict[str, list[tuple[int, int]]] = {}
    records = 0
    for path in paths:
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                center = (start + end) // 2
                windows.setdefault(chrom, []).append((max(0, center - half_width), center + half_width + 1))
                records += 1
    return merge_interval_dict(windows), records


def intersection(
    data_a: dict[str, list[tuple[int, int]]],
    data_b: dict[str, list[tuple[int, int]]],
) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    for chrom in set(data_a).intersection(data_b):
        a = data_a[chrom]
        b = data_b[chrom]
        i = j = 0
        while i < len(a) and j < len(b):
            start = max(a[i][0], b[j][0])
            end = min(a[i][1], b[j][1])
            if start < end:
                out.setdefault(chrom, []).append((start, end))
            if a[i][1] < b[j][1]:
                i += 1
            else:
                j += 1
    return merge_interval_dict(out)


def subtract_intervals(
    data_a: dict[str, list[tuple[int, int]]],
    data_b: dict[str, list[tuple[int, int]]],
) -> dict[str, list[tuple[int, int]]]:
    out: dict[str, list[tuple[int, int]]] = {}
    for chrom, intervals in data_a.items():
        blockers = data_b.get(chrom, [])
        if not blockers:
            out[chrom] = list(intervals)
            continue
        j = 0
        for start, end in intervals:
            pieces = [(start, end)]
            while j < len(blockers) and blockers[j][1] <= start:
                j += 1
            k = j
            while k < len(blockers) and blockers[k][0] < end:
                block_start, block_end = blockers[k]
                next_pieces = []
                for piece_start, piece_end in pieces:
                    if block_end <= piece_start or block_start >= piece_end:
                        next_pieces.append((piece_start, piece_end))
                        continue
                    if piece_start < block_start:
                        next_pieces.append((piece_start, block_start))
                    if block_end < piece_end:
                        next_pieces.append((block_end, piece_end))
                pieces = next_pieces
                if not pieces:
                    break
                k += 1
            out.setdefault(chrom, []).extend(pieces)
    return merge_interval_dict(out)


def collect_denovo_all(validation_dir: Path) -> list[Path]:
    return sorted((validation_dir / "diff_footprints" / "denovo_only").glob("*/beds/*_all.bed"))


def collect_jaspar_bound(validation_dir: Path, condition: str | None = None) -> list[Path]:
    pattern = f"*_{condition}_bound.bed" if condition else "*_bound.bed"
    paths = sorted((validation_dir / "diff_footprints" / "jaspar2026_plus_denovo").glob(f"*/beds/{pattern}"))
    return [path for path in paths if "_denovo_" not in path.name and "_denovo_" not in str(path.parent.parent.name)]


def coverage_table(validation_dir: Path) -> pd.DataFrame:
    denovo_paths = collect_denovo_all(validation_dir)
    denovo_windows, denovo_records = read_center_windows(denovo_paths)
    rows = []
    labels = [
        ("All JASPAR-bound", None),
        ("K562 JASPAR-bound", "K562"),
        ("HepG2 JASPAR-bound", "HepG2"),
    ]
    for label, condition in labels:
        jaspar_paths = collect_jaspar_bound(validation_dir, condition)
        jaspar_windows, jaspar_records = read_center_windows(jaspar_paths)
        covered = intersection(jaspar_windows, denovo_windows)
        uncovered = subtract_intervals(jaspar_windows, denovo_windows)
        jaspar_bp = interval_bp(jaspar_windows)
        rows.append(
            {
                "site_set": label,
                "condition": condition or "either",
                "jaspar_motif_files": len(jaspar_paths),
                "jaspar_site_records": jaspar_records,
                "denovo_motif_files": len(denovo_paths),
                "denovo_site_records": denovo_records,
                "jaspar_merged_windows": interval_count(jaspar_windows),
                "jaspar_bp": jaspar_bp,
                "covered_by_denovo_all_windows": interval_count(covered),
                "covered_by_denovo_all_bp": interval_bp(covered),
                "not_covered_by_denovo_all_windows": interval_count(uncovered),
                "not_covered_by_denovo_all_bp": interval_bp(uncovered),
                "fraction_jaspar_bp_covered": interval_bp(covered) / jaspar_bp if jaspar_bp else 0.0,
                "window_half_width_bp": WINDOW_HALF_WIDTH,
            }
        )
    return pd.DataFrame(rows)


def plot_coverage(table: pd.DataFrame, out_prefix: Path) -> None:
    apply_style(base_size=7)
    fig, ax = plt.subplots(figsize=(4.7, 2.1))
    y = np.arange(len(table))
    covered = table["fraction_jaspar_bp_covered"].to_numpy(dtype=float) * 100
    uncovered = 100 - covered
    ax.barh(y, covered, color="#4f8dcf", edgecolor="white", linewidth=0.5, height=0.62, label="covered by all de novo hits")
    ax.barh(y, uncovered, left=covered, color="#d7dbe2", edgecolor="white", linewidth=0.5, height=0.62, label="JASPAR-bound only")
    for idx, value in enumerate(covered):
        ax.text(min(value - 2, 96), idx, f"{value:.1f}%", ha="right", va="center", color="white", fontsize=6.4)
    ax.set_yticks(y, table["site_set"].tolist())
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel(f"JASPAR-bound site-window bp (+/-{WINDOW_HALF_WIDTH} bp) covered")
    ax.set_title("All de novo motif hits recover most JASPAR-bound sites", loc="left", pad=5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="0.9", linewidth=0.55)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=2, frameon=False, fontsize=6.2)
    bold_all_text(ax)
    fig.tight_layout(pad=0.6)
    fig.subplots_adjust(bottom=0.34)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    svg_path = out_prefix.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight")
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-dir", type=Path, default=DEFAULT_VALIDATION_DIR)
    parser.add_argument("--out-prefix", type=Path, default=DEFAULT_OUT_PREFIX)
    args = parser.parse_args(argv)

    table = coverage_table(args.validation_dir)
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out_prefix.with_suffix(".tsv"), sep="\t", index=False)
    plot_coverage(table, args.out_prefix)
    print(f"Wrote {args.out_prefix.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
