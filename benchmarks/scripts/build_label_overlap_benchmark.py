#!/usr/bin/env python3
"""Build metrics-ready benchmark tables by overlapping predictions with label BEDs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_COLUMNS = ["chrom", "start", "end", "name", "score"]


def read_bed_table(path: str | Path, default_columns: list[str] | None = None) -> pd.DataFrame:
    """Read a BED-like table with optional header or fp-tools '#chrom' header."""

    path = Path(path)
    first_data = ""
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                first_data = line.rstrip("\n")
                break
    if not first_data:
        return pd.DataFrame(columns=default_columns or DEFAULT_COLUMNS)

    if first_data.startswith("#"):
        columns = first_data.lstrip("#").split("\t")
        return pd.read_csv(path, sep="\t", comment="#", names=columns)
    fields = first_data.split("\t")
    has_header = len(fields) >= 3 and not fields[1].isdigit()
    if has_header:
        return pd.read_csv(path, sep="\t")
    columns = list(default_columns or DEFAULT_COLUMNS)
    while len(columns) < len(fields):
        columns.append(f"extra_{len(columns) + 1}")
    return pd.read_csv(path, sep="\t", header=None, names=columns[:len(fields)])


def read_label_intervals(path: str | Path) -> dict[str, list[tuple[int, int]]]:
    labels: dict[str, list[tuple[int, int]]] = {}
    frame = read_bed_table(path, default_columns=["chrom", "start", "end"])
    for _, row in frame.iterrows():
        labels.setdefault(str(row["chrom"]), []).append((int(row["start"]), int(row["end"])))
    return labels


def overlap_bp(chrom: str, start: int, end: int, labels: dict[str, list[tuple[int, int]]]) -> int:
    best = 0
    for label_start, label_end in labels.get(chrom, []):
        best = max(best, max(0, min(end, label_end) - max(start, label_start)))
    return int(best)


def build_label_overlap_table(
    predictions: str | Path,
    labels_bed: str | Path,
    output: str | Path,
    score_col: str = "score",
    min_overlap_bp: int = 1,
    method: str = "fp-tools",
    tf: str = "",
    cell: str = "",
    metadata_cols: list[str] | None = None,
) -> pd.DataFrame:
    """Create a binary label/score table from scored intervals and BED labels."""

    frame = read_bed_table(predictions)
    required = {"chrom", "start", "end", score_col}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Prediction table is missing required columns: {sorted(missing)}")
    labels = read_label_intervals(labels_bed)

    rows = []
    metadata_cols = [column for column in (metadata_cols or []) if column in frame.columns]
    for _, row in frame.iterrows():
        start = int(row["start"])
        end = int(row["end"])
        overlap = overlap_bp(str(row["chrom"]), start, end, labels)
        out_row = {
            "label": int(overlap >= min_overlap_bp),
            "score": float(row[score_col]),
            "method": method,
            "tf": tf,
            "cell": cell,
            "chrom": row["chrom"],
            "start": start,
            "end": end,
            "overlap_bp": overlap,
        }
        for column in metadata_cols:
            out_row[column] = row[column]
        rows.append(out_row)

    out = pd.DataFrame(rows)
    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True, help="BED-like scored prediction intervals.")
    parser.add_argument("--labels-bed", required=True, help="Positive label BED, e.g. TF ChIP/CUT&RUN peaks.")
    parser.add_argument("--out", required=True, help="Output metrics-ready TSV.")
    parser.add_argument("--score-col", default="score")
    parser.add_argument("--min-overlap-bp", type=int, default=1)
    parser.add_argument("--method", default="fp-tools")
    parser.add_argument("--tf", default="")
    parser.add_argument("--cell", default="")
    parser.add_argument("--metadata-cols", nargs="*", default=["name"])
    args = parser.parse_args()

    table = build_label_overlap_table(
        args.predictions,
        args.labels_bed,
        args.out,
        score_col=args.score_col,
        min_overlap_bp=args.min_overlap_bp,
        method=args.method,
        tf=args.tf,
        cell=args.cell,
        metadata_cols=args.metadata_cols,
    )
    print(f"wrote {len(table)} labeled benchmark rows to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
