#!/usr/bin/env python
"""Motif-relaxed and motif-free candidate reranking."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path: str | Path) -> pd.DataFrame:
    """Read TSV/BED-like tables, including fp-tools candidate files with # headers."""

    path = Path(path)
    first = path.read_text(encoding="utf-8").splitlines()[0]
    if first.startswith("#") and "	" in first:
        columns = first.lstrip("#").split("	")
        return pd.read_csv(path, sep="	", comment="#", names=columns)
    return pd.read_csv(path, sep="	")


def read_family_map(path: str | Path | None) -> dict[str, str]:
    if path is None:
        return {}
    mapping = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            mapping[fields[0]] = fields[1]
    return mapping


def _numeric_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(np.zeros(len(frame)), index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").fillna(0.0).astype(float)


def _minmax(series: pd.Series) -> pd.Series:
    if len(series) == 0:
        return series.astype(float)
    min_value = float(series.min())
    max_value = float(series.max())
    if not np.isfinite(min_value) or not np.isfinite(max_value) or max_value == min_value:
        return pd.Series(np.zeros(len(series)), index=series.index, dtype=float)
    return (series - min_value) / (max_value - min_value)


def rerank_sites(
    sites: str | Path,
    output: str | Path,
    score_columns: list[str],
    weights: list[float] | None = None,
    family_map: str | Path | None = None,
    motif_column: str = "motif_id",
    family_bonus: float = 0.0,
    top_per_family: int | None = None,
) -> pd.DataFrame:
    """Rank sites using normalized score columns and optional motif family grouping."""

    frame = read_table(sites).copy()
    if weights is None:
        weights = [1.0] * len(score_columns)
    if len(weights) != len(score_columns):
        raise ValueError("--weights must match --score-columns length")

    rank_score = pd.Series(np.zeros(len(frame)), index=frame.index, dtype=float)
    for column, weight in zip(score_columns, weights):
        norm_column = f"{column}_norm"
        frame[norm_column] = _minmax(_numeric_series(frame, column))
        rank_score += float(weight) * frame[norm_column]

    mapping = read_family_map(family_map)
    if mapping and motif_column in frame.columns:
        frame["motif_family"] = frame[motif_column].map(mapping).fillna(frame[motif_column].astype(str))
        if family_bonus:
            family_size = frame.groupby("motif_family")[motif_column].transform("nunique")
            rank_score += np.log1p(family_size.astype(float)) * float(family_bonus)
    elif "motif_family" not in frame.columns:
        frame["motif_family"] = "."

    frame["rank_score"] = rank_score
    frame = frame.sort_values(["rank_score", "chrom", "start"], ascending=[False, True, True], kind="mergesort")
    if top_per_family is not None and top_per_family > 0:
        frame = frame.groupby("motif_family", group_keys=False).head(int(top_per_family))
        frame = frame.sort_values(["rank_score", "chrom", "start"], ascending=[False, True, True], kind="mergesort")
    if "rank" in frame.columns:
        frame = frame.rename(columns={"rank": "input_rank"})
    frame.insert(0, "rank", range(1, len(frame) + 1))

    out_path = Path(output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out_path, sep="	", index=False)
    return frame


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Rerank motif-relaxed or motif-free candidate tables.")
    parser.add_argument("--sites", required=True, help="Candidate/feature/prediction TSV or BED-like table.")
    parser.add_argument("--out", required=True, help="Output ranked TSV.")
    parser.add_argument("--score-columns", nargs="+", default=["candidate_score"], help="Numeric columns to normalize and combine.")
    parser.add_argument("--weights", nargs="*", type=float, default=None, help="Optional weights matching --score-columns.")
    parser.add_argument("--family-map", default=None, help="Optional TSV mapping motif_id to motif_family.")
    parser.add_argument("--motif-column", default="motif_id", help="Motif identifier column for --family-map.")
    parser.add_argument("--family-bonus", type=float, default=0.0, help="Bonus weight for motif families with multiple motifs.")
    parser.add_argument("--top-per-family", type=int, default=None, help="Keep at most N rows per motif family after ranking.")
    args = parser.parse_args(argv)

    frame = rerank_sites(
        args.sites,
        args.out,
        score_columns=args.score_columns,
        weights=args.weights,
        family_map=args.family_map,
        motif_column=args.motif_column,
        family_bonus=args.family_bonus,
        top_per_family=args.top_per_family,
    )
    print(f"Wrote {len(frame)} ranked sites to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
