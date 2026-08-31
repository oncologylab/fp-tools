#!/usr/bin/env python3
"""Diagnose TF-specific expected-bias removal with fixed geometries."""

from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from search_tf_footprint_models import (
    binary_metrics,
    candidate_from_row,
    extract_profiles,
    score_candidate,
)


REQUIRED_TRACKS = {"raw", "expected", "corrected"}


def safe_correlation(first: np.ndarray, second: np.ndarray) -> float:
    finite = np.isfinite(first) & np.isfinite(second)
    if finite.sum() < 3 or np.std(first[finite]) == 0 or np.std(second[finite]) == 0:
        return np.nan
    return float(np.corrcoef(first[finite], second[finite])[0, 1])


def diagnosis(row: dict[str, float | int | str]) -> str:
    flags = []
    expected_auroc = float(row["expected_auroc"])
    correction_delta = float(row["correction_delta_auroc"])
    residual = abs(float(row["corrected_expected_correlation"]))
    if abs(expected_auroc - 0.5) >= 0.05:
        flags.append("sequence_bias_label_association")
    if residual >= 0.20:
        flags.append("residual_expected_bias")
    if correction_delta <= -0.03:
        flags.append("correction_harms_discrimination")
    if correction_delta >= 0.03:
        flags.append("correction_improves_discrimination")
    return ";".join(flags) if flags else "no_large_bias_effect"


def evaluate(
    sites: pd.DataFrame,
    winners: pd.DataFrame,
    tracks: pd.DataFrame,
    cache_dir: Path,
    flank: int,
    split: str,
) -> pd.DataFrame:
    cache_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for cell, cell_sites in sites.groupby("cell", sort=True):
        cell_sites = cell_sites.reset_index(drop=True)
        for model, model_tracks in tracks[tracks["cell"] == cell].groupby("model", sort=True):
            available = set(model_tracks["track"])
            if not REQUIRED_TRACKS.issubset(available):
                continue
            profiles = {}
            for track in REQUIRED_TRACKS:
                source = model_tracks[model_tracks["track"] == track].iloc[0]
                cache = cache_dir / f"{cell}.{model}.{track}.flank{flank}.npz"
                if cache.is_file():
                    values = np.load(cache)["profiles"]
                else:
                    values, valid = extract_profiles(cell_sites, Path(source.signal), flank)
                    np.savez_compressed(cache, profiles=values, valid=valid)
                if len(values) != len(cell_sites):
                    raise ValueError(f"site/profile row mismatch: {cache}")
                profiles[track] = values
            for winner in winners[winners["cell"] == cell].itertuples(index=False):
                positions = np.flatnonzero(
                    (cell_sites["tf"].to_numpy() == str(winner.tf))
                    & (cell_sites["chromosome_split"].to_numpy() == split)
                )
                labels = cell_sites.iloc[positions]["chip_label"].to_numpy(dtype=int)
                candidate = replace(candidate_from_row(winner), correction=str(model))
                scores = {
                    track: score_candidate(values[positions], candidate)
                    for track, values in profiles.items()
                }
                metrics = {track: binary_metrics(labels, score) for track, score in scores.items()}
                row = {
                    "cell": cell,
                    "tf": str(winner.tf),
                    "model": model,
                    "candidate_geometry": candidate.identifier,
                    "n_sites": metrics["raw"]["n_sites"],
                    "positive_sites": metrics["raw"]["positive_sites"],
                    "raw_auroc": metrics["raw"]["auroc"],
                    "raw_auprc": metrics["raw"]["auprc"],
                    "expected_auroc": metrics["expected"]["auroc"],
                    "expected_auprc": metrics["expected"]["auprc"],
                    "corrected_auroc": metrics["corrected"]["auroc"],
                    "corrected_auprc": metrics["corrected"]["auprc"],
                    "raw_expected_correlation": safe_correlation(scores["raw"], scores["expected"]),
                    "corrected_expected_correlation": safe_correlation(scores["corrected"], scores["expected"]),
                }
                row["correction_delta_auroc"] = float(row["corrected_auroc"]) - float(row["raw_auroc"])
                row["correction_delta_auprc"] = float(row["corrected_auprc"]) - float(row["raw_auprc"])
                row["diagnosis"] = diagnosis(row)
                rows.append(row)
    return pd.DataFrame(rows)


def summarize_diagnostics(rows: pd.DataFrame) -> pd.DataFrame:
    summaries = []
    for (cell, tf), group in rows.groupby(["cell", "tf"], sort=True):
        summaries.append(
            {
                "cell": cell,
                "tf": tf,
                "models": len(group),
                "expected_auroc_min": group["expected_auroc"].min(),
                "expected_auroc_max": group["expected_auroc"].max(),
                "max_abs_expected_label_delta": (group["expected_auroc"] - 0.5).abs().max(),
                "max_abs_raw_expected_correlation": group["raw_expected_correlation"].abs().max(),
                "max_abs_corrected_expected_correlation": group["corrected_expected_correlation"].abs().max(),
                "best_correction_delta_auroc": group["correction_delta_auroc"].max(),
                "worst_correction_delta_auroc": group["correction_delta_auroc"].min(),
                "improving_models": ";".join(group.loc[group["correction_delta_auroc"] >= 0.03, "model"]),
                "harming_models": ";".join(group.loc[group["correction_delta_auroc"] <= -0.03, "model"]),
            }
        )
    return pd.DataFrame(summaries)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--winners", type=Path, required=True)
    parser.add_argument("--tracks", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--split", default="test")
    args = parser.parse_args(argv)
    result = evaluate(
        pd.read_csv(args.sites, sep="\t"),
        pd.read_csv(args.winners, sep="\t"),
        pd.read_csv(args.tracks, sep="\t"),
        args.cache_dir,
        args.flank,
        args.split,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    if args.summary:
        summarize_diagnostics(result).to_csv(args.summary, sep="\t", index=False)
    print(result.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
