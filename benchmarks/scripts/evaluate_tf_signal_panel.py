#!/usr/bin/env python3
"""Evaluate frozen TF geometry across depth and replicate signal panels."""

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


REQUIRED_SIGNALS = {"cell", "sample", "depth", "signal"}


def depth_comparison(
    sample_metrics: pd.DataFrame,
    ensemble_metrics: pd.DataFrame,
    low_depth: str = "10m",
    high_depth: str = "full",
) -> pd.DataFrame:
    low = sample_metrics[sample_metrics["depth"] == low_depth].copy()
    if low.duplicated(["cell", "tf"]).any():
        low = low.groupby(["cell", "tf"], as_index=False)[["auroc", "auprc"]].mean()
    high = ensemble_metrics[ensemble_metrics["depth"] == high_depth].copy()
    merged = low[["cell", "tf", "auroc", "auprc"]].merge(
        high[["cell", "tf", "auroc", "auprc", "replicates"]],
        on=["cell", "tf"],
        suffixes=(f"_{low_depth}", f"_{high_depth}_replicate_mean"),
        validate="one_to_one",
    )
    merged["delta_auroc"] = merged[f"auroc_{high_depth}_replicate_mean"] - merged[f"auroc_{low_depth}"]
    merged["delta_auprc"] = merged[f"auprc_{high_depth}_replicate_mean"] - merged[f"auprc_{low_depth}"]
    return merged


def evaluate_panel(
    sites: pd.DataFrame,
    winners: pd.DataFrame,
    signals: pd.DataFrame,
    cache_dir: Path,
    flank: int,
    split: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    ensemble_rows = []
    for cell, cell_sites in sites.groupby("cell", sort=True):
        cell_sites = cell_sites.reset_index(drop=True)
        sample_scores: dict[tuple[str, str], dict[str, np.ndarray]] = {}
        for signal_row in signals[signals["cell"] == cell].itertuples(index=False):
            cache = cache_dir / f"{signal_row.sample}.flank{flank}.npz"
            if cache.is_file():
                profiles = np.load(cache)["profiles"]
                if len(profiles) != len(cell_sites):
                    raise ValueError(f"profile cache row mismatch: {cache}")
            else:
                profiles, valid = extract_profiles(cell_sites, Path(signal_row.signal), flank)
                np.savez_compressed(cache, profiles=profiles, valid=valid)
            sample_scores[(str(signal_row.depth), str(signal_row.sample))] = {}
            for winner in winners[winners["cell"] == cell].itertuples(index=False):
                positions = np.flatnonzero(
                    (cell_sites["tf"].to_numpy() == str(winner.tf))
                    & (cell_sites["chromosome_split"].to_numpy() == split)
                )
                labels = cell_sites.iloc[positions]["chip_label"].to_numpy(dtype=int)
                candidate = replace(candidate_from_row(winner), correction=str(signal_row.sample))
                scores = score_candidate(profiles[positions], candidate)
                sample_scores[(str(signal_row.depth), str(signal_row.sample))][str(winner.tf)] = scores
                rows.append(
                    {
                        "cell": str(cell), "tf": str(winner.tf), "sample": str(signal_row.sample),
                        "depth": str(signal_row.depth), "chromosome_split": split,
                        **binary_metrics(labels, scores),
                    }
                )
        for depth in sorted({key[0] for key in sample_scores}):
            depth_samples = [key for key in sample_scores if key[0] == depth]
            if len(depth_samples) < 2:
                continue
            for winner in winners[winners["cell"] == cell].itertuples(index=False):
                tf = str(winner.tf)
                positions = np.flatnonzero(
                    (cell_sites["tf"].to_numpy() == tf)
                    & (cell_sites["chromosome_split"].to_numpy() == split)
                )
                labels = cell_sites.iloc[positions]["chip_label"].to_numpy(dtype=int)
                ensemble = np.mean(
                    np.stack([sample_scores[key][tf] for key in depth_samples]), axis=0
                )
                ensemble_rows.append(
                    {
                        "cell": str(cell), "tf": tf, "sample": "replicate_mean",
                        "depth": depth, "chromosome_split": split,
                        "replicates": len(depth_samples),
                        **binary_metrics(labels, ensemble),
                    }
                )
    return pd.DataFrame(rows), pd.DataFrame(ensemble_rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--winners", type=Path, required=True)
    parser.add_argument("--signals", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ensemble-output", type=Path, required=True)
    parser.add_argument("--depth-summary-output", type=Path)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--split", default="test")
    args = parser.parse_args(argv)
    signals = pd.read_csv(args.signals, sep="\t")
    missing = sorted(REQUIRED_SIGNALS.difference(signals.columns))
    if missing:
        parser.error("signal panel is missing columns: " + ", ".join(missing))
    rows, ensembles = evaluate_panel(
        pd.read_csv(args.sites, sep="\t"),
        pd.read_csv(args.winners, sep="\t"),
        signals,
        args.cache_dir,
        args.flank,
        args.split,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(args.output, sep="\t", index=False)
    ensembles.to_csv(args.ensemble_output, sep="\t", index=False)
    if args.depth_summary_output and not ensembles.empty:
        depth_comparison(rows, ensembles).to_csv(
            args.depth_summary_output, sep="\t", index=False
        )
    print(rows.to_string(index=False))
    if not ensembles.empty:
        print("\nReplicate means\n" + ensembles.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
