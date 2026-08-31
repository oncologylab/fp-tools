#!/usr/bin/env python3
"""Compare frozen per-TF shape candidates with legacy scores on identical sites."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from evaluate_bigwig_site_scores import score_centers
from search_tf_footprint_models import binary_metrics, candidate_from_row, score_candidate


def compare(
    sites: pd.DataFrame,
    winners: pd.DataFrame,
    baselines: pd.DataFrame,
    cache_dir: Path,
    flank: int,
    split: str,
) -> pd.DataFrame:
    rows = []
    for cell, cell_sites in sites.groupby("cell", sort=True):
        cell_sites = cell_sites.reset_index(drop=True)
        baseline_rows = baselines[baselines["cell"] == cell]
        if len(baseline_rows) != 1:
            raise ValueError(f"expected exactly one baseline signal for {cell}")
        baseline_row = baseline_rows.iloc[0]
        baseline_scores = score_centers(cell_sites, Path(baseline_row.signal))
        for winner in winners[winners["cell"] == cell].itertuples(index=False):
            candidate = candidate_from_row(winner)
            payload = np.load(cache_dir / f"{cell}.{candidate.correction}.flank{flank}.npz")
            profiles = payload["profiles"]
            if len(profiles) != len(cell_sites):
                raise ValueError(f"profile cache row mismatch for {cell}/{candidate.correction}")
            positions = np.flatnonzero(
                (cell_sites["tf"].to_numpy() == str(winner.tf))
                & (cell_sites["chromosome_split"].to_numpy() == split)
            )
            candidate_scores = score_candidate(profiles[positions], candidate)
            matched = np.isfinite(candidate_scores) & np.isfinite(baseline_scores[positions])
            labels = cell_sites.iloc[positions].loc[matched, "chip_label"].to_numpy(dtype=int)
            baseline_metric = binary_metrics(labels, baseline_scores[positions][matched])
            candidate_metric = binary_metrics(labels, candidate_scores[matched])
            rows.append(
                {
                    "cell": str(cell),
                    "tf": str(winner.tf),
                    "chromosome_split": split,
                    "baseline": str(baseline_row.method),
                    "candidate": candidate.identifier,
                    "n_sites": candidate_metric["n_sites"],
                    "positive_sites": candidate_metric["positive_sites"],
                    "baseline_auroc": baseline_metric["auroc"],
                    "candidate_auroc": candidate_metric["auroc"],
                    "delta_auroc": candidate_metric["auroc"] - baseline_metric["auroc"],
                    "baseline_auprc": baseline_metric["auprc"],
                    "candidate_auprc": candidate_metric["auprc"],
                    "delta_auprc": candidate_metric["auprc"] - baseline_metric["auprc"],
                }
            )
    return pd.DataFrame(rows)


def bootstrap_compare(
    sites: pd.DataFrame,
    winners: pd.DataFrame,
    baselines: pd.DataFrame,
    cache_dir: Path,
    flank: int,
    split: str,
    n_bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    """Bootstrap paired candidate-minus-legacy deltas by chromosome."""

    if n_bootstrap <= 0:
        return pd.DataFrame()
    rng = np.random.default_rng(seed)
    rows = []
    for cell, cell_sites in sites.groupby("cell", sort=True):
        cell_sites = cell_sites.reset_index(drop=True)
        baseline_row = baselines[baselines["cell"] == cell].iloc[0]
        baseline_scores = score_centers(cell_sites, Path(baseline_row.signal))
        for winner in winners[winners["cell"] == cell].itertuples(index=False):
            candidate = candidate_from_row(winner)
            profiles = np.load(cache_dir / f"{cell}.{candidate.correction}.flank{flank}.npz")["profiles"]
            positions = np.flatnonzero(
                (cell_sites["tf"].to_numpy() == str(winner.tf))
                & (cell_sites["chromosome_split"].to_numpy() == split)
            )
            candidate_scores = score_candidate(profiles[positions], candidate)
            finite = np.isfinite(candidate_scores) & np.isfinite(baseline_scores[positions])
            group = cell_sites.iloc[positions].loc[finite]
            labels = group["chip_label"].to_numpy(dtype=int)
            baseline = baseline_scores[positions][finite]
            candidate_values = candidate_scores[finite]
            blocks = [
                np.flatnonzero(group["TFBS_chr"].to_numpy() == chrom)
                for chrom in sorted(group["TFBS_chr"].unique())
            ]
            distributions = {"auroc": [], "auprc": []}
            for _ in range(n_bootstrap):
                selected = rng.integers(0, len(blocks), size=len(blocks))
                sample = np.concatenate([blocks[index] for index in selected])
                sample_labels = labels[sample]
                if len(np.unique(sample_labels)) != 2:
                    continue
                distributions["auroc"].append(
                    roc_auc_score(sample_labels, candidate_values[sample])
                    - roc_auc_score(sample_labels, baseline[sample])
                )
                distributions["auprc"].append(
                    average_precision_score(sample_labels, candidate_values[sample])
                    - average_precision_score(sample_labels, baseline[sample])
                )
            for metric, values in distributions.items():
                array = np.asarray(values, dtype=float)
                rows.append(
                    {
                        "cell": str(cell), "tf": str(winner.tf),
                        "chromosome_split": split, "metric": metric,
                        "delta_mean": float(np.mean(array)) if len(array) else np.nan,
                        "ci_low": float(np.quantile(array, 0.025)) if len(array) else np.nan,
                        "ci_high": float(np.quantile(array, 0.975)) if len(array) else np.nan,
                        "probability_delta_gt_zero": float(np.mean(array > 0)) if len(array) else np.nan,
                        "successful_bootstraps": int(len(array)),
                        "requested_bootstraps": int(n_bootstrap),
                        "n_chromosome_blocks": int(len(blocks)),
                    }
                )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--winners", type=Path, required=True)
    parser.add_argument("--baselines", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--split", default="validation")
    parser.add_argument("--bootstrap", type=int, default=0)
    parser.add_argument("--bootstrap-output", type=Path)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)

    sites = pd.read_csv(args.sites, sep="\t")
    winners = pd.read_csv(args.winners, sep="\t")
    baselines = pd.read_csv(args.baselines, sep="\t")
    required = {"cell", "method", "signal"}
    missing = sorted(required.difference(baselines.columns))
    if missing:
        parser.error("baseline manifest is missing columns: " + ", ".join(missing))
    result = compare(sites, winners, baselines, args.cache_dir, args.flank, args.split)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False)
    bootstrap = bootstrap_compare(
        sites, winners, baselines, args.cache_dir, args.flank, args.split,
        args.bootstrap, args.seed,
    )
    if not bootstrap.empty:
        bootstrap_output = args.bootstrap_output or args.output.with_name(
            args.output.stem + "_bootstrap.tsv"
        )
        bootstrap.to_csv(bootstrap_output, sep="\t", index=False)
    print(result.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
