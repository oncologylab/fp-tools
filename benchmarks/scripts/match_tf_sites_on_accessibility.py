#!/usr/bin/env python3
"""Match TF-positive and TF-negative sites on motif score and ATAC coverage."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from build_footprint_site_labels import propensity_match


def standardized_difference(frame: pd.DataFrame, feature: str) -> float:
    positive = frame.loc[frame["chip_label"] == 1, feature].to_numpy(dtype=float)
    negative = frame.loc[frame["chip_label"] == 0, feature].to_numpy(dtype=float)
    pooled = np.sqrt((np.var(positive) + np.var(negative)) / 2.0)
    return float((np.mean(positive) - np.mean(negative)) / pooled) if pooled > 0 else 0.0


def attach_accessibility(sites: pd.DataFrame, cache_dir: Path, flank: int) -> pd.DataFrame:
    parts = []
    for cell, cell_sites in sites.groupby("cell", sort=True):
        cell_sites = cell_sites.reset_index(drop=True)
        payload = np.load(cache_dir / f"{cell}.raw.flank{flank}.npz")
        profiles = payload["profiles"]
        if len(profiles) != len(cell_sites):
            raise ValueError(f"raw profile cache row mismatch for {cell}")
        output = cell_sites.copy()
        output["accessibility"] = np.sum(np.abs(profiles), axis=1)
        output["central_accessibility"] = np.sum(
            np.abs(profiles[:, flank - 25:flank + 26]), axis=1
        )
        parts.append(output)
    return pd.concat(parts, ignore_index=True)


def match_sites(
    sites: pd.DataFrame,
    features: list[str],
    negative_ratio: int,
    seed: int,
    method: str = "optimal",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    matched_parts = []
    diagnostics = []
    for (cell, tf, split), group in sites.groupby(
        ["cell", "tf", "chromosome_split"], sort=True
    ):
        if group["chip_label"].nunique() != 2:
            continue
        positive = group[group["chip_label"] == 1]
        negative = group[group["chip_label"] == 0]
        maximum_positives = len(negative) // negative_ratio
        if len(positive) > maximum_positives:
            positive = positive.sample(n=maximum_positives, random_state=seed)
            group = pd.concat([positive, negative], ignore_index=True)
        working = group.rename(
            columns={
                "chip_label": "label", "TFBS_chr": "chrom",
                "TFBS_start": "start", "TFBS_end": "end",
            }
        ).copy()
        if method == "propensity":
            matched = propensity_match(working, features, negative_ratio, seed)
        elif method == "optimal":
            matched = optimal_feature_match(working, features, negative_ratio)
        else:
            raise ValueError(f"unknown matching method: {method}")
        matched = matched.rename(
            columns={
                "label": "chip_label", "chrom": "TFBS_chr",
                "start": "TFBS_start", "end": "TFBS_end",
            }
        )
        matched_parts.append(matched)
        record = {
            "cell": cell,
            "tf": tf,
            "chromosome_split": split,
            "positive_sites": int((matched["chip_label"] == 1).sum()),
            "negative_sites": int((matched["chip_label"] == 0).sum()),
        }
        for feature in features:
            record[f"before_smd_{feature}"] = standardized_difference(group, feature)
            record[f"after_smd_{feature}"] = standardized_difference(matched, feature)
        diagnostics.append(record)
    return pd.concat(matched_parts, ignore_index=True), pd.DataFrame(diagnostics)


def optimal_feature_match(
    frame: pd.DataFrame,
    features: list[str],
    negative_ratio: int = 1,
) -> pd.DataFrame:
    """Globally minimize standardized multivariate positive/negative distance."""

    if negative_ratio != 1:
        raise ValueError("optimal feature matching currently requires negative_ratio=1")
    usable = frame[frame["label"].isin([0, 1])].dropna(subset=features).copy()
    positive = usable[usable["label"] == 1]
    negative = usable[usable["label"] == 0]
    if positive.empty or len(negative) < len(positive):
        raise ValueError(
            f"matching requires {len(positive)} negatives for {len(positive)} positives; "
            f"found {len(negative)}"
        )
    values = usable[features].to_numpy(dtype=float)
    center = np.mean(values, axis=0)
    scale = np.std(values, axis=0)
    scale[scale == 0] = 1.0
    positive_values = (positive[features].to_numpy(dtype=float) - center) / scale
    negative_values = (negative[features].to_numpy(dtype=float) - center) / scale
    distance = cdist(positive_values, negative_values, metric="sqeuclidean")
    positive_rows, negative_columns = linear_sum_assignment(distance)
    selected_positive = positive.iloc[positive_rows].copy()
    selected_negative = negative.iloc[negative_columns].copy()
    pair_distance = np.sqrt(distance[positive_rows, negative_columns])
    selected_positive["match_distance"] = pair_distance
    selected_negative["match_distance"] = pair_distance
    return pd.concat([selected_positive, selected_negative], ignore_index=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--diagnostics", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--features", nargs="+", default=["motif_score", "accessibility"])
    parser.add_argument("--negative-ratio", type=int, default=1)
    parser.add_argument("--method", choices=["optimal", "propensity"], default="optimal")
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)

    sites = pd.read_csv(args.sites, sep="\t")
    attached = attach_accessibility(sites, args.cache_dir, args.flank)
    matched, diagnostics = match_sites(
        attached, args.features, args.negative_ratio, args.seed, args.method
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    matched.to_csv(args.output, sep="\t", index=False, compression="gzip")
    diagnostics.to_csv(args.diagnostics, sep="\t", index=False)
    print(diagnostics.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
