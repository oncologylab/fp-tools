#!/usr/bin/env python3
"""Transfer validation-selected DWM profile classifiers across replicates."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from search_tf_footprint_models import binary_metrics
from search_tf_profile_classifiers import (
    FeatureSpec,
    ModelSpec,
    build_model,
    load_cell_profiles,
    make_features,
    model_scores,
)


def select_dwm_winners(validation: pd.DataFrame) -> pd.DataFrame:
    eligible = validation[validation["signals"] == "DWM"].copy()
    eligible = eligible.sort_values(
        ["cell", "tf", "selection_score", "auprc", "auroc"],
        ascending=[True, True, False, False, False],
        kind="mergesort",
    )
    return eligible.drop_duplicates(["cell", "tf"], keep="first").reset_index(drop=True)


def feature_from_row(row) -> FeatureSpec:
    match = re.search(r"\.b(\d+)$", str(row.feature))
    if not match:
        raise ValueError(f"cannot parse feature bin count: {row.feature}")
    return FeatureSpec(
        signals=("DWM",),
        normalization=str(row.normalization),
        folded=bool(row.folded),
        oriented=bool(row.oriented),
        bins=int(match.group(1)),
    )


def resolve_profile_cache(signal, panel_cache: Path, flank: int) -> Path:
    explicit = getattr(signal, "profile_cache", None)
    if explicit is not None and pd.notna(explicit):
        return Path(explicit)
    return panel_cache / f"{signal.sample}.flank{flank}.npz"


def evaluate(
    development_sites: pd.DataFrame,
    test_sites: pd.DataFrame,
    development_cache: Path,
    validation: pd.DataFrame,
    panel: pd.DataFrame,
    panel_cache: Path,
    flank: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    winners = select_dwm_winners(validation)
    sample_rows = []
    ensemble_rows = []
    for cell, cell_test in test_sites.groupby("cell", sort=True):
        cell_test = cell_test.reset_index(drop=True)
        cell_dev = development_sites[development_sites["cell"] == cell].reset_index(drop=True)
        cell_winners = winners[winners["cell"] == cell]
        if cell_winners.empty:
            continue
        development_profiles = load_cell_profiles(development_cache, str(cell), flank)["DWM"]
        models = {}
        specs = {}
        for winner in cell_winners.itertuples(index=False):
            tf = str(winner.tf)
            positions = np.flatnonzero(cell_dev["tf"].to_numpy() == tf)
            labels = cell_dev.iloc[positions]["chip_label"].to_numpy(dtype=int)
            feature = feature_from_row(winner)
            features = make_features(
                {"DWM": development_profiles}, feature, cell_dev["TFBS_strand"].to_numpy()
            )
            model_spec = ModelSpec(str(winner.model_family), float(winner.model_parameter))
            model = build_model(model_spec, seed).fit(features[positions], labels)
            models[tf] = model
            specs[tf] = feature

        scores_by_depth: dict[str, dict[str, list[np.ndarray]]] = {}
        for signal in panel[panel["cell"] == cell].itertuples(index=False):
            cache = resolve_profile_cache(signal, panel_cache, flank)
            profiles = np.load(cache)["profiles"]
            if len(profiles) != len(cell_test):
                raise ValueError(f"site/profile row mismatch: {cache}")
            feature_cache = {
                spec.identifier: make_features(
                    {"DWM": profiles}, spec, cell_test["TFBS_strand"].to_numpy()
                )
                for spec in set(specs.values())
            }
            for tf, model in models.items():
                positions = np.flatnonzero(cell_test["tf"].to_numpy() == tf)
                labels = cell_test.iloc[positions]["chip_label"].to_numpy(dtype=int)
                scores = model_scores(model, feature_cache[specs[tf].identifier][positions])
                metrics = binary_metrics(labels, scores)
                sample_rows.append(
                    {
                        "cell": cell, "tf": tf, "sample": str(signal.sample),
                        "depth": str(signal.depth), "feature": specs[tf].identifier,
                        **metrics,
                    }
                )
                scores_by_depth.setdefault(str(signal.depth), {}).setdefault(tf, []).append(scores)
        for depth, by_tf in scores_by_depth.items():
            for tf, scores in by_tf.items():
                if len(scores) < 2:
                    continue
                positions = np.flatnonzero(cell_test["tf"].to_numpy() == tf)
                labels = cell_test.iloc[positions]["chip_label"].to_numpy(dtype=int)
                metrics = binary_metrics(labels, np.mean(np.stack(scores), axis=0))
                ensemble_rows.append(
                    {
                        "cell": cell, "tf": tf, "sample": "replicate_mean",
                        "depth": depth, "replicates": len(scores), "feature": specs[tf].identifier,
                        **metrics,
                    }
                )
    return pd.DataFrame(sample_rows), pd.DataFrame(ensemble_rows)


def replicate_stability(sample_metrics: pd.DataFrame, depth: str = "full") -> pd.DataFrame:
    selected = sample_metrics[sample_metrics["depth"] == depth]
    return selected.groupby(["cell", "tf"], as_index=False).agg(
        replicates=("sample", "nunique"),
        auroc_mean=("auroc", "mean"),
        auroc_min=("auroc", "min"),
        auroc_max=("auroc", "max"),
        auroc_sd=("auroc", "std"),
        auprc_mean=("auprc", "mean"),
        auprc_min=("auprc", "min"),
        auprc_max=("auprc", "max"),
        auprc_sd=("auprc", "std"),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development-sites", type=Path, required=True)
    parser.add_argument("--test-sites", type=Path, required=True)
    parser.add_argument("--development-cache", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--panel-cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ensemble-output", type=Path, required=True)
    parser.add_argument("--stability-output", type=Path)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    sample, ensemble = evaluate(
        pd.read_csv(args.development_sites, sep="\t"),
        pd.read_csv(args.test_sites, sep="\t"),
        args.development_cache,
        pd.read_csv(args.validation, sep="\t"),
        pd.read_csv(args.panel, sep="\t"),
        args.panel_cache,
        args.flank,
        args.seed,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    sample.to_csv(args.output, sep="\t", index=False)
    ensemble.to_csv(args.ensemble_output, sep="\t", index=False)
    if args.stability_output:
        replicate_stability(sample).to_csv(args.stability_output, sep="\t", index=False)
    print(ensemble.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
