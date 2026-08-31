#!/usr/bin/env python3
"""Test learned TF-profile models with chromosome-separated selection.

This is benchmark research, not a production scorer.  Hyperparameters are
selected on validation chromosomes, then refit on train plus validation and
evaluated once on test chromosomes.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from search_tf_footprint_models import binary_metrics


CORRECTIONS = ("raw", "PWM", "DWM")


@dataclass(frozen=True)
class FeatureSpec:
    signals: tuple[str, ...]
    normalization: str
    folded: bool
    oriented: bool = False
    bins: int = 41

    @property
    def identifier(self) -> str:
        signals = "+".join(self.signals)
        shape = "folded" if self.folded else "full"
        orientation = "oriented" if self.oriented else "genomic"
        return f"{signals}.{self.normalization}.{shape}.{orientation}.b{self.bins}"


@dataclass(frozen=True)
class ModelSpec:
    family: str
    parameter: float

    @property
    def identifier(self) -> str:
        return f"{self.family}.{self.parameter:g}"


class MatchedTemplate:
    def __init__(self, ridge: float):
        self.ridge = ridge

    def fit(self, features: np.ndarray, labels: np.ndarray):
        negative = features[labels == 0]
        positive = features[labels == 1]
        variance = (negative.var(axis=0) + positive.var(axis=0)) / 2.0
        floor = max(float(np.median(variance[variance > 0])) if np.any(variance > 0) else 1.0, 1e-8)
        self.weight_ = (positive.mean(axis=0) - negative.mean(axis=0)) / (
            variance + self.ridge * floor
        )
        self.offset_ = (positive.mean(axis=0) + negative.mean(axis=0)) / 2.0
        norm = np.linalg.norm(self.weight_)
        if norm > 0:
            self.weight_ /= norm
        return self

    def decision_function(self, features: np.ndarray) -> np.ndarray:
        return (features - self.offset_) @ self.weight_


def pool_profiles(profiles: np.ndarray, bins: int) -> np.ndarray:
    """Average fixed profiles into equally spaced bins."""

    edges = np.linspace(0, profiles.shape[1], bins + 1, dtype=int)
    return np.column_stack(
        [profiles[:, edges[index]:edges[index + 1]].mean(axis=1) for index in range(bins)]
    )


def normalize_profiles(profiles: np.ndarray, method: str) -> np.ndarray:
    if method == "none":
        return profiles
    if method != "outer_rms":
        raise ValueError(f"unknown profile normalization: {method}")
    quarter = max(1, profiles.shape[1] // 4)
    outer = np.concatenate([profiles[:, :quarter], profiles[:, -quarter:]], axis=1)
    background = outer.mean(axis=1, keepdims=True)
    scale = np.sqrt(np.mean(np.square(outer - background), axis=1, keepdims=True))
    floor = max(float(np.nanmedian(scale)) * 1e-6, 1e-8)
    return (profiles - background) / (scale + floor)


def fold_profiles(features: np.ndarray) -> np.ndarray:
    """Fold a motif-centred feature matrix without duplicating its centre."""

    width = features.shape[1]
    middle = width // 2
    if width % 2:
        paired = (features[:, :middle][:, ::-1] + features[:, middle + 1:]) / 2.0
        return np.column_stack([features[:, middle], paired])
    return (features[:, :middle][:, ::-1] + features[:, middle:]) / 2.0


def make_features(
    profiles: dict[str, np.ndarray],
    spec: FeatureSpec,
    strands: np.ndarray | None = None,
) -> np.ndarray:
    parts = []
    for signal in spec.signals:
        values = np.asarray(profiles[signal], dtype=float)
        if spec.oriented:
            if strands is None or len(strands) != len(values):
                raise ValueError("oriented features require one motif strand per profile")
            values = values.copy()
            reverse = np.asarray(strands).astype(str) == "-"
            values[reverse] = values[reverse, ::-1]
        values = normalize_profiles(values, spec.normalization)
        values = pool_profiles(values, spec.bins)
        if spec.folded:
            values = fold_profiles(values)
        parts.append(values)
    return np.concatenate(parts, axis=1)


def feature_grid() -> list[FeatureSpec]:
    signal_sets = [(signal,) for signal in CORRECTIONS] + [CORRECTIONS]
    return [
        FeatureSpec(signals, normalization, folded, oriented)
        for signals in signal_sets
        for normalization in ("none", "outer_rms")
        for folded, oriented in ((False, False), (False, True), (True, False))
    ]


def model_grid(feature: FeatureSpec) -> list[ModelSpec]:
    specs = [
        *(ModelSpec("template", value) for value in (0.1, 1.0, 10.0)),
        *(ModelSpec("logistic_l2", value) for value in (0.01, 0.1, 1.0)),
        *(ModelSpec("lda", value) for value in (0.1, 0.5, 0.9)),
    ]
    if feature.normalization == "outer_rms" and (
        feature.signals == CORRECTIONS or feature.signals == ("DWM",)
    ):
        specs.extend(ModelSpec("extra_trees", value) for value in (5.0, 20.0))
    return specs


def build_model(spec: ModelSpec, seed: int):
    if spec.family == "template":
        return MatchedTemplate(spec.parameter)
    if spec.family == "logistic_l2":
        return make_pipeline(
            StandardScaler(),
            LogisticRegression(
                C=spec.parameter,
                class_weight="balanced",
                max_iter=2000,
                solver="liblinear",
                random_state=seed,
            ),
        )
    if spec.family == "lda":
        return LinearDiscriminantAnalysis(solver="lsqr", shrinkage=spec.parameter)
    if spec.family == "extra_trees":
        return ExtraTreesClassifier(
            n_estimators=200,
            min_samples_leaf=int(spec.parameter),
            max_features="sqrt",
            class_weight="balanced",
            n_jobs=1,
            random_state=seed,
        )
    raise ValueError(f"unknown model family: {spec.family}")


def model_scores(model, features: np.ndarray) -> np.ndarray:
    if hasattr(model, "decision_function"):
        return np.asarray(model.decision_function(features), dtype=float)
    return np.asarray(model.predict_proba(features)[:, 1], dtype=float)


def selection_score(metrics: dict[str, float | int]) -> float:
    prevalence = float(metrics["positive_sites"]) / max(int(metrics["n_sites"]), 1)
    adjusted_auprc = (float(metrics["auprc"]) - prevalence) / max(1.0 - prevalence, 1e-6)
    return float(metrics["auroc"]) + adjusted_auprc


def orientation_comparison(rows: pd.DataFrame) -> pd.DataFrame:
    metrics = ["validation_auroc", "validation_auprc", "test_auroc", "test_auprc"]
    pivot = rows.pivot(index=["cell", "tf"], columns="oriented", values=metrics)
    pivot.columns = [
        f"{metric}_{'oriented' if oriented else 'genomic'}"
        for metric, oriented in pivot.columns
    ]
    pivot = pivot.reset_index()
    for metric in metrics:
        pivot[f"delta_{metric}"] = pivot[f"{metric}_oriented"] - pivot[f"{metric}_genomic"]
    return pivot


def load_cell_profiles(cache_dir: Path, cell: str, flank: int) -> dict[str, np.ndarray]:
    output = {}
    for correction in CORRECTIONS:
        path = cache_dir / f"{cell}.{correction}.flank{flank}.npz"
        if not path.is_file():
            raise FileNotFoundError(path)
        output[correction] = np.load(path)["profiles"]
    lengths = {len(values) for values in output.values()}
    if len(lengths) != 1:
        raise ValueError(f"profile cache lengths disagree for {cell}: {sorted(lengths)}")
    return output


def evaluate(
    development_sites: pd.DataFrame,
    test_sites: pd.DataFrame,
    development_cache: Path,
    test_cache: Path,
    flank: int,
    seed: int,
) -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame
]:
    validation_rows = []
    family_test_rows = []
    winner_rows = []
    overall_test_rows = []
    prediction_rows = []
    orientation_test_rows = []
    for cell in sorted(set(development_sites["cell"]).intersection(test_sites["cell"])):
        dev_cell = development_sites[development_sites["cell"] == cell].reset_index(drop=True)
        test_cell = test_sites[test_sites["cell"] == cell].reset_index(drop=True)
        dev_profiles = load_cell_profiles(development_cache, cell, flank)
        heldout_profiles = load_cell_profiles(test_cache, cell, flank)
        if len(dev_cell) != len(next(iter(dev_profiles.values()))) or len(test_cell) != len(next(iter(heldout_profiles.values()))):
            raise ValueError(f"site/profile row mismatch for {cell}")
        dev_features = {
            spec.identifier: make_features(dev_profiles, spec, dev_cell["TFBS_strand"].to_numpy())
            for spec in feature_grid()
        }
        test_features = {
            spec.identifier: make_features(heldout_profiles, spec, test_cell["TFBS_strand"].to_numpy())
            for spec in feature_grid()
        }
        common_tfs = sorted(set(dev_cell["tf"]).intersection(test_cell["tf"]))
        for tf in common_tfs:
            dev_positions = np.flatnonzero(dev_cell["tf"].to_numpy() == tf)
            test_positions = np.flatnonzero(test_cell["tf"].to_numpy() == tf)
            train_positions = dev_positions[dev_cell.iloc[dev_positions]["chromosome_split"].to_numpy() == "train"]
            validation_positions = dev_positions[dev_cell.iloc[dev_positions]["chromosome_split"].to_numpy() == "validation"]
            train_labels = dev_cell.iloc[train_positions]["chip_label"].to_numpy(dtype=int)
            validation_labels = dev_cell.iloc[validation_positions]["chip_label"].to_numpy(dtype=int)
            test_labels = test_cell.iloc[test_positions]["chip_label"].to_numpy(dtype=int)
            if any(len(np.unique(labels)) != 2 for labels in (train_labels, validation_labels, test_labels)):
                continue
            tf_validation = []
            for feature in feature_grid():
                train_x = dev_features[feature.identifier][train_positions]
                validation_x = dev_features[feature.identifier][validation_positions]
                for model_spec in model_grid(feature):
                    model = build_model(model_spec, seed)
                    model.fit(train_x, train_labels)
                    metrics = binary_metrics(validation_labels, model_scores(model, validation_x))
                    row = {
                        "cell": cell,
                        "tf": tf,
                        "feature": feature.identifier,
                        "signals": "+".join(feature.signals),
                        "normalization": feature.normalization,
                        "folded": feature.folded,
                        "oriented": feature.oriented,
                        "model_family": model_spec.family,
                        "model_parameter": model_spec.parameter,
                        **metrics,
                    }
                    row["selection_score"] = selection_score(metrics)
                    validation_rows.append(row)
                    tf_validation.append((row, feature, model_spec))
            for family in sorted({item[2].family for item in tf_validation}):
                eligible = [item for item in tf_validation if item[2].family == family]
                selected = max(eligible, key=lambda item: (item[0]["selection_score"], item[0]["auprc"]))
                row, feature, model_spec = selected
                fit_positions = np.concatenate([train_positions, validation_positions])
                fit_labels = dev_cell.iloc[fit_positions]["chip_label"].to_numpy(dtype=int)
                model = build_model(model_spec, seed)
                model.fit(dev_features[feature.identifier][fit_positions], fit_labels)
                metrics = binary_metrics(test_labels, model_scores(model, test_features[feature.identifier][test_positions]))
                family_test_rows.append(
                    {
                        **{key: row[key] for key in (
                            "cell", "tf", "feature", "signals", "normalization", "folded", "oriented",
                            "model_family", "model_parameter", "selection_score",
                        )},
                        "validation_auroc": row["auroc"],
                        "validation_auprc": row["auprc"],
                        **{f"test_{key}": value for key, value in metrics.items()},
                    }
                )
            for oriented in (False, True):
                eligible = [item for item in tf_validation if item[1].oriented == oriented]
                selected_orientation = max(
                    eligible, key=lambda item: (item[0]["selection_score"], item[0]["auprc"])
                )
                row, feature, model_spec = selected_orientation
                fit_positions = np.concatenate([train_positions, validation_positions])
                fit_labels = dev_cell.iloc[fit_positions]["chip_label"].to_numpy(dtype=int)
                model = build_model(model_spec, seed)
                model.fit(dev_features[feature.identifier][fit_positions], fit_labels)
                metrics = binary_metrics(
                    test_labels,
                    model_scores(model, test_features[feature.identifier][test_positions]),
                )
                orientation_test_rows.append(
                    {
                        **{key: row[key] for key in (
                            "cell", "tf", "feature", "signals", "normalization", "folded",
                            "oriented", "model_family", "model_parameter", "selection_score",
                        )},
                        "validation_auroc": row["auroc"],
                        "validation_auprc": row["auprc"],
                        **{f"test_{key}": value for key, value in metrics.items()},
                    }
                )
            selected = max(tf_validation, key=lambda item: (item[0]["selection_score"], item[0]["auprc"]))
            row, feature, model_spec = selected
            winner_rows.append({**row, "selection_status": "frozen_before_test"})
            fit_positions = np.concatenate([train_positions, validation_positions])
            fit_labels = dev_cell.iloc[fit_positions]["chip_label"].to_numpy(dtype=int)
            model = build_model(model_spec, seed)
            model.fit(dev_features[feature.identifier][fit_positions], fit_labels)
            test_scores = model_scores(model, test_features[feature.identifier][test_positions])
            test_metrics = binary_metrics(test_labels, test_scores)
            overall_test_rows.append(
                {
                    **{key: row[key] for key in (
                        "cell", "tf", "feature", "signals", "normalization", "folded", "oriented",
                        "model_family", "model_parameter", "selection_score",
                    )},
                    "validation_auroc": row["auroc"],
                    "validation_auprc": row["auprc"],
                    **{f"test_{key}": value for key, value in test_metrics.items()},
                }
            )
            for site, label, score in zip(
                test_cell.iloc[test_positions].itertuples(index=False),
                test_labels,
                test_scores,
            ):
                prediction_rows.append(
                    {
                        "cell": cell,
                        "tf": tf,
                        "TFBS_chr": str(site.TFBS_chr),
                        "TFBS_start": int(site.TFBS_start),
                        "TFBS_end": int(site.TFBS_end),
                        "chip_label": int(label),
                        "classifier_score": float(score),
                        "feature": feature.identifier,
                        "model_family": model_spec.family,
                    }
                )
    return (
        pd.DataFrame(validation_rows),
        pd.DataFrame(family_test_rows),
        pd.DataFrame(winner_rows),
        pd.DataFrame(overall_test_rows),
        pd.DataFrame(prediction_rows),
        pd.DataFrame(orientation_test_rows),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--development-sites", type=Path, required=True)
    parser.add_argument("--test-sites", type=Path, required=True)
    parser.add_argument("--development-cache", type=Path, required=True)
    parser.add_argument("--test-cache", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--flank", type=int, default=100)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    validation, family_test, winners, overall_test, predictions, orientation_test = evaluate(
        pd.read_csv(args.development_sites, sep="\t"),
        pd.read_csv(args.test_sites, sep="\t"),
        args.development_cache,
        args.test_cache,
        args.flank,
        args.seed,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    validation.to_csv(args.outdir / "classifier_validation_metrics.tsv.gz", sep="\t", index=False, compression="gzip")
    family_test.to_csv(args.outdir / "classifier_family_test_metrics.tsv", sep="\t", index=False)
    winners.to_csv(args.outdir / "classifier_frozen_winners.tsv", sep="\t", index=False)
    overall_test.to_csv(args.outdir / "classifier_overall_test_metrics.tsv", sep="\t", index=False)
    predictions.to_csv(
        args.outdir / "classifier_overall_test_predictions.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    orientation_test.to_csv(
        args.outdir / "classifier_orientation_test_metrics.tsv", sep="\t", index=False
    )
    orientation_comparison(orientation_test).to_csv(
        args.outdir / "classifier_orientation_comparison.tsv", sep="\t", index=False
    )
    print(overall_test.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
