#!/usr/bin/env python3
"""Measure post-lock supervised information ceilings for holdout profiles.

This is a diagnostic, never a promotion or model-selection input. Functional
PCA bases are learned from label-free train chromosomes. ChIP labels are used
only by fixed elastic-net classifiers evaluated out of fold across locked test
chromosomes. The output distinguishes missing footprint information from a
label-free detector that chose the wrong shape or latent-state orientation.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any, Sequence

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import binary_metrics, stable_seed  # noqa: E402
from evaluate_locked_holdout_policy import (  # noqa: E402
    Artifact,
    load_artifact,
    reference_to_candidate_indexes,
)
from fp_tools.tools.functional_footprints import (  # noqa: E402
    FunctionalPCA,
    normalize_functional_profiles,
)


SCHEMA = "fp-tools-locked-holdout-information-ceiling-v1"
BASELINE_FEATURES = ("motif_score", "log_accessibility")


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_artifact(value: str) -> tuple[str, str, Path]:
    fields = value.split(",", 2)
    if len(fields) != 3 or not all(fields):
        raise argparse.ArgumentTypeError("artifact must use MODEL,CELL,JSON")
    return fields[0], fields[1], Path(fields[2])


def cross_chromosome_predictions(
    features: np.ndarray,
    labels: np.ndarray,
    chromosomes: np.ndarray,
    *,
    seed: int,
    maximum_iterations: int = 5000,
) -> tuple[np.ndarray, int]:
    """Return fixed elastic-net predictions from leave-chromosome-out folds."""

    values = np.asarray(features, dtype=float)
    labels = np.asarray(labels, dtype=int)
    chromosomes = np.asarray(chromosomes).astype(str)
    if values.ndim != 2 or len(values) != len(labels) or len(labels) != len(chromosomes):
        raise ValueError("features, labels, and chromosomes must have equal rows")
    if maximum_iterations < 1:
        raise ValueError("maximum_iterations must be positive")
    finite = np.isfinite(values).all(axis=1)
    predictions = np.full(len(labels), np.nan, dtype=float)
    folds = 0
    for chromosome in sorted(np.unique(chromosomes[finite])):
        test = finite & (chromosomes == chromosome)
        train = finite & (chromosomes != chromosome)
        if min(test.sum(), train.sum()) < 4:
            continue
        if np.unique(labels[test]).size != 2 or np.unique(labels[train]).size != 2:
            continue
        scaler = StandardScaler().fit(values[train])
        model = LogisticRegression(
            solver="saga",
            l1_ratio=0.5,
            C=1.0,
            class_weight="balanced",
            max_iter=maximum_iterations,
            tol=1e-5,
            random_state=stable_seed("ceiling", chromosome, seed=seed),
        )
        model.fit(scaler.transform(values[train]), labels[train])
        predictions[test] = model.predict_proba(scaler.transform(values[test]))[:, 1]
        folds += 1
    return predictions, folds


def scored_metrics(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, Any]:
    valid = np.isfinite(probabilities)
    if valid.sum() < 4 or np.unique(labels[valid]).size != 2:
        return {
            "sites_scored": int(valid.sum()),
            "auroc": np.nan,
            "auprc": np.nan,
            "brier": np.nan,
        }
    return {
        "sites_scored": int(valid.sum()),
        **{
            key: float(value)
            for key, value in binary_metrics(labels[valid], probabilities[valid]).items()
            if key in {"auroc", "auprc", "brier"}
        },
    }


def artifact_channels(artifact: Artifact) -> list[str]:
    channels = ["combined_residual"]
    for channel in ("shared_strand_residual", "antisymmetric_strand_residual"):
        if channel in artifact.profiles:
            channels.append(channel)
    return channels


def information_ceiling_rows(
    scores: pd.DataFrame,
    artifacts: dict[tuple[str, str], Artifact],
    *,
    variance_threshold: float,
    max_components: int,
    maximum_train_per_tf: int,
    seed: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (cell, tf), task in scores.groupby(["cell", "tf"], sort=True):
        task = task.reset_index(drop=True)
        labels = task["label"].to_numpy(dtype=int)
        chromosomes = task["TFBS_chr"].astype(str).to_numpy()
        baseline = task[list(BASELINE_FEATURES)].to_numpy(dtype=float)
        baseline_predictions, baseline_folds = cross_chromosome_predictions(
            baseline,
            labels,
            chromosomes,
            seed=stable_seed(cell, tf, "baseline", seed=seed),
        )
        baseline_metrics = scored_metrics(labels, baseline_predictions)
        reference = artifacts[("DWM", str(cell))]
        reference_indexes = task["artifact_index"].to_numpy(dtype=int)
        model = str(task["bias_configuration"].iloc[0])
        methods = [("DWM", reference, reference_indexes)]
        if model != "DWM":
            candidate = artifacts[(model, str(cell))]
            mapping = reference_to_candidate_indexes(reference, candidate)
            candidate_indexes = mapping[reference_indexes]
            if np.any(candidate_indexes < 0):
                raise ValueError(f"{cell}/{tf} matched sites are absent from {model}")
            methods.append((model, candidate, candidate_indexes))
        frozen = task["candidate_probability"].to_numpy(dtype=float)
        frozen_metrics = scored_metrics(labels, frozen)
        flipped_metrics = scored_metrics(labels, 1.0 - frozen)
        for correction, artifact, evaluation_indexes in methods:
            observed, _expected = artifact.observed_expected()
            train_mask = (
                artifact.valid
                & artifact.sites["tf"].astype(str).eq(str(tf)).to_numpy()
                & artifact.sites["chromosome_split"].astype(str).eq("train").to_numpy()
            )
            train_indexes = np.flatnonzero(train_mask)
            if len(train_indexes) > maximum_train_per_tf:
                rng = np.random.default_rng(
                    stable_seed(cell, tf, correction, "train-cap", seed=seed)
                )
                train_indexes = np.sort(
                    rng.choice(
                        train_indexes,
                        size=maximum_train_per_tf,
                        replace=False,
                    )
                )
            weights = np.sqrt(
                np.maximum(observed[train_indexes].sum(axis=1), 1.0)
            )
            positions = np.arange(
                -(artifact.profiles["combined_residual"].shape[1] // 2),
                artifact.profiles["combined_residual"].shape[1] // 2 + 1,
                dtype=float,
            )
            for channel in artifact_channels(artifact):
                training_profiles = normalize_functional_profiles(
                    artifact.profiles[channel][train_indexes], positions
                )
                evaluation_profiles = normalize_functional_profiles(
                    artifact.profiles[channel][evaluation_indexes], positions
                )
                fpca = FunctionalPCA(
                    variance_threshold=variance_threshold,
                    max_components=max_components,
                    seed=stable_seed(cell, tf, correction, channel, seed=seed),
                ).fit(training_profiles, sample_weight=weights)
                functional_scores = fpca.transform(evaluation_profiles)
                profile_predictions, profile_folds = cross_chromosome_predictions(
                    functional_scores,
                    labels,
                    chromosomes,
                    seed=stable_seed(cell, tf, correction, channel, "profile", seed=seed),
                )
                combined_predictions, combined_folds = cross_chromosome_predictions(
                    np.column_stack((baseline, functional_scores)),
                    labels,
                    chromosomes,
                    seed=stable_seed(cell, tf, correction, channel, "combined", seed=seed),
                )
                profile_metrics = scored_metrics(labels, profile_predictions)
                combined_metrics = scored_metrics(labels, combined_predictions)
                relative_gain = (
                    (combined_metrics["auprc"] - baseline_metrics["auprc"])
                    / max(baseline_metrics["auprc"], 1e-8)
                    if np.isfinite(combined_metrics["auprc"])
                    and np.isfinite(baseline_metrics["auprc"])
                    else np.nan
                )
                rows.append(
                    {
                        "cell": str(cell),
                        "tf": str(tf),
                        "motif_family": str(task["motif_family"].iloc[0]),
                        "correction": correction,
                        "channel": channel,
                        "sites": int(len(task)),
                        "train_profiles": int(len(train_indexes)),
                        "functional_components": int(len(fpca.components_)),
                        "baseline_folds": int(baseline_folds),
                        "profile_folds": int(profile_folds),
                        "combined_folds": int(combined_folds),
                        **{f"baseline_{key}": value for key, value in baseline_metrics.items()},
                        **{f"profile_{key}": value for key, value in profile_metrics.items()},
                        **{f"combined_{key}": value for key, value in combined_metrics.items()},
                        "combined_relative_auprc_gain_over_baseline": float(relative_gain),
                        "shape_information_above_baseline": bool(relative_gain >= 0.10),
                        "frozen_candidate_auroc": frozen_metrics["auroc"],
                        "frozen_candidate_auprc": frozen_metrics["auprc"],
                        "flipped_candidate_auroc": flipped_metrics["auroc"],
                        "flipped_candidate_auprc": flipped_metrics["auprc"],
                        "latent_state_inversion_signal": bool(
                            frozen_metrics["auroc"] < 0.5
                            and flipped_metrics["auroc"] >= 0.6
                        ),
                        "diagnostic_only": True,
                    }
                )
    output = pd.DataFrame(rows)
    dwm = output[output["correction"].eq("DWM") & output["channel"].eq("combined_residual")][
        ["cell", "tf", "combined_auroc", "combined_auprc"]
    ].rename(
        columns={
            "combined_auroc": "dwm_supervised_combined_auroc",
            "combined_auprc": "dwm_supervised_combined_auprc",
        }
    )
    output = output.merge(dwm, on=["cell", "tf"], how="left", validate="many_to_one")
    output["supervised_auroc_gain_over_dwm"] = (
        output["combined_auroc"] - output["dwm_supervised_combined_auroc"]
    )
    output["supervised_relative_auprc_gain_over_dwm"] = (
        output["combined_auprc"] - output["dwm_supervised_combined_auprc"]
    ) / output["dwm_supervised_combined_auprc"].clip(lower=1e-8)
    return output.sort_values(
        ["cell", "tf", "correction", "channel"], kind="mergesort"
    ).reset_index(drop=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument(
        "--artifact",
        action="append",
        type=parse_artifact,
        required=True,
        metavar="MODEL,CELL,JSON",
    )
    parser.add_argument("--variance-threshold", type=float, default=0.95)
    parser.add_argument("--max-components", type=int, default=20)
    parser.add_argument("--maximum-train-per-tf", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    artifact_paths = {(model, cell): path for model, cell, path in args.artifact}
    if len(artifact_paths) != len(args.artifact):
        raise SystemExit("duplicate MODEL,CELL artifact keys")
    artifacts = {
        key: load_artifact(key[0], key[1], path)
        for key, path in sorted(artifact_paths.items())
    }
    scores = pd.read_csv(args.scores, sep="\t")
    required = {
        "cell",
        "tf",
        "motif_family",
        "TFBS_chr",
        "label",
        "motif_score",
        "log_accessibility",
        "artifact_index",
        "candidate_probability",
        "bias_configuration",
    }
    missing = required.difference(scores.columns)
    if missing:
        raise SystemExit("score table lacks columns: " + ", ".join(sorted(missing)))
    with threadpool_limits(limits=1):
        metrics = information_ceiling_rows(
            scores,
            artifacts,
            variance_threshold=args.variance_threshold,
            max_components=args.max_components,
            maximum_train_per_tf=args.maximum_train_per_tf,
            seed=args.seed,
        )
    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics_path = args.outdir / "locked_holdout_information_ceiling.tsv"
    metrics.to_csv(metrics_path, sep="\t", index=False)
    manifest = {
        "schema": SCHEMA,
        "diagnostic_only": True,
        "eligible_for_promotion": False,
        "holdout_labels_used": True,
        "pca_training_labels_used": False,
        "classifier": {
            "model": "fixed elastic-net logistic regression",
            "C": 1.0,
            "l1_ratio": 0.5,
            "outer_validation": "leave-one-test-chromosome-out",
        },
        "scores": {"path": str(args.scores), "sha256": file_sha256(args.scores)},
        "artifacts": {
            f"{model}|{cell}": {"path": str(path), "sha256": file_sha256(path)}
            for (model, cell), path in sorted(artifact_paths.items())
        },
        "variance_threshold": float(args.variance_threshold),
        "max_components": int(args.max_components),
        "maximum_train_per_tf": int(args.maximum_train_per_tf),
        "seed": int(args.seed),
        "rows": int(len(metrics)),
        "output": {"path": str(metrics_path), "sha256": file_sha256(metrics_path)},
    }
    (args.outdir / "locked_holdout_information_ceiling_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    columns = [
        "cell",
        "tf",
        "correction",
        "channel",
        "baseline_auroc",
        "combined_auroc",
        "combined_relative_auprc_gain_over_baseline",
        "supervised_auroc_gain_over_dwm",
        "latent_state_inversion_signal",
    ]
    print(metrics[columns].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
