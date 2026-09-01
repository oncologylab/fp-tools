#!/usr/bin/env python3
"""Measure an evaluation-only supervised ceiling for frozen TF profiles."""

from __future__ import annotations

import argparse
from hashlib import blake2b, sha256
import json
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from sklearn.exceptions import ConvergenceWarning

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from build_strand_functional_profiles import site_hashes  # noqa: E402
from diagnose_locked_holdout_information_ceiling import (  # noqa: E402
    cross_chromosome_predictions,
    scored_metrics,
)
from evaluate_functional_footprints import binary_metrics  # noqa: E402
from evaluate_frozen_functional_policy import (  # noqa: E402
    TEST_RESULT_SCHEMA,
    preflight_test_artifact,
    validate_policy,
)
from evaluate_strand_functional_templates import (  # noqa: E402
    file_sha256,
    load_artifact,
)
from fp_tools.tools.functional_footprints import (  # noqa: E402
    FunctionalPCA,
    normalize_functional_profiles,
)


SCHEMA = "fp-tools-frozen-functional-information-ceiling-v2"
CHANNELS = (
    "combined_residual",
    "shared_strand_residual",
    "antisymmetric_strand_residual",
)


def parse_artifact(value: str) -> tuple[str, str, Path]:
    fields = value.split(",", 2)
    if len(fields) != 3 or not all(fields):
        raise argparse.ArgumentTypeError("artifact must use MODEL,CELL,JSON")
    return fields[0], fields[1], Path(fields[2])


def stable_seed(*values: object, seed: int) -> int:
    digest = blake2b(digest_size=8)
    for value in (seed, *values):
        digest.update(str(value).encode())
        digest.update(b"\0")
    return int.from_bytes(digest.digest(), "little") % (2**32 - 1)


def load_site_scores(test_manifest: Path, policy_id: str) -> tuple[pd.DataFrame, dict]:
    document = json.loads(test_manifest.read_text(encoding="utf-8"))
    if document.get("schema") != TEST_RESULT_SCHEMA:
        raise ValueError("unsupported frozen functional test manifest")
    if document.get("policy_id") != policy_id:
        raise ValueError("test manifest and functional policy differ")
    if document.get("models_refitted_on_test") is not False:
        raise ValueError("test manifest does not certify frozen evaluation")
    record = document.get("outputs", {}).get("site_scores", {})
    path = Path(record.get("path", ""))
    if not path.is_file() or file_sha256(path) != record.get("sha256"):
        raise ValueError("test site-score artifact is absent or changed")
    scores = pd.read_csv(path, sep="\t")
    required = {
        "cell",
        "tf",
        "motif_family",
        "bias_configuration",
        "candidate_id",
        "artifact_index",
        "TFBS_chr",
        "motif_score",
        "log_accessibility",
        "label",
        "candidate_probability",
        "dwm_score",
        "raw_score",
        "direct_score",
    }
    missing = required.difference(scores.columns)
    if missing:
        raise ValueError("test site scores lack columns: " + ", ".join(sorted(missing)))
    if scores.duplicated(["cell", "tf", "artifact_index"]).any():
        raise ValueError("test site scores duplicate TF/artifact indexes")
    return scores, document


def failure_classification(
    *,
    supervised_relative_auprc_gain: float,
    label_free_auroc_gain: float,
    label_free_relative_auprc_gain: float,
    supervised_converged: bool = True,
) -> str:
    if not supervised_converged:
        return "supervised_fit_unstable"
    if not np.isfinite(supervised_relative_auprc_gain):
        return "insufficient_supervised_folds"
    if supervised_relative_auprc_gain < 0.10:
        return "assay_limited"
    if label_free_auroc_gain >= 0.03 and label_free_relative_auprc_gain >= 0.10:
        return "detectable"
    return "shape_model_limited"


def raw_guarded_failure_classification(
    *,
    supervised_relative_auprc_gain_over_raw: float,
    signal_panel_relative_auprc_gain_over_raw: float,
    functional_relative_auprc_gain_over_signal_panel: float,
    label_free_auroc_gain_over_raw: float,
    label_free_relative_auprc_gain_over_raw: float,
    supervised_converged: bool = True,
) -> str:
    """Classify a TF only after preserving raw geometry as a guardrail.

    The classifier is an evaluation-only information ceiling.  Consequently,
    an improvement here diagnoses recoverable information; it never promotes
    the supervised model itself.
    """

    values = (
        supervised_relative_auprc_gain_over_raw,
        signal_panel_relative_auprc_gain_over_raw,
        functional_relative_auprc_gain_over_signal_panel,
        label_free_auroc_gain_over_raw,
        label_free_relative_auprc_gain_over_raw,
    )
    if not supervised_converged:
        return "supervised_fit_unstable"
    if not all(np.isfinite(value) for value in values):
        return "insufficient_supervised_folds"
    if (
        label_free_auroc_gain_over_raw >= 0.03
        and label_free_relative_auprc_gain_over_raw >= 0.10
    ):
        return "detectable_above_raw"
    if supervised_relative_auprc_gain_over_raw < 0.10:
        return "assay_limited_relative_to_raw"
    if functional_relative_auprc_gain_over_signal_panel >= 0.10:
        return "shape_model_limited"
    if signal_panel_relative_auprc_gain_over_raw >= 0.10:
        return "signal_combination_limited"
    return "covariate_or_shape_model_limited"


def supervised_predictions(
    features: np.ndarray,
    labels: np.ndarray,
    chromosomes: np.ndarray,
    *,
    seed: int,
    maximum_iterations: int,
) -> tuple[np.ndarray, int, int]:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        predictions, folds = cross_chromosome_predictions(
            features,
            labels,
            chromosomes,
            seed=seed,
            maximum_iterations=maximum_iterations,
        )
    convergence_warnings = sum(
        issubclass(record.category, ConvergenceWarning) for record in caught
    )
    return predictions, folds, convergence_warnings


def relative_auprc_gain(candidate: dict, baseline: dict) -> float:
    candidate_value = float(candidate["auprc"])
    baseline_value = float(baseline["auprc"])
    if not np.isfinite(candidate_value) or not np.isfinite(baseline_value):
        return np.nan
    return (candidate_value - baseline_value) / max(baseline_value, 1e-8)


def task_rows(
    *,
    record: dict,
    training_sites: pd.DataFrame,
    training_profiles: dict[str, np.ndarray],
    test_arrays: dict[str, np.ndarray],
    scores: pd.DataFrame,
    variance_threshold: float,
    max_components: int,
    maximum_train_per_tf: int,
    classifier_max_iter: int,
    seed: int,
) -> list[dict]:
    cell = str(record["cell"])
    tf = str(record["tf"])
    task = scores[
        scores["cell"].astype(str).eq(cell)
        & scores["tf"].astype(str).eq(tf)
    ].reset_index(drop=True)
    if task.empty:
        return []
    if set(task["candidate_id"].astype(str)) != {str(record["candidate"]["candidate_id"])}:
        raise ValueError(f"{cell}/{tf} site scores use the wrong frozen candidate")
    indexes = task["artifact_index"].to_numpy(dtype=int)
    if np.any(indexes < 0) or np.any(indexes >= len(test_arrays["valid"])):
        raise ValueError(f"{cell}/{tf} site scores contain invalid artifact indexes")
    labels = task["label"].to_numpy(dtype=int)
    chromosomes = task["TFBS_chr"].astype(str).to_numpy()
    baseline_features = task[["motif_score", "log_accessibility"]].to_numpy(
        dtype=float
    )
    baseline_predictions, baseline_folds, baseline_warnings = supervised_predictions(
        baseline_features,
        labels,
        chromosomes,
        seed=stable_seed(cell, tf, "baseline", seed=seed),
        maximum_iterations=classifier_max_iter,
    )
    baseline_metrics = scored_metrics(labels, baseline_predictions)
    candidate_score = task["candidate_probability"].to_numpy(dtype=float)
    dwm_score = task["dwm_score"].to_numpy(dtype=float)
    raw_score = task["raw_score"].to_numpy(dtype=float)
    direct_score = task["direct_score"].to_numpy(dtype=float)
    candidate_metrics = binary_metrics(labels, candidate_score)
    dwm_metrics = binary_metrics(labels, dwm_score)
    raw_metrics = binary_metrics(labels, raw_score)
    direct_metrics = binary_metrics(labels, direct_score)
    label_free_auroc_gain = float(candidate_metrics["auroc"] - dwm_metrics["auroc"])
    label_free_relative_auprc_gain = relative_auprc_gain(
        candidate_metrics,
        dwm_metrics,
    )
    label_free_auroc_gain_over_raw = float(
        candidate_metrics["auroc"] - raw_metrics["auroc"]
    )
    label_free_relative_auprc_gain_over_raw = relative_auprc_gain(
        candidate_metrics,
        raw_metrics,
    )
    raw_covariate_features = np.column_stack((baseline_features, raw_score))
    raw_covariate_predictions, raw_covariate_folds, raw_covariate_warnings = (
        supervised_predictions(
            raw_covariate_features,
            labels,
            chromosomes,
            seed=stable_seed(cell, tf, "raw-covariates", seed=seed),
            maximum_iterations=classifier_max_iter,
        )
    )
    raw_covariate_metrics = scored_metrics(labels, raw_covariate_predictions)
    signal_panel_features = np.column_stack(
        (
            baseline_features,
            raw_score,
            dwm_score,
            direct_score,
            candidate_score,
        )
    )
    signal_panel_predictions, signal_panel_folds, signal_panel_warnings = (
        supervised_predictions(
            signal_panel_features,
            labels,
            chromosomes,
            seed=stable_seed(cell, tf, "signal-panel", seed=seed),
            maximum_iterations=classifier_max_iter,
        )
    )
    signal_panel_metrics = scored_metrics(labels, signal_panel_predictions)

    train_mask = training_sites["tf"].astype(str).eq(tf).to_numpy()
    train_indexes = np.flatnonzero(train_mask)
    training_pool = "tf"
    if len(train_indexes) < 100:
        train_mask = training_sites["motif_family"].astype(str).eq(
            str(record["motif_family"])
        ).to_numpy()
        train_indexes = np.flatnonzero(train_mask)
        training_pool = "family"
    if len(train_indexes) < 100:
        return []
    if len(train_indexes) > maximum_train_per_tf:
        rng = np.random.default_rng(
            stable_seed(cell, tf, "train-cap", seed=seed)
        )
        train_indexes = np.sort(
            rng.choice(train_indexes, maximum_train_per_tf, replace=False)
        )
    positions = np.arange(test_arrays["combined_residual"].shape[1], dtype=float)
    positions -= test_arrays["combined_residual"].shape[1] // 2
    train_observed = (
        training_profiles["plus_observed"][train_indexes]
        + training_profiles["minus_observed"][train_indexes]
    )
    weights = np.sqrt(np.maximum(train_observed.sum(axis=1), 1.0))
    rows = []
    for channel in CHANNELS:
        train_values = normalize_functional_profiles(
            training_profiles[channel][train_indexes],
            positions,
        )
        test_values = normalize_functional_profiles(
            test_arrays[channel][indexes],
            positions,
        )
        fpca = FunctionalPCA(
            variance_threshold=variance_threshold,
            max_components=max_components,
            seed=stable_seed(cell, tf, channel, "fpca", seed=seed),
        ).fit(train_values, sample_weight=weights)
        functional_scores = fpca.transform(test_values)
        profile_predictions, profile_folds, profile_warnings = supervised_predictions(
            functional_scores,
            labels,
            chromosomes,
            seed=stable_seed(cell, tf, channel, "profile", seed=seed),
            maximum_iterations=classifier_max_iter,
        )
        combined_predictions, combined_folds, combined_warnings = supervised_predictions(
            np.column_stack((baseline_features, functional_scores)),
            labels,
            chromosomes,
            seed=stable_seed(cell, tf, channel, "combined", seed=seed),
            maximum_iterations=classifier_max_iter,
        )
        full_predictions, full_folds, full_warnings = supervised_predictions(
            np.column_stack((signal_panel_features, functional_scores)),
            labels,
            chromosomes,
            seed=stable_seed(cell, tf, channel, "full", seed=seed),
            maximum_iterations=classifier_max_iter,
        )
        profile_metrics = scored_metrics(labels, profile_predictions)
        combined_metrics = scored_metrics(labels, combined_predictions)
        full_metrics = scored_metrics(labels, full_predictions)
        relative_gain = relative_auprc_gain(combined_metrics, baseline_metrics)
        full_gain_over_raw = relative_auprc_gain(full_metrics, raw_metrics)
        full_gain_over_raw_covariates = relative_auprc_gain(
            full_metrics,
            raw_covariate_metrics,
        )
        full_gain_over_signal_panel = relative_auprc_gain(
            full_metrics,
            signal_panel_metrics,
        )
        signal_panel_gain_over_raw = relative_auprc_gain(
            signal_panel_metrics,
            raw_metrics,
        )
        rows.append(
            {
                "cell": cell,
                "tf": tf,
                "motif_family": str(record["motif_family"]),
                "candidate_id": str(record["candidate"]["candidate_id"]),
                "channel": channel,
                "test_sites": int(len(task)),
                "positive_sites": int(np.sum(labels == 1)),
                "negative_sites": int(np.sum(labels == 0)),
                "training_pool": training_pool,
                "training_profiles": int(len(train_indexes)),
                "functional_components": int(len(fpca.components_)),
                "baseline_folds": int(baseline_folds),
                "profile_folds": int(profile_folds),
                "combined_folds": int(combined_folds),
                "raw_covariate_folds": int(raw_covariate_folds),
                "signal_panel_folds": int(signal_panel_folds),
                "full_folds": int(full_folds),
                "baseline_convergence_warnings": int(baseline_warnings),
                "profile_convergence_warnings": int(profile_warnings),
                "combined_convergence_warnings": int(combined_warnings),
                "raw_covariate_convergence_warnings": int(raw_covariate_warnings),
                "signal_panel_convergence_warnings": int(signal_panel_warnings),
                "full_convergence_warnings": int(full_warnings),
                **{f"baseline_{key}": value for key, value in baseline_metrics.items()},
                **{f"profile_{key}": value for key, value in profile_metrics.items()},
                **{f"combined_{key}": value for key, value in combined_metrics.items()},
                **{
                    f"raw_covariate_{key}": value
                    for key, value in raw_covariate_metrics.items()
                },
                **{
                    f"signal_panel_{key}": value
                    for key, value in signal_panel_metrics.items()
                },
                **{f"full_{key}": value for key, value in full_metrics.items()},
                "combined_relative_auprc_gain_over_baseline": relative_gain,
                "shape_information_above_baseline": bool(relative_gain >= 0.10),
                "full_relative_auprc_gain_over_raw": full_gain_over_raw,
                "full_relative_auprc_gain_over_raw_covariates": (
                    full_gain_over_raw_covariates
                ),
                "full_relative_auprc_gain_over_signal_panel": (
                    full_gain_over_signal_panel
                ),
                "signal_panel_relative_auprc_gain_over_raw": (
                    signal_panel_gain_over_raw
                ),
                "functional_information_above_signal_panel": bool(
                    full_gain_over_signal_panel >= 0.10
                ),
                "label_free_auroc": float(candidate_metrics["auroc"]),
                "label_free_auprc": float(candidate_metrics["auprc"]),
                "dwm_auroc": float(dwm_metrics["auroc"]),
                "dwm_auprc": float(dwm_metrics["auprc"]),
                "raw_auroc": float(raw_metrics["auroc"]),
                "raw_auprc": float(raw_metrics["auprc"]),
                "direct_auroc": float(direct_metrics["auroc"]),
                "direct_auprc": float(direct_metrics["auprc"]),
                "label_free_auroc_gain_over_dwm": label_free_auroc_gain,
                "label_free_relative_auprc_gain_over_dwm": (
                    label_free_relative_auprc_gain
                ),
                "label_free_auroc_gain_over_raw": label_free_auroc_gain_over_raw,
                "label_free_relative_auprc_gain_over_raw": (
                    label_free_relative_auprc_gain_over_raw
                ),
                "diagnostic_only": True,
            }
        )
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--test-manifest", type=Path, required=True)
    parser.add_argument(
        "--training-artifact",
        action="append",
        type=parse_artifact,
        required=True,
        metavar="MODEL,CELL,JSON",
    )
    parser.add_argument(
        "--test-artifact",
        action="append",
        type=parse_artifact,
        required=True,
        metavar="MODEL,CELL,JSON",
    )
    parser.add_argument("--variance-threshold", type=float, default=0.95)
    parser.add_argument("--max-components", type=int, default=20)
    parser.add_argument("--maximum-train-per-tf", type=int, default=5000)
    parser.add_argument("--classifier-max-iter", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    policy, models = validate_policy(args.policy)
    scores, test_manifest = load_site_scores(
        args.test_manifest,
        str(policy["policy_id"]),
    )
    study_path = Path(policy["study"]["path"])
    study = json.loads(study_path.read_text(encoding="utf-8"))
    training_paths = {
        (model, cell): path for model, cell, path in args.training_artifact
    }
    test_paths = {(model, cell): path for model, cell, path in args.test_artifact}
    policy_keys = {
        (str(record["bias_configuration"]), str(record["cell"]))
        for record, _candidate, _model in models
    }
    if set(training_paths) != policy_keys or set(test_paths) != policy_keys:
        raise ValueError("training/test artifacts do not exactly match policy keys")

    training = {}
    tests = {}
    input_records = []
    for key in sorted(policy_keys):
        training[key] = load_artifact(training_paths[key], key[1], study)[:2]
        test_document, test_arrays = preflight_test_artifact(
            test_paths[key],
            expected_cell=key[1],
        )
        test_sites = pd.read_csv(test_document["sites"], sep="\t").reset_index(drop=True)
        if not np.array_equal(site_hashes(test_sites), test_arrays["site_hash"]):
            raise ValueError(f"test site order mismatch: {test_paths[key]}")
        tests[key] = test_arrays
        input_records.extend(
            [
                {
                    "purpose": "label-free-training",
                    "path": str(training_paths[key]),
                    "sha256": file_sha256(training_paths[key]),
                },
                {
                    "purpose": "test-profiles",
                    "path": str(test_paths[key]),
                    "sha256": file_sha256(test_paths[key]),
                },
            ]
        )

    rows = []
    with threadpool_limits(limits=1):
        for record, _candidate, _model in models:
            key = (str(record["bias_configuration"]), str(record["cell"]))
            training_sites, training_profiles = training[key]
            rows.extend(
                task_rows(
                    record=record,
                    training_sites=training_sites,
                    training_profiles=training_profiles,
                    test_arrays=tests[key],
                    scores=scores,
                    variance_threshold=args.variance_threshold,
                    max_components=args.max_components,
                    maximum_train_per_tf=args.maximum_train_per_tf,
                    classifier_max_iter=args.classifier_max_iter,
                    seed=args.seed,
                )
            )
    results = pd.DataFrame(rows)
    if results.empty:
        raise ValueError("no TF information-ceiling task was evaluable")
    best = (
        results.sort_values(
            ["cell", "tf", "full_auprc", "channel"],
            ascending=[True, True, False, True],
        )
        .groupby(["cell", "tf"], sort=True, as_index=False)
        .first()
    )
    best["failure_classification"] = [
        failure_classification(
            supervised_relative_auprc_gain=float(supervised),
            label_free_auroc_gain=float(auroc),
            label_free_relative_auprc_gain=float(auprc),
            supervised_converged=(int(baseline_warnings) + int(combined_warnings))
            == 0,
        )
        for supervised, auroc, auprc, baseline_warnings, combined_warnings in zip(
            best["combined_relative_auprc_gain_over_baseline"],
            best["label_free_auroc_gain_over_dwm"],
            best["label_free_relative_auprc_gain_over_dwm"],
            best["baseline_convergence_warnings"],
            best["combined_convergence_warnings"],
            strict=True,
        )
    ]
    best["raw_guarded_failure_classification"] = [
        raw_guarded_failure_classification(
            supervised_relative_auprc_gain_over_raw=float(supervised),
            signal_panel_relative_auprc_gain_over_raw=float(signal_panel),
            functional_relative_auprc_gain_over_signal_panel=float(functional),
            label_free_auroc_gain_over_raw=float(auroc),
            label_free_relative_auprc_gain_over_raw=float(auprc),
            supervised_converged=(
                int(raw_covariate_warnings)
                + int(signal_panel_warnings)
                + int(full_warnings)
            )
            == 0,
        )
        for (
            supervised,
            signal_panel,
            functional,
            auroc,
            auprc,
            raw_covariate_warnings,
            signal_panel_warnings,
            full_warnings,
        ) in zip(
            best["full_relative_auprc_gain_over_raw"],
            best["signal_panel_relative_auprc_gain_over_raw"],
            best["full_relative_auprc_gain_over_signal_panel"],
            best["label_free_auroc_gain_over_raw"],
            best["label_free_relative_auprc_gain_over_raw"],
            best["raw_covariate_convergence_warnings"],
            best["signal_panel_convergence_warnings"],
            best["full_convergence_warnings"],
            strict=True,
        )
    ]

    args.outdir.mkdir(parents=True, exist_ok=True)
    results_path = args.outdir / "functional_information_ceiling.tsv"
    best_path = args.outdir / "functional_information_ceiling_best.tsv"
    results.to_csv(results_path, sep="\t", index=False)
    best.to_csv(best_path, sep="\t", index=False)
    manifest = {
        "schema": SCHEMA,
        "policy_id": policy["policy_id"],
        "diagnostic_only": True,
        "used_for_model_selection": False,
        "functional_basis_training_labels_used": False,
        "classifier_test_labels_used": True,
        "raw_guardrail_required": True,
        "supervised_feature_sets": {
            "baseline": ["motif_score", "log_accessibility"],
            "raw_covariates": [
                "motif_score",
                "log_accessibility",
                "raw_score",
            ],
            "signal_panel": [
                "motif_score",
                "log_accessibility",
                "raw_score",
                "dwm_score",
                "direct_score",
                "candidate_probability",
            ],
            "full": ["signal_panel", "label_free_functional_pc_scores"],
        },
        "variance_threshold": args.variance_threshold,
        "max_components": args.max_components,
        "maximum_train_per_tf": args.maximum_train_per_tf,
        "classifier_max_iter": args.classifier_max_iter,
        "seed": args.seed,
        "inputs": [
            {"path": str(args.policy), "sha256": file_sha256(args.policy)},
            {
                "path": str(args.test_manifest),
                "sha256": file_sha256(args.test_manifest),
            },
            *input_records,
        ],
        "source_test_input_id": test_manifest["test_input_id"],
        "outputs": {
            "all_channels": {
                "path": str(results_path),
                "sha256": file_sha256(results_path),
            },
            "best_diagnostic_channel": {
                "path": str(best_path),
                "sha256": file_sha256(best_path),
            },
        },
    }
    canonical = json.dumps(manifest, sort_keys=True, separators=(",", ":"))
    manifest["ceiling_id"] = sha256(canonical.encode()).hexdigest()
    manifest_path = args.outdir / "functional_information_ceiling_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        best[
            [
                "cell",
                "tf",
                "channel",
                "full_relative_auprc_gain_over_raw",
                "full_relative_auprc_gain_over_signal_panel",
                "label_free_relative_auprc_gain_over_raw",
                "raw_guarded_failure_classification",
            ]
        ].to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
