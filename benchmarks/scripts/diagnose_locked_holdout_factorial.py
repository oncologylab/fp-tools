#!/usr/bin/env python3
"""Run an explicitly post-lock detector factorial on opened holdout labels.

This diagnostic answers whether a failed frozen route was caused by its bias
configuration, detector, or latent-state orientation.  It must never be used
as promotion evidence: ChIP labels are already open and every reported winner
is selected post hoc.  A promising configuration has to be frozen and tested
on a new, untouched holdout before package integration.
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
from threadpoolctl import threadpool_limits

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import binary_metrics, stable_seed  # noqa: E402
from evaluate_locked_holdout_policy import (  # noqa: E402
    Artifact,
    combined_candidates,
    fit_combined_route,
    fit_strand_route,
    load_artifact,
    reference_to_candidate_indexes,
    strand_candidates,
)
from fp_tools.tools.functional_footprints import (  # noqa: E402
    profile_descriptors,
    standardized_functional_separation,
)


SCHEMA = "fp-tools-locked-holdout-posthoc-factorial-v1"
SUPPORTED_COMBINED_FAMILIES = {"spline", "gp", "fda", "hybrid"}


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


def probability_metrics(labels: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    values = binary_metrics(
        np.asarray(labels, dtype=int),
        np.asarray(probabilities, dtype=float),
    )
    return {
        key: float(values[key])
        for key in ("auroc", "auprc", "brier")
    }


def orientation_diagnostic(
    labels: np.ndarray,
    probabilities: np.ndarray,
) -> dict[str, float | bool]:
    direct = probability_metrics(labels, probabilities)
    flipped = probability_metrics(labels, 1.0 - np.asarray(probabilities, dtype=float))
    inverted = bool(flipped["auroc"] > direct["auroc"])
    selected = flipped if inverted else direct
    return {
        "flipped_auroc": flipped["auroc"],
        "flipped_auprc": flipped["auprc"],
        "orientation_inverted_posthoc": inverted,
        "oracle_orientation_auroc": selected["auroc"],
        "oracle_orientation_auprc": selected["auprc"],
    }


def supported_candidates(artifact: Artifact) -> list[str]:
    if artifact.schema == "fp-tools-strand-functional-profiles-v1":
        return sorted(strand_candidates())
    return sorted(
        candidate_id
        for candidate_id, candidate in combined_candidates().items()
        if candidate.family in SUPPORTED_COMBINED_FAMILIES
    )


def artifact_positions(artifact: Artifact) -> np.ndarray:
    width = int(artifact.profiles["combined_residual"].shape[1])
    if width % 2 != 1:
        raise ValueError(f"{artifact.path} profile width must be odd")
    flank = width // 2
    return np.arange(-flank, flank + 1, dtype=float)


def _functional_metrics(
    profiles: np.ndarray,
    labels: np.ndarray,
    positions: np.ndarray,
) -> dict[str, float]:
    positive = np.mean(profiles[labels == 1], axis=0)
    negative = np.mean(profiles[labels == 0], axis=0)
    positive_descriptor = profile_descriptors(positive, positions)
    negative_descriptor = profile_descriptors(negative, positions)
    return {
        "functional_separation": float(
            standardized_functional_separation(profiles, labels, positions)
        ),
        "positive_depletion": float(positive_descriptor.depletion),
        "negative_depletion": float(negative_descriptor.depletion),
        "depletion_difference": float(
            positive_descriptor.depletion - negative_descriptor.depletion
        ),
    }


def select_posthoc_winners(metrics: pd.DataFrame) -> pd.DataFrame:
    passing = metrics[
        metrics["status"].eq("ok") & metrics["converged"].astype(bool)
    ].copy()
    if passing.empty:
        return passing
    keys = ["cell", "tf", "correction"]
    return (
        passing.sort_values(
            keys + ["auroc_gain", "relative_auprc_gain", "candidate_id"],
            ascending=[True, True, True, False, False, True],
            kind="mergesort",
        )
        .groupby(keys, as_index=False, sort=True)
        .head(1)
        .assign(
            selection_status="posthoc_diagnostic_only",
            eligible_for_promotion=False,
        )
        .reset_index(drop=True)
    )


def run_factorial(
    scores: pd.DataFrame,
    artifacts: dict[tuple[str, str], Artifact],
    *,
    corrections: set[str] | None,
    maximum_train_per_tf: int,
    seed: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    combined_cache: dict[tuple[str, str, str], tuple[np.ndarray | None, float]] = {}
    strand_cache: dict[tuple[str, str, str], tuple[np.ndarray, float]] = {}
    for (correction, cell), artifact in sorted(artifacts.items()):
        if corrections is not None and correction not in corrections:
            continue
        if correction == "DWM":
            reference = artifact
        else:
            reference = artifacts.get(("DWM", cell))
            if reference is None:
                raise ValueError(f"{correction}/{cell} lacks a DWM reference artifact")
        mapping = reference_to_candidate_indexes(reference, artifact)
        positions = artifact_positions(artifact)
        cell_scores = scores[scores["cell"].astype(str).eq(cell)]
        for tf, task in cell_scores.groupby("tf", sort=True):
            task = task.reset_index(drop=True)
            reference_indexes = task["artifact_index"].to_numpy(dtype=int)
            evaluation_indexes = mapping[reference_indexes]
            if np.any(evaluation_indexes < 0):
                continue
            labels = task["label"].to_numpy(dtype=int)
            family = str(task["motif_family"].iloc[0])
            reference_probabilities = task["reference_probability"].to_numpy(dtype=float)
            reference_metrics = probability_metrics(labels, reference_probabilities)
            for candidate_id in supported_candidates(artifact):
                base = {
                    "cell": cell,
                    "tf": str(tf),
                    "motif_family": family,
                    "correction": correction,
                    "candidate_id": candidate_id,
                    "sites": int(len(task)),
                    "positive_sites": int(np.sum(labels == 1)),
                    "negative_sites": int(np.sum(labels == 0)),
                    "training_labels_used": False,
                    "holdout_labels_used_for_evaluation": True,
                    "posthoc_after_locked_evaluation": True,
                    "eligible_for_promotion": False,
                }
                try:
                    if artifact.schema == "fp-tools-strand-functional-profiles-v1":
                        output = fit_strand_route(
                            artifact,
                            candidate_id,
                            str(tf),
                            family,
                            evaluation_indexes,
                            positions,
                            maximum_train_per_tf,
                            stable_seed(correction, cell, tf, seed=seed),
                            strand_cache,
                        )
                    else:
                        output = fit_combined_route(
                            artifact,
                            candidate_id,
                            str(tf),
                            family,
                            evaluation_indexes,
                            positions,
                            maximum_train_per_tf,
                            stable_seed(correction, cell, tf, seed=seed),
                            combined_cache,
                        )
                    metrics = probability_metrics(labels, output.probabilities)
                    relative_auprc = (
                        metrics["auprc"] - reference_metrics["auprc"]
                    ) / max(reference_metrics["auprc"], 1e-8)
                    descriptors = profile_descriptors(output.footprint_profile, positions)
                    rows.append(
                        {
                            **base,
                            "status": "ok",
                            "converged": bool(output.converged),
                            "iterations": int(output.iterations),
                            "fit_seconds": float(output.fit_seconds),
                            "dispersion": float(output.dispersion),
                            **metrics,
                            "reference_auroc": reference_metrics["auroc"],
                            "reference_auprc": reference_metrics["auprc"],
                            "reference_brier": reference_metrics["brier"],
                            "auroc_gain": metrics["auroc"] - reference_metrics["auroc"],
                            "relative_auprc_gain": relative_auprc,
                            **orientation_diagnostic(labels, output.probabilities),
                            **_functional_metrics(
                                output.evaluation_profiles,
                                labels,
                                positions,
                            ),
                            **{
                                f"model_profile_{key}": float(value)
                                for key, value in vars(descriptors).items()
                            },
                        }
                    )
                except Exception as error:
                    rows.append(
                        {
                            **base,
                            "status": "error",
                            "error": f"{type(error).__name__}: {error}",
                            "converged": False,
                        }
                    )
    return pd.DataFrame(rows).sort_values(
        ["cell", "tf", "correction", "candidate_id"], kind="mergesort"
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
    parser.add_argument("--maximum-train-per-tf", type=int, default=10000)
    parser.add_argument(
        "--evaluate-correction",
        action="append",
        help="Evaluate only this correction; repeat as needed (DWM may still be loaded as reference)",
    )
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    paths = {(model, cell): path for model, cell, path in args.artifact}
    if len(paths) != len(args.artifact):
        raise SystemExit("duplicate MODEL,CELL artifact keys")
    artifacts = {
        key: load_artifact(key[0], key[1], path)
        for key, path in sorted(paths.items())
    }
    scores = pd.read_csv(args.scores, sep="\t")
    required = {
        "cell",
        "tf",
        "motif_family",
        "artifact_index",
        "label",
        "reference_probability",
    }
    missing = required.difference(scores.columns)
    if missing:
        raise SystemExit("score table lacks columns: " + ", ".join(sorted(missing)))
    with threadpool_limits(limits=1):
        metrics = run_factorial(
            scores,
            artifacts,
            corrections=(
                None
                if args.evaluate_correction is None
                else set(args.evaluate_correction)
            ),
            maximum_train_per_tf=args.maximum_train_per_tf,
            seed=args.seed,
        )
    winners = select_posthoc_winners(metrics)
    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics_path = args.outdir / "locked_holdout_posthoc_factorial.tsv.gz"
    winners_path = args.outdir / "locked_holdout_posthoc_winners.tsv"
    metrics.to_csv(metrics_path, sep="\t", index=False, compression="gzip")
    winners.to_csv(winners_path, sep="\t", index=False)
    manifest = {
        "schema": SCHEMA,
        "diagnostic_only": True,
        "posthoc_after_locked_evaluation": True,
        "eligible_for_promotion": False,
        "selection_must_be_retested_on_new_untouched_holdout": True,
        "scores": {"path": str(args.scores), "sha256": file_sha256(args.scores)},
        "artifacts": {
            f"{model}|{cell}": {"path": str(path), "sha256": file_sha256(path)}
            for (model, cell), path in sorted(paths.items())
        },
        "maximum_train_per_tf": int(args.maximum_train_per_tf),
        "seed": int(args.seed),
        "metric_rows": int(len(metrics)),
        "winner_rows": int(len(winners)),
        "outputs": {
            "metrics": {"path": str(metrics_path), "sha256": file_sha256(metrics_path)},
            "winners": {"path": str(winners_path), "sha256": file_sha256(winners_path)},
        },
    }
    (args.outdir / "locked_holdout_posthoc_factorial_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    columns = [
        "cell",
        "tf",
        "correction",
        "candidate_id",
        "auroc",
        "reference_auroc",
        "auroc_gain",
        "relative_auprc_gain",
        "orientation_inverted_posthoc",
    ]
    print(winners[columns].to_string(index=False) if len(winners) else "No successful fits")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
