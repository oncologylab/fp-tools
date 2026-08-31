#!/usr/bin/env python3
"""Select per-TF count mixtures using only held-out unlabeled chromosomes.

Each candidate is fitted on the first portion of the development training
chromosomes and ranked by its marginal likelihood gain on the remaining
training chromosomes.  ChIP labels are opened only after that choice is
frozen, to measure whether the unsupervised selector transfers to occupancy.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
from time import perf_counter

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import binary_metrics  # noqa: E402
from evaluate_functional_template_transfer import selection_score  # noqa: E402
from evaluate_strand_functional_templates import load_artifact  # noqa: E402
from evaluate_strand_label_free_models import (  # noqa: E402
    Candidate,
    _count_profiles,
    candidate_grid,
    parse_artifact,
    validate_unlabeled_training_sites,
)
from fp_tools.tools.functional_footprints import BiasAwareFunctionalMixture  # noqa: E402
from fp_tools.tools.parametric_bias import estimate_nb_dispersion  # noqa: E402


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def log_marginal_gain(
    shape_log_likelihood_ratio: np.ndarray,
    prior_log_odds: np.ndarray,
) -> np.ndarray:
    """Log mixture likelihood minus the unbound-state likelihood."""

    shape = np.asarray(shape_log_likelihood_ratio, dtype=float)
    prior = np.asarray(prior_log_odds, dtype=float)
    if shape.shape != prior.shape:
        raise ValueError("shape and prior log odds must agree")
    log_unbound_weight = -np.logaddexp(0.0, prior)
    log_bound_weight = -np.logaddexp(0.0, -prior)
    return np.logaddexp(log_unbound_weight, shape + log_bound_weight)


def make_model(
    positions: np.ndarray,
    candidate: Candidate,
    dispersion: float,
    *,
    shrinkage: float,
) -> BiasAwareFunctionalMixture:
    return BiasAwareFunctionalMixture(
        positions,
        smoother=candidate.smoother,
        dispersion=dispersion,
        max_iter=60,
        shrinkage=shrinkage,
        accessibility_background=candidate.background,
        prior_constraint="none",
        profile_outer_limit=50.0,
        likelihood_limit=candidate.window,
    )


def chromosome_masks(
    sites: pd.DataFrame,
    study: dict,
    tune_chromosome_count: int,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    chromosomes = list(study["chromosome_split"]["train"])
    if tune_chromosome_count < 1 or tune_chromosome_count >= len(chromosomes):
        raise ValueError("tune chromosome count must leave fit and tune chromosomes")
    tune_chromosomes = chromosomes[-tune_chromosome_count:]
    observed = sites["TFBS_chr"].astype(str)
    fit = observed.isin(chromosomes[:-tune_chromosome_count]).to_numpy()
    tune = observed.isin(tune_chromosomes).to_numpy()
    if np.any(fit & tune):
        raise AssertionError("fit and tune chromosome masks overlap")
    return fit, tune, tune_chromosomes


def fit_family_prior(
    candidate: Candidate,
    observed: np.ndarray,
    expected: np.ndarray,
    indexes: np.ndarray,
    positions: np.ndarray,
) -> dict:
    started = perf_counter()
    dispersion = estimate_nb_dispersion(observed[indexes], expected[indexes])
    model = make_model(positions, candidate, dispersion, shrinkage=0.0)
    result = model.fit(observed[indexes], expected[indexes])
    return {
        "profile": result.footprint_profile,
        "dispersion": dispersion,
        "converged": bool(result.converged),
        "sites": int(len(indexes)),
        "seconds": perf_counter() - started,
    }


def evaluate_cv_candidate(
    candidate: Candidate,
    *,
    bias_configuration: str,
    cell: str,
    tf: str,
    motif_family: str,
    sites: pd.DataFrame,
    observed: np.ndarray,
    expected: np.ndarray,
    fit_mask: np.ndarray,
    tune_mask: np.ndarray,
    family_prior: dict,
    positions: np.ndarray,
    minimum_fit_sites: int,
    minimum_tune_sites: int,
) -> dict:
    started = perf_counter()
    tf_values = sites["tf"].astype(str).to_numpy()
    fit_indexes = np.flatnonzero(fit_mask & (tf_values == tf))
    tune_indexes = np.flatnonzero(tune_mask & (tf_values == tf))
    base = {
        "cell": cell,
        "tf": tf,
        "motif_family": motif_family,
        "bias_configuration": bias_configuration,
        **asdict(candidate),
        "fit_sites": int(len(fit_indexes)),
        "tune_sites": int(len(tune_indexes)),
        "selection_labels_used": False,
        "motif_or_accessibility_features_used": False,
    }
    if len(fit_indexes) < minimum_fit_sites or len(tune_indexes) < minimum_tune_sites:
        return {**base, "status": "insufficient_sites"}
    try:
        model = make_model(
            positions,
            candidate,
            float(family_prior["dispersion"]),
            shrinkage=50.0,
        )
        result = model.fit(
            observed[fit_indexes],
            expected[fit_indexes],
            prior_profile=np.asarray(family_prior["profile"], dtype=float),
        )
        shape, prior = model.predict_log_odds_components(
            observed[tune_indexes],
            expected[tune_indexes],
        )
        gain = log_marginal_gain(shape, prior)
        posterior = model.predict(observed[tune_indexes], expected[tune_indexes])
        descriptors = result.descriptors
        return {
            **base,
            "status": "ok",
            "converged": bool(result.converged),
            "iterations": int(result.iterations),
            "heldout_log_marginal_gain": float(np.sum(gain)),
            "heldout_gain_per_site": float(np.mean(gain)),
            "heldout_gain_per_site_position": float(
                np.mean(gain) / max(int(model.likelihood_mask.sum()), 1)
            ),
            "heldout_bound_fraction": float(np.mean(posterior)),
            "heldout_posterior_entropy": float(
                np.mean(
                    -posterior * np.log(np.maximum(posterior, 1e-12))
                    - (1.0 - posterior) * np.log(np.maximum(1.0 - posterior, 1e-12))
                )
            ),
            "profile_depletion": float(descriptors.depletion),
            "profile_width": float(descriptors.width),
            "profile_asymmetry": float(descriptors.asymmetry),
            "profile_periodicity": float(descriptors.periodicity),
            "profile_plausible": bool(
                descriptors.depletion > 0
                and 1 <= descriptors.width <= 100
                and 0.01 <= np.mean(posterior) <= 0.99
            ),
            "fit_seconds": perf_counter() - started,
        }
    except Exception as error:
        return {
            **base,
            "status": "error",
            "error": f"{type(error).__name__}: {error}",
            "fit_seconds": perf_counter() - started,
        }


def select_candidates(cv_metrics: pd.DataFrame) -> pd.DataFrame:
    passing = cv_metrics[
        cv_metrics["status"].eq("ok")
        & cv_metrics["converged"].astype(bool)
        & cv_metrics["profile_plausible"].astype(bool)
    ].copy()
    keys = ["cell", "tf"]
    if passing.empty:
        return passing
    return (
        passing.sort_values(
            keys
            + [
                "heldout_gain_per_site_position",
                "heldout_gain_per_site",
                "fit_seconds",
                "bias_configuration",
                "candidate_id",
            ],
            ascending=[True, True, False, False, True, True, True],
            kind="mergesort",
        )
        .groupby(keys, as_index=False, sort=True)
        .head(1)
        .reset_index(drop=True)
    )


def refit_and_score(
    selected: pd.Series,
    *,
    training_sites: pd.DataFrame,
    training_profiles: dict[str, np.ndarray],
    evaluation_sites: pd.DataFrame,
    evaluation_profiles: dict[str, np.ndarray],
    positions: np.ndarray,
    outdir: Path,
    minimum_evaluation_sites: int,
) -> tuple[dict, pd.DataFrame | None]:
    candidate = Candidate(
        candidate_id=str(selected["candidate_id"]),
        family="count",
        smoother=str(selected["smoother"]),
        background=str(selected["background"]),
        window=float(selected["window"]),
    )
    tf = str(selected["tf"])
    motif_family = str(selected["motif_family"])
    train_tf = training_sites["tf"].astype(str).eq(tf).to_numpy()
    train_family = training_sites["motif_family"].astype(str).eq(motif_family).to_numpy()
    validation = (
        evaluation_sites["tf"].astype(str).eq(tf)
        & evaluation_sites["chromosome_split"].eq("validation")
    ).to_numpy()
    labels = evaluation_sites.loc[validation, "chip_label"].to_numpy(dtype=int)
    base = {
        **selected.to_dict(),
        "validation_sites": int(np.sum(validation)),
        "validation_positive_sites": int(np.sum(labels == 1)),
        "validation_negative_sites": int(np.sum(labels == 0)),
        "validation_labels_used_for_selection": False,
    }
    if min(np.sum(labels == 0), np.sum(labels == 1)) < minimum_evaluation_sites:
        return {**base, "validation_status": "insufficient_sites"}, None
    observed, expected = _count_profiles(training_profiles)
    evaluation_observed, evaluation_expected = _count_profiles(evaluation_profiles)
    dispersion = estimate_nb_dispersion(observed[train_family], expected[train_family])
    family_model = make_model(positions, candidate, dispersion, shrinkage=0.0)
    family_result = family_model.fit(observed[train_family], expected[train_family])
    model = make_model(positions, candidate, dispersion, shrinkage=50.0)
    result = model.fit(
        observed[train_tf],
        expected[train_tf],
        prior_profile=family_result.footprint_profile,
    )
    probabilities = model.predict(
        evaluation_observed[validation],
        evaluation_expected[validation],
    )
    metrics = binary_metrics(labels, probabilities)
    safe_name = re.sub(r"[^A-Za-z0-9_-]+", "_", f"{selected['cell']}__{tf}__{selected['bias_configuration']}__{candidate.candidate_id}")
    model_prefix = outdir / "models" / safe_name
    npz_path, json_path = model.save(
        model_prefix,
        metadata={
            "selection": "unlabeled held-out chromosome marginal likelihood",
            "cell": str(selected["cell"]),
            "tf": tf,
            "motif_family": motif_family,
            "bias_configuration": str(selected["bias_configuration"]),
        },
    )
    row = {
        **base,
        "validation_status": "ok",
        "refit_converged": bool(result.converged),
        "refit_iterations": int(result.iterations),
        "model_npz": str(npz_path),
        "model_json": str(json_path),
        "selection_score_after_opening_labels": selection_score(metrics),
        **metrics,
    }
    scores = pd.DataFrame(
        {
            "cell": selected["cell"],
            "tf": tf,
            "bias_configuration": selected["bias_configuration"],
            "candidate_id": candidate.candidate_id,
            "site_index": np.flatnonzero(validation),
            "chip_label": labels,
            "binding_probability": probabilities,
        }
    )
    return row, scores


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-artifact", action="append", type=parse_artifact, required=True, metavar="MODEL,CELL,JSON")
    parser.add_argument("--evaluation-artifact", action="append", type=parse_artifact, required=True, metavar="MODEL,CELL,JSON")
    parser.add_argument(
        "--study",
        type=Path,
        default=Path("benchmarks/manifests/footprint_functional_v1.spec.json"),
    )
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tune-chromosome-count", type=int, default=3)
    parser.add_argument("--minimum-fit-sites", type=int, default=75)
    parser.add_argument("--minimum-tune-sites", type=int, default=20)
    parser.add_argument("--minimum-evaluation-sites", type=int, default=25)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args(argv)
    study = json.loads(args.study.read_text(encoding="utf-8"))
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[tasks["split"].eq("development")]
    training_paths = {(model, cell): path for model, cell, path in args.training_artifact}
    evaluation_paths = {(model, cell): path for model, cell, path in args.evaluation_artifact}
    if set(training_paths) != set(evaluation_paths):
        raise SystemExit("training and evaluation artifact keys must match")
    artifacts = {}
    input_rows = []
    for key in sorted(training_paths):
        model, cell = key
        train_sites, train_profiles, train_document = load_artifact(training_paths[key], cell, study)
        validate_unlabeled_training_sites(train_sites, training_paths[key])
        eval_sites, eval_profiles, eval_document = load_artifact(evaluation_paths[key], cell, study)
        artifacts[key] = (train_sites, train_profiles, eval_sites, eval_profiles)
        for purpose, path, document in (
            ("training", training_paths[key], train_document),
            ("evaluation", evaluation_paths[key], eval_document),
        ):
            input_rows.append(
                {
                    "bias_configuration": model,
                    "cell": cell,
                    "purpose": purpose,
                    "path": str(path),
                    "sha256": file_sha256(path),
                    "sites": int(document["sites_total"]),
                }
            )
    candidates = [candidate for candidate in candidate_grid() if candidate.family == "count"]
    positions = np.arange(
        -int(study["profile_flank_bp"]),
        int(study["profile_flank_bp"]) + 1,
        dtype=float,
    )
    family_priors = {}
    family_futures = {}
    with threadpool_limits(limits=1), ThreadPoolExecutor(max_workers=args.workers) as executor:
        for (bias_configuration, cell), artifact in artifacts.items():
            sites, profiles, _eval_sites, _eval_profiles = artifact
            observed, expected = _count_profiles(profiles)
            fit_mask, _tune_mask, _tune_chromosomes = chromosome_masks(
                sites, study, args.tune_chromosome_count
            )
            cell_tasks = tasks[tasks["cell"].astype(str).eq(cell)]
            for motif_family in sorted(cell_tasks["motif_family"].astype(str).unique()):
                family_indexes = np.flatnonzero(
                    fit_mask
                    & sites["motif_family"].astype(str).eq(motif_family).to_numpy()
                )
                for candidate in candidates:
                    key = (bias_configuration, cell, motif_family, candidate.candidate_id)
                    family_futures[
                        executor.submit(
                            fit_family_prior,
                            candidate,
                            observed,
                            expected,
                            family_indexes,
                            positions,
                        )
                    ] = key
        for future in as_completed(family_futures):
            family_priors[family_futures[future]] = future.result()

    rows = []
    futures = {}
    tune_chromosomes = None
    with threadpool_limits(limits=1), ThreadPoolExecutor(max_workers=args.workers) as executor:
        for (bias_configuration, cell), artifact in artifacts.items():
            sites, profiles, _eval_sites, _eval_profiles = artifact
            observed, expected = _count_profiles(profiles)
            fit_mask, tune_mask, tune_chromosomes = chromosome_masks(
                sites, study, args.tune_chromosome_count
            )
            for task in tasks[tasks["cell"].astype(str).eq(cell)].itertuples(index=False):
                for candidate in candidates:
                    family_key = (
                        bias_configuration,
                        cell,
                        str(task.motif_family),
                        candidate.candidate_id,
                    )
                    future = executor.submit(
                        evaluate_cv_candidate,
                        candidate,
                        bias_configuration=bias_configuration,
                        cell=cell,
                        tf=str(task.tf),
                        motif_family=str(task.motif_family),
                        sites=sites,
                        observed=observed,
                        expected=expected,
                        fit_mask=fit_mask,
                        tune_mask=tune_mask,
                        family_prior=family_priors[family_key],
                        positions=positions,
                        minimum_fit_sites=args.minimum_fit_sites,
                        minimum_tune_sites=args.minimum_tune_sites,
                    )
                    futures[future] = (bias_configuration, cell, str(task.tf), candidate.candidate_id)
        for future in as_completed(futures):
            rows.append(future.result())
    cv_metrics = pd.DataFrame(rows)
    selected = select_candidates(cv_metrics)
    args.outdir.mkdir(parents=True, exist_ok=True)
    validation_rows = []
    score_frames = []
    for winner in selected.itertuples(index=False):
        key = (str(winner.bias_configuration), str(winner.cell))
        train_sites, train_profiles, eval_sites, eval_profiles = artifacts[key]
        row, scores = refit_and_score(
            pd.Series(winner._asdict()),
            training_sites=train_sites,
            training_profiles=train_profiles,
            evaluation_sites=eval_sites,
            evaluation_profiles=eval_profiles,
            positions=positions,
            outdir=args.outdir,
            minimum_evaluation_sites=args.minimum_evaluation_sites,
        )
        validation_rows.append(row)
        if scores is not None:
            score_frames.append(scores)
    validation = pd.DataFrame(validation_rows)
    scores = pd.concat(score_frames, ignore_index=True) if score_frames else pd.DataFrame()
    cv_path = args.outdir / "unlabeled_likelihood_cv_metrics.tsv.gz"
    selected_path = args.outdir / "unlabeled_likelihood_selected.tsv"
    validation_path = args.outdir / "unlabeled_likelihood_validation.tsv"
    scores_path = args.outdir / "unlabeled_likelihood_site_scores.tsv.gz"
    cv_metrics.to_csv(cv_path, sep="\t", index=False)
    selected.to_csv(selected_path, sep="\t", index=False)
    validation.to_csv(validation_path, sep="\t", index=False)
    scores.to_csv(scores_path, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-unlabeled-count-model-selection-v1",
        "locked_test_labels_read": False,
        "selection_labels_used": False,
        "validation_labels_opened_only_after_selection": True,
        "motif_or_accessibility_features_used": False,
        "study": str(args.study),
        "study_sha256": file_sha256(args.study),
        "artifacts": input_rows,
        "candidates": [asdict(candidate) for candidate in candidates],
        "tune_chromosomes": tune_chromosomes,
        "minimum_fit_sites": args.minimum_fit_sites,
        "minimum_tune_sites": args.minimum_tune_sites,
        "minimum_evaluation_sites": args.minimum_evaluation_sites,
        "workers": args.workers,
        "outputs": {
            path.name: {"path": str(path), "sha256": file_sha256(path)}
            for path in (cv_path, selected_path, validation_path, scores_path)
        },
    }
    (args.outdir / "unlabeled_likelihood_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    columns = [
        "cell",
        "tf",
        "bias_configuration",
        "candidate_id",
        "heldout_gain_per_site_position",
        "auroc",
        "auprc",
    ]
    print(validation[[column for column in columns if column in validation]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
