#!/usr/bin/env python3
"""Search TF-specific functional models without opening locked holdouts.

This script reuses the checksummed motif-profile cache from a completed
``evaluate_functional_footprints.py`` development run. It evaluates broad
site-accessibility backgrounds, monotone binding priors, spline penalties,
Matérn length scales, shrinkage levels, and FDA dimensions for every eligible
development TF. ChIP labels are used only to rank candidates on validation
chromosomes; all mixture fits use unlabeled train-chromosome motif sites.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter
from typing import Any, Sequence

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import (  # noqa: E402
    _evaluation_profiles,
    binary_metrics,
    chromosome_split,
    stable_seed,
    validate_sites,
)
from fp_tools.tools.functional_footprints import (  # noqa: E402
    BiasAwareFunctionalMixture,
    FdaMixtureModel,
    HybridFdaGpModel,
    deviance_profiles,
    profile_descriptors,
)
from fp_tools.tools.parametric_bias import estimate_nb_dispersion  # noqa: E402


MODEL_FAMILIES = ("spline", "gp", "fda", "hybrid")


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class FunctionalCandidate:
    candidate_id: str
    family: str
    background: str = "none"
    prior_constraint: str = "none"
    shrinkage: float = 50.0
    spline_penalty: float = 10.0
    long_length_scale: float = 50.0
    short_length_scale: float = 10.0
    gp_ridge: float = 1.0
    variance_threshold: float = 0.95
    max_components: int = 20


def _number(value: float) -> str:
    return f"{float(value):g}".replace(".", "p")


def _count_candidate(
    family: str,
    *,
    background: str = "none",
    prior_constraint: str = "none",
    shrinkage: float = 50.0,
    spline_penalty: float = 10.0,
    long_length_scale: float = 50.0,
    short_length_scale: float = 10.0,
    gp_ridge: float = 1.0,
) -> FunctionalCandidate:
    prior = {"none": "free", "motif": "motif", "motif-accessibility": "motif_access"}[
        prior_constraint
    ]
    if family == "spline":
        identity = (
            f"spline.bg_{background}.prior_{prior}."
            f"pen_{_number(spline_penalty)}.shrink_{_number(shrinkage)}"
        )
    else:
        identity = (
            f"gp.bg_{background}.prior_{prior}."
            f"long_{_number(long_length_scale)}.short_{_number(short_length_scale)}."
            f"ridge_{_number(gp_ridge)}.shrink_{_number(shrinkage)}"
        )
    return FunctionalCandidate(
        identity,
        family,
        background=background,
        prior_constraint=prior_constraint,
        shrinkage=shrinkage,
        spline_penalty=spline_penalty,
        long_length_scale=long_length_scale,
        short_length_scale=short_length_scale,
        gp_ridge=gp_ridge,
    )


def _fda_candidate(
    family: str,
    variance_threshold: float,
    max_components: int,
) -> FunctionalCandidate:
    identity = (
        f"{family}.variance_{_number(variance_threshold)}.components_{max_components}"
    )
    return FunctionalCandidate(
        identity,
        family,
        variance_threshold=variance_threshold,
        max_components=max_components,
    )


def candidate_grid(
    grid: str,
    families: Sequence[str] = MODEL_FAMILIES,
) -> list[FunctionalCandidate]:
    """Return a deterministic, de-duplicated staged hyperparameter grid."""

    if grid not in {"compact", "full"}:
        raise ValueError("grid must be compact or full")
    unknown = set(families).difference(MODEL_FAMILIES)
    if unknown:
        raise ValueError("unknown model families: " + ", ".join(sorted(unknown)))
    candidates: dict[str, FunctionalCandidate] = {}

    for family in set(families).intersection({"spline", "gp"}):
        for background in ("none", "linear", "quadratic", "gp-long"):
            for constraint in ("none", "motif-accessibility"):
                candidate = _count_candidate(
                    family,
                    background=background,
                    prior_constraint=constraint,
                )
                candidates[candidate.candidate_id] = candidate

    if grid == "full":
        if "spline" in families:
            for penalty in (1.0, 3.0, 10.0, 30.0, 100.0):
                candidate = _count_candidate(
                    "spline",
                    prior_constraint="motif-accessibility",
                    spline_penalty=penalty,
                )
                candidates[candidate.candidate_id] = candidate
            for shrinkage in (0.0, 10.0, 50.0, 200.0):
                candidate = _count_candidate(
                    "spline",
                    prior_constraint="motif-accessibility",
                    shrinkage=shrinkage,
                )
                candidates[candidate.candidate_id] = candidate
            candidate = _count_candidate("spline", prior_constraint="motif")
            candidates[candidate.candidate_id] = candidate
        if "gp" in families:
            for long_scale in (30.0, 50.0, 80.0):
                for short_scale in (3.0, 6.0, 10.0, 15.0):
                    candidate = _count_candidate(
                        "gp",
                        prior_constraint="motif-accessibility",
                        long_length_scale=long_scale,
                        short_length_scale=short_scale,
                    )
                    candidates[candidate.candidate_id] = candidate
            for ridge in (0.3, 1.0, 3.0):
                candidate = _count_candidate(
                    "gp",
                    prior_constraint="motif-accessibility",
                    gp_ridge=ridge,
                )
                candidates[candidate.candidate_id] = candidate
            for shrinkage in (0.0, 10.0, 50.0, 200.0):
                candidate = _count_candidate(
                    "gp",
                    prior_constraint="motif-accessibility",
                    shrinkage=shrinkage,
                )
                candidates[candidate.candidate_id] = candidate
            candidate = _count_candidate("gp", prior_constraint="motif")
            candidates[candidate.candidate_id] = candidate

    fda_thresholds = (0.95,) if grid == "compact" else (0.80, 0.90, 0.95, 0.99)
    component_limits = (20,) if grid == "compact" else (5, 10, 20)
    for family in set(families).intersection({"fda", "hybrid"}):
        for threshold in fda_thresholds:
            for components in component_limits:
                candidate = _fda_candidate(family, threshold, components)
                candidates[candidate.candidate_id] = candidate
    return [candidates[key] for key in sorted(candidates)]


def _count_model(
    candidate: FunctionalCandidate,
    positions: np.ndarray,
    dispersion: float,
    *,
    shrinkage: float | None = None,
) -> BiasAwareFunctionalMixture:
    return BiasAwareFunctionalMixture(
        positions,
        smoother=candidate.family,
        dispersion=dispersion,
        max_iter=75,
        tolerance=1e-5,
        shrinkage=candidate.shrinkage if shrinkage is None else shrinkage,
        long_length_scale=candidate.long_length_scale,
        short_length_scale=candidate.short_length_scale,
        spline_penalty=candidate.spline_penalty,
        gp_ridge=candidate.gp_ridge,
        accessibility_background=candidate.background,
        prior_constraint=candidate.prior_constraint,
    )


def _selection_score(metrics: dict[str, float | int]) -> float:
    prevalence = float(metrics["prevalence"])
    adjusted_auprc = (float(metrics["auprc"]) - prevalence) / max(
        1.0 - prevalence,
        1e-8,
    )
    return float(metrics["auroc"]) + adjusted_auprc


def evaluate_candidate(
    candidate: FunctionalCandidate,
    *,
    cell: str,
    correction: str,
    tasks: pd.DataFrame,
    unlabeled_sites: pd.DataFrame,
    train_observed: np.ndarray,
    train_expected: np.ndarray,
    development_sites: pd.DataFrame,
    development_observed: np.ndarray,
    development_expected: np.ndarray,
    positions: np.ndarray,
    dispersion: float,
    minimum_evaluation_sites: int,
    seed: int,
) -> tuple[list[dict[str, Any]], list[pd.DataFrame]]:
    """Fit one configuration across every eligible TF in one cell/correction."""

    start_time = perf_counter()
    rows: list[dict[str, Any]] = []
    curves: list[pd.DataFrame] = []
    family_priors: dict[str, np.ndarray] = {}
    global_prior: np.ndarray | None = None
    if candidate.family in {"spline", "gp"} and candidate.shrinkage > 0:
        limit = min(len(unlabeled_sites), 400)
        indexes = np.random.default_rng(
            stable_seed(cell, correction, candidate.candidate_id, "global", seed=seed)
        ).choice(len(unlabeled_sites), size=limit, replace=False)
        global_model = _count_model(candidate, positions, dispersion, shrinkage=0.0)
        global_result = global_model.fit(
            train_observed[indexes],
            train_expected[indexes],
            motif_score=unlabeled_sites.iloc[indexes]["motif_score"].to_numpy(dtype=float),
            accessibility=train_observed[indexes].sum(axis=1),
        )
        global_prior = global_result.footprint_profile
        for family, family_sites in unlabeled_sites.groupby("motif_family", sort=True):
            family_indexes = family_sites.index.to_numpy(dtype=int)
            family_model = _count_model(candidate, positions, dispersion)
            family_result = family_model.fit(
                train_observed[family_indexes],
                train_expected[family_indexes],
                motif_score=unlabeled_sites.iloc[family_indexes]["motif_score"].to_numpy(
                    dtype=float
                ),
                accessibility=train_observed[family_indexes].sum(axis=1),
                prior_profile=global_prior,
            )
            family_priors[str(family)] = family_result.footprint_profile

    residual = None
    if candidate.family in {"fda", "hybrid"}:
        residual = deviance_profiles(train_observed, train_expected, dispersion)

    for task in tasks[tasks["cell"] == cell].itertuples(index=False):
        tf = str(task.tf)
        family = str(task.motif_family)
        train_indexes = np.flatnonzero(unlabeled_sites["tf"].astype(str).to_numpy() == tf)
        validation_indexes = np.flatnonzero(
            (development_sites["tf"].astype(str).to_numpy() == tf)
            & (development_sites["chromosome_split"].astype(str).to_numpy() == "validation")
        )
        if len(train_indexes) < 100 or len(validation_indexes) < minimum_evaluation_sites:
            continue
        fit_start = perf_counter()
        if candidate.family in {"spline", "gp"}:
            model = _count_model(candidate, positions, dispersion)
            result = model.fit(
                train_observed[train_indexes],
                train_expected[train_indexes],
                motif_score=unlabeled_sites.iloc[train_indexes]["motif_score"].to_numpy(
                    dtype=float
                ),
                accessibility=train_observed[train_indexes].sum(axis=1),
                prior_profile=family_priors.get(family, global_prior),
            )
            likelihood_log_odds, prior_log_odds = model.predict_log_odds_components(
                development_observed[validation_indexes],
                development_expected[validation_indexes],
                motif_score=development_sites.iloc[validation_indexes]["motif_score"].to_numpy(
                    dtype=float
                ),
                accessibility=development_observed[validation_indexes].sum(axis=1),
            )
            shape_probabilities = 1.0 / (
                1.0 + np.exp(-np.clip(likelihood_log_odds, -40.0, 40.0))
            )
            prior_probabilities = 1.0 / (
                1.0 + np.exp(-np.clip(prior_log_odds, -40.0, 40.0))
            )
            probabilities = 1.0 / (
                1.0
                + np.exp(-np.clip(likelihood_log_odds + prior_log_odds, -40.0, 40.0))
            )
            profile = result.footprint_profile
            standard_error = result.standard_error
            converged = result.converged
            iterations = result.iterations
            prior_coefficients = result.prior_coefficients
        elif candidate.family == "fda":
            assert residual is not None
            model = FdaMixtureModel(
                variance_threshold=candidate.variance_threshold,
                max_components=candidate.max_components,
                seed=stable_seed(cell, correction, tf, candidate.candidate_id, seed=seed),
            ).fit(
                residual[train_indexes],
                positions=positions,
                sample_weight=np.sqrt(
                    np.maximum(train_observed[train_indexes].sum(axis=1), 1.0)
                ),
            )
            validation_residual = deviance_profiles(
                development_observed[validation_indexes],
                development_expected[validation_indexes],
                dispersion,
            )
            probabilities = model.predict_proba(validation_residual)
            component_profiles = model.component_profiles()
            assert model.binding_component_ is not None
            profile = component_profiles[model.binding_component_] - component_profiles[
                1 - model.binding_component_
            ]
            standard_error = np.full_like(profile, np.nan)
            converged = bool(model.mixture.converged_) if model.mixture is not None else False
            iterations = int(model.mixture.n_iter_) if model.mixture is not None else 0
            prior_coefficients = np.full(3, np.nan)
            shape_probabilities = probabilities
            prior_probabilities = np.full_like(probabilities, 0.5)
        else:
            assert residual is not None
            model = HybridFdaGpModel(
                positions,
                variance_threshold=candidate.variance_threshold,
                max_components=candidate.max_components,
                seed=stable_seed(cell, correction, tf, candidate.candidate_id, seed=seed),
            ).fit(
                residual[train_indexes],
                sample_weight=np.sqrt(
                    np.maximum(train_observed[train_indexes].sum(axis=1), 1.0)
                ),
            )
            validation_residual = deviance_profiles(
                development_observed[validation_indexes],
                development_expected[validation_indexes],
                dispersion,
            )
            probabilities = model.predict_proba(validation_residual)
            assert model.bound_mean_ is not None and model.unbound_mean_ is not None
            profile = model.bound_mean_ - model.unbound_mean_
            standard_error = np.full_like(profile, np.nan)
            converged = bool(
                model.fda.mixture.converged_ if model.fda.mixture is not None else False
            )
            iterations = int(
                model.fda.mixture.n_iter_ if model.fda.mixture is not None else 0
            )
            prior_coefficients = np.full(3, np.nan)
            shape_probabilities = probabilities
            prior_probabilities = np.full_like(probabilities, 0.5)

        labels = development_sites.iloc[validation_indexes]["chip_label"].to_numpy(dtype=int)
        metrics = binary_metrics(labels, probabilities)
        shape_metrics = binary_metrics(labels, shape_probabilities)
        prior_metrics = binary_metrics(labels, prior_probabilities)
        descriptors = asdict(profile_descriptors(profile, positions))
        rows.append(
            {
                "cell": cell,
                "tf": tf,
                "motif_family": family,
                "correction": correction,
                **asdict(candidate),
                "training_labels_used": False,
                "training_sites": int(len(train_indexes)),
                "fit_seconds": perf_counter() - fit_start,
                "converged": converged,
                "iterations": iterations,
                "prior_intercept": float(prior_coefficients[0]),
                "prior_motif": float(prior_coefficients[1]),
                "prior_accessibility": float(prior_coefficients[2]),
                "dispersion": dispersion,
                **descriptors,
                **metrics,
                "shape_auroc": float(shape_metrics["auroc"]),
                "shape_auprc": float(shape_metrics["auprc"]),
                "shape_brier": float(shape_metrics["brier"]),
                "prior_auroc": float(prior_metrics["auroc"]),
                "prior_auprc": float(prior_metrics["auprc"]),
                "prior_brier": float(prior_metrics["brier"]),
                "profile_incremental_auprc": float(metrics["auprc"])
                - float(prior_metrics["auprc"]),
                "selection_score": _selection_score(metrics),
                "status": "ok",
            }
        )
        curves.append(
            pd.DataFrame(
                {
                    "cell": cell,
                    "tf": tf,
                    "motif_family": family,
                    "correction": correction,
                    "candidate_id": candidate.candidate_id,
                    "position": positions.astype(int),
                    "footprint_profile": profile,
                    "standard_error": standard_error,
                }
            )
        )
    elapsed = perf_counter() - start_time
    for row in rows:
        row["candidate_seconds_all_tfs"] = elapsed
    return rows, curves


def select_winners(metrics: pd.DataFrame) -> pd.DataFrame:
    passing = metrics[(metrics["status"] == "ok") & metrics["converged"].astype(bool)].copy()
    if passing.empty:
        return passing
    return (
        passing.sort_values(
            ["cell", "tf", "selection_score", "auprc", "auroc", "candidate_id"],
            ascending=[True, True, False, False, False, True],
            kind="mergesort",
        )
        .groupby(["cell", "tf"], as_index=False)
        .first()
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-run", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--corrections", nargs="+")
    parser.add_argument("--families", nargs="+", choices=MODEL_FAMILIES, default=list(MODEL_FAMILIES))
    parser.add_argument("--grid", choices=("compact", "full"), default="compact")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--minimum-evaluation-sites", type=int)
    parser.add_argument("--seed", type=int)
    args = parser.parse_args(argv)
    if args.workers < 1:
        raise SystemExit("--workers must be positive")

    base_manifest_path = args.base_run / "functional_benchmark_manifest.json"
    base_manifest = json.loads(base_manifest_path.read_text(encoding="utf-8"))
    if base_manifest.get("test_unlocked"):
        raise SystemExit("hyperparameter search refuses a base run that opened test labels")
    study_path = Path(base_manifest["study"])
    development_path = Path(base_manifest["development_sites"])
    tracks_path = Path(base_manifest["tracks"])
    genome = Path(base_manifest["genome"]) if base_manifest.get("genome") else None
    study = json.loads(study_path.read_text(encoding="utf-8"))
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[tasks["split"] == "development"].copy()
    development = validate_sites(pd.read_csv(development_path, sep="\t"), development_path)
    development["chromosome_split"] = development["TFBS_chr"].map(
        lambda value: chromosome_split(str(value), study)
    )
    tracks = pd.read_csv(tracks_path, sep="\t")
    corrections = tuple(args.corrections or base_manifest["corrections"])
    unknown_corrections = set(corrections).difference(tracks["model"].astype(str))
    if unknown_corrections:
        raise SystemExit("unknown corrections: " + ", ".join(sorted(unknown_corrections)))
    minimum_sites = int(
        args.minimum_evaluation_sites or base_manifest["minimum_evaluation_sites"]
    )
    seed = int(args.seed if args.seed is not None else base_manifest["seed"])
    flank = int(base_manifest["flank"])
    positions = np.arange(-flank, flank + 1, dtype=float)
    candidates = candidate_grid(args.grid, args.families)
    args.outdir.mkdir(parents=True, exist_ok=True)

    futures = {}
    rows: list[dict[str, Any]] = []
    curve_frames: list[pd.DataFrame] = []
    with threadpool_limits(limits=1), ThreadPoolExecutor(max_workers=args.workers) as executor:
        for cell in sorted(tasks["cell"].unique()):
            unlabeled_path = args.base_run / f"{cell}.unlabeled_training_sites.tsv.gz"
            unlabeled = pd.read_csv(unlabeled_path, sep="\t")
            if any("chip" in column.lower() or "label" in column.lower() for column in unlabeled):
                raise SystemExit(f"unlabeled training table contains a label column: {unlabeled_path}")
            cell_development = development[development["cell"] == cell].reset_index(drop=True)
            for correction in corrections:
                train_observed, train_expected = _evaluation_profiles(
                    unlabeled,
                    tracks,
                    cell,
                    correction,
                    args.base_run / "profile_cache",
                    "unlabeled",
                    flank,
                    genome,
                )
                development_observed, development_expected = _evaluation_profiles(
                    cell_development,
                    tracks,
                    cell,
                    correction,
                    args.base_run / "profile_cache",
                    "development",
                    flank,
                    genome,
                )
                local_development = cell_development.copy()
                local_development["accessibility"] = development_observed.sum(axis=1)
                dispersion = estimate_nb_dispersion(train_observed, train_expected)
                for candidate in candidates:
                    future = executor.submit(
                        evaluate_candidate,
                        candidate,
                        cell=cell,
                        correction=correction,
                        tasks=tasks,
                        unlabeled_sites=unlabeled,
                        train_observed=train_observed,
                        train_expected=train_expected,
                        development_sites=local_development,
                        development_observed=development_observed,
                        development_expected=development_expected,
                        positions=positions,
                        dispersion=dispersion,
                        minimum_evaluation_sites=minimum_sites,
                        seed=seed,
                    )
                    futures[future] = (cell, correction, candidate)
        for future in as_completed(futures):
            cell, correction, candidate = futures[future]
            try:
                candidate_rows, candidate_curves = future.result()
                rows.extend(candidate_rows)
                curve_frames.extend(candidate_curves)
            except Exception as error:  # retain the rest of a large factorial run
                rows.append(
                    {
                        "cell": cell,
                        "tf": "*",
                        "motif_family": "*",
                        "correction": correction,
                        **asdict(candidate),
                        "status": "error",
                        "error": f"{type(error).__name__}: {error}",
                    }
                )

    metrics = pd.DataFrame(rows)
    if len(metrics):
        metrics = metrics.sort_values(
            ["cell", "tf", "correction", "candidate_id"],
            kind="mergesort",
        )
    curves = pd.concat(curve_frames, ignore_index=True) if curve_frames else pd.DataFrame()
    winners = select_winners(metrics)
    metrics_path = args.outdir / "functional_hyperparameter_metrics.tsv.gz"
    curves_path = args.outdir / "functional_hyperparameter_profiles.tsv.gz"
    winners_path = args.outdir / "functional_hyperparameter_winners.tsv"
    metrics.to_csv(metrics_path, sep="\t", index=False)
    curves.to_csv(curves_path, sep="\t", index=False)
    winners.to_csv(winners_path, sep="\t", index=False)
    document = {
        "schema": "fp-tools-functional-hyperparameter-search-v1",
        "base_run": str(args.base_run),
        "base_manifest": str(base_manifest_path),
        "base_manifest_sha256": file_sha256(base_manifest_path),
        "locked_test_labels_read": False,
        "selection_split": "validation",
        "corrections": list(corrections),
        "families": list(args.families),
        "grid": args.grid,
        "candidate_count": len(candidates),
        "candidates": [asdict(candidate) for candidate in candidates],
        "workers": args.workers,
        "minimum_evaluation_sites": minimum_sites,
        "seed": seed,
        "outputs": {
            "metrics": {"path": str(metrics_path), "sha256": file_sha256(metrics_path)},
            "profiles": {"path": str(curves_path), "sha256": file_sha256(curves_path)},
            "winners": {"path": str(winners_path), "sha256": file_sha256(winners_path)},
        },
    }
    (args.outdir / "functional_hyperparameter_manifest.json").write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(winners[["cell", "tf", "correction", "candidate_id", "auroc", "auprc"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
