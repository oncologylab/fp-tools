#!/usr/bin/env python3
"""Fit strand-aware spline/FDA/GP mixtures without occupancy labels.

Training artifacts must contain no ChIP/label columns. Labels are read only
from separate evaluation artifacts on validation chromosomes after fitting.
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

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import binary_metrics, stable_seed  # noqa: E402
from evaluate_functional_template_transfer import selection_score  # noqa: E402
from evaluate_strand_functional_templates import load_artifact  # noqa: E402
from fp_tools.tools.functional_footprints import (  # noqa: E402
    BiasAwareFunctionalMixture,
    ConditionalMultinomialMixture,
    CovariateAnchoredFdaModel,
    CovariateResidualizedFdaModel,
    FdaMixtureModel,
    HybridFdaGpModel,
    normalize_functional_profiles,
    profile_descriptors,
    standardized_functional_separation,
)
from fp_tools.tools.parametric_bias import estimate_nb_dispersion  # noqa: E402


CHANNELS = (
    "combined_residual",
    "shared_strand_residual",
    "antisymmetric_strand_residual",
)


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    family: str
    smoother: str = ""
    background: str = "none"
    window: float = 50.0
    channel: str = "combined_residual"
    training_pool: str = "tf"
    anchor_strength: float = 0.0
    covariate_ridge: float = 0.0


def candidate_grid() -> list[Candidate]:
    candidates = []
    for family in ("count", "conditional"):
        for smoother in ("spline", "gp"):
            for background in ("none", "linear", "gp-long"):
                for window in (30.0, 50.0, 80.0):
                    candidates.append(
                        Candidate(
                            f"{family}_{smoother}.bg_{background}.window_{int(window)}",
                            family,
                            smoother=smoother,
                            background=background,
                            window=window,
                        )
                    )
    for smoother in ("spline", "gp"):
        for background in ("none", "linear", "gp-long"):
            for window in (30.0, 50.0, 80.0):
                candidates.append(
                    Candidate(
                        "conditional-protected_"
                        f"{smoother}.bg_{background}.window_{int(window)}",
                        "conditional",
                        smoother=smoother,
                        background=background,
                        window=window,
                    )
                )
    for family in ("fda", "hybrid"):
        for channel in CHANNELS:
            for training_pool in ("tf", "family"):
                candidates.append(
                    Candidate(
                        f"{family}.{channel}.pool_{training_pool}",
                        family,
                        channel=channel,
                        training_pool=training_pool,
                    )
                )
    for channel in CHANNELS:
        for training_pool in ("tf", "family"):
            for anchor_strength in (0.5, 1.0, 2.0):
                anchor_label = str(anchor_strength).replace(".", "p")
                candidates.append(
                    Candidate(
                        f"anchored-fda.{channel}.pool_{training_pool}.anchor_{anchor_label}",
                        "anchored-fda",
                        channel=channel,
                        training_pool=training_pool,
                        anchor_strength=anchor_strength,
                    )
                )
            for covariate_ridge in (1.0, 10.0, 100.0):
                ridge_label = str(covariate_ridge).replace(".", "p")
                candidates.append(
                    Candidate(
                        f"residualized-fda.{channel}.pool_{training_pool}.ridge_{ridge_label}",
                        "residualized-fda",
                        channel=channel,
                        training_pool=training_pool,
                        covariate_ridge=covariate_ridge,
                    )
                )
    return candidates


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


def validate_unlabeled_training_sites(sites: pd.DataFrame, path: str | Path) -> None:
    forbidden = [column for column in sites if "label" in column.lower() or "chip" in column.lower()]
    if forbidden:
        raise ValueError(
            f"unlabeled training artifact {path} contains forbidden columns: {', '.join(forbidden)}"
        )
    if "chromosome_split" not in sites:
        raise ValueError(f"unlabeled training artifact {path} lacks chromosome_split")
    observed_splits = set(sites["chromosome_split"].dropna().astype(str))
    if observed_splits != {"train"}:
        raise ValueError(
            f"unlabeled training artifact {path} must contain only train chromosomes; "
            f"observed {sorted(observed_splits)}"
        )


def _count_profiles(profiles: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    observed = profiles["plus_observed"] + profiles["minus_observed"]
    expected = profiles["plus_expected"] + profiles["minus_expected"]
    return observed, expected


def _evaluate_candidate(
    candidate: Candidate,
    *,
    bias_configuration: str,
    cell: str,
    tf: str,
    motif_family: str,
    train_sites: pd.DataFrame,
    train_profiles: dict[str, np.ndarray],
    evaluation_sites: pd.DataFrame,
    evaluation_profiles: dict[str, np.ndarray],
    positions: np.ndarray,
    minimum_evaluation_sites: int,
    seed: int,
) -> tuple[dict, pd.DataFrame | None]:
    started = perf_counter()
    tf_train = np.flatnonzero(train_sites["tf"].astype(str).to_numpy() == tf)
    family_train = np.flatnonzero(
        train_sites["motif_family"].astype(str).to_numpy() == motif_family
    )
    validation = np.flatnonzero(
        (evaluation_sites["tf"].astype(str).to_numpy() == tf)
        & (evaluation_sites["chromosome_split"].astype(str).to_numpy() == "validation")
    )
    labels = evaluation_sites.iloc[validation]["chip_label"].to_numpy(dtype=int)
    base = {
        "cell": cell,
        "tf": tf,
        "motif_family": motif_family,
        "bias_configuration": bias_configuration,
        **asdict(candidate),
        "training_labels_used": False,
        "motif_or_accessibility_features_used": candidate.family
        in {"anchored-fda", "residualized-fda"},
        "evaluation_motif_or_accessibility_features_used": candidate.family
        == "residualized-fda",
        "tf_training_sites": int(len(tf_train)),
        "family_training_sites": int(len(family_train)),
        "validation_sites": int(len(validation)),
        "validation_positive_sites": int(np.sum(labels == 1)),
        "validation_negative_sites": int(np.sum(labels == 0)),
    }
    if (
        len(tf_train) < 100
        or min(np.sum(labels == 0), np.sum(labels == 1)) < minimum_evaluation_sites
    ):
        return {**base, "status": "insufficient_sites"}, None
    try:
        if candidate.family in {"count", "conditional"}:
            train_observed, train_expected = _count_profiles(train_profiles)
            evaluation_observed, evaluation_expected = _count_profiles(evaluation_profiles)
            dispersion = estimate_nb_dispersion(
                train_observed[family_train], train_expected[family_train]
            )
            model_class = (
                BiasAwareFunctionalMixture
                if candidate.family == "count"
                else ConditionalMultinomialMixture
            )
            conditional_kwargs = (
                {
                    "profile_constraint": "canonical-protection"
                    if candidate.candidate_id.startswith("conditional-protected_")
                    else "none"
                }
                if candidate.family == "conditional"
                else {}
            )
            family_model = model_class(
                positions,
                smoother=candidate.smoother,
                dispersion=dispersion,
                max_iter=50,
                shrinkage=0.0,
                accessibility_background=candidate.background,
                prior_constraint="none",
                profile_outer_limit=50.0,
                likelihood_limit=candidate.window,
                **conditional_kwargs,
            )
            family_result = family_model.fit(
                train_observed[family_train],
                train_expected[family_train],
            )
            model = model_class(
                positions,
                smoother=candidate.smoother,
                dispersion=dispersion,
                max_iter=60,
                shrinkage=50.0,
                accessibility_background=candidate.background,
                prior_constraint="none",
                profile_outer_limit=50.0,
                likelihood_limit=candidate.window,
                **conditional_kwargs,
            )
            result = model.fit(
                train_observed[tf_train],
                train_expected[tf_train],
                prior_profile=family_result.footprint_profile,
            )
            shape_log_odds, _prior = model.predict_log_odds_components(
                evaluation_observed[validation],
                evaluation_expected[validation],
            )
            probabilities = 1.0 / (1.0 + np.exp(-np.clip(shape_log_odds, -30.0, 30.0)))
            profile = result.footprint_profile
            standard_error = result.standard_error
            validation_shape = normalize_functional_profiles(
                evaluation_profiles["combined_residual"][validation],
                positions,
            )
            converged = bool(result.converged)
            iterations = int(result.iterations)
        else:
            indexes = tf_train if candidate.training_pool == "tf" else family_train
            train_shape = train_profiles[candidate.channel][indexes]
            validation_shape = evaluation_profiles[candidate.channel][validation]
            weights = np.sqrt(
                np.maximum(
                    (
                        train_profiles["plus_observed"][indexes]
                        + train_profiles["minus_observed"][indexes]
                    ).sum(axis=1),
                    1.0,
                )
            )
            model_seed = stable_seed(
                bias_configuration,
                cell,
                tf,
                candidate.candidate_id,
                seed=seed,
            )
            if candidate.family == "anchored-fda":
                coverage = (
                    train_profiles["plus_observed"][indexes]
                    + train_profiles["minus_observed"][indexes]
                ).sum(axis=1)
                model = CovariateAnchoredFdaModel(
                    max_components=20,
                    anchor_strength=candidate.anchor_strength,
                    seed=model_seed,
                ).fit(
                    train_shape,
                    motif_score=train_sites.iloc[indexes]["motif_score"].to_numpy(
                        dtype=float
                    ),
                    accessibility=coverage,
                    positions=positions,
                    sample_weight=weights,
                )
                shape_log_odds, _anchor = model.predict_log_odds_components(
                    validation_shape
                )
                probabilities = 1.0 / (
                    1.0 + np.exp(-np.clip(shape_log_odds, -40.0, 40.0))
                )
                profile = model.profile_difference()
                converged = bool(model.converged_)
                iterations = int(model.iterations_)
            elif candidate.family == "residualized-fda":
                coverage = (
                    train_profiles["plus_observed"][indexes]
                    + train_profiles["minus_observed"][indexes]
                ).sum(axis=1)
                evaluation_coverage = (
                    evaluation_profiles["plus_observed"][validation]
                    + evaluation_profiles["minus_observed"][validation]
                ).sum(axis=1)
                model = CovariateResidualizedFdaModel(
                    max_components=20,
                    covariate_ridge=candidate.covariate_ridge,
                    seed=model_seed,
                ).fit(
                    train_shape,
                    motif_score=train_sites.iloc[indexes]["motif_score"].to_numpy(
                        dtype=float
                    ),
                    accessibility=coverage,
                    positions=positions,
                    sample_weight=weights,
                )
                probabilities = model.predict_proba(
                    validation_shape,
                    motif_score=evaluation_sites.iloc[validation]["motif_score"].to_numpy(
                        dtype=float
                    ),
                    accessibility=evaluation_coverage,
                )
                profile = model.profile_difference()
                converged = True
                iterations = int(model.mixture.n_iter_) if model.mixture is not None else 0
            elif candidate.family == "fda":
                model = FdaMixtureModel(max_components=20, seed=model_seed).fit(
                    train_shape,
                    positions=positions,
                    sample_weight=weights,
                )
                probabilities = model.predict_proba(validation_shape)
                components = model.component_profiles()
                assert model.binding_component_ is not None
                profile = components[model.binding_component_] - components[1 - model.binding_component_]
                converged = bool(model.mixture.converged_) if model.mixture is not None else False
                iterations = int(model.mixture.n_iter_) if model.mixture is not None else 0
            else:
                model = HybridFdaGpModel(
                    positions,
                    max_components=20,
                    seed=model_seed,
                ).fit(train_shape, sample_weight=weights)
                probabilities = model.predict_proba(validation_shape)
                assert model.bound_mean_ is not None and model.unbound_mean_ is not None
                profile = model.bound_mean_ - model.unbound_mean_
                converged = bool(model.fda.mixture.converged_) if model.fda.mixture is not None else False
                iterations = int(model.fda.mixture.n_iter_) if model.fda.mixture is not None else 0
            standard_error = np.full_like(profile, np.nan)
        metrics = binary_metrics(labels, probabilities)
        descriptors = asdict(profile_descriptors(profile, positions))
        normalized_validation = normalize_functional_profiles(validation_shape, positions)
        positive_mean = np.mean(normalized_validation[labels == 1], axis=0)
        negative_mean = np.mean(normalized_validation[labels == 0], axis=0)
        row = {
            **base,
            "status": "ok",
            "fit_seconds": perf_counter() - started,
            "converged": converged,
            "iterations": iterations,
            "selection_score": selection_score(metrics),
            "functional_separation": standardized_functional_separation(
                normalized_validation,
                labels,
                positions,
            ),
            **descriptors,
            **metrics,
        }
        curve = pd.DataFrame(
            {
                "cell": cell,
                "tf": tf,
                "motif_family": motif_family,
                "bias_configuration": bias_configuration,
                "candidate_id": candidate.candidate_id,
                "position": positions.astype(int),
                "footprint_profile": profile,
                "standard_error": standard_error,
                "positive_mean": positive_mean,
                "negative_mean": negative_mean,
                "positive_minus_negative": positive_mean - negative_mean,
            }
        )
        return row, curve
    except Exception as error:
        return {
            **base,
            "status": "error",
            "error": f"{type(error).__name__}: {error}",
            "fit_seconds": perf_counter() - started,
        }, None


def select_winners(metrics: pd.DataFrame) -> pd.DataFrame:
    passing = metrics[(metrics["status"] == "ok") & metrics["converged"].astype(bool)].copy()
    keys = ["cell", "tf", "bias_configuration"]
    if passing.empty:
        return passing
    return (
        passing.sort_values(
            keys + ["selection_score", "auprc", "auroc", "candidate_id"],
            ascending=[True, True, True, False, False, False, True],
            kind="mergesort",
        )
        .groupby(keys, as_index=False, sort=True)
        .head(1)
        .reset_index(drop=True)
    )


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
    parser.add_argument("--minimum-evaluation-sites", type=int, default=100)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument(
        "--candidate-prefix",
        action="append",
        help="Evaluate only candidate IDs beginning with this prefix; repeat as needed",
    )
    args = parser.parse_args(argv)
    study = json.loads(args.study.read_text(encoding="utf-8"))
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[tasks["split"] == "development"]
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
                    "model": model,
                    "cell": cell,
                    "purpose": purpose,
                    "path": str(path),
                    "sha256": file_sha256(path),
                    "sites": int(document["sites_total"]),
                }
            )
    candidates = candidate_grid()
    if args.candidate_prefix is not None:
        candidates = [
            candidate
            for candidate in candidates
            if any(
                candidate.candidate_id.startswith(prefix)
                for prefix in args.candidate_prefix
            )
        ]
        if not candidates:
            raise SystemExit("no candidates match --candidate-prefix")
    positions = np.arange(-int(study["profile_flank_bp"]), int(study["profile_flank_bp"]) + 1, dtype=float)
    futures = {}
    rows = []
    curves = []
    with threadpool_limits(limits=1), ThreadPoolExecutor(max_workers=args.workers) as executor:
        for (bias_configuration, cell), artifact in artifacts.items():
            train_sites, train_profiles, eval_sites, eval_profiles = artifact
            for task in tasks[tasks["cell"].astype(str) == cell].itertuples(index=False):
                for candidate in candidates:
                    future = executor.submit(
                        _evaluate_candidate,
                        candidate,
                        bias_configuration=bias_configuration,
                        cell=cell,
                        tf=str(task.tf),
                        motif_family=str(task.motif_family),
                        train_sites=train_sites,
                        train_profiles=train_profiles,
                        evaluation_sites=eval_sites,
                        evaluation_profiles=eval_profiles,
                        positions=positions,
                        minimum_evaluation_sites=args.minimum_evaluation_sites,
                        seed=args.seed,
                    )
                    futures[future] = (bias_configuration, cell, str(task.tf), candidate.candidate_id)
        for future in as_completed(futures):
            row, curve = future.result()
            rows.append(row)
            if curve is not None:
                curves.append(curve)
    metrics = pd.DataFrame(rows)
    profiles = pd.concat(curves, ignore_index=True) if curves else pd.DataFrame()
    winners = select_winners(metrics)
    args.outdir.mkdir(parents=True, exist_ok=True)
    metrics_path = args.outdir / "strand_label_free_metrics.tsv.gz"
    profiles_path = args.outdir / "strand_label_free_profiles.tsv.gz"
    winners_path = args.outdir / "strand_label_free_winners.tsv"
    metrics.to_csv(metrics_path, sep="\t", index=False)
    profiles.to_csv(profiles_path, sep="\t", index=False)
    winners.to_csv(winners_path, sep="\t", index=False)
    manifest = {
        "schema": "fp-tools-strand-label-free-evaluation-v1",
        "locked_test_labels_read": False,
        "training_labels_used": False,
        "motif_or_accessibility_features_used": any(
            candidate.family in {"anchored-fda", "residualized-fda"}
            for candidate in candidates
        ),
        "evaluation_motif_or_accessibility_features_used": any(
            candidate.family == "residualized-fda" for candidate in candidates
        ),
        "study": str(args.study),
        "study_sha256": file_sha256(args.study),
        "artifacts": input_rows,
        "candidate_count": len(candidates),
        "candidates": [asdict(candidate) for candidate in candidates],
        "minimum_evaluation_sites": args.minimum_evaluation_sites,
        "workers": args.workers,
        "seed": args.seed,
        "metrics_rows": int(len(metrics)),
        "winner_rows": int(len(winners)),
        "outputs": {
            "metrics": {"path": str(metrics_path), "sha256": file_sha256(metrics_path)},
            "profiles": {"path": str(profiles_path), "sha256": file_sha256(profiles_path)},
            "winners": {"path": str(winners_path), "sha256": file_sha256(winners_path)},
        },
    }
    (args.outdir / "strand_label_free_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    columns = ["cell", "tf", "bias_configuration", "candidate_id", "auroc", "auprc"]
    print(winners[columns].to_string(index=False) if len(winners) else "no eligible winners")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
