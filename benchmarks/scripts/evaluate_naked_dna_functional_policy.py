#!/usr/bin/env python3
"""Apply the frozen functional-detector policy to naked-DNA controls.

Models are refitted only from the frozen, label-free K562/HepG2 training
artifacts.  Naked-DNA profiles never influence model choice or fitting.  Every
promoted family is compared with the exact DWM reference recorded at policy
freeze, on identical motif sites and with a fixed posterior threshold.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any, Sequence

import numpy as np
import pandas as pd
from scipy.special import expit

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from build_combined_functional_profiles import SCHEMA as COMBINED_SCHEMA  # noqa: E402
from build_strand_functional_profiles import site_hashes  # noqa: E402
from evaluate_functional_footprints import stable_seed  # noqa: E402
from evaluate_strand_label_free_models import (  # noqa: E402
    Candidate as StrandCandidate,
    candidate_grid as strand_candidate_grid,
)
from search_functional_model_grid import (  # noqa: E402
    FunctionalCandidate,
    _count_model,
    candidate_grid as dwm_candidate_grid,
)
from select_count_models_by_unlabeled_likelihood import make_model  # noqa: E402
from fp_tools.tools.functional_footprints import (  # noqa: E402
    BiasAwareFunctionalMixture,
    CovariateAnchoredFdaModel,
    CovariateResidualizedFdaModel,
    FdaMixtureModel,
    HybridFdaGpModel,
    deviance_profiles,
    normalize_functional_profiles,
    orient_profiles,
    profile_descriptors,
)
from fp_tools.tools.parametric_bias import estimate_nb_dispersion  # noqa: E402


STRAND_SCHEMA = "fp-tools-strand-functional-profiles-v1"
PROMOTION_METHODS = {
    "frozen_policy_candidate": "frozen_policy",
    "frozen_dwm_reference": "DWM_reference",
}
STRAND_ARRAYS = (
    "plus_observed",
    "minus_observed",
    "plus_expected",
    "minus_expected",
    "combined_residual",
    "shared_strand_residual",
    "antisymmetric_strand_residual",
)


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_training_artifact(value: str) -> tuple[str, str, Path]:
    fields = value.split(",", 2)
    if len(fields) != 3 or not all(fields):
        raise argparse.ArgumentTypeError("training artifact must use MODEL,CELL,JSON")
    return fields[0], fields[1], Path(fields[2])


def parse_naked_artifact(value: str) -> tuple[str, str, str, Path]:
    fields = value.split(",", 3)
    if len(fields) != 4 or not all(fields):
        raise argparse.ArgumentTypeError(
            "naked strand artifact must use MODEL,CELL,REPLICATE,JSON"
        )
    return fields[0], fields[1], fields[2], Path(fields[3])


def parse_dwm_artifact(value: str) -> tuple[str, str, Path]:
    fields = value.split(",", 2)
    if len(fields) != 3 or not all(fields):
        raise argparse.ArgumentTypeError(
            "naked DWM artifact must use CELL,REPLICATE,JSON"
        )
    return fields[0], fields[1], Path(fields[2])


def parse_dwm_training_artifact(value: str) -> tuple[str, Path]:
    fields = value.split(",", 1)
    if len(fields) != 2 or not all(fields):
        raise argparse.ArgumentTypeError(
            "DWM training artifact must use CELL,JSON"
        )
    return fields[0], Path(fields[1])


def _validate_label_free_sites(sites: pd.DataFrame, path: Path) -> None:
    forbidden = [
        column
        for column in sites.columns
        if "label" in column.lower() or "chip" in column.lower()
    ]
    if forbidden:
        raise ValueError(f"{path} contains forbidden columns: {', '.join(forbidden)}")


def _load_profile_artifact(
    path: Path,
    *,
    expected_schema: str,
    expected_cell: str,
    arrays_required: Sequence[str],
) -> tuple[pd.DataFrame, dict[str, np.ndarray], np.ndarray, np.ndarray, dict]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != expected_schema:
        raise ValueError(f"unsupported artifact schema: {path}")
    if document.get("metadata", {}).get("labels_used") is not False:
        raise ValueError(f"artifact does not certify label-free construction: {path}")
    profiles_path = Path(document["profiles_npz"])
    sites_path = Path(document["sites"])
    if file_sha256(profiles_path) != document["profiles_sha256"]:
        raise ValueError(f"profile checksum mismatch: {profiles_path}")
    if file_sha256(sites_path) != document["sites_sha256"]:
        raise ValueError(f"site checksum mismatch: {sites_path}")
    sites = pd.read_csv(sites_path, sep="\t").reset_index(drop=True)
    _validate_label_free_sites(sites, sites_path)
    if "cell" not in sites or set(sites["cell"].astype(str)) != {expected_cell}:
        raise ValueError(f"artifact sites do not exclusively contain {expected_cell}: {path}")
    hashes = site_hashes(sites)
    with np.load(profiles_path, allow_pickle=False) as arrays:
        missing = sorted(set(arrays_required).difference(arrays.files))
        if missing:
            raise ValueError(f"artifact is missing arrays {missing}: {path}")
        if not np.array_equal(np.asarray(arrays["site_hash"], dtype=np.uint64), hashes):
            raise ValueError(f"site order hash mismatch: {path}")
        profiles = {
            name: np.asarray(arrays[name], dtype=np.float64)
            for name in arrays_required
        }
        valid = np.asarray(arrays["valid"], dtype=bool)
    if valid.shape != (len(sites),):
        raise ValueError(f"artifact valid mask has wrong shape: {path}")
    return sites, profiles, valid, hashes, document


def load_strand_artifact(path: Path, cell: str):
    return _load_profile_artifact(
        path,
        expected_schema=STRAND_SCHEMA,
        expected_cell=cell,
        arrays_required=STRAND_ARRAYS,
    )


def load_combined_artifact(path: Path, cell: str):
    return _load_profile_artifact(
        path,
        expected_schema=COMBINED_SCHEMA,
        expected_cell=cell,
        arrays_required=("observed", "expected", "combined_residual"),
    )


def load_dwm_training(base_run: Path, cell: str, flank: int):
    sites_path = base_run / f"{cell}.unlabeled_training_sites.tsv.gz"
    sites = pd.read_csv(sites_path, sep="\t").reset_index(drop=True)
    _validate_label_free_sites(sites, sites_path)
    hashes = site_hashes(sites)
    profiles: dict[str, np.ndarray] = {}
    valid = np.ones(len(sites), dtype=bool)
    inputs = []
    for track in ("raw", "expected"):
        path = base_run / "profile_cache" / f"unlabeled.{cell}.DWM.{track}.flank{flank}.npz"
        with np.load(path, allow_pickle=False) as arrays:
            if not np.array_equal(np.asarray(arrays["site_hash"], dtype=np.uint64), hashes):
                raise ValueError(f"DWM training cache site order mismatch: {path}")
            values = np.asarray(arrays["profiles"], dtype=np.float64)
            profiles["observed" if track == "raw" else "expected"] = orient_profiles(
                values, sites["TFBS_strand"].astype(str)
            )
            valid &= np.asarray(arrays["valid"], dtype=bool)
        inputs.append({"path": str(path), "sha256": file_sha256(path)})
    profiles["combined_residual"] = deviance_profiles(
        profiles["observed"], profiles["expected"], 0.0
    )
    return sites, profiles, valid, hashes, inputs


def load_dwm_training_source(
    cell: str,
    *,
    artifact_paths: dict[str, Path],
    base_run: Path | None,
    flank: int,
):
    """Load an explicit combined artifact, with legacy cache fallback."""

    if cell in artifact_paths:
        path = artifact_paths[cell]
        sites, profiles, valid, hashes, document = load_combined_artifact(
            path, cell
        )
        return sites, profiles, valid, hashes, [
            {
                "path": str(path),
                "sha256": file_sha256(path),
                "profiles_sha256": str(document["profiles_sha256"]),
                "sites_sha256": str(document["sites_sha256"]),
            }
        ]
    if base_run is not None:
        return load_dwm_training(base_run, cell, flank)
    raise ValueError(f"no DWM training source was provided for {cell}")


@dataclass
class FittedDetector:
    model_family: str
    model: Any
    dispersion: float
    channel: str = "combined_residual"
    source: str = "strand"


def _count_arrays(profiles: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    if "observed" in profiles:
        return profiles["observed"], profiles["expected"]
    return (
        profiles["plus_observed"] + profiles["minus_observed"],
        profiles["plus_expected"] + profiles["minus_expected"],
    )


def fit_strand_detector(
    candidate: StrandCandidate,
    sites: pd.DataFrame,
    profiles: dict[str, np.ndarray],
    *,
    tf: str,
    motif_family: str,
    positions: np.ndarray,
    seed: int,
) -> FittedDetector:
    tf_indexes = np.flatnonzero(sites["tf"].astype(str).to_numpy() == tf)
    family_indexes = np.flatnonzero(
        sites["motif_family"].astype(str).to_numpy() == motif_family
    )
    if len(tf_indexes) < 100 or len(family_indexes) < 100:
        raise ValueError(f"insufficient strand training sites for {tf}/{motif_family}")
    if candidate.family == "count":
        observed, expected = _count_arrays(profiles)
        dispersion = estimate_nb_dispersion(
            observed[family_indexes], expected[family_indexes]
        )
        family_model = make_model(
            positions, candidate, dispersion, shrinkage=0.0
        )
        family_result = family_model.fit(
            observed[family_indexes], expected[family_indexes]
        )
        model = make_model(positions, candidate, dispersion, shrinkage=50.0)
        model.fit(
            observed[tf_indexes],
            expected[tf_indexes],
            prior_profile=family_result.footprint_profile,
        )
        return FittedDetector("count", model, dispersion, source="strand")
    indexes = tf_indexes if candidate.training_pool == "tf" else family_indexes
    values = profiles[candidate.channel][indexes]
    observed, _expected = _count_arrays(profiles)
    weights = np.sqrt(np.maximum(observed[indexes].sum(axis=1), 1.0))
    model_seed = stable_seed(tf, motif_family, candidate.candidate_id, seed=seed)
    if candidate.family == "anchored-fda":
        coverage = observed[indexes].sum(axis=1)
        model = CovariateAnchoredFdaModel(
            max_components=20,
            anchor_strength=candidate.anchor_strength,
            seed=model_seed,
        ).fit(
            values,
            motif_score=sites.iloc[indexes]["motif_score"].to_numpy(dtype=float),
            accessibility=coverage,
            positions=positions,
            sample_weight=weights,
        )
    elif candidate.family == "residualized-fda":
        coverage = observed[indexes].sum(axis=1)
        model = CovariateResidualizedFdaModel(
            max_components=20,
            covariate_ridge=candidate.covariate_ridge,
            seed=model_seed,
        ).fit(
            values,
            motif_score=sites.iloc[indexes]["motif_score"].to_numpy(dtype=float),
            accessibility=coverage,
            positions=positions,
            sample_weight=weights,
        )
    elif candidate.family == "fda":
        model = FdaMixtureModel(max_components=20, seed=model_seed).fit(
            values, positions=positions, sample_weight=weights
        )
    elif candidate.family == "hybrid":
        model = HybridFdaGpModel(
            positions, max_components=20, seed=model_seed
        ).fit(values, sample_weight=weights)
    else:
        raise ValueError(f"unsupported strand candidate family: {candidate.family}")
    return FittedDetector(
        candidate.family,
        model,
        0.0,
        channel=candidate.channel,
        source="strand",
    )


def fit_dwm_detector(
    candidate: FunctionalCandidate,
    sites: pd.DataFrame,
    profiles: dict[str, np.ndarray],
    *,
    cell: str,
    tf: str,
    motif_family: str,
    positions: np.ndarray,
    seed: int,
) -> FittedDetector:
    observed, expected = _count_arrays(profiles)
    dispersion = estimate_nb_dispersion(observed, expected)
    tf_indexes = np.flatnonzero(sites["tf"].astype(str).to_numpy() == tf)
    family_indexes = np.flatnonzero(
        sites["motif_family"].astype(str).to_numpy() == motif_family
    )
    if len(tf_indexes) < 100 or len(family_indexes) < 100:
        raise ValueError(f"insufficient DWM training sites for {tf}/{motif_family}")
    if candidate.family in {"spline", "gp"}:
        global_prior = None
        family_prior = None
        if candidate.shrinkage > 0:
            limit = min(len(sites), 400)
            indexes = np.random.default_rng(
                stable_seed(cell, candidate.candidate_id, "global", seed=seed)
            ).choice(len(sites), size=limit, replace=False)
            global_model = _count_model(
                candidate, positions, dispersion, shrinkage=0.0
            )
            global_result = global_model.fit(
                observed[indexes],
                expected[indexes],
                motif_score=sites.iloc[indexes]["motif_score"].to_numpy(dtype=float),
                accessibility=observed[indexes].sum(axis=1),
            )
            global_prior = global_result.footprint_profile
            family_model = _count_model(candidate, positions, dispersion)
            family_result = family_model.fit(
                observed[family_indexes],
                expected[family_indexes],
                motif_score=sites.iloc[family_indexes]["motif_score"].to_numpy(
                    dtype=float
                ),
                accessibility=observed[family_indexes].sum(axis=1),
                prior_profile=global_prior,
            )
            family_prior = family_result.footprint_profile
        model = _count_model(candidate, positions, dispersion)
        model.fit(
            observed[tf_indexes],
            expected[tf_indexes],
            motif_score=sites.iloc[tf_indexes]["motif_score"].to_numpy(dtype=float),
            accessibility=observed[tf_indexes].sum(axis=1),
            prior_profile=family_prior if family_prior is not None else global_prior,
        )
        return FittedDetector("count", model, dispersion, source="dwm")
    residual = deviance_profiles(observed, expected, dispersion)
    weights = np.sqrt(np.maximum(observed[tf_indexes].sum(axis=1), 1.0))
    model_seed = stable_seed(cell, tf, candidate.candidate_id, seed=seed)
    if candidate.family == "fda":
        model = FdaMixtureModel(
            variance_threshold=candidate.variance_threshold,
            max_components=candidate.max_components,
            seed=model_seed,
        ).fit(
            residual[tf_indexes], positions=positions, sample_weight=weights
        )
    elif candidate.family == "hybrid":
        model = HybridFdaGpModel(
            positions,
            variance_threshold=candidate.variance_threshold,
            max_components=candidate.max_components,
            seed=model_seed,
        ).fit(residual[tf_indexes], sample_weight=weights)
    else:
        raise ValueError(f"unsupported DWM reference family: {candidate.family}")
    return FittedDetector(candidate.family, model, dispersion, source="dwm")


def predict_detector(
    fitted: FittedDetector,
    profiles: dict[str, np.ndarray],
    sites: pd.DataFrame | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    observed, expected = _count_arrays(profiles)
    total_signal = observed.sum(axis=1)
    if fitted.model_family == "count":
        assert isinstance(fitted.model, BiasAwareFunctionalMixture)
        shape_log_odds, _prior = fitted.model.predict_log_odds_components(
            observed, expected
        )
        probabilities = expit(np.clip(shape_log_odds, -40.0, 40.0))
        residual = deviance_profiles(observed, expected, fitted.dispersion)
    else:
        if fitted.source == "strand":
            residual = profiles[fitted.channel]
        else:
            residual = deviance_profiles(observed, expected, fitted.dispersion)
        if fitted.model_family == "anchored-fda":
            shape_log_odds, _anchor = fitted.model.predict_log_odds_components(
                residual
            )
            probabilities = expit(np.clip(shape_log_odds, -40.0, 40.0))
        elif fitted.model_family == "residualized-fda":
            if sites is None:
                raise ValueError("residualized FDA prediction requires motif sites")
            probabilities = fitted.model.predict_proba(
                residual,
                motif_score=sites["motif_score"].to_numpy(dtype=float),
                accessibility=total_signal,
            )
        else:
            probabilities = fitted.model.predict_proba(residual)
    return np.asarray(probabilities, dtype=float), total_signal, residual


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return float("nan"), float("nan")
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = (proportion + z * z / (2.0 * total)) / denominator
    radius = z * np.sqrt(
        proportion * (1.0 - proportion) / total + z * z / (4.0 * total * total)
    ) / denominator
    return float(max(0.0, center - radius)), float(min(1.0, center + radius))


def summarize_false_positives(
    probabilities: np.ndarray,
    total_signal: np.ndarray,
    valid: np.ndarray,
    *,
    threshold: float,
) -> tuple[dict[str, float | int], np.ndarray, np.ndarray]:
    probabilities = np.asarray(probabilities, dtype=float)
    total_signal = np.asarray(total_signal, dtype=float)
    valid = np.asarray(valid, dtype=bool) & np.isfinite(probabilities)
    informative = valid & (total_signal > 0)
    calls = informative & (probabilities >= threshold)
    n_valid = int(valid.sum())
    n_informative = int(informative.sum())
    n_calls = int(calls.sum())
    low, high = wilson_interval(n_calls, n_valid)
    informative_low, informative_high = wilson_interval(n_calls, n_informative)
    return (
        {
            "sites_valid": n_valid,
            "sites_informative": n_informative,
            "false_positive_calls": n_calls,
            "false_positive_rate": n_calls / n_valid if n_valid else np.nan,
            "false_positive_rate_lower_95": low,
            "false_positive_rate_upper_95": high,
            "informative_false_positive_rate": (
                n_calls / n_informative if n_informative else np.nan
            ),
            "informative_false_positive_rate_lower_95": informative_low,
            "informative_false_positive_rate_upper_95": informative_high,
            "mean_probability_valid": (
                float(np.mean(probabilities[valid])) if n_valid else np.nan
            ),
            "mean_probability_informative": (
                float(np.mean(probabilities[informative]))
                if n_informative
                else np.nan
            ),
        },
        informative,
        calls,
    )


def aggregate_profile_rows(
    residual: np.ndarray,
    informative: np.ndarray,
    calls: np.ndarray,
    positions: np.ndarray,
    metadata: dict[str, Any],
) -> list[dict[str, Any]]:
    normalized = normalize_functional_profiles(residual, positions)
    rows: list[dict[str, Any]] = []
    for name, mask in (("informative", informative), ("false_positive_call", calls)):
        if not np.any(mask):
            continue
        mean = np.mean(normalized[mask], axis=0)
        descriptors = profile_descriptors(mean, positions)
        for position, value in zip(positions, mean):
            rows.append(
                {
                    **metadata,
                    "group": name,
                    "sites": int(mask.sum()),
                    "position": int(position),
                    "mean_normalized_residual": float(value),
                    "aggregate_depletion": float(descriptors.depletion),
                    "aggregate_width": float(descriptors.width),
                    "aggregate_asymmetry": float(descriptors.asymmetry),
                }
            )
    return rows


def _candidate_lookup() -> tuple[dict[str, StrandCandidate], dict[str, FunctionalCandidate]]:
    strand = {candidate.candidate_id: candidate for candidate in strand_candidate_grid()}
    dwm = {candidate.candidate_id: candidate for candidate in dwm_candidate_grid("compact")}
    return strand, dwm


def evaluate(
    *,
    study: dict,
    policy: pd.DataFrame,
    strand_training_paths: dict[tuple[str, str], Path],
    naked_strand_paths: dict[tuple[str, str, str], Path],
    naked_dwm_paths: dict[tuple[str, str], Path],
    dwm_training_paths: dict[str, Path],
    dwm_base_run: Path | None,
    threshold: float,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    promoted = policy[policy["passes_development_gates"].astype(bool)].copy()
    routes = promoted.set_index("motif_family")
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[
        tasks["split"].eq("development")
        & tasks["motif_family"].astype(str).isin(routes.index.astype(str))
    ].copy()
    positions = np.arange(
        -int(study["profile_flank_bp"]),
        int(study["profile_flank_bp"]) + 1,
        dtype=float,
    )
    strand_candidates, dwm_candidates = _candidate_lookup()
    score_frames: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []
    aggregate_rows: list[dict[str, Any]] = []
    input_records: list[dict[str, str]] = []
    training_cache: dict[tuple[str, str], tuple] = {}
    dwm_training_cache: dict[str, tuple] = {}

    for task in tasks.sort_values(["cell", "tf"]).itertuples(index=False):
        cell = str(task.cell)
        tf = str(task.tf)
        motif_family = str(task.motif_family)
        route = routes.loc[motif_family]
        bias_configuration = str(route["candidate_bias_configuration"])
        candidate_id = str(route["candidate_id"])
        reference_id = str(route["reference_candidate_id"])
        if candidate_id not in strand_candidates:
            raise ValueError(f"unknown frozen strand candidate: {candidate_id}")
        if reference_id not in dwm_candidates:
            raise ValueError(f"unknown frozen DWM candidate: {reference_id}")
        training_key = (bias_configuration, cell)
        if training_key not in training_cache:
            training_path = strand_training_paths[training_key]
            training_cache[training_key] = load_strand_artifact(training_path, cell)
            input_records.append(
                {"purpose": "strand_training", "path": str(training_path), "sha256": file_sha256(training_path)}
            )
        train_sites, train_profiles, train_valid, train_hashes, _document = training_cache[
            training_key
        ]
        if not np.all(train_valid):
            train_sites = train_sites.loc[train_valid].reset_index(drop=True)
            train_profiles = {name: values[train_valid] for name, values in train_profiles.items()}
            train_hashes = train_hashes[train_valid]
        candidate_model = fit_strand_detector(
            strand_candidates[candidate_id],
            train_sites,
            train_profiles,
            tf=tf,
            motif_family=motif_family,
            positions=positions,
            seed=stable_seed(cell, tf, candidate_id, seed=seed),
        )
        if cell not in dwm_training_cache:
            dwm_training_cache[cell] = load_dwm_training_source(
                cell,
                artifact_paths=dwm_training_paths,
                base_run=dwm_base_run,
                flank=int(study["profile_flank_bp"]),
            )
            for record in dwm_training_cache[cell][4]:
                input_records.append({"purpose": "dwm_training", **record})
        dwm_sites, dwm_profiles, dwm_valid, dwm_hashes, _inputs = dwm_training_cache[cell]
        if not np.array_equal(train_hashes, dwm_hashes[train_valid]):
            raise ValueError(f"strand and DWM training sites differ for {cell}")
        reference_model = fit_dwm_detector(
            dwm_candidates[reference_id],
            dwm_sites.loc[dwm_valid].reset_index(drop=True),
            {name: values[dwm_valid] for name, values in dwm_profiles.items()},
            cell=cell,
            tf=tf,
            motif_family=motif_family,
            positions=positions,
            seed=stable_seed(cell, tf, reference_id, seed=seed),
        )

        replicates = sorted(
            replicate
            for model, artifact_cell, replicate in naked_strand_paths
            if model == bias_configuration and artifact_cell == cell
        )
        if not replicates:
            raise ValueError(f"no naked-DNA strand artifact for {bias_configuration}/{cell}")
        for replicate in replicates:
            strand_path = naked_strand_paths[(bias_configuration, cell, replicate)]
            dwm_path = naked_dwm_paths[(cell, replicate)]
            naked_sites, naked_strand, strand_valid, strand_hashes, _strand_doc = (
                load_strand_artifact(strand_path, cell)
            )
            dwm_naked_sites, naked_dwm, dwm_naked_valid, dwm_naked_hashes, _dwm_doc = (
                load_combined_artifact(dwm_path, cell)
            )
            if not np.array_equal(strand_hashes, dwm_naked_hashes):
                raise ValueError(
                    f"candidate and DWM naked-DNA sites differ for {cell}/{replicate}"
                )
            input_records.extend(
                [
                    {"purpose": "naked_strand", "path": str(strand_path), "sha256": file_sha256(strand_path)},
                    {"purpose": "naked_dwm", "path": str(dwm_path), "sha256": file_sha256(dwm_path)},
                ]
            )
            tf_mask = naked_sites["tf"].astype(str).eq(tf).to_numpy()
            common_valid = tf_mask & strand_valid & dwm_naked_valid
            for method, fitted, profiles, method_candidate_id in (
                ("frozen_policy_candidate", candidate_model, naked_strand, candidate_id),
                ("frozen_dwm_reference", reference_model, naked_dwm, reference_id),
            ):
                probabilities, total_signal, residual = predict_detector(
                    fitted, profiles, naked_sites
                )
                metrics, informative, calls = summarize_false_positives(
                    probabilities,
                    total_signal,
                    common_valid,
                    threshold=threshold,
                )
                metadata = {
                    "cell": cell,
                    "tf": tf,
                    "motif_family": motif_family,
                    "replicate": replicate,
                    "method": method,
                    "candidate_id": method_candidate_id,
                    "bias_configuration": (
                        bias_configuration if method == "frozen_policy_candidate" else "DWM"
                    ),
                }
                summary_rows.append({**metadata, **metrics})
                aggregate_rows.extend(
                    aggregate_profile_rows(
                        residual,
                        informative,
                        calls,
                        positions,
                        metadata,
                    )
                )
                selected = np.flatnonzero(tf_mask)
                score_frames.append(
                    pd.DataFrame(
                        {
                            **{key: value for key, value in metadata.items()},
                            "site_hash": strand_hashes[selected],
                            "TFBS_chr": naked_sites.loc[tf_mask, "TFBS_chr"].to_numpy(),
                            "TFBS_start": naked_sites.loc[tf_mask, "TFBS_start"].to_numpy(),
                            "TFBS_end": naked_sites.loc[tf_mask, "TFBS_end"].to_numpy(),
                            "total_signal": total_signal[selected],
                            "binding_probability": probabilities[selected],
                            "valid": common_valid[selected],
                            "informative": informative[selected],
                            "false_positive_call": calls[selected],
                        }
                    )
                )

    summaries = pd.DataFrame(summary_rows)
    scores = pd.concat(score_frames, ignore_index=True) if score_frames else pd.DataFrame()
    aggregates = pd.DataFrame(aggregate_rows)
    paired = summaries.pivot_table(
        index=["cell", "tf", "motif_family", "replicate"],
        columns="method",
        values=["false_positive_rate", "informative_false_positive_rate"],
    )
    paired.columns = ["__".join(column) for column in paired.columns]
    paired = paired.reset_index()
    primary_candidate = "false_positive_rate__frozen_policy_candidate"
    primary_reference = "false_positive_rate__frozen_dwm_reference"
    informative_candidate = "informative_false_positive_rate__frozen_policy_candidate"
    informative_reference = "informative_false_positive_rate__frozen_dwm_reference"
    paired["false_positive_rate_increase"] = paired[primary_candidate] - paired[primary_reference]
    paired["informative_false_positive_rate_increase"] = (
        paired[informative_candidate] - paired[informative_reference]
    )
    gates = study["promotion_gates"]
    maximum_rate = float(gates["maximum_naked_dna_false_positive_rate"])
    maximum_increase = float(gates["maximum_naked_dna_false_positive_rate_increase"])
    gate = {
        "schema": "fp-tools-naked-dna-functional-policy-gate-v1",
        "posterior_threshold": float(threshold),
        "minimum_signal_for_call": "strictly greater than zero",
        "all_site_gate_required": True,
        "informative_site_gate_required": True,
        "maximum_false_positive_rate": maximum_rate,
        "maximum_false_positive_rate_increase": maximum_increase,
        "pairs": int(len(paired)),
        "maximum_candidate_false_positive_rate": float(paired[primary_candidate].max()),
        "maximum_candidate_informative_false_positive_rate": float(
            paired[informative_candidate].max()
        ),
        "maximum_false_positive_rate_increase_observed": float(
            paired["false_positive_rate_increase"].max()
        ),
        "maximum_informative_false_positive_rate_increase_observed": float(
            paired["informative_false_positive_rate_increase"].max()
        ),
    }
    gate["passes"] = bool(
        gate["maximum_candidate_false_positive_rate"] <= maximum_rate
        and gate["maximum_candidate_informative_false_positive_rate"] <= maximum_rate
        and gate["maximum_false_positive_rate_increase_observed"] <= maximum_increase
        and gate["maximum_informative_false_positive_rate_increase_observed"]
        <= maximum_increase
    )
    gate["inputs"] = list(
        {record["path"]: record for record in input_records}.values()
    )
    return summaries, scores, aggregates, {"gate": gate, "paired": paired}


def render_profiles(aggregates: pd.DataFrame, output: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    output.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(output) as pdf:
        for (cell, tf, replicate), group in aggregates.groupby(
            ["cell", "tf", "replicate"], sort=True
        ):
            figure, axis = plt.subplots(figsize=(7.4, 4.4))
            for method, color in (
                ("frozen_dwm_reference", "#555555"),
                ("frozen_policy_candidate", "#C23B33"),
            ):
                selected = group[
                    group["method"].eq(method) & group["group"].eq("informative")
                ].sort_values("position")
                if selected.empty:
                    continue
                axis.plot(
                    selected["position"],
                    selected["mean_normalized_residual"],
                    linewidth=1.8,
                    color=color,
                    label=method.replace("frozen_", "").replace("_", " "),
                )
            axis.axhline(0, color="#888888", linewidth=0.7)
            axis.axvline(0, color="#888888", linewidth=0.7, linestyle="--")
            axis.set_title(f"Naked DNA: {cell} {tf} — {replicate}")
            axis.set_xlabel("Position from motif center (bp)")
            axis.set_ylabel("Mean normalized residual")
            axis.legend(frameon=False)
            figure.tight_layout()
            pdf.savefig(figure)
            plt.close(figure)


def promotion_false_positive_table(summaries: pd.DataFrame) -> pd.DataFrame:
    """Map detailed naked-DNA rows onto the final promotion-auditor schema."""

    required = {
        "cell",
        "tf",
        "motif_family",
        "replicate",
        "method",
        "false_positive_rate",
        "informative_false_positive_rate",
    }
    missing = required.difference(summaries.columns)
    if missing:
        raise ValueError(
            "naked-DNA summary lacks promotion columns: "
            + ", ".join(sorted(missing))
        )
    observed = set(summaries["method"].astype(str))
    if observed != set(PROMOTION_METHODS):
        raise ValueError(
            "naked-DNA summary methods differ from the frozen paired policy: "
            f"{sorted(observed)}"
        )
    output = summaries[
        [
            "cell",
            "tf",
            "motif_family",
            "replicate",
            "method",
            "false_positive_rate",
            "informative_false_positive_rate",
            "sites_valid",
            "sites_informative",
        ]
    ].copy()
    output.insert(
        4,
        "candidate_id",
        output["method"].astype(str).map(PROMOTION_METHODS),
    )
    output = output.rename(columns={"method": "naked_dna_method"})
    return output.sort_values(
        ["cell", "tf", "replicate", "candidate_id"], kind="mergesort"
    ).reset_index(drop=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--study",
        type=Path,
        default=Path("benchmarks/manifests/footprint_functional_v1.spec.json"),
    )
    parser.add_argument(
        "--policy",
        type=Path,
        default=Path("benchmarks/manifests/compact/functional_detector_policy_v1.tsv"),
    )
    parser.add_argument(
        "--dwm-base-run",
        type=Path,
        help="Legacy development run containing DWM profile caches.",
    )
    parser.add_argument(
        "--dwm-training-artifact",
        action="append",
        type=parse_dwm_training_artifact,
        default=[],
        metavar="CELL,JSON",
        help="Explicit label-free combined DWM training artifact; repeat by cell.",
    )
    parser.add_argument(
        "--strand-training-artifact",
        action="append",
        type=parse_training_artifact,
        required=True,
        metavar="MODEL,CELL,JSON",
    )
    parser.add_argument(
        "--naked-strand-artifact",
        action="append",
        type=parse_naked_artifact,
        required=True,
        metavar="MODEL,CELL,REPLICATE,JSON",
    )
    parser.add_argument(
        "--naked-dwm-artifact",
        action="append",
        type=parse_dwm_artifact,
        required=True,
        metavar="CELL,REPLICATE,JSON",
    )
    parser.add_argument("--posterior-threshold", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    if not 0 < args.posterior_threshold < 1:
        raise SystemExit("--posterior-threshold must be between zero and one")
    study = json.loads(args.study.read_text(encoding="utf-8"))
    policy = pd.read_csv(args.policy, sep="\t")
    training_paths = {
        (model, cell): path for model, cell, path in args.strand_training_artifact
    }
    naked_strand_paths = {
        (model, cell, replicate): path
        for model, cell, replicate, path in args.naked_strand_artifact
    }
    naked_dwm_paths = {
        (cell, replicate): path for cell, replicate, path in args.naked_dwm_artifact
    }
    dwm_training_paths = dict(args.dwm_training_artifact)
    if len(dwm_training_paths) != len(args.dwm_training_artifact):
        raise SystemExit("duplicate --dwm-training-artifact cells")
    if args.dwm_base_run is None and not dwm_training_paths:
        raise SystemExit("provide --dwm-base-run or --dwm-training-artifact")
    summaries, scores, aggregates, result = evaluate(
        study=study,
        policy=policy,
        strand_training_paths=training_paths,
        naked_strand_paths=naked_strand_paths,
        naked_dwm_paths=naked_dwm_paths,
        dwm_training_paths=dwm_training_paths,
        dwm_base_run=args.dwm_base_run,
        threshold=args.posterior_threshold,
        seed=args.seed,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary_path = args.outdir / "naked_dna_false_positive_summary.tsv"
    scores_path = args.outdir / "naked_dna_site_scores.tsv.gz"
    aggregate_path = args.outdir / "naked_dna_aggregate_profiles.tsv.gz"
    paired_path = args.outdir / "naked_dna_policy_vs_dwm.tsv"
    promotion_path = args.outdir / "naked_dna_promotion_false_positive_rates.tsv"
    plot_path = args.outdir / "naked_dna_aggregate_profiles.pdf"
    promotion_rates = promotion_false_positive_table(summaries)
    summaries.to_csv(summary_path, sep="\t", index=False)
    scores.to_csv(scores_path, sep="\t", index=False)
    aggregates.to_csv(aggregate_path, sep="\t", index=False)
    result["paired"].to_csv(paired_path, sep="\t", index=False)
    promotion_rates.to_csv(promotion_path, sep="\t", index=False)
    render_profiles(aggregates, plot_path)
    gate = result["gate"]
    gate.update(
        {
            "study": str(args.study),
            "study_sha256": file_sha256(args.study),
            "policy": str(args.policy),
            "policy_sha256": file_sha256(args.policy),
            "dwm_base_run": str(args.dwm_base_run),
            "seed": int(args.seed),
            "outputs": {
                path.name: {"path": str(path), "sha256": file_sha256(path)}
                for path in (
                    summary_path,
                    scores_path,
                    aggregate_path,
                    paired_path,
                    promotion_path,
                    plot_path,
                )
            },
        }
    )
    gate_path = args.outdir / "naked_dna_functional_policy_gate.json"
    gate_path.write_text(
        json.dumps(gate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(result["paired"].to_string(index=False))
    print(json.dumps({key: value for key, value in gate.items() if key != "inputs"}, indent=2, sort_keys=True))
    return 0 if gate["passes"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
