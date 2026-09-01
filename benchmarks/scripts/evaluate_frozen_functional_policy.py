#!/usr/bin/env python3
"""Apply a frozen label-free functional policy to a no-refit test split."""

from __future__ import annotations

import argparse
from dataclasses import asdict
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

from build_strand_functional_profiles import site_hashes  # noqa: E402
from evaluate_functional_footprints import (  # noqa: E402
    binary_metrics,
    chromosome_split,
)
from evaluate_parametric_factorization import (  # noqa: E402
    align_baseline,
    block_bootstrap_delta,
    load_dwm_baseline,
    load_safe_configuration,
    orient_aligned_baseline,
    parse_name_path,
    residual_score,
)
from evaluate_strand_functional_templates import PROFILE_ARRAYS  # noqa: E402
from evaluate_strand_label_free_models import (  # noqa: E402
    Candidate,
    file_sha256,
    parse_artifact,
)
from fp_tools.tools.functional_footprints import (  # noqa: E402
    BiasAwareFunctionalMixture,
    CovariateAnchoredFdaModel,
    CovariateResidualizedFdaModel,
    FdaMixtureModel,
    HybridFdaGpModel,
    normalize_functional_profiles,
    profile_descriptors,
    standardized_functional_separation,
)
from fp_tools.tools.parametric_factorization import (  # noqa: E402
    FrozenParametricFactorization,
)
from freeze_label_free_functional_models import (  # noqa: E402
    POLICY_SCHEMA,
    immutable_write_json,
    load_frozen_model,
)


TEST_FREEZE_SCHEMA = "fp-tools-frozen-functional-test-inputs-v1"
TEST_RESULT_SCHEMA = "fp-tools-frozen-functional-test-results-v1"


def _canonical_id(document: dict, identifier: str) -> str:
    content = dict(document)
    observed = str(content.pop(identifier))
    canonical = json.dumps(content, sort_keys=True, separators=(",", ":"))
    expected = sha256(canonical.encode()).hexdigest()
    if observed != expected:
        raise ValueError(f"{identifier} does not match the frozen document")
    return observed


def validate_policy(path: Path) -> tuple[dict, list[tuple[dict, Candidate, object]]]:
    policy = json.loads(path.read_text(encoding="utf-8"))
    if policy.get("schema") != POLICY_SCHEMA:
        raise ValueError("unsupported label-free functional policy")
    _canonical_id(policy, "policy_id")
    if policy.get("training_labels_used") is not False:
        raise ValueError("functional policy does not certify label-free fitting")
    if policy.get("test_labels_used_for_policy_selection") is not False:
        raise ValueError("functional policy used test labels for selection")
    records = [
        policy["study"],
        policy["selection_manifest"],
        policy["winners"],
    ]
    for artifact in policy["training_artifacts"]:
        records.extend(
            [artifact["artifact"], artifact["profiles"], artifact["sites"]]
        )
    for record in records:
        if file_sha256(record["path"]) != record["sha256"]:
            raise ValueError(f"frozen policy input changed: {record['path']}")

    models = []
    for record in policy["models"]:
        for key in ("model_npz", "model_json"):
            artifact = record[key]
            if file_sha256(artifact["path"]) != artifact["sha256"]:
                raise ValueError(f"frozen model changed: {artifact['path']}")
        candidate = Candidate(**record["candidate"])
        model = load_frozen_model(
            candidate.family,
            Path(record["model_npz"]["path"]),
        )
        models.append((record, candidate, model))
    return policy, models


def preflight_test_artifact(
    path: Path,
    *,
    expected_cell: str,
) -> tuple[dict, dict[str, np.ndarray]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != "fp-tools-strand-functional-profiles-v1":
        raise ValueError(f"unsupported test profile artifact: {path}")
    metadata = document.get("metadata", {})
    if metadata.get("labels_used") is not False:
        raise ValueError(f"test profiles do not certify label-free construction: {path}")
    metadata_cell = metadata.get("cell")
    if metadata_cell is not None and str(metadata_cell) != expected_cell:
        raise ValueError(f"test profile cell does not match {expected_cell}: {path}")
    profiles_path = Path(document["profiles_npz"])
    sites_path = Path(document["sites"])
    if file_sha256(profiles_path) != document["profiles_sha256"]:
        raise ValueError(f"test profile checksum mismatch: {profiles_path}")
    if file_sha256(sites_path) != document["sites_sha256"]:
        raise ValueError(f"test site checksum mismatch: {sites_path}")
    with np.load(profiles_path, allow_pickle=False) as source:
        missing = sorted(set(PROFILE_ARRAYS + ("valid", "site_hash")).difference(source.files))
        if missing:
            raise ValueError(
                f"test artifact is missing arrays {', '.join(missing)}: {profiles_path}"
            )
        arrays = {
            name: np.asarray(source[name])
            for name in PROFILE_ARRAYS + ("valid", "site_hash")
        }
    if any(len(values) != int(document["sites_total"]) for values in arrays.values()):
        raise ValueError(f"test profile arrays have inconsistent lengths: {path}")
    return document, arrays


def write_test_input_freeze(
    path: Path,
    *,
    policy_path: Path,
    policy: dict,
    reference_configuration_path: Path,
    reference_configuration: dict,
    prior_test_manifest: Path,
    test_artifacts: dict[tuple[str, str], Path],
    test_documents: dict[tuple[str, str], dict],
    baselines: dict[str, Path],
) -> dict:
    if not prior_test_manifest.is_file():
        raise ValueError("prior test manifest is required for access-state provenance")
    inputs = []
    for key in sorted(test_artifacts):
        document = test_documents[key]
        inputs.extend(
            [
                {
                    "purpose": "test-profile-manifest",
                    "path": str(test_artifacts[key]),
                    "sha256": file_sha256(test_artifacts[key]),
                },
                {
                    "purpose": "test-profile-arrays",
                    "path": str(document["profiles_npz"]),
                    "sha256": document["profiles_sha256"],
                },
                {
                    "purpose": "test-site-labels",
                    "path": str(document["sites"]),
                    "sha256": document["sites_sha256"],
                },
            ]
        )
    for cell, baseline in sorted(baselines.items()):
        inputs.append(
            {
                "purpose": f"{cell}-conventional-DWM",
                "path": str(baseline),
                "sha256": file_sha256(baseline),
            }
        )
    document = {
        "schema": TEST_FREEZE_SCHEMA,
        "policy_id": policy["policy_id"],
        "policy": {"path": str(policy_path), "sha256": file_sha256(policy_path)},
        "reference_configuration_id": reference_configuration["configuration_id"],
        "reference_configuration": {
            "path": str(reference_configuration_path),
            "sha256": file_sha256(reference_configuration_path),
        },
        "prior_test_manifest": {
            "path": str(prior_test_manifest),
            "sha256": file_sha256(prior_test_manifest),
        },
        "test_labels_previously_opened_for_other_model": True,
        "test_labels_used_for_functional_policy_selection": False,
        "functional_models_refitted_on_test": False,
        "inputs": inputs,
    }
    canonical = json.dumps(document, sort_keys=True, separators=(",", ":"))
    document["test_input_id"] = sha256(canonical.encode()).hexdigest()
    immutable_write_json(path, document)
    return document


def load_test_sites(document: dict, arrays: dict[str, np.ndarray], study: dict) -> pd.DataFrame:
    sites = pd.read_csv(document["sites"], sep="\t").reset_index(drop=True)
    if not np.array_equal(site_hashes(sites), arrays["site_hash"]):
        raise ValueError("test site order does not match profile arrays")
    observed_split = sites["TFBS_chr"].map(
        lambda chromosome: chromosome_split(str(chromosome), study)
    )
    if set(observed_split.astype(str)) != {"test"}:
        raise ValueError("test artifacts contain non-test chromosomes")
    sites["chromosome_split"] = observed_split
    return sites


def candidate_score_and_profile(
    candidate: Candidate,
    model: object,
    arrays: dict[str, np.ndarray],
    indexes: np.ndarray,
    positions: np.ndarray,
    motif_score: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if candidate.family == "count":
        if not isinstance(model, BiasAwareFunctionalMixture):
            raise TypeError("count policy has the wrong serialized model type")
        observed = arrays["plus_observed"][indexes] + arrays["minus_observed"][indexes]
        expected = arrays["plus_expected"][indexes] + arrays["minus_expected"][indexes]
        log_odds, _prior = model.predict_log_odds_components(observed, expected)
        score = 1.0 / (1.0 + np.exp(-np.clip(log_odds, -30.0, 30.0)))
        aggregate = normalize_functional_profiles(
            arrays["combined_residual"][indexes], positions
        )
        assert model.result_ is not None
        fitted_profile = model.result_.footprint_profile
        return score, aggregate, fitted_profile
    aggregate = normalize_functional_profiles(
        arrays[candidate.channel][indexes], positions
    )
    if candidate.family == "fda":
        if not isinstance(model, FdaMixtureModel):
            raise TypeError("FDA policy has the wrong serialized model type")
        score = model.predict_proba(arrays[candidate.channel][indexes])
        components = model.component_profiles()
        assert model.binding_component_ is not None
        fitted_profile = (
            components[model.binding_component_]
            - components[1 - model.binding_component_]
        )
        return score, aggregate, fitted_profile
    if candidate.family == "anchored-fda":
        if not isinstance(model, CovariateAnchoredFdaModel):
            raise TypeError("anchored FDA policy has the wrong serialized model type")
        shape_log_odds, _anchor = model.predict_log_odds_components(
            arrays[candidate.channel][indexes]
        )
        score = 1.0 / (1.0 + np.exp(-np.clip(shape_log_odds, -40.0, 40.0)))
        return score, aggregate, model.profile_difference()
    if candidate.family == "residualized-fda":
        if not isinstance(model, CovariateResidualizedFdaModel):
            raise TypeError(
                "residualized FDA policy has the wrong serialized model type"
            )
        coverage = (
            arrays["plus_observed"][indexes]
            + arrays["minus_observed"][indexes]
        ).sum(axis=1)
        if motif_score is None:
            raise ValueError("residualized FDA requires motif-score covariates")
        score = model.predict_proba(
            arrays[candidate.channel][indexes],
            motif_score=motif_score,
            accessibility=coverage,
        )
        return score, aggregate, model.profile_difference()
    if candidate.family == "hybrid":
        if not isinstance(model, HybridFdaGpModel):
            raise TypeError("hybrid policy has the wrong serialized model type")
        score = model.predict_proba(arrays[candidate.channel][indexes])
        assert model.bound_mean_ is not None and model.unbound_mean_ is not None
        return score, aggregate, model.bound_mean_ - model.unbound_mean_
    raise ValueError(f"unsupported frozen candidate family: {candidate.family}")


def aggregate_curve(
    profiles: np.ndarray,
    labels: np.ndarray,
    chromosomes: np.ndarray,
    *,
    iterations: int,
    seed: int,
) -> dict[str, np.ndarray]:
    profiles = np.asarray(profiles, dtype=float)
    labels = np.asarray(labels, dtype=int)
    chromosomes = np.asarray(chromosomes, dtype=str)
    positive = profiles[labels == 1]
    negative = profiles[labels == 0]
    positive_mean = np.nanmean(positive, axis=0)
    negative_mean = np.nanmean(negative, axis=0)
    difference = positive_mean - negative_mean
    unique = np.asarray(sorted(set(chromosomes)))
    sufficient = {}
    for chromosome in unique:
        sufficient[chromosome] = {}
        member = chromosomes == chromosome
        for label in (0, 1):
            selected = profiles[member & (labels == label)]
            sufficient[chromosome][label] = (
                np.nansum(selected, axis=0),
                int(len(selected)),
            )
    rng = np.random.default_rng(seed)
    bootstrap = []
    for _ in range(iterations):
        draw = rng.choice(unique, size=len(unique), replace=True)
        means = []
        usable = True
        for label in (0, 1):
            total = np.zeros(profiles.shape[1], dtype=float)
            count = 0
            for chromosome in draw:
                value, number = sufficient[str(chromosome)][label]
                total += value
                count += number
            if count == 0:
                usable = False
                break
            means.append(total / count)
        if usable:
            bootstrap.append(means[1] - means[0])
    bootstrap_values = np.asarray(bootstrap, dtype=float)
    if len(bootstrap_values):
        lower, upper = np.nanquantile(bootstrap_values, [0.025, 0.975], axis=0)
    else:
        lower = np.full_like(difference, np.nan)
        upper = np.full_like(difference, np.nan)
    return {
        "positive_mean": positive_mean,
        "negative_mean": negative_mean,
        "difference": difference,
        "lower_95": lower,
        "upper_95": upper,
    }


def metric_record(
    *,
    record: dict,
    candidate: Candidate,
    method: str,
    labels: np.ndarray,
    score: np.ndarray,
    profiles: np.ndarray,
    positions: np.ndarray,
    prediction_seconds: float,
    minimum_sites_per_class: int,
    split_name: str = "test",
) -> dict:
    metrics = binary_metrics(labels, score)
    positive = int(np.sum(labels == 1))
    negative = int(np.sum(labels == 0))
    output = {
        "cell": record["cell"],
        "tf": record["tf"],
        "motif_family": record["motif_family"],
        "bias_configuration": record["bias_configuration"],
        "candidate_id": candidate.candidate_id,
        "candidate_family": candidate.family,
        "method": method,
        "split": split_name,
        "status": (
            "eligible"
            if min(positive, negative) >= minimum_sites_per_class
            else "underpowered"
        ),
        "n_sites": len(labels),
        "n_positive": positive,
        "n_negative": negative,
        "prediction_seconds": prediction_seconds,
        "functional_separation": standardized_functional_separation(
            profiles, labels, positions
        ),
        **{
            key: metrics[key]
            for key in ("auroc", "auprc", "brier", "prevalence")
        },
    }
    if np.all((score >= 0.0) & (score <= 1.0)):
        output["calibration_error"] = expected_calibration_error(labels, score)
    else:
        output["calibration_error"] = np.nan
    return output


def expected_calibration_error(
    labels: np.ndarray,
    probability: np.ndarray,
    bins: int = 10,
) -> float:
    edges = np.linspace(0.0, 1.0, bins + 1)
    value = 0.0
    for index in range(bins):
        selected = (probability >= edges[index]) & (
            probability <= edges[index + 1]
            if index == bins - 1
            else probability < edges[index + 1]
        )
        if np.any(selected):
            value += float(np.mean(selected)) * abs(
                float(np.mean(labels[selected])) - float(np.mean(probability[selected]))
            )
    return value


def site_score_frame(
    *,
    record: dict,
    candidate: Candidate,
    sites: pd.DataFrame,
    indexes: np.ndarray,
    candidate_score: np.ndarray,
    dwm_score: np.ndarray,
    direct_score: np.ndarray,
) -> pd.DataFrame:
    """Return one auditable score row per evaluated motif occurrence."""

    selected = sites.iloc[indexes].reset_index(drop=True).copy()
    if not (
        len(selected)
        == len(candidate_score)
        == len(dwm_score)
        == len(direct_score)
    ):
        raise ValueError("site metadata and frozen score arrays have different lengths")
    output = pd.DataFrame(
        {
            "cell": str(record["cell"]),
            "tf": str(record["tf"]),
            "motif_family": str(record["motif_family"]),
            "bias_configuration": str(record["bias_configuration"]),
            "candidate_id": candidate.candidate_id,
            "artifact_index": np.asarray(indexes, dtype=int),
            "TFBS_chr": selected["TFBS_chr"].astype(str),
            "TFBS_start": selected["TFBS_start"].to_numpy(dtype=int),
            "TFBS_end": selected["TFBS_end"].to_numpy(dtype=int),
            "TFBS_strand": selected["TFBS_strand"].astype(str),
            "motif_score": selected["motif_score"].to_numpy(dtype=float),
            "accessibility": selected["accessibility"].to_numpy(dtype=float),
            "label": selected["chip_label"].to_numpy(dtype=int),
            "candidate_probability": np.asarray(candidate_score, dtype=float),
            "dwm_score": np.asarray(dwm_score, dtype=float),
            "direct_score": np.asarray(direct_score, dtype=float),
        }
    )
    output["log_accessibility"] = np.log1p(
        np.maximum(output["accessibility"].to_numpy(dtype=float), 0.0)
    )
    if "motif" in selected:
        output["motif_id"] = selected["motif"].astype(str)
    elif "motif_id" in selected:
        output["motif_id"] = selected["motif_id"].astype(str)
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument(
        "--test-artifact",
        action="append",
        type=parse_artifact,
        required=True,
        metavar="MODEL,CELL,JSON",
    )
    parser.add_argument(
        "--dwm-baseline",
        action="append",
        type=parse_name_path,
        required=True,
        metavar="CELL=NPZ",
    )
    parser.add_argument("--reference-configuration", type=Path, required=True)
    parser.add_argument("--prior-test-manifest", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--minimum-sites-per-class", type=int, default=200)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--aggregate-bootstrap-iterations", type=int, default=500)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)

    policy, models = validate_policy(args.policy)
    study = json.loads(
        Path(policy["study"]["path"]).read_text(encoding="utf-8")
    )
    reference_configuration = load_safe_configuration(args.reference_configuration)
    reference_model = FrozenParametricFactorization.load(
        reference_configuration["factorization_model"]["path"]
    )
    dispersion = float(reference_model.total_dispersion_)

    test_paths = {(model, cell): path for model, cell, path in args.test_artifact}
    baseline_paths = dict(args.dwm_baseline)
    policy_keys = {
        (str(record["bias_configuration"]), str(record["cell"]))
        for record, _candidate, _model in models
    }
    if set(test_paths) != policy_keys:
        raise ValueError("test artifacts do not exactly match policy cells/models")
    if set(baseline_paths) != {cell for _model, cell in policy_keys}:
        raise ValueError("DWM baselines do not exactly match policy cells")

    test_documents = {}
    test_arrays = {}
    for key, path in sorted(test_paths.items()):
        document, arrays = preflight_test_artifact(path, expected_cell=key[1])
        test_documents[key] = document
        test_arrays[key] = arrays
    baselines = {}
    for cell, path in sorted(baseline_paths.items()):
        baselines[cell], _inputs = load_dwm_baseline(path)

    args.outdir.mkdir(parents=True, exist_ok=True)
    test_input_path = args.outdir / "frozen_functional_test_inputs.freeze.json"
    test_input_freeze = write_test_input_freeze(
        test_input_path,
        policy_path=args.policy,
        policy=policy,
        reference_configuration_path=args.reference_configuration,
        reference_configuration=reference_configuration,
        prior_test_manifest=args.prior_test_manifest,
        test_artifacts=test_paths,
        test_documents=test_documents,
        baselines=baseline_paths,
    )

    sites_by_key = {
        key: load_test_sites(test_documents[key], test_arrays[key], study)
        for key in sorted(test_paths)
    }
    aligned_dwm = {}
    valid_by_key = {}
    for key in sorted(test_paths):
        arrays = test_arrays[key]
        sites = sites_by_key[key]
        expected, baseline_valid = align_baseline(arrays, baselines[key[1]])
        aligned_dwm[key] = orient_aligned_baseline(expected, baselines[key[1]], sites)
        valid_by_key[key] = (
            arrays["valid"].astype(bool)
            & baseline_valid
            & np.isfinite(aligned_dwm[key]).all(axis=1)
            & np.isfinite(arrays["combined_residual"]).all(axis=1)
        )

    metrics_rows = []
    bootstrap_rows = []
    curve_rows = []
    site_score_frames = []
    with threadpool_limits(limits=1):
        for record, candidate, model in models:
            key = (str(record["bias_configuration"]), str(record["cell"]))
            sites = sites_by_key[key]
            arrays = test_arrays[key]
            selected = (
                sites["tf"].astype(str).eq(str(record["tf"])).to_numpy()
                & valid_by_key[key]
            )
            indexes = np.flatnonzero(selected)
            if not len(indexes):
                continue
            labels = sites.iloc[indexes]["chip_label"].to_numpy(dtype=int)
            if np.unique(labels).size != 2:
                continue
            positions = np.arange(arrays["plus_observed"].shape[1], dtype=float)
            positions -= arrays["plus_observed"].shape[1] // 2
            observed = arrays["plus_observed"][indexes] + arrays["minus_observed"][indexes]
            started = perf_counter()
            candidate_score, candidate_profiles, fitted_profile = candidate_score_and_profile(
                candidate,
                model,
                arrays,
                indexes,
                positions,
                motif_score=sites.iloc[indexes]["motif_score"].to_numpy(dtype=float),
            )
            candidate_seconds = perf_counter() - started
            dwm_score, dwm_profiles = residual_score(
                observed,
                aligned_dwm[key][indexes],
                positions,
                "deviance",
                dispersion,
            )
            direct_expected = arrays["plus_expected"][indexes] + arrays["minus_expected"][indexes]
            direct_score, direct_profiles = residual_score(
                observed,
                direct_expected,
                positions,
                "deviance",
                dispersion,
            )
            site_score_frames.append(
                site_score_frame(
                    record=record,
                    candidate=candidate,
                    sites=sites,
                    indexes=indexes,
                    candidate_score=candidate_score,
                    dwm_score=dwm_score,
                    direct_score=direct_score,
                )
            )
            dwm_aggregate = normalize_functional_profiles(dwm_profiles, positions)
            direct_aggregate = normalize_functional_profiles(
                direct_profiles, positions
            )
            method_values = (
                (
                    "DWM_conventional_geometry",
                    dwm_score,
                    dwm_aggregate,
                    0.0,
                    None,
                ),
                (
                    "LOG21_direct_geometry",
                    direct_score,
                    direct_aggregate,
                    0.0,
                    None,
                ),
                (
                    f"frozen_{candidate.candidate_id}",
                    candidate_score,
                    candidate_profiles,
                    candidate_seconds,
                    fitted_profile,
                ),
            )
            for method, score, profiles, seconds, model_profile in method_values:
                metrics_rows.append(
                    metric_record(
                        record=record,
                        candidate=candidate,
                        method=method,
                        labels=labels,
                        score=score,
                        profiles=profiles,
                        positions=positions,
                        prediction_seconds=seconds,
                        minimum_sites_per_class=args.minimum_sites_per_class,
                    )
                )
                curve = aggregate_curve(
                    profiles,
                    labels,
                    sites.iloc[indexes]["TFBS_chr"].astype(str).to_numpy(),
                    iterations=args.aggregate_bootstrap_iterations,
                    seed=args.seed,
                )
                descriptors = asdict(profile_descriptors(curve["difference"], positions))
                for offset, position in enumerate(positions.astype(int)):
                    curve_rows.append(
                        {
                            "cell": record["cell"],
                            "tf": record["tf"],
                            "motif_family": record["motif_family"],
                            "candidate_id": candidate.candidate_id,
                            "method": method,
                            "position": position,
                            "positive_mean": curve["positive_mean"][offset],
                            "negative_mean": curve["negative_mean"][offset],
                            "positive_minus_negative": curve["difference"][offset],
                            "lower_95": curve["lower_95"][offset],
                            "upper_95": curve["upper_95"][offset],
                            "frozen_model_profile": (
                                np.nan if model_profile is None else model_profile[offset]
                            ),
                            **descriptors,
                        }
                    )
            positive = int(np.sum(labels == 1))
            negative = int(np.sum(labels == 0))
            if min(positive, negative) >= args.minimum_sites_per_class:
                task_sites = sites.iloc[indexes].reset_index(drop=True)
                for method, score in (
                    ("LOG21_direct_geometry", direct_score),
                    (f"frozen_{candidate.candidate_id}", candidate_score),
                ):
                    bootstrap_rows.append(
                        {
                            "cell": record["cell"],
                            "tf": record["tf"],
                            "motif_family": record["motif_family"],
                            "candidate_id": candidate.candidate_id,
                            "method": method,
                            **block_bootstrap_delta(
                                task_sites,
                                score,
                                dwm_score,
                                iterations=args.bootstrap_iterations,
                                seed=args.seed,
                            ),
                        }
                    )

    metrics = pd.DataFrame(metrics_rows)
    baseline = metrics[metrics["method"] == "DWM_conventional_geometry"][
        ["cell", "tf", "auroc", "auprc", "functional_separation"]
    ].rename(
        columns={
            "auroc": "dwm_auroc",
            "auprc": "dwm_auprc",
            "functional_separation": "dwm_functional_separation",
        }
    )
    metrics = metrics.merge(baseline, on=["cell", "tf"], validate="many_to_one")
    metrics["auroc_gain_over_dwm"] = metrics["auroc"] - metrics["dwm_auroc"]
    metrics["relative_auprc_gain_over_dwm"] = (
        metrics["auprc"] - metrics["dwm_auprc"]
    ) / metrics["dwm_auprc"].clip(lower=1e-8)
    metrics["functional_separation_relative_change_over_dwm"] = (
        metrics["functional_separation"]
        / metrics["dwm_functional_separation"].clip(lower=1e-8)
        - 1.0
    )
    metrics_path = args.outdir / "frozen_functional_test_metrics.tsv"
    bootstrap_path = args.outdir / "frozen_functional_test_bootstrap.tsv"
    curves_path = args.outdir / "frozen_functional_test_profiles.tsv.gz"
    site_scores_path = args.outdir / "frozen_functional_test_site_scores.tsv.gz"
    metrics.to_csv(metrics_path, sep="\t", index=False)
    pd.DataFrame(bootstrap_rows).to_csv(bootstrap_path, sep="\t", index=False)
    pd.DataFrame(curve_rows).to_csv(curves_path, sep="\t", index=False)
    pd.concat(site_score_frames, ignore_index=True).to_csv(
        site_scores_path,
        sep="\t",
        index=False,
    )
    manifest = {
        "schema": TEST_RESULT_SCHEMA,
        "policy_id": policy["policy_id"],
        "test_input_id": test_input_freeze["test_input_id"],
        "test_input_freeze": {
            "path": str(test_input_path),
            "sha256": file_sha256(test_input_path),
        },
        "models_refitted_on_test": False,
        "test_labels_used_for_policy_selection": False,
        "test_labels_previously_opened_for_other_model": True,
        "minimum_sites_per_class": args.minimum_sites_per_class,
        "dispersion": dispersion,
        "outputs": {
            "metrics": {"path": str(metrics_path), "sha256": file_sha256(metrics_path)},
            "bootstrap": {
                "path": str(bootstrap_path),
                "sha256": file_sha256(bootstrap_path),
            },
            "profiles": {"path": str(curves_path), "sha256": file_sha256(curves_path)},
            "site_scores": {
                "path": str(site_scores_path),
                "sha256": file_sha256(site_scores_path),
            },
        },
    }
    manifest_path = args.outdir / "frozen_functional_test_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    candidate_metrics = metrics[metrics["method"].str.startswith("frozen_")]
    columns = [
        "cell",
        "tf",
        "status",
        "candidate_id",
        "auroc_gain_over_dwm",
        "relative_auprc_gain_over_dwm",
        "functional_separation_relative_change_over_dwm",
    ]
    print(candidate_metrics[columns].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
