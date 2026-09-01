#!/usr/bin/env python3
"""Refit validation-selected functional detectors and freeze them before test.

The winner table may contain validation metrics, but model fitting reads only
the separate label-free chr1--15 training artifacts. The resulting policy is
immutable and records checksums for every fitted model and training input.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import stable_seed  # noqa: E402
from evaluate_strand_functional_templates import load_artifact  # noqa: E402
from evaluate_strand_label_free_models import (  # noqa: E402
    Candidate,
    candidate_grid,
    file_sha256,
    parse_artifact,
    validate_unlabeled_training_sites,
)
from fp_tools.tools.functional_footprints import (  # noqa: E402
    BiasAwareFunctionalMixture,
    ConditionalMultinomialMixture,
    CovariateAnchoredFdaModel,
    CovariateResidualizedFdaModel,
    FdaMixtureModel,
    HybridFdaGpModel,
)
from fp_tools.tools.parametric_bias import estimate_nb_dispersion  # noqa: E402


POLICY_SCHEMA = "fp-tools-label-free-functional-policy-v1"


def safe_token(value: object) -> str:
    token = re.sub(r"[^A-Za-z0-9_-]+", "_", str(value)).strip("_-")
    if not token:
        raise ValueError("empty model filename token")
    return token


def immutable_write_json(path: Path, document: dict) -> None:
    rendered = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if path.exists() and path.read_text(encoding="utf-8") != rendered:
        raise ValueError(f"refusing to replace a different immutable policy: {path}")
    path.write_text(rendered, encoding="utf-8")


def validate_selection(
    winners_path: Path,
    selection_manifest_path: Path,
) -> tuple[pd.DataFrame, dict]:
    manifest = json.loads(selection_manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "fp-tools-strand-label-free-evaluation-v1":
        raise ValueError("unsupported functional selection manifest")
    if manifest.get("locked_test_labels_read") is not False:
        raise ValueError("selection manifest does not certify unopened test labels")
    if manifest.get("training_labels_used") is not False:
        raise ValueError("selection manifest used labels for model fitting")
    output = manifest.get("outputs", {}).get("winners", {})
    if output.get("sha256") != file_sha256(winners_path):
        raise ValueError("winner table checksum does not match selection manifest")
    winners = pd.read_csv(winners_path, sep="\t")
    required = {
        "cell",
        "tf",
        "motif_family",
        "bias_configuration",
        "candidate_id",
        "family",
        "status",
        "converged",
    }
    missing = required.difference(winners.columns)
    if missing:
        raise ValueError("winner table is missing columns: " + ", ".join(sorted(missing)))
    if winners.empty or not winners["status"].astype(str).eq("ok").all():
        raise ValueError("winner table contains no usable winners")
    if not winners["converged"].astype(bool).all():
        raise ValueError("winner table contains a non-converged model")
    duplicate = winners.duplicated(["cell", "tf", "bias_configuration"])
    if duplicate.any():
        raise ValueError("winner table contains duplicate TF policies")
    return winners, manifest


def candidate_from_row(row: pd.Series) -> Candidate:
    candidates = {candidate.candidate_id: candidate for candidate in candidate_grid()}
    candidate_id = str(row["candidate_id"])
    if candidate_id not in candidates:
        raise ValueError(f"unknown functional candidate: {candidate_id}")
    candidate = candidates[candidate_id]
    for field, expected in asdict(candidate).items():
        if field == "candidate_id":
            continue
        observed = row.get(field)
        if isinstance(expected, float):
            if not np.isclose(float(observed), expected):
                raise ValueError(
                    f"winner {candidate_id} has inconsistent {field}: {observed}"
                )
        elif ("" if pd.isna(observed) else str(observed)) != str(expected):
            raise ValueError(
                f"winner {candidate_id} has inconsistent {field}: {observed}"
            )
    return candidate


def fit_model(
    candidate: Candidate,
    *,
    bias_configuration: str,
    cell: str,
    tf: str,
    motif_family: str,
    sites: pd.DataFrame,
    profiles: dict[str, np.ndarray],
    positions: np.ndarray,
    seed: int,
):
    tf_indexes = np.flatnonzero(sites["tf"].astype(str).to_numpy() == tf)
    family_indexes = np.flatnonzero(
        sites["motif_family"].astype(str).to_numpy() == motif_family
    )
    if len(tf_indexes) < 100 or len(family_indexes) < 100:
        raise ValueError(f"insufficient label-free training sites for {cell} {tf}")
    if candidate.family in {"count", "conditional"}:
        observed = profiles["plus_observed"] + profiles["minus_observed"]
        expected = profiles["plus_expected"] + profiles["minus_expected"]
        dispersion = estimate_nb_dispersion(
            observed[family_indexes], expected[family_indexes]
        )
        model_class = (
            BiasAwareFunctionalMixture
            if candidate.family == "count"
            else ConditionalMultinomialMixture
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
        )
        family_result = family_model.fit(
            observed[family_indexes],
            expected[family_indexes],
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
        )
        model.fit(
            observed[tf_indexes],
            expected[tf_indexes],
            prior_profile=family_result.footprint_profile,
        )
        return model, len(tf_indexes), len(family_indexes)

    indexes = tf_indexes if candidate.training_pool == "tf" else family_indexes
    training = profiles[candidate.channel][indexes]
    weights = np.sqrt(
        np.maximum(
            (
                profiles["plus_observed"][indexes]
                + profiles["minus_observed"][indexes]
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
    coverage = (
        profiles["plus_observed"][indexes]
        + profiles["minus_observed"][indexes]
    ).sum(axis=1)
    if candidate.family == "anchored-fda":
        model = CovariateAnchoredFdaModel(
            max_components=20,
            anchor_strength=candidate.anchor_strength,
            seed=model_seed,
        ).fit(
            training,
            motif_score=sites.iloc[indexes]["motif_score"].to_numpy(dtype=float),
            accessibility=coverage,
            positions=positions,
            sample_weight=weights,
        )
    elif candidate.family == "residualized-fda":
        model = CovariateResidualizedFdaModel(
            max_components=20,
            covariate_ridge=candidate.covariate_ridge,
            seed=model_seed,
        ).fit(
            training,
            motif_score=sites.iloc[indexes]["motif_score"].to_numpy(dtype=float),
            accessibility=coverage,
            positions=positions,
            sample_weight=weights,
        )
    elif candidate.family == "fda":
        model = FdaMixtureModel(max_components=20, seed=model_seed).fit(
            training,
            positions=positions,
            sample_weight=weights,
        )
    elif candidate.family == "hybrid":
        model = HybridFdaGpModel(
            positions,
            max_components=20,
            seed=model_seed,
        ).fit(training, sample_weight=weights)
    else:
        raise ValueError(f"unsupported policy candidate: {candidate.family}")
    return model, len(tf_indexes), len(family_indexes)


def load_frozen_model(family: str, path: Path):
    if family == "count":
        return BiasAwareFunctionalMixture.load(path)
    if family == "conditional":
        return ConditionalMultinomialMixture.load(path)
    if family == "fda":
        return FdaMixtureModel.load(path)
    if family == "hybrid":
        return HybridFdaGpModel.load(path)
    if family == "anchored-fda":
        return CovariateAnchoredFdaModel.load(path)
    if family == "residualized-fda":
        return CovariateResidualizedFdaModel.load(path)
    raise ValueError(f"unsupported frozen model family: {family}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--training-artifact",
        action="append",
        type=parse_artifact,
        required=True,
        metavar="MODEL,CELL,JSON",
    )
    parser.add_argument("--winners", type=Path, required=True)
    parser.add_argument("--selection-manifest", type=Path, required=True)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)

    winners, selection_manifest = validate_selection(
        args.winners,
        args.selection_manifest,
    )
    study = json.loads(args.study.read_text(encoding="utf-8"))
    if file_sha256(args.study) != selection_manifest.get("study_sha256"):
        raise ValueError("study checksum does not match selection manifest")
    if int(selection_manifest.get("seed", -1)) != args.seed:
        raise ValueError("seed does not match selection manifest")

    training_paths = {
        (model, cell): path for model, cell, path in args.training_artifact
    }
    required_artifacts = set(
        zip(
            winners["bias_configuration"].astype(str),
            winners["cell"].astype(str),
        )
    )
    if set(training_paths) != required_artifacts:
        raise ValueError("training artifacts do not exactly match winner cells/models")

    training = {}
    input_records = []
    for key in sorted(training_paths):
        path = training_paths[key]
        sites, profiles, document = load_artifact(path, key[1], study)
        validate_unlabeled_training_sites(sites, path)
        training[key] = (sites, profiles)
        input_records.append(
            {
                "bias_configuration": key[0],
                "cell": key[1],
                "artifact": {"path": str(path), "sha256": file_sha256(path)},
                "profiles": {
                    "path": str(document["profiles_npz"]),
                    "sha256": document["profiles_sha256"],
                },
                "sites": {
                    "path": str(document["sites"]),
                    "sha256": document["sites_sha256"],
                },
            }
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    model_dir = args.outdir / "models"
    model_dir.mkdir(parents=True, exist_ok=True)
    positions = np.arange(
        -int(study["profile_flank_bp"]),
        int(study["profile_flank_bp"]) + 1,
        dtype=float,
    )
    model_records = []
    with threadpool_limits(limits=1):
        for _, row in winners.sort_values(
            ["bias_configuration", "cell", "tf"], kind="mergesort"
        ).iterrows():
            candidate = candidate_from_row(row)
            cell = str(row["cell"])
            tf = str(row["tf"])
            motif_family = str(row["motif_family"])
            bias_configuration = str(row["bias_configuration"])
            sites, profiles = training[(bias_configuration, cell)]
            prefix = model_dir / "__".join(
                safe_token(value)
                for value in (bias_configuration, cell, tf, candidate.candidate_id)
            )
            npz_path = prefix.with_suffix(".npz")
            json_path = prefix.with_suffix(".json")
            if npz_path.exists() != json_path.exists():
                raise ValueError(f"incomplete resumable functional model: {prefix}")
            if npz_path.exists():
                load_frozen_model(candidate.family, npz_path)
                metadata = json.loads(json_path.read_text(encoding="utf-8")).get(
                    "metadata", {}
                )
                expected_metadata = {
                    "bias_configuration": bias_configuration,
                    "candidate_id": candidate.candidate_id,
                    "cell": cell,
                    "tf": tf,
                    "motif_family": motif_family,
                    "training_labels_used": False,
                    "training_split": "train",
                    "selection_split": "validation",
                }
                if metadata != expected_metadata:
                    raise ValueError(f"resumable model metadata differs: {prefix}")
                tf_sites = int(np.sum(sites["tf"].astype(str).eq(tf)))
                family_sites = int(
                    np.sum(sites["motif_family"].astype(str).eq(motif_family))
                )
            else:
                model, tf_sites, family_sites = fit_model(
                    candidate,
                    bias_configuration=bias_configuration,
                    cell=cell,
                    tf=tf,
                    motif_family=motif_family,
                    sites=sites,
                    profiles=profiles,
                    positions=positions,
                    seed=args.seed,
                )
                model.save(
                    prefix,
                    metadata={
                        "bias_configuration": bias_configuration,
                        "candidate_id": candidate.candidate_id,
                        "cell": cell,
                        "tf": tf,
                        "motif_family": motif_family,
                        "training_labels_used": False,
                        "training_split": "train",
                        "selection_split": "validation",
                    },
                )
            model_records.append(
                {
                    "bias_configuration": bias_configuration,
                    "cell": cell,
                    "tf": tf,
                    "motif_family": motif_family,
                    "candidate": asdict(candidate),
                    "model_type": candidate.family,
                    "model_npz": {
                        "path": str(npz_path),
                        "sha256": file_sha256(npz_path),
                    },
                    "model_json": {
                        "path": str(json_path),
                        "sha256": file_sha256(json_path),
                    },
                    "tf_training_sites": tf_sites,
                    "family_training_sites": family_sites,
                    "validation_metrics": {
                        key: float(row[key])
                        for key in (
                            "auroc",
                            "auprc",
                            "brier",
                            "functional_separation",
                            "selection_score",
                        )
                    },
                }
            )

    policy = {
        "schema": POLICY_SCHEMA,
        "selection_split": "validation",
        "training_split": "train",
        "test_labels_opened_before_policy": True,
        "test_labels_used_for_policy_selection": False,
        "global_test_label_access_reason": "prior parametric-factorization test",
        "training_labels_used": False,
        "models_refitted_on_test": False,
        "study": {"path": str(args.study), "sha256": file_sha256(args.study)},
        "selection_manifest": {
            "path": str(args.selection_manifest),
            "sha256": file_sha256(args.selection_manifest),
        },
        "winners": {"path": str(args.winners), "sha256": file_sha256(args.winners)},
        "seed": args.seed,
        "training_artifacts": input_records,
        "models": model_records,
    }
    canonical = json.dumps(policy, sort_keys=True, separators=(",", ":"))
    policy["policy_id"] = sha256(canonical.encode()).hexdigest()
    policy_path = args.outdir / "functional_policy.freeze.json"
    immutable_write_json(policy_path, policy)
    print(f"froze {len(model_records)} models: {policy_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
