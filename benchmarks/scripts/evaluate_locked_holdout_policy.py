#!/usr/bin/env python3
"""Evaluate the frozen TF-family policy on locked ENCODE holdouts.

This program has two deliberately separate modes.  Freeze mode validates only
label-free motif/profile artifacts and records their hashes.  Evaluation mode
requires both that freeze and an explicit ``--unlock-holdout`` flag before it
opens any ChIP peak file.  Mixture models are fitted on chr1--15 without ChIP
labels; labels are used only for matched scoring on chr19--22 and chrX.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import md5, sha256  # noqa: S324 - MD5 verifies published ENCODE metadata
import json
from pathlib import Path
import sys
from time import perf_counter
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from build_footprint_site_labels import (  # noqa: E402
    label_sites,
    propensity_match,
    read_peaks,
)
from build_strand_functional_profiles import site_hashes  # noqa: E402
from evaluate_functional_footprints import binary_metrics, stable_seed  # noqa: E402
from evaluate_strand_label_free_models import (  # noqa: E402
    Candidate as StrandCandidate,
    candidate_grid as strand_candidate_grid,
)
from fp_tools.tools.functional_footprints import (  # noqa: E402
    BiasAwareFunctionalMixture,
    FdaMixtureModel,
    HybridFdaGpModel,
    deviance_profiles,
    normalize_functional_profiles,
    profile_descriptors,
    standardized_functional_separation,
)
from fp_tools.tools.parametric_bias import estimate_nb_dispersion  # noqa: E402
from fp_tools.utils.fasta import open_fasta  # noqa: E402
from search_functional_model_grid import (  # noqa: E402
    FunctionalCandidate,
    _count_model,
    candidate_grid as combined_candidate_grid,
)


FREEZE_SCHEMA = "fp-tools-locked-holdout-evaluation-freeze-v1"
RESULT_SCHEMA = "fp-tools-locked-holdout-policy-evaluation-v1"
PROMOTION_CANDIDATE = "frozen_policy"
PROMOTION_REFERENCE = "DWM_reference"
ARTIFACT_SCHEMAS = {
    "fp-tools-combined-functional-profiles-v1",
    "fp-tools-strand-functional-profiles-v1",
}
MATCH_FEATURES = (
    "motif_score",
    "log_accessibility",
    "gc_fraction",
    "peak_position_signed",
    "peak_position_abs",
)
TEST_SPLIT = "test"
TRAIN_SPLIT = "train"


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_md5(path: str | Path) -> str:
    digest = md5()  # noqa: S324 - checksum comparison, not cryptographic security
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_hash(document: dict[str, Any]) -> str:
    payload = json.dumps(document, sort_keys=True, separators=(",", ":")).encode()
    return sha256(payload).hexdigest()


def parse_key_path(value: str) -> tuple[str, str, Path]:
    fields = value.split(",", 2)
    if len(fields) != 3 or not all(fields):
        raise argparse.ArgumentTypeError("artifact must use MODEL,CELL,JSON")
    return fields[0], fields[1], Path(fields[2])


def parse_replicate_key_path(value: str) -> tuple[str, str, str, Path]:
    fields = value.split(",", 3)
    if len(fields) != 4 or not all(fields):
        raise argparse.ArgumentTypeError(
            "replicate artifact must use REPLICATE,MODEL,CELL,JSON"
        )
    return fields[0], fields[1], fields[2], Path(fields[3])


def parse_cell_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("site source must use CELL=TSV")
    cell, path = value.split("=", 1)
    if not cell or not path:
        raise argparse.ArgumentTypeError("site source must use CELL=TSV")
    return cell, Path(path)


def forbidden_label_columns(frame: pd.DataFrame) -> list[str]:
    return [
        column
        for column in frame.columns
        if "chip" in column.lower() or "label" in column.lower()
    ]


@dataclass
class Artifact:
    model: str
    cell: str
    path: Path
    document: dict[str, Any]
    sites: pd.DataFrame
    profiles: dict[str, np.ndarray]

    @property
    def schema(self) -> str:
        return str(self.document["schema"])

    @property
    def valid(self) -> np.ndarray:
        return np.asarray(self.profiles["valid"], dtype=bool)

    def observed_expected(self) -> tuple[np.ndarray, np.ndarray]:
        if self.schema == "fp-tools-combined-functional-profiles-v1":
            return self.profiles["observed"], self.profiles["expected"]
        return (
            self.profiles["plus_observed"] + self.profiles["minus_observed"],
            self.profiles["plus_expected"] + self.profiles["minus_expected"],
        )


@dataclass
class ModelOutput:
    probabilities: np.ndarray
    evaluation_profiles: np.ndarray
    footprint_profile: np.ndarray
    converged: bool
    iterations: int
    fit_seconds: float
    dispersion: float


def _resolve_recorded_path(document_path: Path, value: str) -> Path:
    path = Path(value)
    if path.is_file():
        return path
    relative = document_path.parent / path
    if relative.is_file():
        return relative
    raise FileNotFoundError(path)


def load_artifact(model: str, cell: str, path: Path) -> Artifact:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") not in ARTIFACT_SCHEMAS:
        raise ValueError(f"{path} has unsupported artifact schema {document.get('schema')!r}")
    sites_path = _resolve_recorded_path(path, str(document["sites"]))
    profiles_path = _resolve_recorded_path(path, str(document["profiles_npz"]))
    if file_sha256(sites_path) != str(document["sites_sha256"]):
        raise ValueError(f"site checksum mismatch for {path}")
    if file_sha256(profiles_path) != str(document["profiles_sha256"]):
        raise ValueError(f"profile checksum mismatch for {path}")
    sites = pd.read_csv(sites_path, sep="\t").reset_index(drop=True)
    forbidden = forbidden_label_columns(sites)
    if forbidden:
        raise ValueError(f"{path} contains forbidden label columns: {', '.join(forbidden)}")
    required = {
        "cell",
        "tf",
        "motif_id",
        "motif_family",
        "TFBS_chr",
        "TFBS_start",
        "TFBS_end",
        "TFBS_strand",
        "motif_score",
        "peak_start",
        "peak_end",
        "chromosome_split",
    }
    missing = required.difference(sites.columns)
    if missing:
        raise ValueError(f"{path} is missing site columns: {', '.join(sorted(missing))}")
    if set(sites["cell"].astype(str)) != {cell}:
        raise ValueError(f"{path} does not contain exactly cell {cell}")
    with np.load(profiles_path, allow_pickle=False) as payload:
        profiles = {key: payload[key] for key in payload.files}
    required_arrays = {"valid", "site_hash"}
    if document["schema"] == "fp-tools-combined-functional-profiles-v1":
        required_arrays.update({"observed", "expected", "combined_residual"})
    else:
        required_arrays.update(
            {
                "plus_observed",
                "minus_observed",
                "plus_expected",
                "minus_expected",
                "combined_residual",
                "shared_strand_residual",
                "antisymmetric_strand_residual",
            }
        )
    missing_arrays = required_arrays.difference(profiles)
    if missing_arrays:
        raise ValueError(f"{path} is missing arrays: {', '.join(sorted(missing_arrays))}")
    if any(len(array) != len(sites) for array in profiles.values()):
        raise ValueError(f"{path} profile row counts do not match its site table")
    if not np.array_equal(np.asarray(profiles["site_hash"]), site_hashes(sites)):
        raise ValueError(f"{path} site hashes do not match its site table")
    if int(document["sites_total"]) != len(sites):
        raise ValueError(f"{path} sites_total does not match its site table")
    if bool(document.get("metadata", {}).get("labels_used", False)):
        raise ValueError(f"{path} reports that labels were used")
    return Artifact(model, cell, path, document, sites, profiles)


def validate_site_source(artifact: Artifact, source: Path) -> None:
    sites = pd.read_csv(source, sep="\t").reset_index(drop=True)
    forbidden = forbidden_label_columns(sites)
    if forbidden:
        raise ValueError(f"{source} contains forbidden label columns: {', '.join(forbidden)}")
    if len(sites) != len(artifact.sites) or not np.array_equal(
        site_hashes(sites), site_hashes(artifact.sites)
    ):
        raise ValueError(f"{artifact.path} does not preserve the frozen site source {source}")


def validate_exact_site_alignment(artifacts: dict[tuple[str, str], Artifact]) -> None:
    """Require candidate artifacts to be exact ordered subsets of cell DWM sites."""

    for cell in sorted({key[1] for key in artifacts}):
        reference = artifacts[("DWM", cell)]
        reference_hashes = np.asarray(reference.profiles["site_hash"], dtype=np.uint64)
        if len(np.unique(reference_hashes)) != len(reference_hashes):
            raise ValueError(f"DWM artifact for {cell} contains duplicate site hashes")
        lookup = {int(value): index for index, value in enumerate(reference_hashes)}
        for (model, item_cell), artifact in artifacts.items():
            if item_cell != cell or model == "DWM":
                continue
            candidate_hashes = np.asarray(artifact.profiles["site_hash"], dtype=np.uint64)
            if len(np.unique(candidate_hashes)) != len(candidate_hashes):
                raise ValueError(f"{artifact.path} contains duplicate site hashes")
            if any(int(value) not in lookup for value in candidate_hashes):
                raise ValueError(
                    f"{artifact.path} contains sites absent from DWM artifact {reference.path}"
                )
            reference_indexes = np.asarray(
                [lookup[int(value)] for value in candidate_hashes], dtype=int
            )
            key_columns = [
                "cell", "tf", "motif_id", "motif_family", "TFBS_chr",
                "TFBS_start", "TFBS_end", "TFBS_strand", "chromosome_split",
            ]
            left = reference.sites.iloc[reference_indexes][key_columns].reset_index(drop=True)
            right = artifact.sites[key_columns].reset_index(drop=True)
            if not left.equals(right):
                raise ValueError(
                    f"site rows differ between {reference.path} and {artifact.path}"
                )


def validate_replicate_artifacts(
    artifacts: dict[tuple[str, str], Artifact],
    replicate_artifacts: dict[tuple[str, str, str], Artifact],
    routes: pd.DataFrame,
) -> None:
    """Validate complete per-cell replicate groups against pooled site sources."""

    if not replicate_artifacts:
        raise ValueError("at least two biological replicate artifact groups are required")
    for cell in sorted(routes["cell"].astype(str).unique()):
        required_models = {"DWM"}.union(
            routes.loc[
                routes["cell"].astype(str).eq(cell), "bias_configuration"
            ].astype(str)
        )
        replicates = sorted(
            {
                replicate
                for replicate, _model, item_cell in replicate_artifacts
                if item_cell == cell
            }
        )
        if len(replicates) < 2:
            raise ValueError(f"{cell} requires at least two biological replicate groups")
        pooled_reference = artifacts[("DWM", cell)]
        for replicate in replicates:
            observed_models = {
                model
                for item_replicate, model, item_cell in replicate_artifacts
                if item_replicate == replicate and item_cell == cell
            }
            if observed_models != required_models:
                raise ValueError(
                    f"{cell}/{replicate} replicate models differ from frozen routes; "
                    f"expected {sorted(required_models)}, observed {sorted(observed_models)}"
                )
            group = {
                (model, cell): replicate_artifacts[(replicate, model, cell)]
                for model in required_models
            }
            validate_exact_site_alignment(group)
            replicate_reference = group[("DWM", cell)]
            pooled_hashes = np.asarray(pooled_reference.profiles["site_hash"], dtype=np.uint64)
            replicate_hashes = np.asarray(replicate_reference.profiles["site_hash"], dtype=np.uint64)
            if not np.array_equal(pooled_hashes, replicate_hashes):
                raise ValueError(
                    f"{cell}/{replicate} DWM sites differ from the pooled DWM artifact"
                )


def reference_to_candidate_indexes(
    reference: Artifact, candidate: Artifact
) -> np.ndarray:
    """Map each reference row to its candidate-artifact row, or -1 if absent."""

    output = np.full(len(reference.sites), -1, dtype=int)
    lookup = {
        int(value): index
        for index, value in enumerate(np.asarray(reference.profiles["site_hash"], dtype=np.uint64))
    }
    for candidate_index, value in enumerate(
        np.asarray(candidate.profiles["site_hash"], dtype=np.uint64)
    ):
        output[lookup[int(value)]] = candidate_index
    return output


def map_indexes_by_hash(
    source: Artifact, target: Artifact, source_indexes: np.ndarray
) -> np.ndarray:
    """Map selected source rows into a target artifact by immutable site hash."""

    target_lookup = {
        int(value): index
        for index, value in enumerate(
            np.asarray(target.profiles["site_hash"], dtype=np.uint64)
        )
    }
    source_hashes = np.asarray(source.profiles["site_hash"], dtype=np.uint64)
    mapped = np.asarray(
        [target_lookup.get(int(source_hashes[index]), -1) for index in source_indexes],
        dtype=int,
    )
    if np.any(mapped < 0):
        raise ValueError("selected pooled sites are absent from a replicate artifact")
    return mapped


def valid_on_reference(reference: Artifact, target: Artifact) -> np.ndarray:
    mapping = reference_to_candidate_indexes(reference, target)
    available = mapping >= 0
    output = np.zeros(len(reference.sites), dtype=bool)
    output[available] = target.valid[mapping[available]]
    return output


def validate_tables(
    study: dict[str, Any],
    routes: pd.DataFrame,
    policy: pd.DataFrame,
    chip: pd.DataFrame,
) -> pd.DataFrame:
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[tasks["split"].astype(str).eq("locked_holdout")].copy()
    task_keys = ["cell", "tf", "motif_id", "motif_family", "role"]
    if tasks.duplicated(["cell", "tf"]).any():
        raise ValueError("locked study tasks contain duplicate cell/TF rows")
    for name, table in (("routes", routes), ("ChIP manifest", chip)):
        if table.duplicated(["cell", "tf"]).any():
            raise ValueError(f"{name} contains duplicate cell/TF rows")
        merged = tasks[task_keys].merge(table[task_keys], on=task_keys, how="outer", indicator=True)
        if not merged["_merge"].eq("both").all():
            raise ValueError(f"{name} does not exactly match the locked study tasks")
    accessions = tasks[["cell", "tf", "chip_accession"]].merge(
        chip[["cell", "tf", "file_accession"]], on=["cell", "tf"], validate="one_to_one"
    )
    if not accessions["chip_accession"].astype(str).eq(accessions["file_accession"].astype(str)).all():
        raise ValueError("ChIP manifest accessions differ from the preregistered study")
    if policy.duplicated("motif_family").any():
        raise ValueError("detector policy contains duplicate motif families")
    known = set(policy["motif_family"].astype(str))
    promoted = routes[~routes["bias_configuration"].astype(str).eq("DWM")]
    missing_policy = set(promoted["motif_family"].astype(str)).difference(known)
    if missing_policy:
        raise ValueError("promoted routes lack policy rows: " + ", ".join(sorted(missing_policy)))
    return tasks


def required_artifact_keys(routes: pd.DataFrame) -> set[tuple[str, str]]:
    keys = {("DWM", str(cell)) for cell in routes["cell"].unique()}
    keys.update(
        (str(row.bias_configuration), str(row.cell))
        for row in routes.itertuples(index=False)
        if str(row.bias_configuration) != "DWM"
    )
    return keys


def hash_record(path: Path) -> dict[str, Any]:
    return {"path": str(path), "sha256": file_sha256(path), "bytes": path.stat().st_size}


def freeze_options(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "positive_summit_distance": int(args.positive_summit_distance),
        "negative_peak_distance": int(args.negative_peak_distance),
        "minimum_sites_per_class": int(args.minimum_sites_per_class),
        "maximum_train_per_tf": int(args.maximum_train_per_tf),
        "matching_features": list(MATCH_FEATURES),
        "maximum_matching_smd": float(args.maximum_matching_smd),
        "bootstrap": int(args.bootstrap),
        "profile_bootstrap": int(args.profile_bootstrap),
        "seed": int(args.seed),
        "training_split": TRAIN_SPLIT,
        "evaluation_split": TEST_SPLIT,
    }


def build_freeze_document(
    args: argparse.Namespace,
    artifacts: dict[tuple[str, str], Artifact],
    replicate_artifacts: dict[tuple[str, str, str], Artifact],
    site_sources: dict[str, Path],
    routes: pd.DataFrame,
) -> dict[str, Any]:
    script = Path(__file__).resolve()
    document: dict[str, Any] = {
        "schema": FREEZE_SCHEMA,
        "policy_frozen_before_holdout": True,
        "locked_holdout_labels_read": False,
        "evaluator": hash_record(script),
        "study": hash_record(args.study),
        "routes": hash_record(args.routes),
        "policy": hash_record(args.policy),
        "chip_manifest_metadata_only": hash_record(args.chip_manifest),
        "genome": hash_record(args.genome),
        "site_sources": {
            cell: hash_record(path) for cell, path in sorted(site_sources.items())
        },
        "artifacts": {
            f"{model}|{cell}": {
                "document": hash_record(artifact.path),
                "profiles": hash_record(_resolve_recorded_path(artifact.path, artifact.document["profiles_npz"])),
                "sites": hash_record(_resolve_recorded_path(artifact.path, artifact.document["sites"])),
                "schema": artifact.schema,
                "sites_total": len(artifact.sites),
                "sites_valid": int(artifact.valid.sum()),
            }
            for (model, cell), artifact in sorted(artifacts.items())
        },
        "replicate_artifacts": {
            f"{replicate}|{model}|{cell}": {
                "document": hash_record(artifact.path),
                "profiles": hash_record(_resolve_recorded_path(artifact.path, artifact.document["profiles_npz"])),
                "sites": hash_record(_resolve_recorded_path(artifact.path, artifact.document["sites"])),
                "schema": artifact.schema,
                "sites_total": len(artifact.sites),
                "sites_valid": int(artifact.valid.sum()),
            }
            for (replicate, model, cell), artifact in sorted(replicate_artifacts.items())
        },
        "options": freeze_options(args),
        "frozen_routes": routes.sort_values(["cell", "tf"], kind="mergesort").to_dict("records"),
    }
    document["freeze_id"] = canonical_hash(document)
    return document


def verify_hash_record(record: dict[str, Any]) -> None:
    path = Path(record["path"])
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.stat().st_size != int(record["bytes"]):
        raise ValueError(f"frozen byte size changed for {path}")
    if file_sha256(path) != str(record["sha256"]):
        raise ValueError(f"frozen SHA-256 changed for {path}")


def verify_freeze(
    freeze_path: Path,
    args: argparse.Namespace,
    artifacts: dict[tuple[str, str], Artifact],
    replicate_artifacts: dict[tuple[str, str, str], Artifact],
    site_sources: dict[str, Path],
    routes: pd.DataFrame,
) -> dict[str, Any]:
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if freeze.get("schema") != FREEZE_SCHEMA:
        raise ValueError(f"{freeze_path} has an unsupported freeze schema")
    freeze_id = freeze.pop("freeze_id", None)
    if freeze_id != canonical_hash(freeze):
        raise ValueError(f"{freeze_path} freeze_id is invalid")
    freeze["freeze_id"] = freeze_id
    expected = build_freeze_document(
        args, artifacts, replicate_artifacts, site_sources, routes
    )
    if expected != freeze:
        raise ValueError("current evaluator inputs/options differ from the locked freeze")
    for key in ("evaluator", "study", "routes", "policy", "chip_manifest_metadata_only", "genome"):
        verify_hash_record(freeze[key])
    for record in freeze["site_sources"].values():
        verify_hash_record(record)
    for artifact in freeze["artifacts"].values():
        for key in ("document", "profiles", "sites"):
            verify_hash_record(artifact[key])
    for artifact in freeze["replicate_artifacts"].values():
        for key in ("document", "profiles", "sites"):
            verify_hash_record(artifact[key])
    return freeze


def validate_chip_files(chip: pd.DataFrame) -> None:
    for row in chip.itertuples(index=False):
        path = Path(row.local_path)
        if not path.is_file():
            raise FileNotFoundError(
                f"locked ChIP file is absent: {path}; download only after creating the freeze"
            )
        if path.stat().st_size != int(row.expected_bytes):
            raise ValueError(f"ENCODE size mismatch for {row.file_accession}")
        if file_md5(path) != str(row.checksum):
            raise ValueError(f"ENCODE MD5 mismatch for {row.file_accession}")


def _subsample(indexes: np.ndarray, limit: int, *seed_parts: object, seed: int) -> np.ndarray:
    indexes = np.asarray(indexes, dtype=int)
    if len(indexes) <= limit:
        return indexes
    rng = np.random.default_rng(stable_seed(*seed_parts, seed=seed))
    return np.sort(rng.choice(indexes, size=limit, replace=False))


def training_indexes(
    artifact: Artifact,
    *,
    tf: str | None = None,
    family: str | None = None,
    maximum_per_tf: int,
    seed: int,
) -> np.ndarray:
    sites = artifact.sites
    base = artifact.valid & sites["chromosome_split"].astype(str).eq(TRAIN_SPLIT).to_numpy()
    if tf is not None:
        indexes = np.flatnonzero(base & sites["tf"].astype(str).eq(tf).to_numpy())
        return _subsample(indexes, maximum_per_tf, artifact.model, artifact.cell, tf, seed=seed)
    if family is not None:
        parts = []
        selected = sites[base & sites["motif_family"].astype(str).eq(family)]
        for item_tf, group in selected.groupby("tf", sort=True):
            parts.append(
                _subsample(
                    group.index.to_numpy(dtype=int),
                    maximum_per_tf,
                    artifact.model,
                    artifact.cell,
                    family,
                    item_tf,
                    seed=seed,
                )
            )
        return np.concatenate(parts) if parts else np.asarray([], dtype=int)
    indexes = np.flatnonzero(base)
    return _subsample(indexes, 400, artifact.model, artifact.cell, "global", seed=seed)


def combined_candidates() -> dict[str, FunctionalCandidate]:
    candidates = combined_candidate_grid("full")
    return {candidate.candidate_id: candidate for candidate in candidates}


def strand_candidates() -> dict[str, StrandCandidate]:
    return {candidate.candidate_id: candidate for candidate in strand_candidate_grid()}


def _probability(log_odds: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(log_odds, -40.0, 40.0)))


def fit_combined_route(
    artifact: Artifact,
    candidate_id: str,
    tf: str,
    family: str,
    evaluation_indexes: np.ndarray,
    positions: np.ndarray,
    maximum_train_per_tf: int,
    seed: int,
    prior_cache: dict[tuple[str, str, str], tuple[np.ndarray | None, float]],
) -> ModelOutput:
    started = perf_counter()
    candidates = combined_candidates()
    if candidate_id not in candidates:
        raise ValueError(f"unknown combined candidate {candidate_id}")
    candidate = candidates[candidate_id]
    observed, expected = artifact.observed_expected()
    tf_train = training_indexes(
        artifact,
        tf=tf,
        maximum_per_tf=maximum_train_per_tf,
        seed=seed,
    )
    if len(tf_train) < 100:
        raise ValueError(f"{artifact.cell}/{tf} has only {len(tf_train)} training sites")
    family_train = training_indexes(
        artifact,
        family=family,
        maximum_per_tf=maximum_train_per_tf,
        seed=seed,
    )
    dispersion = estimate_nb_dispersion(observed[family_train], expected[family_train])
    if candidate.family in {"spline", "gp"}:
        cache_key = (str(artifact.path), candidate_id, family)
        if cache_key not in prior_cache:
            global_prior = None
            if candidate.shrinkage > 0:
                global_indexes = training_indexes(
                    artifact,
                    maximum_per_tf=maximum_train_per_tf,
                    seed=seed,
                )
                global_model = _count_model(candidate, positions, dispersion, shrinkage=0.0)
                global_result = global_model.fit(
                    observed[global_indexes],
                    expected[global_indexes],
                    motif_score=artifact.sites.iloc[global_indexes]["motif_score"].to_numpy(dtype=float),
                    accessibility=observed[global_indexes].sum(axis=1),
                )
                family_model = _count_model(candidate, positions, dispersion)
                family_result = family_model.fit(
                    observed[family_train],
                    expected[family_train],
                    motif_score=artifact.sites.iloc[family_train]["motif_score"].to_numpy(dtype=float),
                    accessibility=observed[family_train].sum(axis=1),
                    prior_profile=global_result.footprint_profile,
                )
                global_prior = family_result.footprint_profile
            prior_cache[cache_key] = (global_prior, float(dispersion))
        prior_profile, dispersion = prior_cache[cache_key]
        model = _count_model(candidate, positions, dispersion)
        result = model.fit(
            observed[tf_train],
            expected[tf_train],
            motif_score=artifact.sites.iloc[tf_train]["motif_score"].to_numpy(dtype=float),
            accessibility=observed[tf_train].sum(axis=1),
            prior_profile=prior_profile,
        )
        likelihood, prior = model.predict_log_odds_components(
            observed[evaluation_indexes],
            expected[evaluation_indexes],
            motif_score=artifact.sites.iloc[evaluation_indexes]["motif_score"].to_numpy(dtype=float),
            accessibility=observed[evaluation_indexes].sum(axis=1),
        )
        probabilities = _probability(likelihood + prior)
        adjusted_expected = model._background(
            observed[evaluation_indexes], expected[evaluation_indexes]
        )
        profiles = normalize_functional_profiles(
            deviance_profiles(observed[evaluation_indexes], adjusted_expected, dispersion),
            positions,
        )
        footprint = result.footprint_profile
        converged = bool(result.converged)
        iterations = int(result.iterations)
    else:
        residual = artifact.profiles["combined_residual"]
        weights = np.sqrt(np.maximum(observed[tf_train].sum(axis=1), 1.0))
        model_seed = stable_seed(artifact.cell, tf, candidate_id, seed=seed)
        if candidate.family == "fda":
            model = FdaMixtureModel(
                variance_threshold=candidate.variance_threshold,
                max_components=candidate.max_components,
                seed=model_seed,
            ).fit(residual[tf_train], positions=positions, sample_weight=weights)
            probabilities = model.predict_proba(residual[evaluation_indexes])
            components = model.component_profiles()
            assert model.binding_component_ is not None
            footprint = components[model.binding_component_] - components[1 - model.binding_component_]
            converged = bool(model.mixture.converged_) if model.mixture is not None else False
            iterations = int(model.mixture.n_iter_) if model.mixture is not None else 0
        elif candidate.family == "hybrid":
            model = HybridFdaGpModel(
                positions,
                variance_threshold=candidate.variance_threshold,
                max_components=candidate.max_components,
                seed=model_seed,
            ).fit(residual[tf_train], sample_weight=weights)
            probabilities = model.predict_proba(residual[evaluation_indexes])
            assert model.bound_mean_ is not None and model.unbound_mean_ is not None
            footprint = model.bound_mean_ - model.unbound_mean_
            converged = bool(model.fda.mixture.converged_) if model.fda.mixture is not None else False
            iterations = int(model.fda.mixture.n_iter_) if model.fda.mixture is not None else 0
        else:
            raise ValueError(f"locked evaluator does not support {candidate.family}")
        profiles = normalize_functional_profiles(residual[evaluation_indexes], positions)
    return ModelOutput(
        np.asarray(probabilities),
        np.asarray(profiles),
        np.asarray(footprint),
        converged,
        iterations,
        perf_counter() - started,
        float(dispersion),
    )


def fit_strand_route(
    artifact: Artifact,
    candidate_id: str,
    tf: str,
    family: str,
    evaluation_indexes: np.ndarray,
    positions: np.ndarray,
    maximum_train_per_tf: int,
    seed: int,
    prior_cache: dict[tuple[str, str, str], tuple[np.ndarray, float]],
) -> ModelOutput:
    started = perf_counter()
    candidates = strand_candidates()
    if candidate_id not in candidates:
        raise ValueError(f"unknown strand candidate {candidate_id}")
    candidate = candidates[candidate_id]
    observed, expected = artifact.observed_expected()
    tf_train = training_indexes(
        artifact, tf=tf, maximum_per_tf=maximum_train_per_tf, seed=seed
    )
    family_train = training_indexes(
        artifact, family=family, maximum_per_tf=maximum_train_per_tf, seed=seed
    )
    if len(tf_train) < 100 or len(family_train) < 100:
        raise ValueError(f"{artifact.cell}/{tf} has insufficient unlabeled training sites")
    dispersion = estimate_nb_dispersion(observed[family_train], expected[family_train])
    if candidate.family == "count":
        cache_key = (str(artifact.path), candidate_id, family)
        if cache_key not in prior_cache:
            family_model = BiasAwareFunctionalMixture(
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
            family_result = family_model.fit(observed[family_train], expected[family_train])
            prior_cache[cache_key] = (family_result.footprint_profile, float(dispersion))
        prior_profile, dispersion = prior_cache[cache_key]
        model = BiasAwareFunctionalMixture(
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
        result = model.fit(
            observed[tf_train], expected[tf_train], prior_profile=prior_profile
        )
        likelihood, _prior = model.predict_log_odds_components(
            observed[evaluation_indexes], expected[evaluation_indexes]
        )
        probabilities = _probability(likelihood)
        profiles = normalize_functional_profiles(
            artifact.profiles["combined_residual"][evaluation_indexes], positions
        )
        footprint = result.footprint_profile
        converged = bool(result.converged)
        iterations = int(result.iterations)
    else:
        indexes = tf_train if candidate.training_pool == "tf" else family_train
        train_profiles = artifact.profiles[candidate.channel][indexes]
        evaluation_profiles = artifact.profiles[candidate.channel][evaluation_indexes]
        weights = np.sqrt(np.maximum(observed[indexes].sum(axis=1), 1.0))
        model_seed = stable_seed(artifact.model, artifact.cell, tf, candidate_id, seed=seed)
        if candidate.family == "fda":
            model = FdaMixtureModel(max_components=20, seed=model_seed).fit(
                train_profiles, positions=positions, sample_weight=weights
            )
            probabilities = model.predict_proba(evaluation_profiles)
            components = model.component_profiles()
            assert model.binding_component_ is not None
            footprint = components[model.binding_component_] - components[1 - model.binding_component_]
            converged = bool(model.mixture.converged_) if model.mixture is not None else False
            iterations = int(model.mixture.n_iter_) if model.mixture is not None else 0
        elif candidate.family == "hybrid":
            model = HybridFdaGpModel(positions, max_components=20, seed=model_seed).fit(
                train_profiles, sample_weight=weights
            )
            probabilities = model.predict_proba(evaluation_profiles)
            assert model.bound_mean_ is not None and model.unbound_mean_ is not None
            footprint = model.bound_mean_ - model.unbound_mean_
            converged = bool(model.fda.mixture.converged_) if model.fda.mixture is not None else False
            iterations = int(model.fda.mixture.n_iter_) if model.fda.mixture is not None else 0
        else:
            raise ValueError(f"locked evaluator does not support {candidate.family}")
        profiles = normalize_functional_profiles(evaluation_profiles, positions)
    return ModelOutput(
        np.asarray(probabilities),
        np.asarray(profiles),
        np.asarray(footprint),
        converged,
        iterations,
        perf_counter() - started,
        float(dispersion),
    )


def gc_fractions(sites: pd.DataFrame, genome: Path, flank: int = 100) -> np.ndarray:
    output = np.full(len(sites), np.nan, dtype=float)
    with open_fasta(genome) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        for index, row in enumerate(sites.itertuples(index=False)):
            chrom = str(row.TFBS_chr)
            center = (int(row.TFBS_start) + int(row.TFBS_end)) // 2
            start, end = center - flank, center + flank + 1
            if chrom not in lengths or start < 0 or end > int(lengths[chrom]):
                continue
            sequence = fasta.fetch(chrom, start, end).upper()
            bases = sum(sequence.count(base) for base in "ACGT")
            output[index] = (sequence.count("G") + sequence.count("C")) / bases if bases else np.nan
    return output


def standardized_difference(frame: pd.DataFrame, feature: str) -> float:
    positive = frame.loc[frame["label"].eq(1), feature].to_numpy(dtype=float)
    negative = frame.loc[frame["label"].eq(0), feature].to_numpy(dtype=float)
    pooled = np.sqrt((np.var(positive) + np.var(negative)) / 2.0)
    return float((np.mean(positive) - np.mean(negative)) / pooled) if pooled > 0 else 0.0


def build_matched_test_sites(
    sites: pd.DataFrame,
    valid: np.ndarray,
    observed: np.ndarray,
    chip_path: Path,
    genome: Path,
    *,
    positive_summit_distance: int,
    negative_peak_distance: int,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    test_indexes = np.flatnonzero(
        valid & sites["chromosome_split"].astype(str).eq(TEST_SPLIT).to_numpy()
    )
    test = sites.iloc[test_indexes].copy().reset_index().rename(columns={"index": "artifact_index"})
    label_input = pd.DataFrame(
        {
            "chrom": test["TFBS_chr"],
            "start": test["TFBS_start"],
            "end": test["TFBS_end"],
            "strand": test["TFBS_strand"],
            "site_id": test["artifact_index"].astype(str),
            "motif_score": test["motif_score"],
        }
    )
    labelled = label_sites(
        label_input,
        read_peaks(chip_path),
        positive_summit_distance=positive_summit_distance,
        negative_peak_distance=negative_peak_distance,
    )
    for column in ("label", "label_reason", "nearest_peak_distance", "nearest_summit_distance"):
        test[column] = labelled[column].to_numpy()
    test["accessibility"] = observed[test["artifact_index"].to_numpy(dtype=int)].sum(axis=1)
    test["log_accessibility"] = np.log1p(np.maximum(test["accessibility"], 0.0))
    test["gc_fraction"] = gc_fractions(test, genome)
    center = (test["TFBS_start"].to_numpy(dtype=float) + test["TFBS_end"].to_numpy(dtype=float)) / 2.0
    peak_center = (test["peak_start"].to_numpy(dtype=float) + test["peak_end"].to_numpy(dtype=float)) / 2.0
    half_width = np.maximum(
        (test["peak_end"].to_numpy(dtype=float) - test["peak_start"].to_numpy(dtype=float)) / 2.0,
        1.0,
    )
    test["peak_position_signed"] = (center - peak_center) / half_width
    test["peak_position_abs"] = np.abs(test["peak_position_signed"])
    natural = test[test["label"].isin([0, 1])].dropna(subset=list(MATCH_FEATURES)).copy()
    # ``propensity_match`` uses the generic BED-style keys for deterministic
    # tie-breaking.  Keep the original TFBS columns as the authoritative site
    # coordinates and add aliases only for the matcher.
    natural["chrom"] = natural["TFBS_chr"]
    natural["start"] = natural["TFBS_start"]
    natural["end"] = natural["TFBS_end"]
    positives = natural[natural["label"].eq(1)]
    negatives = natural[natural["label"].eq(0)]
    if positives.empty or negatives.empty:
        raise ValueError("test sites lack a natural positive or high-confidence negative class")
    if len(positives) > len(negatives):
        positives = positives.sample(
            n=len(negatives), random_state=stable_seed("positive-cap", seed=seed)
        )
        natural = pd.concat([positives, negatives], ignore_index=True)
    before = {feature: standardized_difference(natural, feature) for feature in MATCH_FEATURES}
    matched = propensity_match(natural, MATCH_FEATURES, negative_ratio=1, seed=seed)
    after = {feature: standardized_difference(matched, feature) for feature in MATCH_FEATURES}
    diagnostics = {
        "test_sites_valid": int(len(test)),
        "positive_summit_supported": int(test["label"].eq(1).sum()),
        "negative_far_from_peak": int(test["label"].eq(0).sum()),
        "indeterminate": int(test["label"].eq(-1).sum()),
        "matched_per_class": int(matched["label"].eq(1).sum()),
        "before_smd": before,
        "after_smd": after,
        "maximum_absolute_after_smd": float(max(abs(value) for value in after.values())),
    }
    return matched.sort_values(["TFBS_chr", "TFBS_start", "label"], kind="mergesort").reset_index(drop=True), diagnostics


def expected_calibration_error(labels: np.ndarray, probabilities: np.ndarray, bins: int = 10) -> float:
    labels = np.asarray(labels, dtype=int)
    probabilities = np.clip(np.asarray(probabilities, dtype=float), 0.0, 1.0)
    edges = np.linspace(0.0, 1.0, bins + 1)
    total = len(labels)
    if total == 0:
        return np.nan
    value = 0.0
    for index in range(bins):
        if index + 1 == bins:
            member = (probabilities >= edges[index]) & (probabilities <= edges[index + 1])
        else:
            member = (probabilities >= edges[index]) & (probabilities < edges[index + 1])
        if member.any():
            value += member.mean() * abs(labels[member].mean() - probabilities[member].mean())
    return float(value)


def bootstrap_delta(
    frame: pd.DataFrame,
    candidate: np.ndarray,
    reference: np.ndarray,
    *,
    iterations: int,
    seed: int,
) -> dict[str, float | int]:
    chromosomes = np.asarray(sorted(frame["TFBS_chr"].astype(str).unique()))
    rng = np.random.default_rng(seed)
    auroc, relative_auprc = [], []
    labels = frame["label"].to_numpy(dtype=int)
    chrom = frame["TFBS_chr"].astype(str).to_numpy()
    for _ in range(iterations):
        sampled = rng.choice(chromosomes, size=len(chromosomes), replace=True)
        indexes = np.concatenate([np.flatnonzero(chrom == item) for item in sampled])
        if np.unique(labels[indexes]).size != 2:
            continue
        candidate_metrics = binary_metrics(labels[indexes], candidate[indexes])
        reference_metrics = binary_metrics(labels[indexes], reference[indexes])
        auroc.append(float(candidate_metrics["auroc"]) - float(reference_metrics["auroc"]))
        relative_auprc.append(
            (float(candidate_metrics["auprc"]) - float(reference_metrics["auprc"]))
            / max(float(reference_metrics["auprc"]), 1e-8)
        )
    def summary(values: list[float], prefix: str) -> dict[str, float]:
        array = np.asarray(values, dtype=float)
        return {
            f"{prefix}_bootstrap_mean": float(np.mean(array)) if len(array) else np.nan,
            f"{prefix}_bootstrap_lower_95": float(np.quantile(array, 0.025)) if len(array) else np.nan,
            f"{prefix}_bootstrap_upper_95": float(np.quantile(array, 0.975)) if len(array) else np.nan,
            f"{prefix}_bootstrap_probability_positive": float(np.mean(array > 0)) if len(array) else np.nan,
        }
    return {
        "bootstrap_requested": int(iterations),
        "bootstrap_successful": int(len(auroc)),
        **summary(auroc, "auroc_gain"),
        **summary(relative_auprc, "relative_auprc_gain"),
    }


def mean_band(
    profiles: np.ndarray,
    *,
    iterations: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    profiles = np.asarray(profiles, dtype=float)
    mean = np.mean(profiles, axis=0)
    if len(profiles) < 2 or iterations < 2:
        return mean, np.full_like(mean, np.nan), np.full_like(mean, np.nan)
    rng = np.random.default_rng(seed)
    sampled = np.empty((iterations, profiles.shape[1]), dtype=np.float32)
    for index in range(iterations):
        rows = rng.integers(0, len(profiles), size=len(profiles))
        sampled[index] = np.mean(profiles[rows], axis=0)
    return mean, np.quantile(sampled, 0.025, axis=0), np.quantile(sampled, 0.975, axis=0)


def aggregate_rows(
    *,
    cell: str,
    tf: str,
    method: str,
    profiles: np.ndarray,
    labels: np.ndarray,
    positions: np.ndarray,
    iterations: int,
    seed: int,
) -> list[dict[str, Any]]:
    rows = []
    for label, group in ((0, "matched_negative"), (1, "chip_positive")):
        mean, lower, upper = mean_band(
            profiles[labels == label],
            iterations=iterations,
            seed=stable_seed(cell, tf, method, group, seed=seed),
        )
        rows.extend(
            {
                "cell": cell,
                "tf": tf,
                "method": method,
                "group": group,
                "position": int(position),
                "mean": float(value),
                "lower_95": float(low),
                "upper_95": float(high),
            }
            for position, value, low, high in zip(positions, mean, lower, upper)
        )
    return rows


def reference_candidate(policy: pd.DataFrame, route: Any) -> str:
    family = str(route.motif_family)
    match = policy[policy["motif_family"].astype(str).eq(family)]
    if len(match):
        return str(match.iloc[0]["reference_candidate_id"])
    if str(route.bias_configuration) == "DWM":
        return str(route.candidate_id)
    raise ValueError(f"promoted family {family} lacks a reference route")


def _profile_metrics(profiles: np.ndarray, labels: np.ndarray, positions: np.ndarray) -> dict[str, float]:
    positive = profiles[labels == 1].mean(axis=0)
    negative = profiles[labels == 0].mean(axis=0)
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


def render_pdf(metrics: pd.DataFrame, profiles: pd.DataFrame, output: Path) -> None:
    colors = {"chip_positive": "#C23B33", "matched_negative": "#2A6FBB"}
    with PdfPages(output) as pdf:
        for row in metrics.itertuples(index=False):
            if str(row.status) != "ok":
                continue
            figure, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
            subset = profiles[
                profiles["cell"].astype(str).eq(str(row.cell))
                & profiles["tf"].astype(str).eq(str(row.tf))
            ]
            for axis, method, title in (
                (axes[0], "candidate", f"Frozen route\n{row.bias_configuration} / {row.candidate_id}"),
                (axes[1], "reference", f"DWM reference\n{row.reference_candidate_id}"),
            ):
                for group in ("chip_positive", "matched_negative"):
                    curve = subset[
                        subset["method"].eq(method) & subset["group"].eq(group)
                    ].sort_values("position")
                    axis.plot(curve["position"], curve["mean"], color=colors[group], label=group)
                    axis.fill_between(
                        curve["position"].to_numpy(dtype=float),
                        curve["lower_95"].to_numpy(dtype=float),
                        curve["upper_95"].to_numpy(dtype=float),
                        color=colors[group],
                        alpha=0.18,
                    )
                axis.axvline(0, color="#777777", lw=0.8, ls="--")
                axis.set_title(title, fontsize=9)
                axis.set_xlabel("Position relative to motif (bp)")
                axis.set_ylabel("Normalized signed-deviance residual")
                axis.legend(frameon=False, fontsize=8)
            axes[2].axis("off")
            text = (
                f"{row.cell} / {row.tf} ({row.motif_family})\n\n"
                f"Matched sites/class: {row.positive_sites:,}\n"
                f"Candidate AUROC: {row.candidate_auroc:.3f}\n"
                f"Reference AUROC: {row.reference_auroc:.3f}\n"
                f"AUROC gain: {row.auroc_gain:+.3f}\n\n"
                f"Candidate AUPRC: {row.candidate_auprc:.3f}\n"
                f"Reference AUPRC: {row.reference_auprc:.3f}\n"
                f"Relative AUPRC gain: {row.relative_auprc_gain:+.1%}\n\n"
                f"Functional-separation gain: {row.functional_separation_gain:+.1%}\n"
                f"Max matched SMD: {row.maximum_absolute_after_smd:.3f}\n"
                f"Bootstrap P(AUROC gain > 0): "
                f"{row.auroc_gain_bootstrap_probability_positive:.3f}"
            )
            axes[2].text(0.02, 0.98, text, ha="left", va="top", family="monospace", fontsize=9)
            figure.suptitle("Locked, label-free-trained footprint transfer", fontsize=12)
            pdf.savefig(figure)
            plt.close(figure)


def promotion_summary(metrics: pd.DataFrame, study: dict[str, Any]) -> dict[str, Any]:
    ok = metrics[metrics["status"].eq("ok")].copy()
    promoted = ok[~ok["bias_configuration"].astype(str).eq("DWM")].copy()
    positive = ok[ok["role"].astype(str).eq("positive_control")]
    gates = study["promotion_gates"]
    families_improved = int(
        promoted.groupby("motif_family")["auroc_gain"].mean().gt(0).sum()
    ) if len(promoted) else 0
    cells_improved = int(
        promoted.groupby("cell")["auroc_gain"].mean().gt(0).sum()
    ) if len(promoted) else 0
    mean_auroc = float(promoted["auroc_gain"].mean()) if len(promoted) else np.nan
    relative_auprc = float(promoted["relative_auprc_gain"].mean()) if len(promoted) else np.nan
    maximum_control_loss = float(np.maximum(-positive["auroc_gain"], 0).max()) if len(positive) else np.nan
    checks = {
        "minimum_mean_auroc_gain": bool(
            len(promoted) and mean_auroc >= float(gates["minimum_mean_auroc_gain"])
        ),
        "minimum_relative_auprc_gain": bool(
            len(promoted) and relative_auprc >= float(gates["minimum_relative_auprc_gain"])
        ),
        "minimum_difficult_tf_families_improved": bool(
            families_improved >= int(gates["minimum_difficult_tf_families_improved"])
        ),
        "minimum_holdout_cells_improved": bool(
            cells_improved >= int(gates["minimum_holdout_cells_improved"])
        ),
        "maximum_positive_control_auroc_loss": bool(
            np.isfinite(maximum_control_loss)
            and maximum_control_loss <= float(gates["maximum_positive_control_auroc_loss"])
        ),
        "matching_balance": bool(
            len(ok) and ok["matching_pass"].astype(bool).all()
        ),
        "biological_replicate_direction_stable": bool(
            len(promoted)
            and "replicate_direction_stable" in promoted
            and promoted["replicate_direction_stable"].fillna(False).astype(bool).all()
        ),
    }
    return {
        "promoted_task_rows": int(len(promoted)),
        "mean_auroc_gain": mean_auroc,
        "mean_relative_auprc_gain": relative_auprc,
        "families_improved": families_improved,
        "holdout_cells_improved": cells_improved,
        "maximum_positive_control_auroc_loss": maximum_control_loss,
        "checks": checks,
        "locked_holdout_performance_gate_passed": bool(all(checks.values())),
    }


def promotion_metric_table(metrics: pd.DataFrame) -> pd.DataFrame:
    """Convert successful paired rows into the fail-closed auditor schema."""

    rows: list[dict[str, Any]] = []
    for item in metrics[metrics["status"].eq("ok")].itertuples(index=False):
        common = {
            "cell": str(item.cell),
            "tf": str(item.tf),
            "motif_id": str(item.motif_id),
            "motif_family": str(item.motif_family),
            "role": str(item.role),
            "split": "locked_holdout",
            "bias_configuration": str(item.bias_configuration),
            "route_candidate_id": str(item.candidate_id),
            "reference_candidate_id": str(item.reference_candidate_id),
        }
        rows.extend(
            (
                {
                    **common,
                    "method": PROMOTION_CANDIDATE,
                    "auroc": float(item.candidate_auroc),
                    "auprc": float(item.candidate_auprc),
                    "brier": float(item.candidate_brier),
                    "ece": float(item.candidate_ece),
                },
                {
                    **common,
                    "method": PROMOTION_REFERENCE,
                    "auroc": float(item.reference_auroc),
                    "auprc": float(item.reference_auprc),
                    "brier": float(item.reference_brier),
                    "ece": float(item.reference_ece),
                },
            )
        )
    return pd.DataFrame(rows)


def promotion_descriptor_table(metrics: pd.DataFrame) -> pd.DataFrame:
    """Emit group-level depletion descriptors for the final gate."""

    rows: list[dict[str, Any]] = []
    for item in metrics[metrics["status"].eq("ok")].itertuples(index=False):
        common = {
            "cell": str(item.cell),
            "tf": str(item.tf),
            "motif_family": str(item.motif_family),
        }
        for correction, prefix in (
            (PROMOTION_CANDIDATE, "candidate"),
            (PROMOTION_REFERENCE, "reference"),
        ):
            rows.extend(
                (
                    {
                        **common,
                        "correction": correction,
                        "group": "chip_positive",
                        "depletion": float(getattr(item, f"{prefix}_positive_depletion")),
                    },
                    {
                        **common,
                        "correction": correction,
                        "group": "matched_negative",
                        "depletion": float(getattr(item, f"{prefix}_negative_depletion")),
                    },
                )
            )
    return pd.DataFrame(rows)


def promotion_stability_table(metrics: pd.DataFrame) -> pd.DataFrame:
    """Expose biological-replicate direction as explicit gate evidence."""

    ok = metrics[metrics["status"].eq("ok")].copy()
    if ok.empty:
        return pd.DataFrame()
    return pd.DataFrame(
        {
            "cell": ok["cell"].astype(str),
            "tf": ok["tf"].astype(str),
            "candidate_id": PROMOTION_CANDIDATE,
            "direction_consistent": ok["replicate_direction_stable"]
            .fillna(False)
            .astype(bool),
            "biological_replicates": ok["biological_replicates"].fillna(0).astype(int),
        }
    )


def run_evaluation(
    args: argparse.Namespace,
    study: dict[str, Any],
    routes: pd.DataFrame,
    policy: pd.DataFrame,
    chip: pd.DataFrame,
    artifacts: dict[tuple[str, str], Artifact],
    replicate_artifacts: dict[tuple[str, str, str], Artifact],
    freeze: dict[str, Any],
) -> int:
    positions = np.arange(
        -int(study["profile_flank_bp"]), int(study["profile_flank_bp"]) + 1, dtype=float
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    metric_rows: list[dict[str, Any]] = []
    score_frames: list[pd.DataFrame] = []
    profile_rows: list[dict[str, Any]] = []
    matching_rows: list[dict[str, Any]] = []
    replicate_metric_rows: list[dict[str, Any]] = []
    combined_prior_cache: dict[tuple[str, str, str], tuple[np.ndarray | None, float]] = {}
    strand_prior_cache: dict[tuple[str, str, str], tuple[np.ndarray, float]] = {}
    for route in routes.sort_values(["cell", "tf"], kind="mergesort").itertuples(index=False):
        cell, tf, family = str(route.cell), str(route.tf), str(route.motif_family)
        candidate_artifact = artifacts[(str(route.bias_configuration), cell)]
        reference_artifact = artifacts[("DWM", cell)]
        observed, _expected = reference_artifact.observed_expected()
        candidate_map = reference_to_candidate_indexes(reference_artifact, candidate_artifact)
        candidate_available = candidate_map >= 0
        candidate_valid_on_reference = np.zeros(len(reference_artifact.sites), dtype=bool)
        candidate_valid_on_reference[candidate_available] = candidate_artifact.valid[
            candidate_map[candidate_available]
        ]
        joint_valid = reference_artifact.valid & candidate_valid_on_reference
        for replicate in sorted(
            {
                item_replicate
                for item_replicate, _model, item_cell in replicate_artifacts
                if item_cell == cell
            }
        ):
            replicate_reference = replicate_artifacts[(replicate, "DWM", cell)]
            replicate_candidate = replicate_artifacts[
                (replicate, str(route.bias_configuration), cell)
            ]
            joint_valid &= valid_on_reference(reference_artifact, replicate_reference)
            joint_valid &= valid_on_reference(reference_artifact, replicate_candidate)
        chip_row = chip[
            chip["cell"].astype(str).eq(cell) & chip["tf"].astype(str).eq(tf)
        ].iloc[0]
        base: dict[str, Any] = {
            "cell": cell,
            "tf": tf,
            "motif_id": str(route.motif_id),
            "motif_family": family,
            "role": str(route.role),
            "route_source": str(route.route_source),
            "bias_configuration": str(route.bias_configuration),
            "candidate_id": str(route.candidate_id),
            "reference_candidate_id": reference_candidate(policy, route),
            "training_labels_used": False,
            "evaluation_split": TEST_SPLIT,
            "chip_accession": str(chip_row.file_accession),
        }
        try:
            task_mask = reference_artifact.sites["tf"].astype(str).eq(tf).to_numpy()
            matched, matching = build_matched_test_sites(
                reference_artifact.sites,
                joint_valid & task_mask,
                observed,
                Path(chip_row.local_path),
                args.genome,
                positive_summit_distance=args.positive_summit_distance,
                negative_peak_distance=args.negative_peak_distance,
                seed=stable_seed(cell, tf, "matching", seed=args.seed),
            )
            matching_rows.append(
                {
                    "cell": cell,
                    "tf": tf,
                    **{key: value for key, value in matching.items() if not isinstance(value, dict)},
                    **{f"before_smd_{key}": value for key, value in matching["before_smd"].items()},
                    **{f"after_smd_{key}": value for key, value in matching["after_smd"].items()},
                }
            )
            if matching["matched_per_class"] < args.minimum_sites_per_class:
                metric_rows.append({**base, **matching, "status": "insufficient_matched_sites"})
                continue
            indexes = matched["artifact_index"].to_numpy(dtype=int)
            candidate_indexes = candidate_map[indexes]
            if np.any(candidate_indexes < 0):
                raise ValueError("matched reference sites are absent from the candidate artifact")
            labels = matched["label"].to_numpy(dtype=int)
            if str(route.bias_configuration) == "DWM":
                candidate = fit_combined_route(
                    reference_artifact,
                    str(route.candidate_id),
                    tf,
                    family,
                    indexes,
                    positions,
                    args.maximum_train_per_tf,
                    args.seed,
                    combined_prior_cache,
                )
            else:
                candidate = fit_strand_route(
                    candidate_artifact,
                    str(route.candidate_id),
                    tf,
                    family,
                    candidate_indexes,
                    positions,
                    args.maximum_train_per_tf,
                    args.seed,
                    strand_prior_cache,
                )
            reference = fit_combined_route(
                reference_artifact,
                base["reference_candidate_id"],
                tf,
                family,
                indexes,
                positions,
                args.maximum_train_per_tf,
                args.seed,
                combined_prior_cache,
            )
            candidate_metrics = binary_metrics(labels, candidate.probabilities)
            reference_metrics = binary_metrics(labels, reference.probabilities)
            candidate_profiles = _profile_metrics(candidate.evaluation_profiles, labels, positions)
            reference_profiles = _profile_metrics(reference.evaluation_profiles, labels, positions)
            relative_auprc = (
                (float(candidate_metrics["auprc"]) - float(reference_metrics["auprc"]))
                / max(float(reference_metrics["auprc"]), 1e-8)
            )
            separation_gain = (
                candidate_profiles["functional_separation"]
                - reference_profiles["functional_separation"]
            ) / max(abs(reference_profiles["functional_separation"]), 1e-8)
            bootstrap = bootstrap_delta(
                matched,
                candidate.probabilities,
                reference.probabilities,
                iterations=args.bootstrap,
                seed=stable_seed(cell, tf, "bootstrap", seed=args.seed),
            )
            matching_pass = bool(
                matching["maximum_absolute_after_smd"] <= args.maximum_matching_smd
            )
            row = {
                **base,
                "status": "ok",
                "matching_pass": matching_pass,
                **{key: value for key, value in matching.items() if not isinstance(value, dict)},
                "positive_sites": int(candidate_metrics["positive_sites"]),
                "negative_sites": int(candidate_metrics["negative_sites"]),
                "candidate_auroc": float(candidate_metrics["auroc"]),
                "candidate_auprc": float(candidate_metrics["auprc"]),
                "candidate_brier": float(candidate_metrics["brier"]),
                "candidate_ece": expected_calibration_error(labels, candidate.probabilities),
                "reference_auroc": float(reference_metrics["auroc"]),
                "reference_auprc": float(reference_metrics["auprc"]),
                "reference_brier": float(reference_metrics["brier"]),
                "reference_ece": expected_calibration_error(labels, reference.probabilities),
                "auroc_gain": float(candidate_metrics["auroc"]) - float(reference_metrics["auroc"]),
                "relative_auprc_gain": float(relative_auprc),
                "candidate_functional_separation": candidate_profiles["functional_separation"],
                "reference_functional_separation": reference_profiles["functional_separation"],
                "functional_separation_gain": float(separation_gain),
                "candidate_positive_depletion": candidate_profiles["positive_depletion"],
                "candidate_negative_depletion": candidate_profiles["negative_depletion"],
                "reference_positive_depletion": reference_profiles["positive_depletion"],
                "reference_negative_depletion": reference_profiles["negative_depletion"],
                "candidate_depletion_difference": candidate_profiles["depletion_difference"],
                "reference_depletion_difference": reference_profiles["depletion_difference"],
                "candidate_converged": candidate.converged,
                "reference_converged": reference.converged,
                "candidate_iterations": candidate.iterations,
                "reference_iterations": reference.iterations,
                "candidate_fit_seconds": candidate.fit_seconds,
                "reference_fit_seconds": reference.fit_seconds,
                "candidate_dispersion": candidate.dispersion,
                "reference_dispersion": reference.dispersion,
                **bootstrap,
            }
            metric_rows.append(row)
            cell_replicates = sorted(
                {
                    replicate
                    for replicate, _model, item_cell in replicate_artifacts
                    if item_cell == cell
                }
            )
            for replicate in cell_replicates:
                replicate_reference_artifact = replicate_artifacts[
                    (replicate, "DWM", cell)
                ]
                replicate_candidate_artifact = replicate_artifacts[
                    (replicate, str(route.bias_configuration), cell)
                ]
                replicate_reference_indexes = map_indexes_by_hash(
                    reference_artifact, replicate_reference_artifact, indexes
                )
                replicate_candidate_indexes = map_indexes_by_hash(
                    reference_artifact, replicate_candidate_artifact, indexes
                )
                replicate_base = {
                    "cell": cell,
                    "tf": tf,
                    "motif_family": family,
                    "role": str(route.role),
                    "replicate": replicate,
                    "bias_configuration": str(route.bias_configuration),
                    "candidate_id": str(route.candidate_id),
                    "reference_candidate_id": base["reference_candidate_id"],
                    "training_labels_used": False,
                }
                try:
                    if str(route.bias_configuration) == "DWM":
                        replicate_candidate = fit_combined_route(
                            replicate_reference_artifact,
                            str(route.candidate_id),
                            tf,
                            family,
                            replicate_reference_indexes,
                            positions,
                            args.maximum_train_per_tf,
                            args.seed,
                            combined_prior_cache,
                        )
                    else:
                        replicate_candidate = fit_strand_route(
                            replicate_candidate_artifact,
                            str(route.candidate_id),
                            tf,
                            family,
                            replicate_candidate_indexes,
                            positions,
                            args.maximum_train_per_tf,
                            args.seed,
                            strand_prior_cache,
                        )
                    replicate_reference = fit_combined_route(
                        replicate_reference_artifact,
                        base["reference_candidate_id"],
                        tf,
                        family,
                        replicate_reference_indexes,
                        positions,
                        args.maximum_train_per_tf,
                        args.seed,
                        combined_prior_cache,
                    )
                    candidate_rep_metrics = binary_metrics(
                        labels, replicate_candidate.probabilities
                    )
                    reference_rep_metrics = binary_metrics(
                        labels, replicate_reference.probabilities
                    )
                    replicate_metric_rows.append(
                        {
                            **replicate_base,
                            "status": "ok",
                            "positive_sites": int(candidate_rep_metrics["positive_sites"]),
                            "negative_sites": int(candidate_rep_metrics["negative_sites"]),
                            "candidate_auroc": float(candidate_rep_metrics["auroc"]),
                            "candidate_auprc": float(candidate_rep_metrics["auprc"]),
                            "reference_auroc": float(reference_rep_metrics["auroc"]),
                            "reference_auprc": float(reference_rep_metrics["auprc"]),
                            "auroc_gain": float(candidate_rep_metrics["auroc"])
                            - float(reference_rep_metrics["auroc"]),
                            "relative_auprc_gain": (
                                float(candidate_rep_metrics["auprc"])
                                - float(reference_rep_metrics["auprc"])
                            ) / max(float(reference_rep_metrics["auprc"]), 1e-8),
                            "candidate_converged": replicate_candidate.converged,
                            "reference_converged": replicate_reference.converged,
                        }
                    )
                except Exception as error:
                    replicate_metric_rows.append(
                        {
                            **replicate_base,
                            "status": "error",
                            "error": f"{type(error).__name__}: {error}",
                        }
                    )
            score = matched[
                [
                    "cell", "tf", "motif_id", "motif_family", "TFBS_chr",
                    "TFBS_start", "TFBS_end", "TFBS_strand", "motif_score",
                    "peak_start", "peak_end", "label", "label_reason",
                    "nearest_peak_distance", "nearest_summit_distance",
                    "accessibility", "log_accessibility", "gc_fraction",
                    "peak_position_signed", "peak_position_abs", "artifact_index",
                ]
            ].copy()
            score["candidate_probability"] = candidate.probabilities
            score["reference_probability"] = reference.probabilities
            score["bias_configuration"] = str(route.bias_configuration)
            score["candidate_id"] = str(route.candidate_id)
            score["reference_candidate_id"] = base["reference_candidate_id"]
            score_frames.append(score)
            profile_rows.extend(
                aggregate_rows(
                    cell=cell,
                    tf=tf,
                    method="candidate",
                    profiles=candidate.evaluation_profiles,
                    labels=labels,
                    positions=positions,
                    iterations=args.profile_bootstrap,
                    seed=args.seed,
                )
            )
            profile_rows.extend(
                aggregate_rows(
                    cell=cell,
                    tf=tf,
                    method="reference",
                    profiles=reference.evaluation_profiles,
                    labels=labels,
                    positions=positions,
                    iterations=args.profile_bootstrap,
                    seed=args.seed,
                )
            )
        except Exception as error:
            metric_rows.append(
                {**base, "status": "error", "error": f"{type(error).__name__}: {error}"}
            )
    metrics = pd.DataFrame(metric_rows).sort_values(["cell", "tf"], kind="mergesort")
    replicate_metrics = pd.DataFrame(replicate_metric_rows).sort_values(
        ["cell", "tf", "replicate"], kind="mergesort"
    )
    replicate_ok = replicate_metrics[replicate_metrics["status"].eq("ok")].copy()
    if len(replicate_ok):
        stability = (
            replicate_ok.groupby(["cell", "tf"], sort=True)
            .agg(
                biological_replicates=("replicate", "nunique"),
                replicate_min_auroc_gain=("auroc_gain", "min"),
                replicate_max_auroc_gain=("auroc_gain", "max"),
                replicate_mean_auroc_gain=("auroc_gain", "mean"),
                replicate_min_relative_auprc_gain=("relative_auprc_gain", "min"),
            )
            .reset_index()
        )
        stability["replicate_direction_stable"] = (
            stability["replicate_min_auroc_gain"] >= 0
        )
        metrics = metrics.merge(stability, on=["cell", "tf"], how="left", validate="one_to_one")
    scores = pd.concat(score_frames, ignore_index=True) if score_frames else pd.DataFrame()
    profiles = pd.DataFrame(profile_rows)
    matching_frame = pd.DataFrame(matching_rows).sort_values(["cell", "tf"], kind="mergesort")
    metrics_path = args.outdir / "locked_holdout_metrics.tsv"
    scores_path = args.outdir / "locked_holdout_site_scores.tsv.gz"
    profiles_path = args.outdir / "locked_holdout_aggregate_profiles.tsv.gz"
    matching_path = args.outdir / "locked_holdout_matching_diagnostics.tsv"
    replicate_metrics_path = args.outdir / "locked_holdout_replicate_metrics.tsv"
    promotion_metrics_path = args.outdir / "locked_holdout_promotion_metrics.tsv"
    promotion_descriptors_path = args.outdir / "locked_holdout_promotion_descriptors.tsv"
    promotion_stability_path = args.outdir / "locked_holdout_promotion_stability.tsv"
    pdf_path = args.outdir / "locked_holdout_aggregate_panels.pdf"
    promotion_metrics = promotion_metric_table(metrics)
    promotion_descriptors = promotion_descriptor_table(metrics)
    promotion_stability = promotion_stability_table(metrics)
    metrics.to_csv(metrics_path, sep="\t", index=False)
    scores.to_csv(scores_path, sep="\t", index=False)
    profiles.to_csv(profiles_path, sep="\t", index=False)
    matching_frame.to_csv(matching_path, sep="\t", index=False)
    replicate_metrics.to_csv(replicate_metrics_path, sep="\t", index=False)
    promotion_metrics.to_csv(promotion_metrics_path, sep="\t", index=False)
    promotion_descriptors.to_csv(promotion_descriptors_path, sep="\t", index=False)
    promotion_stability.to_csv(promotion_stability_path, sep="\t", index=False)
    render_pdf(metrics, profiles, pdf_path)
    promotion = promotion_summary(metrics, study)
    manifest = {
        "schema": RESULT_SCHEMA,
        "locked_holdout_labels_read": True,
        "training_labels_used": False,
        "model_selection_after_holdout": False,
        "freeze": hash_record(args.freeze),
        "freeze_id": freeze["freeze_id"],
        "study": hash_record(args.study),
        "routes": hash_record(args.routes),
        "policy": hash_record(args.policy),
        "chip_manifest": hash_record(args.chip_manifest),
        "options": freeze_options(args),
        "promotion": promotion,
        "status_counts": metrics["status"].value_counts().sort_index().to_dict(),
        "outputs": {
            "metrics": hash_record(metrics_path),
            "site_scores": hash_record(scores_path),
            "aggregate_profiles": hash_record(profiles_path),
            "matching_diagnostics": hash_record(matching_path),
            "replicate_metrics": hash_record(replicate_metrics_path),
            "promotion_metrics": hash_record(promotion_metrics_path),
            "promotion_descriptors": hash_record(promotion_descriptors_path),
            "promotion_stability": hash_record(promotion_stability_path),
            "aggregate_pdf": hash_record(pdf_path),
        },
    }
    manifest_path = args.outdir / "locked_holdout_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(metrics.to_string(index=False))
    print(json.dumps(promotion, indent=2, sort_keys=True))
    return 0 if not metrics["status"].eq("error").any() else 2


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", action="append", type=parse_key_path, required=True, metavar="MODEL,CELL,JSON")
    parser.add_argument(
        "--replicate-artifact",
        action="append",
        type=parse_replicate_key_path,
        required=True,
        metavar="REPLICATE,MODEL,CELL,JSON",
    )
    parser.add_argument("--site-source", action="append", type=parse_cell_path, required=True, metavar="CELL=TSV")
    parser.add_argument("--study", type=Path, default=Path("benchmarks/manifests/footprint_functional_v1.spec.json"))
    parser.add_argument("--routes", type=Path, default=Path("benchmarks/manifests/compact/functional_holdout_routes_v1.tsv"))
    parser.add_argument("--policy", type=Path, default=Path("benchmarks/manifests/compact/functional_detector_policy_v1.tsv"))
    parser.add_argument("--chip-manifest", type=Path, default=Path("benchmarks/manifests/compact/functional_holdout_chip_peaks.tsv"))
    parser.add_argument("--genome", type=Path, default=Path("data/public/raw/genome/hg38.fa"))
    parser.add_argument("--freeze-out", type=Path)
    parser.add_argument("--freeze", type=Path)
    parser.add_argument("--unlock-holdout", action="store_true")
    parser.add_argument("--outdir", type=Path)
    parser.add_argument("--positive-summit-distance", type=int, default=100)
    parser.add_argument("--negative-peak-distance", type=int, default=500)
    parser.add_argument("--minimum-sites-per-class", type=int, default=200)
    parser.add_argument("--maximum-train-per-tf", type=int, default=10000)
    parser.add_argument("--maximum-matching-smd", type=float, default=0.25)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--profile-bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args(argv)
    if args.minimum_sites_per_class < 2 or args.maximum_train_per_tf < 100:
        raise SystemExit("site limits are too small for the preregistered evaluation")
    if args.bootstrap < 2 or args.profile_bootstrap < 2:
        raise SystemExit("bootstrap counts must be at least two")
    if args.positive_summit_distance < 0 or args.negative_peak_distance < 0:
        raise SystemExit("label distances must be non-negative")
    study = json.loads(args.study.read_text(encoding="utf-8"))
    routes = pd.read_csv(args.routes, sep="\t")
    policy = pd.read_csv(args.policy, sep="\t")
    chip = pd.read_csv(args.chip_manifest, sep="\t", dtype={"checksum": str})
    validate_tables(study, routes, policy, chip)
    artifact_paths = {(model, cell): path for model, cell, path in args.artifact}
    if len(artifact_paths) != len(args.artifact):
        raise SystemExit("duplicate MODEL,CELL artifact keys")
    required = required_artifact_keys(routes)
    if set(artifact_paths) != required:
        missing = sorted(required.difference(artifact_paths))
        extra = sorted(set(artifact_paths).difference(required))
        raise SystemExit(f"artifact keys differ from frozen routes; missing={missing}, extra={extra}")
    artifacts = {
        key: load_artifact(key[0], key[1], path) for key, path in sorted(artifact_paths.items())
    }
    validate_exact_site_alignment(artifacts)
    replicate_paths = {
        (replicate, model, cell): path
        for replicate, model, cell, path in args.replicate_artifact
    }
    if len(replicate_paths) != len(args.replicate_artifact):
        raise SystemExit("duplicate REPLICATE,MODEL,CELL artifact keys")
    replicate_artifacts = {
        key: load_artifact(key[1], key[2], path)
        for key, path in sorted(replicate_paths.items())
    }
    validate_replicate_artifacts(artifacts, replicate_artifacts, routes)
    site_sources = dict(args.site_source)
    if set(site_sources) != set(routes["cell"].astype(str)):
        raise SystemExit("site-source cells must exactly match route cells")
    for cell, source in site_sources.items():
        validate_site_source(artifacts[("DWM", cell)], source)
    if not args.unlock_holdout:
        if args.freeze_out is None:
            raise SystemExit("labels remain locked; provide --freeze-out to preregister inputs")
        if args.freeze is not None or args.outdir is not None:
            raise SystemExit("freeze mode does not accept --freeze or --outdir")
        document = build_freeze_document(
            args, artifacts, replicate_artifacts, site_sources, routes
        )
        args.freeze_out.parent.mkdir(parents=True, exist_ok=True)
        args.freeze_out.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"wrote locked holdout freeze {document['freeze_id']} to {args.freeze_out}")
        print("ChIP labels were not opened")
        return 0
    if args.freeze is None or args.outdir is None:
        raise SystemExit("--unlock-holdout requires --freeze and --outdir")
    if args.freeze_out is not None:
        raise SystemExit("evaluation mode does not accept --freeze-out")
    freeze = verify_freeze(
        args.freeze, args, artifacts, replicate_artifacts, site_sources, routes
    )
    validate_chip_files(chip)
    with threadpool_limits(limits=1):
        return run_evaluation(
            args,
            study,
            routes,
            policy,
            chip,
            artifacts,
            replicate_artifacts,
            freeze,
        )


if __name__ == "__main__":
    raise SystemExit(main())
