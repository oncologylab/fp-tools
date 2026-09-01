#!/usr/bin/env python3
"""Freeze per-TF call thresholds from development matched negatives."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from build_strand_functional_profiles import site_hashes  # noqa: E402
from evaluate_functional_footprints import load_or_extract_profiles  # noqa: E402
from evaluate_frozen_functional_policy import (  # noqa: E402
    candidate_score_and_profile,
    preflight_test_artifact,
    validate_policy,
)
from evaluate_parametric_factorization import (  # noqa: E402
    align_baseline,
    load_dwm_baseline,
    load_safe_configuration,
    orient_aligned_baseline,
    parse_name_path,
    residual_score,
)
from evaluate_strand_label_free_models import (  # noqa: E402
    file_sha256,
    parse_artifact,
)
from fp_tools.tools.parametric_factorization import (  # noqa: E402
    FrozenParametricFactorization,
)
from freeze_label_free_functional_models import immutable_write_json  # noqa: E402


SCHEMA = "fp-tools-functional-call-thresholds-v1"


def upper_tail_threshold(scores: np.ndarray, target_rate: float) -> tuple[float, int]:
    values = np.asarray(scores, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0 or not 0.0 < target_rate < 1.0:
        raise ValueError("scores must be finite and target_rate must be in (0, 1)")
    allowed = int(np.floor(target_rate * len(values)))
    if allowed == 0:
        return float(np.nextafter(np.max(values), np.inf)), 0
    for threshold in np.unique(values):
        calls = int(np.sum(values >= threshold))
        if calls <= allowed:
            return float(threshold), calls
    return float(np.nextafter(np.max(values), np.inf)), 0


def validation_site_score_frame(
    *,
    record: dict,
    candidate_id: str,
    sites: pd.DataFrame,
    indexes: np.ndarray,
    site_hash: np.ndarray,
    candidate_score: np.ndarray,
    dwm_score: np.ndarray,
) -> pd.DataFrame:
    selected = sites.iloc[indexes].reset_index(drop=True)
    if not (
        len(selected)
        == len(site_hash)
        == len(candidate_score)
        == len(dwm_score)
    ):
        raise ValueError("validation sites and score arrays have different lengths")
    return pd.DataFrame(
        {
            "cell": str(record["cell"]),
            "tf": str(record["tf"]),
            "motif_family": str(record["motif_family"]),
            "bias_configuration": str(record["bias_configuration"]),
            "candidate_id": candidate_id,
            "artifact_index": np.asarray(indexes, dtype=int),
            "site_hash": np.asarray(site_hash, dtype=np.uint64),
            "TFBS_chr": selected["TFBS_chr"].astype(str),
            "label": selected["chip_label"].to_numpy(dtype=int),
            "candidate_score": np.asarray(candidate_score, dtype=float),
            "dwm_score": np.asarray(dwm_score, dtype=float),
        }
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument(
        "--validation-artifact",
        action="append",
        type=parse_artifact,
        required=True,
        metavar="MODEL,CELL,JSON",
    )
    parser.add_argument(
        "--dwm-baseline",
        action="append",
        type=parse_name_path,
        metavar="CELL=NPZ",
    )
    parser.add_argument(
        "--dwm-signal",
        action="append",
        type=parse_name_path,
        metavar="CELL=EXPECTED_BIGWIG",
    )
    parser.add_argument("--reference-configuration", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--target-negative-call-rate", type=float, default=0.025)
    args = parser.parse_args(argv)

    policy, models = validate_policy(args.policy)
    configuration = load_safe_configuration(args.reference_configuration)
    factorization = FrozenParametricFactorization.load(
        configuration["factorization_model"]["path"]
    )
    dispersion = float(factorization.total_dispersion_)
    artifacts = {
        (model, cell): path for model, cell, path in args.validation_artifact
    }
    baselines_paths = dict(args.dwm_baseline or [])
    baseline_signals = dict(args.dwm_signal or [])
    policy_keys = {
        (str(record["bias_configuration"]), str(record["cell"]))
        for record, _candidate, _model in models
    }
    if set(artifacts) != policy_keys:
        raise ValueError("validation artifacts do not exactly match policy cells/models")
    expected_cells = {cell for _model, cell in policy_keys}
    if bool(baselines_paths) == bool(baseline_signals):
        raise ValueError("provide exactly one of --dwm-baseline or --dwm-signal")
    provided_cells = set(baselines_paths or baseline_signals)
    if provided_cells != expected_cells:
        raise ValueError("DWM inputs do not exactly match policy cells")

    datasets = {}
    baselines = {}
    for key, path in sorted(artifacts.items()):
        document, arrays = preflight_test_artifact(path, expected_cell=key[1])
        sites = pd.read_csv(document["sites"], sep="\t").reset_index(drop=True)
        if not np.array_equal(site_hashes(sites), arrays["site_hash"]):
            raise ValueError(f"validation site order mismatch: {path}")
        if "cell" not in sites or set(sites["cell"].astype(str)) != {key[1]}:
            raise ValueError(f"validation sites do not exclusively contain {key[1]}")
        datasets[key] = (sites, arrays, document)
    args.outdir.mkdir(parents=True, exist_ok=True)
    if baselines_paths:
        for cell, path in sorted(baselines_paths.items()):
            baselines[cell], _inputs = load_dwm_baseline(path)
    else:
        for key in sorted(datasets):
            cell = key[1]
            sites, _arrays, _document = datasets[key]
            cache = args.outdir / "dwm_profile_cache" / f"{cell}.npz"
            profiles = load_or_extract_profiles(
                sites,
                baseline_signals[cell],
                cache,
                100,
            )
            with np.load(cache, allow_pickle=False) as source:
                baselines[cell] = {
                    "expected": profiles,
                    "valid": np.asarray(source["valid"], dtype=bool),
                    "site_hash": np.asarray(source["site_hash"], dtype=np.uint64),
                    "orientation_aligned": np.asarray(False),
                }

    rows = []
    site_score_frames = []
    with threadpool_limits(limits=1):
        for record, candidate, model in models:
            key = (str(record["bias_configuration"]), str(record["cell"]))
            sites, arrays, _document = datasets[key]
            baseline_expected, baseline_valid = align_baseline(
                arrays, baselines[key[1]]
            )
            baseline_expected = orient_aligned_baseline(
                baseline_expected,
                baselines[key[1]],
                sites,
            )
            selected = (
                sites["tf"].astype(str).eq(str(record["tf"])).to_numpy()
                & sites["chromosome_split"].astype(str).eq("validation").to_numpy()
                & arrays["valid"].astype(bool)
                & baseline_valid
            )
            indexes = np.flatnonzero(selected)
            labels = sites.iloc[indexes]["chip_label"].to_numpy(dtype=int)
            negative = labels == 0
            if not np.any(negative):
                raise ValueError(
                    f"no validation negatives for {record['cell']} {record['tf']}"
                )
            candidate_score, _profiles, _fitted = candidate_score_and_profile(
                candidate,
                model,
                arrays,
                indexes,
                np.arange(arrays["plus_observed"].shape[1], dtype=float)
                - arrays["plus_observed"].shape[1] // 2,
                motif_score=sites.iloc[indexes]["motif_score"].to_numpy(dtype=float),
            )
            observed = arrays["plus_observed"][indexes] + arrays["minus_observed"][indexes]
            dwm_score, _dwm_profiles = residual_score(
                observed,
                baseline_expected[indexes],
                np.arange(observed.shape[1], dtype=float) - observed.shape[1] // 2,
                "deviance",
                dispersion,
            )
            site_score_frames.append(
                validation_site_score_frame(
                    record=record,
                    candidate_id=candidate.candidate_id,
                    sites=sites,
                    indexes=indexes,
                    site_hash=arrays["site_hash"][indexes],
                    candidate_score=candidate_score,
                    dwm_score=dwm_score,
                )
            )
            for method, score in (
                ("candidate", candidate_score),
                ("DWM", dwm_score),
            ):
                threshold, calls = upper_tail_threshold(
                    score[negative], args.target_negative_call_rate
                )
                rows.append(
                    {
                        "cell": record["cell"],
                        "tf": record["tf"],
                        "motif_family": record["motif_family"],
                        "bias_configuration": record["bias_configuration"],
                        "candidate_id": candidate.candidate_id,
                        "method": method,
                        "threshold": threshold,
                        "validation_negative_sites": int(np.sum(negative)),
                        "validation_negative_calls": calls,
                        "validation_negative_call_rate": calls / np.sum(negative),
                        "target_negative_call_rate": args.target_negative_call_rate,
                    }
                )

    thresholds = pd.DataFrame(rows)
    thresholds_path = args.outdir / "functional_call_thresholds.tsv"
    site_scores_path = args.outdir / "functional_validation_site_scores.tsv.gz"
    thresholds.to_csv(thresholds_path, sep="\t", index=False)
    pd.concat(site_score_frames, ignore_index=True).to_csv(
        site_scores_path,
        sep="\t",
        index=False,
    )
    document = {
        "schema": SCHEMA,
        "policy_id": policy["policy_id"],
        "policy": {"path": str(args.policy), "sha256": file_sha256(args.policy)},
        "reference_configuration": {
            "path": str(args.reference_configuration),
            "sha256": file_sha256(args.reference_configuration),
        },
        "validation_artifacts": [
            {"model": key[0], "cell": key[1], "path": str(path), "sha256": file_sha256(path)}
            for key, path in sorted(artifacts.items())
        ],
        "dwm_baselines": [
            {"cell": cell, "path": str(path), "sha256": file_sha256(path)}
            for cell, path in sorted((baselines_paths or baseline_signals).items())
        ],
        "target_negative_call_rate": args.target_negative_call_rate,
        "thresholds": {
            "path": str(thresholds_path),
            "sha256": file_sha256(thresholds_path),
        },
        "validation_site_scores": {
            "path": str(site_scores_path),
            "sha256": file_sha256(site_scores_path),
        },
        "training_labels_used": False,
        "validation_labels_used_for_thresholds": True,
        "naked_dna_read": False,
    }
    canonical = json.dumps(document, sort_keys=True, separators=(",", ":"))
    document["threshold_id"] = sha256(canonical.encode()).hexdigest()
    output = args.outdir / "functional_call_thresholds.freeze.json"
    immutable_write_json(output, document)
    print(thresholds.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
