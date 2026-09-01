#!/usr/bin/env python3
"""Screen a serialized functional policy on independent naked-DNA profiles."""

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
from evaluate_frozen_functional_policy import (  # noqa: E402
    candidate_score_and_profile,
    validate_policy,
)
from evaluate_naked_dna_functional_policy import wilson_interval  # noqa: E402
from evaluate_parametric_factorization import (  # noqa: E402
    load_safe_configuration,
    residual_score,
)
from evaluate_strand_label_free_models import file_sha256  # noqa: E402
from fp_tools.tools.functional_footprints import normalize_functional_profiles  # noqa: E402
from fp_tools.tools.parametric_factorization import (  # noqa: E402
    FrozenParametricFactorization,
)
from freeze_functional_call_thresholds import SCHEMA as THRESHOLD_SCHEMA  # noqa: E402
from freeze_label_free_functional_models import immutable_write_json  # noqa: E402


SCHEMA = "fp-tools-frozen-functional-naked-dna-v1"


def parse_candidate(value: str) -> tuple[str, str, str, Path]:
    fields = value.split(",", 3)
    if len(fields) != 4 or not all(fields):
        raise argparse.ArgumentTypeError(
            "candidate artifact must use MODEL,CELL,REPLICATE,JSON"
        )
    return fields[0], fields[1], fields[2], Path(fields[3])


def parse_dwm(value: str) -> tuple[str, str, Path]:
    fields = value.split(",", 2)
    if len(fields) != 3 or not all(fields):
        raise argparse.ArgumentTypeError("DWM artifact must use CELL,REPLICATE,JSON")
    return fields[0], fields[1], Path(fields[2])


def validate_threshold_freeze(path: Path, policy_id: str) -> tuple[dict, pd.DataFrame]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != THRESHOLD_SCHEMA:
        raise ValueError("unsupported functional threshold freeze")
    if document.get("policy_id") != policy_id:
        raise ValueError("threshold freeze does not match functional policy")
    content = dict(document)
    observed = str(content.pop("threshold_id"))
    canonical = json.dumps(content, sort_keys=True, separators=(",", ":"))
    if observed != sha256(canonical.encode()).hexdigest():
        raise ValueError("functional threshold ID does not match its document")
    record = document["thresholds"]
    if file_sha256(record["path"]) != record["sha256"]:
        raise ValueError("functional threshold table checksum mismatch")
    return document, pd.read_csv(record["path"], sep="\t")


def preflight_artifact(
    path: Path,
    *,
    schema: str,
) -> tuple[dict, dict[str, np.ndarray]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema") != schema:
        raise ValueError(f"unsupported naked-DNA artifact: {path}")
    if document.get("metadata", {}).get("labels_used") is not False:
        raise ValueError(f"naked-DNA artifact does not certify label-free construction: {path}")
    if file_sha256(document["profiles_npz"]) != document["profiles_sha256"]:
        raise ValueError(f"naked-DNA profile checksum mismatch: {path}")
    if file_sha256(document["sites"]) != document["sites_sha256"]:
        raise ValueError(f"naked-DNA site checksum mismatch: {path}")
    with np.load(document["profiles_npz"], allow_pickle=False) as source:
        arrays = {name: np.asarray(source[name]) for name in source.files}
    return document, arrays


def load_sites(document: dict, arrays: dict[str, np.ndarray], cell: str) -> pd.DataFrame:
    sites = pd.read_csv(document["sites"], sep="\t").reset_index(drop=True)
    forbidden = [
        column
        for column in sites.columns
        if "label" in column.lower() or "chip" in column.lower()
    ]
    if forbidden:
        raise ValueError(
            "naked-DNA sites contain forbidden columns: " + ", ".join(forbidden)
        )
    if "cell" not in sites or set(sites["cell"].astype(str)) != {cell}:
        raise ValueError(f"naked-DNA sites do not exclusively contain {cell}")
    if not np.array_equal(site_hashes(sites), arrays["site_hash"]):
        raise ValueError("naked-DNA site order does not match profile arrays")
    return sites


def rate_record(
    score: np.ndarray,
    valid: np.ndarray,
    informative: np.ndarray,
    threshold: float,
) -> tuple[dict, np.ndarray]:
    finite = np.asarray(valid, dtype=bool) & np.isfinite(score)
    informative = finite & np.asarray(informative, dtype=bool)
    calls = informative & (score >= threshold)
    valid_calls = int(np.sum(calls))
    valid_total = int(np.sum(finite))
    informative_total = int(np.sum(informative))
    valid_low, valid_high = wilson_interval(valid_calls, valid_total)
    info_low, info_high = wilson_interval(valid_calls, informative_total)
    return (
        {
            "valid_sites": valid_total,
            "informative_sites": informative_total,
            "calls": valid_calls,
            "false_positive_rate": valid_calls / valid_total if valid_total else np.nan,
            "false_positive_rate_lower_95": valid_low,
            "false_positive_rate_upper_95": valid_high,
            "informative_false_positive_rate": (
                valid_calls / informative_total if informative_total else np.nan
            ),
            "informative_false_positive_rate_lower_95": info_low,
            "informative_false_positive_rate_upper_95": info_high,
        },
        calls,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--threshold-freeze", type=Path, required=True)
    parser.add_argument(
        "--candidate-artifact",
        action="append",
        type=parse_candidate,
        required=True,
    )
    parser.add_argument("--dwm-artifact", action="append", type=parse_dwm, required=True)
    parser.add_argument("--reference-configuration", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)

    policy, models = validate_policy(args.policy)
    threshold_document, thresholds = validate_threshold_freeze(
        args.threshold_freeze,
        policy["policy_id"],
    )
    configuration = load_safe_configuration(args.reference_configuration)
    factorization = FrozenParametricFactorization.load(
        configuration["factorization_model"]["path"]
    )
    dispersion = float(factorization.total_dispersion_)
    candidate_paths = {
        (model, cell, replicate): path
        for model, cell, replicate, path in args.candidate_artifact
    }
    dwm_paths = {
        (cell, replicate): path for cell, replicate, path in args.dwm_artifact
    }
    policy_keys = {
        (str(record["bias_configuration"]), str(record["cell"]))
        for record, _candidate, _model in models
    }
    candidate_keys = {(model, cell) for model, cell, _replicate in candidate_paths}
    if candidate_keys != policy_keys:
        raise ValueError("candidate naked-DNA artifacts do not match policy cells/models")
    if {(cell, replicate) for _model, cell, replicate in candidate_paths} != set(dwm_paths):
        raise ValueError("candidate and DWM naked-DNA replicates do not match")

    candidates = {}
    dwm = {}
    input_records = []
    for key, path in sorted(candidate_paths.items()):
        document, arrays = preflight_artifact(
            path,
            schema="fp-tools-strand-functional-profiles-v1",
        )
        candidates[key] = (document, arrays)
        input_records.append(
            {"purpose": "candidate", "path": str(path), "sha256": file_sha256(path)}
        )
    for key, path in sorted(dwm_paths.items()):
        document, arrays = preflight_artifact(
            path,
            schema="fp-tools-combined-functional-profiles-v1",
        )
        dwm[key] = (document, arrays)
        input_records.append(
            {"purpose": "DWM", "path": str(path), "sha256": file_sha256(path)}
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    input_freeze = {
        "schema": "fp-tools-frozen-functional-naked-dna-inputs-v1",
        "policy_id": policy["policy_id"],
        "threshold_id": threshold_document["threshold_id"],
        "policy": {"path": str(args.policy), "sha256": file_sha256(args.policy)},
        "threshold_freeze": {
            "path": str(args.threshold_freeze),
            "sha256": file_sha256(args.threshold_freeze),
        },
        "inputs": input_records,
        "models_refitted_on_naked_dna": False,
        "thresholds_changed_on_naked_dna": False,
    }
    canonical = json.dumps(input_freeze, sort_keys=True, separators=(",", ":"))
    input_freeze["naked_input_id"] = sha256(canonical.encode()).hexdigest()
    input_freeze_path = args.outdir / "naked_dna_inputs.freeze.json"
    immutable_write_json(input_freeze_path, input_freeze)

    rate_rows = []
    site_rows = []
    profile_rows = []
    with threadpool_limits(limits=1):
        for record, candidate, model in models:
            model_key = (str(record["bias_configuration"]), str(record["cell"]))
            matching = [
                key for key in candidate_paths if key[:2] == model_key
            ]
            for key in matching:
                replicate = key[2]
                candidate_document, candidate_arrays = candidates[key]
                dwm_document, dwm_arrays = dwm[(key[1], replicate)]
                sites = load_sites(candidate_document, candidate_arrays, key[1])
                dwm_sites = load_sites(dwm_document, dwm_arrays, key[1])
                if not np.array_equal(
                    candidate_arrays["site_hash"], dwm_arrays["site_hash"]
                ) or not sites.equals(dwm_sites):
                    raise ValueError("candidate and DWM naked-DNA sites differ")
                selected = sites["tf"].astype(str).eq(str(record["tf"])).to_numpy()
                indexes = np.flatnonzero(selected)
                if not len(indexes):
                    continue
                positions = np.arange(candidate_arrays["plus_observed"].shape[1], dtype=float)
                positions -= candidate_arrays["plus_observed"].shape[1] // 2
                candidate_score, candidate_profiles, _fitted = candidate_score_and_profile(
                    candidate,
                    model,
                    candidate_arrays,
                    indexes,
                    positions,
                    motif_score=sites.iloc[indexes]["motif_score"].to_numpy(dtype=float),
                )
                dwm_score, dwm_profiles = residual_score(
                    dwm_arrays["observed"][indexes],
                    dwm_arrays["expected"][indexes],
                    positions,
                    "deviance",
                    dispersion,
                )
                dwm_profiles = normalize_functional_profiles(dwm_profiles, positions)
                candidate_valid = candidate_arrays["valid"][indexes].astype(bool)
                dwm_valid = dwm_arrays["valid"][indexes].astype(bool)
                candidate_info = (
                    candidate_arrays["plus_observed"][indexes]
                    + candidate_arrays["minus_observed"][indexes]
                ).sum(axis=1) > 0
                dwm_info = dwm_arrays["observed"][indexes].sum(axis=1) > 0
                threshold_rows = thresholds[
                    thresholds["cell"].astype(str).eq(str(record["cell"]))
                    & thresholds["tf"].astype(str).eq(str(record["tf"]))
                ]
                if set(threshold_rows["method"].astype(str)) != {"candidate", "DWM"}:
                    raise ValueError("threshold table lacks candidate/DWM pair")
                scores = {
                    "candidate": (
                        candidate_score,
                        candidate_profiles,
                        candidate_valid,
                        candidate_info,
                    ),
                    "DWM": (dwm_score, dwm_profiles, dwm_valid, dwm_info),
                }
                calls_by_method = {}
                rates_by_method = {}
                finite_by_method = {}
                for method, (score, profiles, valid, informative) in scores.items():
                    threshold = float(
                        threshold_rows.loc[
                            threshold_rows["method"].astype(str).eq(method),
                            "threshold",
                        ].iloc[0]
                    )
                    rate, calls = rate_record(score, valid, informative, threshold)
                    calls_by_method[method] = calls
                    rates_by_method[method] = rate
                    finite_by_method[method] = np.asarray(valid, dtype=bool) & np.isfinite(
                        score
                    )
                    rate_rows.append(
                        {
                            "cell": record["cell"],
                            "tf": record["tf"],
                            "motif_family": record["motif_family"],
                            "candidate_id": candidate.candidate_id,
                            "replicate": replicate,
                            "method": method,
                            "threshold": threshold,
                            **rate,
                        }
                    )
                    normalized = normalize_functional_profiles(profiles, positions)
                    for group, mask in (
                        ("informative", valid & informative),
                        ("called", calls),
                    ):
                        if not np.any(mask):
                            continue
                        mean = np.nanmean(normalized[mask], axis=0)
                        for offset, position in enumerate(positions.astype(int)):
                            profile_rows.append(
                                {
                                    "cell": record["cell"],
                                    "tf": record["tf"],
                                    "motif_family": record["motif_family"],
                                    "candidate_id": candidate.candidate_id,
                                    "replicate": replicate,
                                    "method": method,
                                    "group": group,
                                    "position": position,
                                    "mean": mean[offset],
                                }
                            )
                    for local, site_index in enumerate(indexes):
                        site_rows.append(
                            {
                                "cell": record["cell"],
                                "tf": record["tf"],
                                "replicate": replicate,
                                "site_hash": int(candidate_arrays["site_hash"][site_index]),
                                "method": method,
                                "score": score[local],
                                "threshold": threshold,
                                "valid": bool(valid[local]),
                                "informative": bool(informative[local]),
                                "called": bool(calls[local]),
                            }
                        )
                # A zero-cut motif site is an informative negative result for the
                # point-rate denominator, not a missing observation.  Pair on
                # finite scores and report nonzero-cut support separately.
                paired = finite_by_method["candidate"] & finite_by_method["DWM"]
                candidate_calls = calls_by_method["candidate"] & paired
                dwm_calls = calls_by_method["DWM"] & paired
                total = int(np.sum(paired))
                candidate_count = int(np.sum(candidate_calls))
                dwm_count = int(np.sum(dwm_calls))
                _low, upper = wilson_interval(candidate_count, total)
                candidate_rate = candidate_count / total if total else np.nan
                dwm_rate = dwm_count / total if total else np.nan
                candidate_independent = rates_by_method["candidate"]
                candidate_point = float(candidate_independent["false_positive_rate"])
                candidate_upper = float(
                    candidate_independent["false_positive_rate_upper_95"]
                )
                increase = candidate_rate - dwm_rate
                rate_rows.append(
                    {
                        "cell": record["cell"],
                        "tf": record["tf"],
                        "motif_family": record["motif_family"],
                        "candidate_id": candidate.candidate_id,
                        "replicate": replicate,
                        "method": "paired_safety",
                        "paired_sites": total,
                        "paired_informative_sites": int(
                            np.sum(paired & candidate_info & dwm_info)
                        ),
                        "candidate_calls": candidate_count,
                        "dwm_calls": dwm_count,
                        "candidate_finite_sites": int(
                            candidate_independent["valid_sites"]
                        ),
                        "candidate_informative_sites": int(
                            candidate_independent["informative_sites"]
                        ),
                        "candidate_false_positive_rate": candidate_point,
                        "candidate_wilson_upper_95": candidate_upper,
                        "paired_candidate_false_positive_rate": candidate_rate,
                        "paired_dwm_false_positive_rate": dwm_rate,
                        "candidate_minus_dwm": increase,
                        "paired_candidate_wilson_upper_95": upper,
                        "passes_point_rate": bool(candidate_point <= 0.05),
                        "passes_wilson": bool(candidate_upper <= 0.05),
                        "passes_increase": bool(increase <= 0.01),
                        "passes_safety": bool(
                            candidate_point <= 0.05
                            and candidate_upper <= 0.05
                            and increase <= 0.01
                        ),
                    }
                )

    rates = pd.DataFrame(rate_rows)
    rates_path = args.outdir / "naked_dna_false_positive_rates.tsv"
    scores_path = args.outdir / "naked_dna_site_scores.tsv.gz"
    profiles_path = args.outdir / "naked_dna_profiles.tsv.gz"
    rates.to_csv(rates_path, sep="\t", index=False)
    pd.DataFrame(site_rows).to_csv(scores_path, sep="\t", index=False)
    pd.DataFrame(profile_rows).to_csv(profiles_path, sep="\t", index=False)
    paired = rates[rates["method"].astype(str).eq("paired_safety")]
    manifest = {
        "schema": SCHEMA,
        "policy_id": policy["policy_id"],
        "threshold_id": threshold_document["threshold_id"],
        "naked_input_id": input_freeze["naked_input_id"],
        "models_refitted_on_naked_dna": False,
        "thresholds_changed_on_naked_dna": False,
        "tasks": int(len(paired)),
        "tasks_passing_safety": int(paired["passes_safety"].eq(True).sum()),
        "all_tasks_pass_safety": bool(paired["passes_safety"].eq(True).all()),
        "outputs": {
            "rates": {"path": str(rates_path), "sha256": file_sha256(rates_path)},
            "scores": {"path": str(scores_path), "sha256": file_sha256(scores_path)},
            "profiles": {"path": str(profiles_path), "sha256": file_sha256(profiles_path)},
        },
    }
    (args.outdir / "naked_dna_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    columns = [
        "cell",
        "tf",
        "paired_sites",
        "candidate_false_positive_rate",
        "paired_dwm_false_positive_rate",
        "candidate_wilson_upper_95",
        "candidate_minus_dwm",
        "passes_safety",
    ]
    print(paired[columns].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
