#!/usr/bin/env python3
"""Evaluate raw-preserving partial subtraction of frozen sequence bias.

This research-only ablation asks whether full residualization discards useful
footprint geometry.  Tune mode selects one global parametric shrinkage and one
TF-specific shrinkage on validation chromosomes.  Test mode validates the
immutable policy and applies it without refitting on the already-opened
development test chromosomes.  The production DWM correction is unchanged.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(REPOSITORY / "src"))
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_frozen_functional_policy import aggregate_curve  # noqa: E402
from evaluate_functional_footprints import binary_metrics  # noqa: E402
from evaluate_parametric_factorization import (  # noqa: E402
    DIFFICULT_ROLES,
    align_baseline,
    block_bootstrap_delta,
    geometry_score,
    load_dwm_baseline,
    load_profiles,
    load_safe_configuration,
    orient_aligned_baseline,
    residual_score,
    sha256_file,
)
from fp_tools.tools.functional_footprints import (  # noqa: E402
    standardized_functional_separation,
)
from fp_tools.tools.parametric_factorization import (  # noqa: E402
    FrozenBiasStrengthCalibrator,
    expected_profile_counts,
)


POLICY_SCHEMA = "fp-tools-frozen-bias-shrinkage-policy-v1"
TEST_SCHEMA = "fp-tools-frozen-bias-shrinkage-test-v1"
INPUT_SCHEMA = "fp-tools-frozen-bias-shrinkage-test-inputs-v1"
PARAMETRIC_SOURCES = ("parametric_direct", "parametric_lambda")
DEFAULT_ALPHAS = (
    0.025,
    0.05,
    0.10,
    0.15,
    0.20,
    0.30,
    0.40,
    0.50,
    0.65,
    0.80,
    1.00,
    1.25,
    1.50,
    2.00,
)


def parse_name_path(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("values must use CELL=PATH")
    name, raw_path = value.split("=", 1)
    if not name or not raw_path:
        raise argparse.ArgumentTypeError("values must use CELL=PATH")
    return name, Path(raw_path)


def parse_alphas(value: str) -> tuple[float, ...]:
    try:
        alphas = tuple(float(item) for item in value.split(",") if item)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("alphas must be comma-separated numbers") from exc
    if not alphas or any(not np.isfinite(item) or item <= 0 or item > 2 for item in alphas):
        raise argparse.ArgumentTypeError("alphas must be finite values in (0, 2]")
    if len(set(alphas)) != len(alphas):
        raise argparse.ArgumentTypeError("alphas must not contain duplicates")
    return tuple(sorted(alphas))


def artifact_paths(prefix: Path) -> tuple[Path, Path, Path]:
    base = prefix.with_suffix("") if prefix.suffix in {".json", ".npz"} else prefix
    return (
        Path(str(base) + ".npz"),
        Path(str(base) + ".json"),
        Path(str(base) + ".sites.tsv.gz"),
    )


def input_record(path: Path, purpose: str) -> dict[str, str]:
    return {"path": str(path), "sha256": sha256_file(path), "purpose": purpose}


def canonical_id(document: dict) -> str:
    payload = json.dumps(document, sort_keys=True, separators=(",", ":"))
    return sha256(payload.encode()).hexdigest()


def write_immutable_json(path: Path, document: dict) -> None:
    rendered = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if path.exists() and path.read_text(encoding="utf-8") != rendered:
        raise ValueError(f"immutable artifact differs: {path}")
    path.write_text(rendered, encoding="utf-8")


def verify_records(records: Iterable[dict]) -> None:
    for record in records:
        path = Path(record["path"])
        if not path.is_file() or sha256_file(path) != record["sha256"]:
            raise ValueError(f"frozen input is absent or changed: {path}")


def validate_policy(path: Path) -> dict:
    policy = json.loads(path.read_text(encoding="utf-8"))
    if policy.get("schema") != POLICY_SCHEMA:
        raise ValueError("unsupported frozen shrinkage policy")
    claimed = str(policy.get("policy_id", ""))
    unsigned = dict(policy)
    unsigned.pop("policy_id", None)
    if canonical_id(unsigned) != claimed:
        raise ValueError("frozen shrinkage policy ID is invalid")
    if policy.get("selection_split") != "validation":
        raise ValueError("shrinkage policy was not selected on validation chromosomes")
    if policy.get("bias_coefficients_refitted") is not False:
        raise ValueError("shrinkage policy does not certify frozen bias coefficients")
    verify_records(policy.get("inputs", []))
    return policy


def load_calibrator(configuration: dict) -> FrozenBiasStrengthCalibrator:
    record = configuration["bias_calibration"]
    if sha256_file(record["path"]) != record["sha256"]:
        raise ValueError("frozen bias calibration changed")
    return FrozenBiasStrengthCalibrator.load(record["path"])


def load_panel(
    *,
    cell: str,
    artifact: Path,
    dwm_path: Path,
    split: str,
    split_chromosomes: set[str],
    calibrator: FrozenBiasStrengthCalibrator,
) -> tuple[dict, list[dict[str, str]]]:
    arrays, sites, document = load_profiles(artifact, require_log_bias=True)
    if document.get("metadata", {}).get("labels_used") is not False:
        raise ValueError(f"profile artifact does not certify label-free construction: {artifact}")
    if set(sites["cell"].astype(str)) != {cell}:
        raise ValueError(f"profile artifact is not exclusive to {cell}: {artifact}")
    selected_split = sites["chromosome_split"].astype(str).eq(split).to_numpy()
    observed_chromosomes = set(sites.loc[selected_split, "TFBS_chr"].astype(str))
    if not observed_chromosomes.issubset(split_chromosomes):
        raise ValueError(f"{cell} {split} artifact contains chromosomes outside the frozen split")
    baseline, baseline_paths = load_dwm_baseline(dwm_path)
    dwm_expected, dwm_valid = align_baseline(arrays, baseline)
    dwm_expected = orient_aligned_baseline(dwm_expected, baseline, sites)
    observed = arrays["plus_observed"] + arrays["minus_observed"]
    direct_expected = arrays["plus_expected"] + arrays["minus_expected"]
    valid = (
        selected_split
        & arrays["valid"].astype(bool)
        & dwm_valid
        & np.isfinite(observed).all(axis=1)
        & np.isfinite(direct_expected).all(axis=1)
        & np.isfinite(dwm_expected).all(axis=1)
        & np.isfinite(arrays["combined_log_bias"]).all(axis=1)
        & (observed.sum(axis=1) > 0)
    )
    indexes = np.flatnonzero(valid)
    if not len(indexes):
        raise ValueError(f"{cell} has no valid {split} profiles")
    observed = observed[indexes].astype(np.float64)
    log_bias = arrays["combined_log_bias"][indexes].astype(np.float64)
    lambda_expected = expected_profile_counts(
        observed,
        calibrator.strength(cell) * log_bias,
    )
    selected_sites = sites.iloc[indexes].reset_index(drop=True)
    positions = np.arange(observed.shape[1], dtype=float) - observed.shape[1] // 2
    npz_path, json_path, sites_path = artifact_paths(artifact)
    inputs = [
        input_record(npz_path, f"{cell}-{split}-profiles"),
        input_record(json_path, f"{cell}-{split}-profile-manifest"),
        input_record(sites_path, f"{cell}-{split}-sites"),
        *[input_record(path, f"{cell}-{split}-DWM") for path in baseline_paths],
    ]
    return (
        {
            "cell": cell,
            "sites": selected_sites,
            "observed": observed,
            "expected": {
                "parametric_direct": direct_expected[indexes].astype(np.float64),
                "parametric_lambda": lambda_expected,
            },
            "dwm_expected": dwm_expected[indexes].astype(np.float64),
            "positions": positions,
            "bias_strength": calibrator.strength(cell),
        },
        inputs,
    )


def load_panels(
    artifacts: dict[str, Path],
    baselines: dict[str, Path],
    *,
    split: str,
    study: dict,
    calibrator: FrozenBiasStrengthCalibrator,
) -> tuple[dict[str, dict], list[dict[str, str]]]:
    if set(artifacts) != set(baselines):
        raise ValueError("artifact and DWM-baseline cells must match exactly")
    split_chromosomes = set(study["chromosome_split"][split])
    panels = {}
    records = []
    for cell in sorted(artifacts):
        panel, inputs = load_panel(
            cell=cell,
            artifact=artifacts[cell],
            dwm_path=baselines[cell],
            split=split,
            split_chromosomes=split_chromosomes,
            calibrator=calibrator,
        )
        panels[cell] = panel
        records.extend(inputs)
    return panels, records


def method_metrics(
    labels: np.ndarray,
    scores: np.ndarray,
    profiles: np.ndarray,
    positions: np.ndarray,
) -> dict:
    return {
        **binary_metrics(labels, scores),
        "functional_separation": standardized_functional_separation(
            profiles,
            labels,
            positions,
        ),
    }


def validation_grid(
    panels: dict[str, dict],
    alphas: tuple[float, ...],
    *,
    dispersion: float,
    minimum_sites_per_class: int,
) -> pd.DataFrame:
    rows = []
    for cell, panel in sorted(panels.items()):
        sites = panel["sites"]
        observed = panel["observed"]
        positions = panel["positions"]
        raw_scores = geometry_score(observed, positions)
        dwm_scores, dwm_profiles = residual_score(
            observed,
            panel["dwm_expected"],
            positions,
            "deviance",
            dispersion,
        )
        expected_geometry = {
            source: geometry_score(expected, positions)
            for source, expected in panel["expected"].items()
        }
        for tf, group in sites.groupby("tf", sort=True):
            indexes = group.index.to_numpy(dtype=int)
            labels = group["chip_label"].to_numpy(dtype=int)
            role = str(group["role"].iloc[0])
            family = str(group["motif_family"].iloc[0])
            positive = int(np.sum(labels == 1))
            negative = int(np.sum(labels == 0))
            status = (
                "eligible"
                if min(positive, negative) >= minimum_sites_per_class
                else "underpowered"
            )
            raw = method_metrics(
                labels,
                raw_scores[indexes],
                observed[indexes],
                positions,
            )
            dwm = method_metrics(
                labels,
                dwm_scores[indexes],
                dwm_profiles[indexes],
                positions,
            )
            base = {
                "cell": cell,
                "tf": str(tf),
                "motif_family": family,
                "role": role,
                "status": status,
                "positive_sites": positive,
                "negative_sites": negative,
                "raw_auroc": float(raw["auroc"]),
                "raw_auprc": float(raw["auprc"]),
                "raw_functional_separation": float(raw["functional_separation"]),
                "dwm_auroc": float(dwm["auroc"]),
                "dwm_auprc": float(dwm["auprc"]),
            }
            rows.append(
                {
                    **base,
                    "method": "raw",
                    "source": "raw",
                    "alpha": 0.0,
                    **raw,
                    "auroc_gain_over_raw": 0.0,
                    "relative_auprc_gain_over_raw": 0.0,
                }
            )
            rows.append(
                {
                    **base,
                    "method": "DWM_conventional_deviance",
                    "source": "DWM",
                    "alpha": 1.0,
                    **dwm,
                    "auroc_gain_over_raw": float(dwm["auroc"] - raw["auroc"]),
                    "relative_auprc_gain_over_raw": float(
                        (dwm["auprc"] - raw["auprc"])
                        / max(float(raw["auprc"]), 1e-8)
                    ),
                }
            )
            for source in PARAMETRIC_SOURCES:
                expected = panel["expected"][source]
                bias_scores = expected_geometry[source]
                for alpha in alphas:
                    scores = raw_scores[indexes] - alpha * bias_scores[indexes]
                    profiles = observed[indexes] - alpha * expected[indexes]
                    metrics = method_metrics(labels, scores, profiles, positions)
                    rows.append(
                        {
                            **base,
                            "method": "partial_bias_subtraction",
                            "source": source,
                            "alpha": alpha,
                            **metrics,
                            "auroc_gain_over_raw": float(
                                metrics["auroc"] - raw["auroc"]
                            ),
                            "relative_auprc_gain_over_raw": float(
                                (metrics["auprc"] - raw["auprc"])
                                / max(float(raw["auprc"]), 1e-8)
                            ),
                        }
                    )
    return pd.DataFrame(rows)


def select_policy_rows(
    grid: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict, dict[str, dict]]:
    candidates = grid[
        grid["status"].eq("eligible")
        & grid["method"].eq("partial_bias_subtraction")
    ].copy()
    ranking = candidates[candidates["role"].isin(DIFFICULT_ROLES)].copy()
    if ranking.empty:
        ranking = candidates.copy()
    summary = ranking.groupby(["source", "alpha"], as_index=False).agg(
        mean_auroc_gain_over_raw=("auroc_gain_over_raw", "mean"),
        mean_relative_auprc_gain_over_raw=("relative_auprc_gain_over_raw", "mean"),
        task_sd_relative_auprc_gain_over_raw=("relative_auprc_gain_over_raw", "std"),
        task_count=("tf", "size"),
        family_count=("motif_family", "nunique"),
        cell_count=("cell", "nunique"),
    )
    ctcf = candidates[candidates["tf"].eq("CTCF")].groupby(["source", "alpha"])[
        "auroc_gain_over_raw"
    ].min()
    summary["minimum_ctcf_auroc_gain_over_raw"] = [
        float(ctcf.get((source, alpha), np.nan))
        for source, alpha in zip(summary["source"], summary["alpha"], strict=True)
    ]
    summary["passes_ctcf_gate"] = (
        summary["minimum_ctcf_auroc_gain_over_raw"] >= -0.02
    )
    summary["standard_error"] = summary[
        "task_sd_relative_auprc_gain_over_raw"
    ].fillna(0.0) / np.sqrt(summary["task_count"].clip(lower=1))
    eligible = summary[summary["passes_ctcf_gate"]].copy()
    if eligible.empty:
        raise RuntimeError("all partial-subtraction settings failed the CTCF gate")
    eligible = eligible.sort_values(
        ["mean_relative_auprc_gain_over_raw", "mean_auroc_gain_over_raw", "alpha", "source"],
        ascending=[False, False, True, True],
        kind="mergesort",
    )
    winner = eligible.iloc[0]
    global_choice = {
        "source": str(winner["source"]),
        "alpha": float(winner["alpha"]),
        "selection_scope": "global_difficult_tasks",
        "mean_auroc_gain_over_raw": float(winner["mean_auroc_gain_over_raw"]),
        "mean_relative_auprc_gain_over_raw": float(
            winner["mean_relative_auprc_gain_over_raw"]
        ),
        "minimum_ctcf_auroc_gain_over_raw": float(
            winner["minimum_ctcf_auroc_gain_over_raw"]
        ),
    }
    per_tf_rows = []
    per_tf_choices = {}
    raw_options = grid[
        grid["status"].eq("eligible") & grid["method"].eq("raw")
    ].copy()
    for tf in sorted(set(candidates["tf"])):
        task_candidates = candidates[candidates["tf"].eq(tf)]
        task_raw = raw_options[raw_options["tf"].eq(tf)]
        aggregates = task_candidates.groupby(["source", "alpha"], as_index=False).agg(
            mean_auroc_gain_over_raw=("auroc_gain_over_raw", "mean"),
            mean_relative_auprc_gain_over_raw=("relative_auprc_gain_over_raw", "mean"),
            minimum_auroc_gain_over_raw=("auroc_gain_over_raw", "min"),
            cell_count=("cell", "nunique"),
        )
        raw_row = pd.DataFrame(
            [
                {
                    "source": "raw",
                    "alpha": 0.0,
                    "mean_auroc_gain_over_raw": 0.0,
                    "mean_relative_auprc_gain_over_raw": 0.0,
                    "minimum_auroc_gain_over_raw": 0.0,
                    "cell_count": int(task_raw["cell"].nunique()),
                }
            ]
        )
        aggregates = pd.concat([aggregates, raw_row], ignore_index=True)
        aggregates["passes_cell_nonregression"] = (
            aggregates["minimum_auroc_gain_over_raw"] >= -0.02
        )
        ranked = aggregates[aggregates["passes_cell_nonregression"]].sort_values(
            ["mean_relative_auprc_gain_over_raw", "mean_auroc_gain_over_raw", "alpha", "source"],
            ascending=[False, False, True, True],
            kind="mergesort",
        )
        selected = ranked.iloc[0]
        choice = {
            "source": str(selected["source"]),
            "alpha": float(selected["alpha"]),
            "selection_scope": "tf_across_cells",
            "validation_cells": int(selected["cell_count"]),
            "mean_auroc_gain_over_raw": float(selected["mean_auroc_gain_over_raw"]),
            "mean_relative_auprc_gain_over_raw": float(
                selected["mean_relative_auprc_gain_over_raw"]
            ),
        }
        per_tf_choices[str(tf)] = choice
        aggregates.insert(0, "tf", str(tf))
        aggregates["selected"] = (
            aggregates["source"].eq(choice["source"])
            & np.isclose(aggregates["alpha"], choice["alpha"])
        )
        per_tf_rows.append(aggregates)
    summary["selected"] = (
        summary["source"].eq(global_choice["source"])
        & np.isclose(summary["alpha"], global_choice["alpha"])
    )
    return (
        summary.sort_values(["selected", "mean_relative_auprc_gain_over_raw"], ascending=[False, False]),
        pd.concat(per_tf_rows, ignore_index=True),
        global_choice,
        per_tf_choices,
    )


def choice_profile(panel: dict, indexes: np.ndarray, choice: dict) -> tuple[np.ndarray, np.ndarray]:
    observed = panel["observed"][indexes]
    if choice["source"] == "raw":
        profiles = observed
    else:
        profiles = observed - float(choice["alpha"]) * panel["expected"][
            str(choice["source"])
        ][indexes]
    return geometry_score(profiles, panel["positions"]), profiles


def evaluate_test(
    panels: dict[str, dict],
    policy: dict,
    *,
    dispersion: float,
    minimum_sites_per_class: int,
    bootstrap_iterations: int,
    aggregate_bootstrap_iterations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metrics_rows = []
    bootstrap_rows = []
    curve_rows = []
    score_rows = []
    for cell, panel in sorted(panels.items()):
        sites = panel["sites"]
        observed = panel["observed"]
        positions = panel["positions"]
        raw_scores = geometry_score(observed, positions)
        dwm_scores, dwm_profiles = residual_score(
            observed,
            panel["dwm_expected"],
            positions,
            "deviance",
            dispersion,
        )
        for tf, group in sites.groupby("tf", sort=True):
            indexes = group.index.to_numpy(dtype=int)
            labels = group["chip_label"].to_numpy(dtype=int)
            positive = int(np.sum(labels == 1))
            negative = int(np.sum(labels == 0))
            status = (
                "eligible"
                if min(positive, negative) >= minimum_sites_per_class
                else "underpowered"
            )
            global_scores, global_profiles = choice_profile(
                panel,
                indexes,
                policy["global_choice"],
            )
            tf_choice = policy["per_tf_choices"].get(str(tf), {"source": "raw", "alpha": 0.0})
            tf_scores, tf_profiles = choice_profile(panel, indexes, tf_choice)
            methods = {
                "raw": (raw_scores[indexes], observed[indexes]),
                "DWM_conventional_deviance": (dwm_scores[indexes], dwm_profiles[indexes]),
                "frozen_global_shrinkage": (global_scores, global_profiles),
                "frozen_tf_specific_shrinkage": (tf_scores, tf_profiles),
            }
            raw_metric = method_metrics(labels, *methods["raw"], positions)
            dwm_metric = method_metrics(labels, *methods["DWM_conventional_deviance"], positions)
            for method, (scores, profiles) in methods.items():
                result = method_metrics(labels, scores, profiles, positions)
                metrics_rows.append(
                    {
                        "cell": cell,
                        "tf": str(tf),
                        "motif_family": str(group["motif_family"].iloc[0]),
                        "role": str(group["role"].iloc[0]),
                        "method": method,
                        "status": status,
                        **result,
                        "raw_auroc": float(raw_metric["auroc"]),
                        "raw_auprc": float(raw_metric["auprc"]),
                        "dwm_auroc": float(dwm_metric["auroc"]),
                        "dwm_auprc": float(dwm_metric["auprc"]),
                        "auroc_gain_over_raw": float(result["auroc"] - raw_metric["auroc"]),
                        "relative_auprc_gain_over_raw": float(
                            (result["auprc"] - raw_metric["auprc"])
                            / max(float(raw_metric["auprc"]), 1e-8)
                        ),
                        "auroc_gain_over_dwm": float(result["auroc"] - dwm_metric["auroc"]),
                        "relative_auprc_gain_over_dwm": float(
                            (result["auprc"] - dwm_metric["auprc"])
                            / max(float(dwm_metric["auprc"]), 1e-8)
                        ),
                    }
                )
                curve = aggregate_curve(
                    profiles,
                    labels,
                    group["TFBS_chr"].astype(str).to_numpy(),
                    iterations=aggregate_bootstrap_iterations,
                    seed=seed,
                )
                for offset, position in enumerate(positions.astype(int)):
                    curve_rows.append(
                        {
                            "cell": cell,
                            "tf": str(tf),
                            "method": method,
                            "position": position,
                            "positive_mean": curve["positive_mean"][offset],
                            "negative_mean": curve["negative_mean"][offset],
                            "positive_minus_negative": curve["difference"][offset],
                            "lower_95": curve["lower_95"][offset],
                            "upper_95": curve["upper_95"][offset],
                        }
                    )
            frame = group[
                ["TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_strand", "chip_label"]
            ].reset_index(drop=True)
            frame.insert(0, "tf", str(tf))
            frame.insert(0, "cell", cell)
            frame["raw_score"] = methods["raw"][0]
            frame["dwm_score"] = methods["DWM_conventional_deviance"][0]
            frame["global_shrinkage_score"] = global_scores
            frame["tf_specific_shrinkage_score"] = tf_scores
            score_rows.extend(frame.to_dict("records"))
            if status == "eligible":
                task_sites = group.reset_index(drop=True)
                for method in ("frozen_global_shrinkage", "frozen_tf_specific_shrinkage"):
                    for baseline_name in ("raw", "DWM_conventional_deviance"):
                        bootstrap_rows.append(
                            {
                                "cell": cell,
                                "tf": str(tf),
                                "method": method,
                                "baseline": baseline_name,
                                **block_bootstrap_delta(
                                    task_sites,
                                    methods[method][0],
                                    methods[baseline_name][0],
                                    iterations=bootstrap_iterations,
                                    seed=seed,
                                ),
                            }
                        )
    return (
        pd.DataFrame(metrics_rows),
        pd.DataFrame(bootstrap_rows),
        pd.DataFrame(curve_rows),
        pd.DataFrame(score_rows),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("tune", "test"), required=True)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--reference-configuration", type=Path, required=True)
    parser.add_argument("--artifact", action="append", type=parse_name_path, required=True)
    parser.add_argument("--dwm-baseline", action="append", type=parse_name_path, required=True)
    parser.add_argument("--policy", type=Path)
    parser.add_argument(
        "--alphas",
        type=parse_alphas,
        default=DEFAULT_ALPHAS,
        help="Comma-separated partial-subtraction strengths used in tune mode.",
    )
    parser.add_argument("--dispersion", type=float, default=1.6144318781609228)
    parser.add_argument("--minimum-sites-per-class", type=int, default=200)
    parser.add_argument("--bootstrap-iterations", type=int, default=1000)
    parser.add_argument("--aggregate-bootstrap-iterations", type=int, default=500)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.dispersion < 0 or args.minimum_sites_per_class < 1:
        raise ValueError("dispersion must be non-negative and minimum support positive")
    configuration = load_safe_configuration(args.reference_configuration)
    study = json.loads(args.study.read_text(encoding="utf-8"))
    calibrator = load_calibrator(configuration)
    artifacts = dict(args.artifact)
    baselines = dict(args.dwm_baseline)
    if len(artifacts) != len(args.artifact) or len(baselines) != len(args.dwm_baseline):
        raise ValueError("cell inputs must not be duplicated")
    args.outdir.mkdir(parents=True, exist_ok=True)

    if args.mode == "tune":
        if args.policy is not None:
            raise ValueError("tune mode writes a policy and does not accept --policy")
        panels, profile_inputs = load_panels(
            artifacts,
            baselines,
            split="validation",
            study=study,
            calibrator=calibrator,
        )
        grid = validation_grid(
            panels,
            args.alphas,
            dispersion=args.dispersion,
            minimum_sites_per_class=args.minimum_sites_per_class,
        )
        global_summary, tf_summary, global_choice, per_tf_choices = select_policy_rows(grid)
        grid_path = args.outdir / "bias_shrinkage_validation_grid.tsv.gz"
        global_path = args.outdir / "bias_shrinkage_global_summary.tsv"
        tf_path = args.outdir / "bias_shrinkage_tf_summary.tsv"
        grid.to_csv(grid_path, sep="\t", index=False, compression="gzip")
        global_summary.to_csv(global_path, sep="\t", index=False)
        tf_summary.to_csv(tf_path, sep="\t", index=False)
        inputs = [
            input_record(args.study, "study"),
            input_record(args.reference_configuration, "safe-reference-configuration"),
            *profile_inputs,
        ]
        policy = {
            "schema": POLICY_SCHEMA,
            "selection_split": "validation",
            "selection_labels_used": True,
            "bias_coefficients_refitted": False,
            "raw_guardrail": True,
            "minimum_sites_per_class": args.minimum_sites_per_class,
            "dispersion": args.dispersion,
            "alphas": list(args.alphas),
            "bias_strengths": {
                cell: float(panel["bias_strength"]) for cell, panel in sorted(panels.items())
            },
            "global_choice": global_choice,
            "per_tf_choices": per_tf_choices,
            "inputs": inputs,
            "outputs": {
                "validation_grid": input_record(grid_path, "validation-grid"),
                "global_summary": input_record(global_path, "global-summary"),
                "tf_summary": input_record(tf_path, "tf-summary"),
            },
        }
        policy["policy_id"] = canonical_id(policy)
        policy_path = args.outdir / "bias_shrinkage_policy.freeze.json"
        write_immutable_json(policy_path, policy)
        print(json.dumps({"policy": str(policy_path), "global_choice": global_choice}, indent=2))
        return 0

    if args.policy is None:
        raise ValueError("test mode requires --policy")
    policy = validate_policy(args.policy)
    configuration_record = next(
        (
            record
            for record in policy["inputs"]
            if record.get("purpose") == "safe-reference-configuration"
        ),
        None,
    )
    if (
        configuration_record is None
        or Path(configuration_record["path"]) != args.reference_configuration
        or configuration_record["sha256"] != sha256_file(args.reference_configuration)
    ):
        raise ValueError("test reference configuration differs from the frozen policy")
    panels, profile_inputs = load_panels(
        artifacts,
        baselines,
        split="test",
        study=study,
        calibrator=calibrator,
    )
    for cell, panel in panels.items():
        expected_strength = policy["bias_strengths"].get(cell)
        if expected_strength is None or not np.isclose(
            float(expected_strength), panel["bias_strength"], rtol=0, atol=1e-12
        ):
            raise ValueError(f"{cell} bias strength differs from the frozen policy")
    input_freeze = {
        "schema": INPUT_SCHEMA,
        "policy": input_record(args.policy, "frozen-shrinkage-policy"),
        "reference_configuration": input_record(
            args.reference_configuration, "safe-reference-configuration"
        ),
        "inputs": [input_record(args.study, "study"), *profile_inputs],
        "test_labels_opened": True,
        "bias_coefficients_refitted": False,
        "shrinkage_refitted": False,
    }
    input_freeze["test_input_id"] = canonical_id(input_freeze)
    input_path = args.outdir / "bias_shrinkage_test_inputs.freeze.json"
    write_immutable_json(input_path, input_freeze)
    metrics, bootstrap, curves, scores = evaluate_test(
        panels,
        policy,
        dispersion=float(policy["dispersion"]),
        minimum_sites_per_class=int(policy["minimum_sites_per_class"]),
        bootstrap_iterations=args.bootstrap_iterations,
        aggregate_bootstrap_iterations=args.aggregate_bootstrap_iterations,
        seed=args.seed,
    )
    metrics_path = args.outdir / "bias_shrinkage_test_metrics.tsv"
    bootstrap_path = args.outdir / "bias_shrinkage_test_bootstrap.tsv"
    curves_path = args.outdir / "bias_shrinkage_test_profiles.tsv.gz"
    scores_path = args.outdir / "bias_shrinkage_test_site_scores.tsv.gz"
    metrics.to_csv(metrics_path, sep="\t", index=False)
    bootstrap.to_csv(bootstrap_path, sep="\t", index=False)
    curves.to_csv(curves_path, sep="\t", index=False, compression="gzip")
    scores.to_csv(scores_path, sep="\t", index=False, compression="gzip")
    manifest = {
        "schema": TEST_SCHEMA,
        "policy_id": policy["policy_id"],
        "models_refitted_on_test": False,
        "raw_guardrail": True,
        "test_input_freeze": input_record(input_path, "test-input-freeze"),
        "outputs": {
            "metrics": input_record(metrics_path, "test-metrics"),
            "bootstrap": input_record(bootstrap_path, "test-bootstrap"),
            "profiles": input_record(curves_path, "test-profiles"),
            "site_scores": input_record(scores_path, "test-site-scores"),
        },
    }
    manifest["test_result_id"] = canonical_id(manifest)
    manifest_path = args.outdir / "bias_shrinkage_test_manifest.json"
    write_immutable_json(manifest_path, manifest)
    eligible = metrics[
        metrics["status"].eq("eligible")
        & metrics["method"].isin(
            ["frozen_global_shrinkage", "frozen_tf_specific_shrinkage"]
        )
    ]
    print(
        eligible.groupby("method", as_index=False).agg(
            tasks=("tf", "size"),
            mean_auroc_gain_over_raw=("auroc_gain_over_raw", "mean"),
            mean_relative_auprc_gain_over_raw=("relative_auprc_gain_over_raw", "mean"),
            mean_auroc_gain_over_dwm=("auroc_gain_over_dwm", "mean"),
            mean_relative_auprc_gain_over_dwm=("relative_auprc_gain_over_dwm", "mean"),
        ).to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
