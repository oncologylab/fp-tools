#!/usr/bin/env python3
"""Apply every locked functional-footprinting promotion gate, failing closed."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_optional(path: Path | None) -> pd.DataFrame | None:
    return pd.read_csv(path, sep="\t") if path is not None else None


def _method_column(metrics: pd.DataFrame, candidate: str, baseline: str) -> str:
    if "candidate_id" in metrics and {candidate, baseline}.issubset(set(metrics["candidate_id"].astype(str))):
        return "candidate_id"
    return "method"


def prepare_pairs(
    metrics: pd.DataFrame,
    study: dict,
    candidate: str,
    baseline: str,
    task_split: str,
) -> pd.DataFrame:
    method_column = _method_column(metrics, candidate, baseline)
    required = {"cell", "tf", method_column, "auroc", "auprc"}
    missing = required.difference(metrics.columns)
    if missing:
        raise ValueError("metrics are missing columns: " + ", ".join(sorted(missing)))
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[tasks["split"] == task_split].copy()
    selected = metrics[metrics[method_column].astype(str).isin([candidate, baseline])].copy()
    if "split" in selected:
        expected_metric_split = "validation" if task_split == "development" else "locked_holdout"
        if expected_metric_split in set(selected["split"].astype(str)):
            selected = selected[selected["split"].astype(str) == expected_metric_split]
    keys = ["cell", "tf"]
    if "motif_id" in selected and "motif_id" in tasks:
        keys.append("motif_id")
    expected_tasks = tasks[keys].drop_duplicates()
    if len(expected_tasks) != len(tasks):
        raise ValueError("study contains duplicate promotion task keys")
    if selected.duplicated(keys + [method_column]).any():
        raise ValueError("metrics contain duplicate task/candidate rows")
    wide = selected.pivot(index=keys, columns=method_column, values=["auroc", "auprc"])
    required_columns = [(metric, method) for metric in ("auroc", "auprc") for method in (baseline, candidate)]
    if any(column not in wide for column in required_columns):
        raise ValueError("candidate or baseline metrics are absent")
    wide = wide.dropna(subset=required_columns).reset_index()
    wide.columns = [
        "_".join(str(value) for value in column if value != "")
        if isinstance(column, tuple)
        else column
        for column in wide.columns
    ]
    observed_tasks = wide[keys].drop_duplicates()
    coverage = expected_tasks.merge(
        observed_tasks,
        on=keys,
        how="outer",
        indicator=True,
    )
    if not coverage["_merge"].eq("both").all():
        missing = int(coverage["_merge"].eq("left_only").sum())
        unexpected = int(coverage["_merge"].eq("right_only").sum())
        raise ValueError(
            "promotion metrics do not exactly cover the frozen task set "
            f"(missing={missing}, unexpected={unexpected})"
        )
    pairs = tasks.merge(wide, on=keys, how="inner", validate="one_to_one")
    if pairs.empty:
        raise ValueError("no complete candidate/baseline task pairs are available")
    pairs["delta_auroc"] = pairs[f"auroc_{candidate}"] - pairs[f"auroc_{baseline}"]
    denominator = pairs[f"auprc_{baseline}"].abs().clip(lower=np.finfo(float).eps)
    pairs["relative_delta_auprc"] = (
        pairs[f"auprc_{candidate}"] - pairs[f"auprc_{baseline}"]
    ) / denominator
    return pairs


def clustered_bootstrap_intervals(
    pairs: pd.DataFrame,
    *,
    iterations: int,
    seed: int,
) -> dict[str, float]:
    groups = [group.index.to_numpy() for _, group in pairs.groupby("motif_family", sort=True)]
    rng = np.random.default_rng(seed)
    auroc = np.empty(iterations, dtype=float)
    auprc = np.empty(iterations, dtype=float)
    for index in range(iterations):
        selected = rng.integers(0, len(groups), size=len(groups))
        rows = np.concatenate([groups[value] for value in selected])
        auroc[index] = pairs.loc[rows, "delta_auroc"].mean()
        auprc[index] = pairs.loc[rows, "relative_delta_auprc"].mean()
    return {
        "mean_auroc_lower_95": float(np.quantile(auroc, 0.025)),
        "mean_auroc_upper_95": float(np.quantile(auroc, 0.975)),
        "relative_auprc_lower_95": float(np.quantile(auprc, 0.025)),
        "relative_auprc_upper_95": float(np.quantile(auprc, 0.975)),
        "probability_both_positive": float(np.mean((auroc > 0) & (auprc > 0))),
    }


def functional_separation_evidence(
    descriptors: pd.DataFrame | None,
    pairs: pd.DataFrame,
    candidate_correction: str,
    baseline_correction: str,
    minimum_gain: float,
) -> tuple[pd.DataFrame, int]:
    if descriptors is None or descriptors.empty:
        return pd.DataFrame(), 0
    required = {"cell", "tf", "motif_family", "correction", "group", "depletion"}
    if not required.issubset(descriptors.columns):
        raise ValueError("profile descriptors lack functional-separation columns")
    selected = descriptors[
        descriptors["correction"].astype(str).isin([candidate_correction, baseline_correction])
    ]
    pivot = selected.pivot_table(
        index=["cell", "tf", "motif_family", "correction"],
        columns="group",
        values="depletion",
        aggfunc="first",
    ).reset_index()
    if not {"chip_positive", "matched_negative"}.issubset(pivot.columns):
        return pd.DataFrame(), 0
    pivot["functional_separation"] = np.abs(
        pivot["chip_positive"] - pivot["matched_negative"]
    )
    wide = pivot.pivot_table(
        index=["cell", "tf", "motif_family"],
        columns="correction",
        values="functional_separation",
        aggfunc="first",
    ).reset_index()
    if candidate_correction not in wide or baseline_correction not in wide:
        return pd.DataFrame(), 0
    wide["relative_functional_separation_gain"] = (
        wide[candidate_correction] - wide[baseline_correction]
    ) / wide[baseline_correction].abs().clip(lower=np.finfo(float).eps)
    roles = pairs[["cell", "tf", "role"]].drop_duplicates()
    wide = wide.merge(roles, on=["cell", "tf"], how="inner")
    difficult = wide[wide["role"] == "difficult"]
    improved_families = int(
        difficult.loc[
            difficult["relative_functional_separation_gain"] >= minimum_gain,
            "motif_family",
        ].nunique()
    )
    return wide, improved_families


def negative_control_evidence(
    frame: pd.DataFrame | None,
    candidate: str,
    baseline: str,
    gates: dict,
) -> tuple[bool, dict[str, float]]:
    evidence = {"candidate_fpr": np.nan, "baseline_fpr": np.nan, "fpr_increase": np.nan}
    if frame is None or frame.empty:
        return False, evidence
    method_column = "candidate_id" if "candidate_id" in frame else "method"
    if not {method_column, "false_positive_rate"}.issubset(frame.columns):
        raise ValueError("negative controls require method/candidate_id and false_positive_rate")
    candidate_rows = frame[frame[method_column].astype(str) == candidate]
    baseline_rows = frame[frame[method_column].astype(str) == baseline]
    if candidate_rows.empty or baseline_rows.empty:
        return False, evidence
    candidate_fpr = float(candidate_rows["false_positive_rate"].mean())
    baseline_fpr = float(baseline_rows["false_positive_rate"].mean())
    increase = candidate_fpr - baseline_fpr
    evidence = {"candidate_fpr": candidate_fpr, "baseline_fpr": baseline_fpr, "fpr_increase": increase}
    return bool(
        candidate_fpr <= float(gates["maximum_naked_dna_false_positive_rate"])
        and increase <= float(gates["maximum_naked_dna_false_positive_rate_increase"])
    ), evidence


def resource_evidence(
    frame: pd.DataFrame | None,
    candidate: str,
    baseline: str,
    gates: dict,
) -> tuple[bool, dict[str, float]]:
    evidence = {"runtime_ratio": np.nan, "memory_ratio": np.nan, "model_size_mb": np.nan}
    if frame is None or frame.empty:
        return False, evidence
    method_column = "candidate_id" if "candidate_id" in frame else "method"
    required = {method_column, "runtime_seconds", "peak_memory_mb", "model_size_mb"}
    if not required.issubset(frame.columns):
        raise ValueError("resource metrics are missing required columns")
    candidate_rows = frame[frame[method_column].astype(str) == candidate]
    baseline_rows = frame[frame[method_column].astype(str) == baseline]
    if candidate_rows.empty or baseline_rows.empty:
        return False, evidence
    runtime_ratio = float(candidate_rows["runtime_seconds"].mean() / max(baseline_rows["runtime_seconds"].mean(), 1e-12))
    memory_ratio = float(candidate_rows["peak_memory_mb"].max() / max(baseline_rows["peak_memory_mb"].max(), 1e-12))
    model_size = float(candidate_rows["model_size_mb"].max())
    evidence = {"runtime_ratio": runtime_ratio, "memory_ratio": memory_ratio, "model_size_mb": model_size}
    return bool(
        runtime_ratio <= float(gates["maximum_runtime_ratio_to_dwm"])
        and memory_ratio <= float(gates["maximum_memory_ratio_to_dwm"])
        and model_size <= float(gates["maximum_model_size_mb"])
    ), evidence


def simple_candidate_check(
    frame: pd.DataFrame | None,
    candidate: str,
    value_column: str,
    predicate,
) -> tuple[bool, list[float]]:
    if frame is None or frame.empty:
        return False, []
    method_column = "candidate_id" if "candidate_id" in frame else "method"
    if method_column not in frame or value_column not in frame:
        raise ValueError(f"evidence requires {method_column} and {value_column}")
    values = frame.loc[frame[method_column].astype(str) == candidate, value_column].tolist()
    return bool(values and all(predicate(value) for value in values)), values


def complexity_evidence(
    frame: pd.DataFrame | None,
    candidate: str,
    gates: dict,
) -> tuple[bool, dict[str, object]]:
    method = candidate.split(":", 1)[-1].lower()
    requires_gp_gate = method in {"gp", "hybrid", "functional-gp", "functional_gp"}
    evidence: dict[str, object] = {
        "required": requires_gp_gate,
        "mean_relative_auprc_gain_over_spline": np.nan,
        "uncertainty_calibration_improved": False,
    }
    if not requires_gp_gate:
        return True, evidence
    if frame is None or frame.empty:
        return False, evidence
    method_column = "candidate_id" if "candidate_id" in frame else "method"
    required = {
        method_column,
        "relative_auprc_gain_over_spline",
        "uncertainty_calibration_improved",
    }
    if not required.issubset(frame.columns):
        raise ValueError("GP/spline complexity evidence is missing required columns")
    selected = frame[frame[method_column].astype(str) == candidate]
    if selected.empty:
        return False, evidence
    mean_gain = float(selected["relative_auprc_gain_over_spline"].mean())
    calibration = bool(selected["uncertainty_calibration_improved"].astype(bool).all())
    evidence.update(
        {
            "mean_relative_auprc_gain_over_spline": mean_gain,
            "uncertainty_calibration_improved": calibration,
        }
    )
    return bool(
        mean_gain >= float(gates["minimum_gp_relative_auprc_gain_over_spline"])
        or calibration
    ), evidence


def evaluate_promotion(
    pairs: pd.DataFrame,
    study: dict,
    *,
    candidate: str,
    baseline: str,
    descriptors: pd.DataFrame | None = None,
    negative_controls: pd.DataFrame | None = None,
    resources: pd.DataFrame | None = None,
    uncertainty: pd.DataFrame | None = None,
    stability: pd.DataFrame | None = None,
    leakage: pd.DataFrame | None = None,
    complexity: pd.DataFrame | None = None,
    locked_holdout_scored: bool = False,
    bootstrap: int = 5000,
    seed: int = 2026,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    gates = study["promotion_gates"]
    pairs = pairs.copy()
    pairs["improved"] = (
        (pairs["delta_auroc"] >= float(gates["minimum_mean_auroc_gain"]))
        & (pairs["relative_delta_auprc"] >= float(gates["minimum_relative_auprc_gain"]))
    )
    mean_auroc = float(pairs["delta_auroc"].mean())
    mean_auprc = float(pairs["relative_delta_auprc"].mean())
    controls = pairs[pairs["role"] == "positive_control"]
    worst_control = float(controls["delta_auroc"].min()) if len(controls) else np.nan
    improved_families = int(pairs.loc[pairs["improved"], "motif_family"].nunique())
    improved_cells = int(pairs.loc[pairs["improved"], "cell"].nunique())
    bootstrap_evidence = clustered_bootstrap_intervals(pairs, iterations=bootstrap, seed=seed)

    candidate_correction = candidate.split(":", 1)[0]
    baseline_correction = baseline.split(":", 1)[0]
    separation, separation_families = functional_separation_evidence(
        descriptors,
        pairs,
        candidate_correction,
        baseline_correction,
        float(gates["minimum_functional_separation_gain"]),
    )
    negative_pass, negative_evidence = negative_control_evidence(
        negative_controls, candidate, baseline, gates
    )
    resource_pass, resource_values = resource_evidence(resources, candidate, baseline, gates)
    uncertainty_pass, uncertainty_values = simple_candidate_check(
        uncertainty,
        candidate,
        "empirical_coverage",
        lambda value: float(gates["minimum_uncertainty_coverage"])
        <= float(value)
        <= float(gates["maximum_uncertainty_coverage"]),
    )
    stability_pass, stability_values = simple_candidate_check(
        stability,
        candidate,
        "direction_consistent",
        lambda value: bool(value),
    )
    leakage_pass, leakage_values = simple_candidate_check(
        leakage,
        candidate,
        "potential_motif_response_requires_review",
        lambda value: not bool(value),
    )
    complexity_pass, complexity_values = complexity_evidence(
        complexity, candidate, gates
    )

    checks = {
        "mean_auroc_gain": mean_auroc >= float(gates["minimum_mean_auroc_gain"]),
        "mean_relative_auprc_gain": mean_auprc >= float(gates["minimum_relative_auprc_gain"]),
        "functional_separation": separation_families >= 2,
        "difficult_tf_families_improved": improved_families >= int(gates["minimum_difficult_tf_families_improved"]),
        "cell_contexts_improved": improved_cells >= int(gates["minimum_holdout_cells_improved"]),
        "positive_control_non_regression": bool(
            len(controls)
            and worst_control >= -float(gates["maximum_positive_control_auroc_loss"])
        ),
        "principal_bootstrap_ci_excludes_zero": bool(
            bootstrap_evidence["mean_auroc_lower_95"] > 0
            and bootstrap_evidence["relative_auprc_lower_95"] > 0
            and bootstrap_evidence["probability_both_positive"]
            >= float(gates["minimum_bootstrap_probability_positive"])
        ),
        "naked_dna_false_positive_control": negative_pass,
        "no_detectable_tf_motif_bias_response": leakage_pass,
        "replicate_seed_depth_stability": stability_pass,
        "cpu_memory_model_size": resource_pass,
        "functional_uncertainty_coverage": uncertainty_pass,
        "gp_or_hybrid_justifies_complexity": complexity_pass,
        "locked_external_validation_scored": bool(locked_holdout_scored),
    }
    summary = {
        "passed": bool(all(checks.values())),
        "candidate": candidate,
        "baseline": baseline,
        "task_count": int(len(pairs)),
        "mean_auroc_gain": mean_auroc,
        "mean_relative_auprc_gain": mean_auprc,
        "worst_positive_control_auroc_gain": worst_control,
        "improved_tf_families": improved_families,
        "improved_cells": improved_cells,
        "functional_separation_families": separation_families,
        "bootstrap": bootstrap_evidence,
        "negative_control": negative_evidence,
        "resources": resource_values,
        "uncertainty_coverages": uncertainty_values,
        "stability_values": stability_values,
        "bias_motif_review_flags": leakage_values,
        "complexity": complexity_values,
        "checks": checks,
        "missing_evidence_fails_closed": True,
    }
    return pairs, separation, summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--task-split", choices=("development", "locked_holdout"), default="development")
    parser.add_argument("--unlock-holdout", action="store_true")
    parser.add_argument("--descriptors", type=Path)
    parser.add_argument("--negative-controls", type=Path)
    parser.add_argument("--resources", type=Path)
    parser.add_argument("--uncertainty", type=Path)
    parser.add_argument("--stability", type=Path)
    parser.add_argument("--leakage", type=Path)
    parser.add_argument("--complexity", type=Path)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.task_split == "locked_holdout" and not args.unlock_holdout:
        raise SystemExit("locked holdout evaluation requires --unlock-holdout after candidate freeze")
    study = json.loads(args.study.read_text(encoding="utf-8"))
    pairs = prepare_pairs(
        pd.read_csv(args.metrics, sep="\t"),
        study,
        args.candidate,
        args.baseline,
        args.task_split,
    )
    pairs, separation, summary = evaluate_promotion(
        pairs,
        study,
        candidate=args.candidate,
        baseline=args.baseline,
        descriptors=read_optional(args.descriptors),
        negative_controls=read_optional(args.negative_controls),
        resources=read_optional(args.resources),
        uncertainty=read_optional(args.uncertainty),
        stability=read_optional(args.stability),
        leakage=read_optional(args.leakage),
        complexity=read_optional(args.complexity),
        locked_holdout_scored=args.task_split == "locked_holdout",
        bootstrap=args.bootstrap,
        seed=args.seed,
    )
    summary.update(
        {
            "task_split": args.task_split,
            "study": str(args.study),
            "study_sha256": file_sha256(args.study),
            "metrics": str(args.metrics),
            "metrics_sha256": file_sha256(args.metrics),
        }
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.outdir / "functional_promotion_task_differences.tsv", sep="\t", index=False)
    separation.to_csv(args.outdir / "functional_promotion_separation.tsv", sep="\t", index=False)
    (args.outdir / "functional_promotion_gate.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
