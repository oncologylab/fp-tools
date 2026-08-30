#!/usr/bin/env python3
"""Apply prespecified promotion gates to a frozen footprint candidate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


TASK_COLUMNS = ["cell", "tf", "motif_id"]


def prepare_pairs(
    metrics: pd.DataFrame,
    study: dict,
    candidate: str,
    baseline: str,
    split: str,
) -> pd.DataFrame:
    required = TASK_COLUMNS + ["method", "auroc", "auprc"]
    missing = [column for column in required if column not in metrics]
    if missing:
        raise ValueError(f"metrics are missing columns: {', '.join(missing)}")
    tasks = pd.DataFrame(study["tasks"])
    tasks = tasks[tasks["split"] == split]
    selected = metrics[metrics["method"].isin([baseline, candidate])]
    duplicated = selected.duplicated(TASK_COLUMNS + ["method"])
    if duplicated.any():
        raise ValueError("metrics contain duplicate task/method rows")
    wide = selected.pivot(index=TASK_COLUMNS, columns="method", values=["auroc", "auprc"])
    required_columns = [(metric, method) for metric in ("auroc", "auprc") for method in (baseline, candidate)]
    absent = [column for column in required_columns if column not in wide]
    if absent:
        raise ValueError("candidate or baseline metrics are absent")
    wide = wide.dropna(subset=required_columns).reset_index()
    wide.columns = [
        "_".join(str(item) for item in column if item != "") if isinstance(column, tuple) else column
        for column in wide.columns
    ]
    pairs = tasks.merge(wide, on=TASK_COLUMNS, how="left", validate="one_to_one")
    metric_columns = [f"{metric}_{method}" for metric in ("auroc", "auprc") for method in (baseline, candidate)]
    pairs = pairs.dropna(subset=metric_columns).copy()
    if pairs.empty:
        raise ValueError(f"no complete {split} candidate/baseline task pairs are available")
    pairs["delta_auroc"] = pairs[f"auroc_{candidate}"] - pairs[f"auroc_{baseline}"]
    denominator = pairs[f"auprc_{baseline}"].abs().clip(lower=np.finfo(float).eps)
    pairs["relative_delta_auprc"] = (
        pairs[f"auprc_{candidate}"] - pairs[f"auprc_{baseline}"]
    ) / denominator
    return pairs


def clustered_pass_probability(
    pairs: pd.DataFrame,
    minimum_auroc_gain: float,
    minimum_relative_auprc_gain: float,
    n_bootstrap: int,
    seed: int,
) -> float:
    families = [group.index.to_numpy() for _, group in pairs.groupby("motif_family", sort=True)]
    rng = np.random.default_rng(seed)
    passed = 0
    for _ in range(n_bootstrap):
        selected = rng.integers(0, len(families), size=len(families))
        indexes = np.concatenate([families[index] for index in selected])
        sample = pairs.loc[indexes]
        passed += int(
            sample["delta_auroc"].median() >= minimum_auroc_gain
            and sample["relative_delta_auprc"].median() >= minimum_relative_auprc_gain
        )
    return passed / n_bootstrap


def evaluate_gates(
    pairs: pd.DataFrame,
    gates: dict,
    n_bootstrap: int = 5000,
    seed: int = 2026,
    negative_controls: pd.DataFrame | None = None,
    candidate: str = "candidate",
    baseline: str = "baseline",
) -> tuple[pd.DataFrame, dict[str, object]]:
    auroc_threshold = float(gates["minimum_auroc_gain"])
    auprc_threshold = float(gates["minimum_relative_auprc_gain"])
    pairs = pairs.copy()
    pairs["improved"] = (
        (pairs["delta_auroc"] >= auroc_threshold)
        & (pairs["relative_delta_auprc"] >= auprc_threshold)
    )
    controls = pairs[pairs["role"] == "positive_control"]
    median_auroc = float(pairs["delta_auroc"].median())
    median_auprc = float(pairs["relative_delta_auprc"].median())
    worst_control = float(controls["delta_auroc"].min()) if len(controls) else np.nan
    improved_families = int(pairs.loc[pairs["improved"], "motif_family"].nunique())
    improved_cells = int(pairs.loc[pairs["improved"], "cell"].nunique())
    probability = clustered_pass_probability(
        pairs,
        auroc_threshold,
        auprc_threshold,
        n_bootstrap,
        seed,
    )
    naked_candidate_fpr = np.nan
    naked_fpr_increase = np.nan
    naked_control_pass = True
    if "maximum_naked_dna_false_positive_rate" in gates:
        naked_control_pass = False
        if negative_controls is not None:
            required = {"method", "false_positive_rate"}
            missing = required.difference(negative_controls.columns)
            if missing:
                raise ValueError(
                    "negative controls are missing columns: " + ", ".join(sorted(missing))
                )
            candidate_rows = negative_controls[negative_controls["method"] == candidate]
            baseline_rows = negative_controls[negative_controls["method"] == baseline]
            if len(candidate_rows) and len(baseline_rows):
                naked_candidate_fpr = float(candidate_rows["false_positive_rate"].median())
                naked_baseline_fpr = float(baseline_rows["false_positive_rate"].median())
                naked_fpr_increase = naked_candidate_fpr - naked_baseline_fpr
                naked_control_pass = bool(
                    naked_candidate_fpr <= float(gates["maximum_naked_dna_false_positive_rate"])
                    and naked_fpr_increase
                    <= float(gates["maximum_naked_dna_false_positive_rate_increase"])
                )
    checks = {
        "median_auroc_gain": median_auroc >= auroc_threshold,
        "median_relative_auprc_gain": median_auprc >= auprc_threshold,
        "strong_control_non_regression": bool(
            len(controls) and worst_control >= -float(gates["maximum_strong_control_auroc_loss"])
        ),
        "difficult_tf_families_improved": improved_families >= int(gates["minimum_difficult_tf_families_improved"]),
        "holdout_cells_improved": improved_cells >= int(gates["minimum_holdout_cells_improved"]),
        "bootstrap_detectability": probability >= float(gates["minimum_detectability_probability"]),
        "naked_dna_false_positive_control": naked_control_pass,
    }
    summary: dict[str, object] = {
        "passed": all(checks.values()),
        "task_count": int(len(pairs)),
        "median_auroc_gain": median_auroc,
        "median_relative_auprc_gain": median_auprc,
        "worst_positive_control_auroc_gain": worst_control,
        "improved_tf_families": improved_families,
        "improved_cells": improved_cells,
        "bootstrap_pass_probability": probability,
        "bootstrap_iterations": int(n_bootstrap),
        "naked_dna_candidate_false_positive_rate": naked_candidate_fpr,
        "naked_dna_false_positive_rate_increase": naked_fpr_increase,
        "seed": int(seed),
        "checks": checks,
    }
    return pairs, summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--study", type=Path, default=Path("benchmarks/manifests/footprint_detectability_v1.spec.json"))
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--split", choices=["development", "locked_holdout"], default="development")
    parser.add_argument("--unlock-holdout", action="store_true")
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--negative-controls",
        type=Path,
        help="TSV with method and naked-DNA false_positive_rate columns.",
    )
    args = parser.parse_args(argv)
    if args.split == "locked_holdout" and not args.unlock_holdout:
        raise SystemExit("locked_holdout evaluation requires --unlock-holdout after the candidate is frozen")

    study = json.loads(args.study.read_text(encoding="utf-8"))
    pairs = prepare_pairs(
        pd.read_csv(args.metrics, sep="\t"),
        study,
        args.candidate,
        args.baseline,
        args.split,
    )
    pairs, summary = evaluate_gates(
        pairs,
        study["promotion_gates"],
        n_bootstrap=args.bootstrap,
        seed=args.seed,
        negative_controls=(
            pd.read_csv(args.negative_controls, sep="\t")
            if args.negative_controls
            else None
        ),
        candidate=args.candidate,
        baseline=args.baseline,
    )
    summary.update({"candidate": args.candidate, "baseline": args.baseline, "split": args.split})
    args.outdir.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.outdir / "promotion_task_differences.tsv", sep="\t", index=False)
    (args.outdir / "promotion_gate.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
