#!/usr/bin/env python3
"""Collate integrity-checked frozen parametric evidence for every development TF.

This is a research report, not a promotion decision.  It keeps raw signal as
an explicit guardrail so recovery from an overcorrected DWM baseline is not
misreported as newly recovered TF information.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_strand_label_free_models import file_sha256  # noqa: E402
from freeze_label_free_functional_models import immutable_write_json  # noqa: E402


SCHEMA = "fp-tools-frozen-parametric-tf-evidence-v1"
EXPECTED_SCHEMAS = {
    "test": "fp-tools-frozen-bias-shrinkage-test-v1",
    "replicate": "fp-tools-frozen-functional-depth-matrix-v2",
    "depth": (
        "fp-tools-frozen-functional-depth-matrix-v2",
        "fp-tools-combined-frozen-functional-depth-matrix-v2",
    ),
    "naked": "fp-tools-frozen-bias-shrinkage-naked-dna-v1",
    "ceiling": "fp-tools-frozen-functional-information-ceiling-v2",
}


def canonical_id(document: dict) -> str:
    encoded = json.dumps(document, sort_keys=True, separators=(",", ":"))
    return sha256(encoded.encode()).hexdigest()


def _resolve_record_path(value: str, manifest: Path) -> Path:
    path = Path(value)
    if path.is_file():
        return path
    relative = manifest.parent / path
    if relative.is_file():
        return relative
    raise FileNotFoundError(path)


def load_output(
    manifest: Path,
    *,
    expected_schema: str | tuple[str, ...],
    output: str,
) -> tuple[dict, pd.DataFrame, Path]:
    document = json.loads(manifest.read_text(encoding="utf-8"))
    schemas = (expected_schema,) if isinstance(expected_schema, str) else expected_schema
    if document.get("schema") not in schemas:
        raise ValueError(f"unsupported manifest schema in {manifest}")
    record = document.get("outputs", {}).get(output)
    if not isinstance(record, dict) or "path" not in record or "sha256" not in record:
        raise ValueError(f"manifest lacks checksummed {output} output: {manifest}")
    path = _resolve_record_path(str(record["path"]), manifest)
    if file_sha256(path) != str(record["sha256"]):
        raise ValueError(f"checksum mismatch for {output}: {path}")
    return document, pd.read_csv(path, sep="\t"), path


def resolve_base_study(study_path: Path, document: dict) -> Path:
    value = Path(str(document["base_study"]))
    candidates = (value, study_path.parent / value, REPOSITORY / value)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(value)


def development_tasks(study_path: Path) -> tuple[pd.DataFrame, Path]:
    study = json.loads(study_path.read_text(encoding="utf-8"))
    base_path = resolve_base_study(study_path, study)
    base = json.loads(base_path.read_text(encoding="utf-8"))
    tasks = pd.DataFrame(base["tasks"])
    tasks = tasks[
        tasks["cell"].isin(study["development_cells"])
        & tasks["split"].eq("development")
    ].copy()
    keys = ["cell", "tf"]
    if tasks.duplicated(keys).any():
        raise ValueError("development study has duplicate cell/TF tasks")
    expected = int(study["development_task_count"])
    if len(tasks) != expected:
        raise ValueError(
            f"development task count differs from freeze: {len(tasks)} != {expected}"
        )
    keep = ["cell", "tf", "motif_id", "motif_family", "role"]
    return tasks[keep].sort_values(keys).reset_index(drop=True), base_path


def method_metrics(metrics: pd.DataFrame, method: str, prefix: str) -> pd.DataFrame:
    selected = metrics[metrics["method"].eq(method)].copy()
    if selected.duplicated(["cell", "tf"]).any():
        raise ValueError(f"duplicate {method} test metric")
    columns = {
        "status": f"{prefix}_status",
        "n_sites": f"{prefix}_sites",
        "positive_sites": f"{prefix}_positive_sites",
        "negative_sites": f"{prefix}_negative_sites",
        "auroc": f"{prefix}_auroc",
        "auprc": f"{prefix}_auprc",
        "brier": f"{prefix}_brier",
        "functional_separation": f"{prefix}_functional_separation",
    }
    available = [column for column in columns if column in selected]
    return selected[["cell", "tf", *available]].rename(columns=columns)


def bootstrap_metrics(
    bootstrap: pd.DataFrame,
    *,
    method: str,
    baseline: str,
    prefix: str,
) -> pd.DataFrame:
    selected = bootstrap[
        bootstrap["method"].eq(method) & bootstrap["baseline"].eq(baseline)
    ].copy()
    if selected.duplicated(["cell", "tf"]).any():
        raise ValueError(f"duplicate {method}/{baseline} bootstrap row")
    rename = {
        column: f"{prefix}_{column}"
        for column in selected
        if column not in {"cell", "tf", "method", "baseline"}
    }
    return selected.drop(columns=["method", "baseline"]).rename(columns=rename)


def replicate_evidence(summary: pd.DataFrame, method: str) -> pd.DataFrame:
    selected = summary[summary["method"].eq(method)].copy()
    if selected.duplicated(["cell", "tf"]).any():
        raise ValueError(f"duplicate {method} replicate summary")
    keep = [
        "cell",
        "tf",
        "samples",
        "observations",
        "auroc_gain_over_dwm_mean",
        "relative_auprc_gain_over_dwm_mean",
        "auroc_gain_over_raw_mean",
        "relative_auprc_gain_over_raw_mean",
        "auroc_gain_positive_fraction",
        "auprc_gain_positive_fraction",
        "auroc_gain_over_raw_positive_fraction",
        "auprc_gain_over_raw_positive_fraction",
    ]
    selected = selected[[column for column in keep if column in selected]]
    return selected.rename(
        columns={
            column: f"replicate_{column}"
            for column in selected
            if column not in {"cell", "tf"}
        }
    )


def depth_evidence(metrics: pd.DataFrame, method: str) -> pd.DataFrame:
    selected = metrics[metrics["method"].eq(method)].copy()
    rows = []
    for (cell, tf), group in selected.groupby(["cell", "tf"], sort=True):
        both_raw = (group["auroc_gain_over_raw"] > 0) & (
            group["relative_auprc_gain_over_raw"] > 0
        )
        both_dwm = (group["auroc_gain_over_dwm"] > 0) & (
            group["relative_auprc_gain_over_dwm"] > 0
        )
        preferred_depth = next(
            (
                depth
                for depth in ("full", "50m", "25m", "10m")
                if depth in set(group["depth"].astype(str))
            ),
            None,
        )
        high = group[group["depth"].astype(str).eq(str(preferred_depth))]
        rows.append(
            {
                "cell": cell,
                "tf": tf,
                "depth_levels": ",".join(sorted(group["depth"].astype(str).unique())),
                "depth_seed_observations": int(len(group)),
                "depth_seeds": int(group["seed"].nunique()),
                "depth_both_gain_over_raw_fraction": float(both_raw.mean()),
                "depth_both_gain_over_dwm_fraction": float(both_dwm.mean()),
                "depth_all_gain_over_raw": bool(both_raw.all()),
                "depth_high_endpoint": str(preferred_depth),
                "depth_high_auroc_gain_over_raw_mean": float(
                    high["auroc_gain_over_raw"].mean()
                ),
                "depth_high_relative_auprc_gain_over_raw_mean": float(
                    high["relative_auprc_gain_over_raw"].mean()
                ),
                "depth_high_auroc_gain_over_dwm_mean": float(
                    high["auroc_gain_over_dwm"].mean()
                ),
                "depth_high_relative_auprc_gain_over_dwm_mean": float(
                    high["relative_auprc_gain_over_dwm"].mean()
                ),
            }
        )
    return pd.DataFrame(rows)


def naked_evidence(rates: pd.DataFrame, method: str) -> pd.DataFrame:
    selected = rates[rates["method"].eq(method)].copy()
    if selected.duplicated(["cell", "tf"]).any():
        raise ValueError(f"duplicate {method} naked-DNA row")
    keep = [
        "cell",
        "tf",
        "finite_sites",
        "calls",
        "false_positive_rate",
        "false_positive_rate_upper_95",
        "false_positive_rate_increase_over_dwm",
        "passes_safety",
    ]
    selected = selected[[column for column in keep if column in selected]]
    return selected.rename(
        columns={
            column: f"naked_{column}"
            for column in selected
            if column not in {"cell", "tf"}
        }
    )


def ceiling_evidence(ceiling: pd.DataFrame) -> pd.DataFrame:
    if ceiling.duplicated(["cell", "tf"]).any():
        raise ValueError("duplicate information-ceiling task")
    keep = [
        "cell",
        "tf",
        "raw_guarded_failure_classification",
        "full_relative_auprc_gain_over_raw",
        "full_relative_auprc_gain_over_signal_panel",
        "label_free_auroc_gain_over_raw",
        "label_free_relative_auprc_gain_over_raw",
    ]
    return ceiling[[column for column in keep if column in ceiling]].rename(
        columns={
            column: f"ceiling_{column}"
            for column in ceiling
            if column not in {"cell", "tf"}
        }
    )


def classify_task(row: pd.Series) -> tuple[str, str]:
    if (
        pd.isna(row.get("candidate_auroc"))
        or str(row.get("candidate_status", "underpowered")) != "eligible"
    ):
        return "underpowered", "collect_more_sites_or_labels"
    if not bool(row.get("naked_passes_safety", False)):
        return "safety_limited", "reject_candidate"

    test_raw_positive = (
        float(row.get("candidate_auroc_gain_over_raw", -np.inf)) > 0
        and float(row.get("candidate_relative_auprc_gain_over_raw", -np.inf)) > 0
    )
    replicate_raw_fraction = min(
        float(row.get("replicate_auroc_gain_over_raw_positive_fraction", 0.0)),
        float(row.get("replicate_auprc_gain_over_raw_positive_fraction", 0.0)),
    )
    depth_raw_fraction = float(row.get("depth_both_gain_over_raw_fraction", 0.0))
    raw_auroc_lower = float(row.get("raw_auroc_gain_lower_95", -np.inf))
    raw_auprc_lower = float(row.get("raw_relative_auprc_gain_lower_95", -np.inf))
    if (
        test_raw_positive
        and replicate_raw_fraction == 1.0
        and depth_raw_fraction == 1.0
    ):
        if raw_auroc_lower > 0 and raw_auprc_lower > 0:
            return "robust_gain_over_raw", "frozen_tf_specific_shrinkage"
        return (
            "replicate_stable_partial_correction_gain",
            "frozen_tf_specific_shrinkage_research_only",
        )

    raw_auroc = float(row.get("raw_auroc", np.nan))
    dwm_auroc = float(row.get("dwm_auroc", np.nan))
    candidate_raw_gain = float(row.get("candidate_auroc_gain_over_raw", np.nan))
    if (
        np.isfinite(raw_auroc)
        and np.isfinite(dwm_auroc)
        and raw_auroc - dwm_auroc >= 0.03
        and abs(candidate_raw_gain) <= 0.01
    ):
        return "dwm_overcorrection", "raw_geometry_research_baseline"

    ceiling = row.get("ceiling_raw_guarded_failure_classification")
    if isinstance(ceiling, str) and ceiling:
        if "shape_model" in ceiling or "signal_combination" in ceiling:
            return ceiling, "functional_shape_research"
        if "assay_limited" in ceiling:
            return ceiling, "no_footprint_call_recommended"
        return ceiling, "retain_raw_guardrail"
    return "shape_or_assay_limited", "retain_raw_guardrail"


def add_candidate_deltas(report: pd.DataFrame) -> pd.DataFrame:
    for baseline in ("raw", "dwm"):
        report[f"candidate_auroc_gain_over_{baseline}"] = (
            report["candidate_auroc"] - report[f"{baseline}_auroc"]
        )
        report[f"candidate_relative_auprc_gain_over_{baseline}"] = (
            report["candidate_auprc"] - report[f"{baseline}_auprc"]
        ) / report[f"{baseline}_auprc"].clip(lower=1e-8)
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--test-manifest", type=Path, required=True)
    parser.add_argument("--replicate-manifest", type=Path, required=True)
    parser.add_argument("--depth-manifest", type=Path, required=True)
    parser.add_argument("--naked-manifest", type=Path, required=True)
    parser.add_argument("--ceiling-manifest", type=Path, required=True)
    parser.add_argument(
        "--candidate-method", default="frozen_tf_specific_shrinkage"
    )
    parser.add_argument("--outdir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    tasks, base_study = development_tasks(args.study)
    test_doc, test, test_path = load_output(
        args.test_manifest,
        expected_schema=EXPECTED_SCHEMAS["test"],
        output="metrics",
    )
    _test_doc, bootstrap, bootstrap_path = load_output(
        args.test_manifest,
        expected_schema=EXPECTED_SCHEMAS["test"],
        output="bootstrap",
    )
    replicate_doc, replicates, replicate_path = load_output(
        args.replicate_manifest,
        expected_schema=EXPECTED_SCHEMAS["replicate"],
        output="replicate_summary",
    )
    depth_doc, depth, depth_path = load_output(
        args.depth_manifest,
        expected_schema=EXPECTED_SCHEMAS["depth"],
        output="metrics",
    )
    naked_doc, naked, naked_path = load_output(
        args.naked_manifest,
        expected_schema=EXPECTED_SCHEMAS["naked"],
        output="rates",
    )
    ceiling_doc, ceiling, ceiling_path = load_output(
        args.ceiling_manifest,
        expected_schema=EXPECTED_SCHEMAS["ceiling"],
        output="best_diagnostic_channel",
    )

    shrinkage_ids = {
        str(test_doc["policy_id"]),
        str(replicate_doc["bias_shrinkage_policy_id"]),
        str(depth_doc["bias_shrinkage_policy_id"]),
        str(naked_doc["policy_id"]),
    }
    if len(shrinkage_ids) != 1:
        raise ValueError("evidence manifests use different shrinkage policies")
    if replicate_doc["policy_id"] != ceiling_doc["policy_id"]:
        raise ValueError("replicate and information-ceiling policies differ")

    report = tasks.copy()
    for frame in (
        method_metrics(test, "raw", "raw"),
        method_metrics(test, "DWM_conventional_deviance", "dwm"),
        method_metrics(test, args.candidate_method, "candidate"),
        bootstrap_metrics(
            bootstrap,
            method=args.candidate_method,
            baseline="raw",
            prefix="raw",
        ),
        bootstrap_metrics(
            bootstrap,
            method=args.candidate_method,
            baseline="DWM_conventional_deviance",
            prefix="dwm",
        ),
        replicate_evidence(replicates, args.candidate_method),
        depth_evidence(depth, args.candidate_method),
        naked_evidence(naked, args.candidate_method),
        ceiling_evidence(ceiling),
    ):
        report = report.merge(frame, on=["cell", "tf"], how="left", validate="one_to_one")
    report = add_candidate_deltas(report)
    classifications = report.apply(classify_task, axis=1, result_type="expand")
    report["failure_classification"] = classifications[0]
    report["recommended_configuration"] = classifications[1]
    report["research_only"] = True
    report["shrinkage_policy_id"] = next(iter(shrinkage_ids))

    summary = (
        report.groupby(
            ["failure_classification", "recommended_configuration"],
            dropna=False,
            sort=True,
        )
        .size()
        .rename("tasks")
        .reset_index()
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    report_path = args.outdir / "frozen_parametric_per_tf_evidence.tsv"
    summary_path = args.outdir / "frozen_parametric_evidence_summary.tsv"
    report.to_csv(report_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    inputs = [
        args.study,
        base_study,
        args.test_manifest,
        test_path,
        bootstrap_path,
        args.replicate_manifest,
        replicate_path,
        args.depth_manifest,
        depth_path,
        args.naked_manifest,
        naked_path,
        args.ceiling_manifest,
        ceiling_path,
    ]
    manifest = {
        "schema": SCHEMA,
        "candidate_method": args.candidate_method,
        "shrinkage_policy_id": next(iter(shrinkage_ids)),
        "functional_policy_id": replicate_doc["policy_id"],
        "development_tasks": int(len(report)),
        "eligible_tasks": int(report["candidate_auroc"].notna().sum()),
        "raw_signal_guardrail": True,
        "promotion_decision": False,
        "inputs": [
            {"path": str(path), "sha256": file_sha256(path)} for path in inputs
        ],
        "outputs": {
            "per_tf_evidence": {
                "path": str(report_path),
                "sha256": file_sha256(report_path),
            },
            "summary": {
                "path": str(summary_path),
                "sha256": file_sha256(summary_path),
            },
        },
    }
    manifest["evidence_id"] = canonical_id(manifest)
    immutable_write_json(args.outdir / "frozen_parametric_evidence_manifest.json", manifest)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
