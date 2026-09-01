#!/usr/bin/env python3
"""Collate frozen TF-detector evidence with raw and DWM guardrails.

The detector policy is selected on validation chromosomes.  This report joins
its locked test result, chromosome-block bootstrap, biological-replicate
stability, depth/seed stability, and independent naked-DNA safety.  It is a
research classification and never constitutes package promotion.
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
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_strand_label_free_models import file_sha256  # noqa: E402
from freeze_label_free_functional_models import immutable_write_json  # noqa: E402
from summarize_frozen_parametric_tf_evidence import (  # noqa: E402
    development_tasks,
    load_output,
)


SCHEMA = "fp-tools-frozen-detector-evidence-v1"
TEST_SCHEMA = "fp-tools-frozen-functional-consensus-test-v1"
NAKED_SCHEMA = "fp-tools-frozen-functional-consensus-naked-dna-v1"
DEPTH_SCHEMAS = (
    "fp-tools-frozen-functional-depth-matrix-v2",
    "fp-tools-combined-frozen-functional-depth-matrix-v2",
)


def canonical_id(document: dict) -> str:
    encoded = json.dumps(document, sort_keys=True, separators=(",", ":"))
    return sha256(encoded.encode()).hexdigest()


def checked_reference(record: dict, owner: Path) -> tuple[dict, Path]:
    if not isinstance(record, dict) or "path" not in record or "sha256" not in record:
        raise ValueError(f"invalid checksummed reference in {owner}")
    path = Path(str(record["path"]))
    if not path.is_file():
        path = owner.parent / path
    if not path.is_file():
        raise FileNotFoundError(record["path"])
    if file_sha256(path) != str(record["sha256"]):
        raise ValueError(f"checksum mismatch: {path}")
    return json.loads(path.read_text(encoding="utf-8")), path


def candidate_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Return the single frozen detector selected for each observation."""

    selected = frame[frame["method"].astype(str).str.startswith("frozen_")].copy()
    identity = [
        column
        for column in ("cell", "sample", "tf", "depth", "seed")
        if column in selected
    ]
    if selected.duplicated(identity).any():
        duplicate = selected.loc[selected.duplicated(identity, keep=False), identity]
        raise ValueError(f"multiple frozen candidates for an observation:\n{duplicate}")
    return selected


def test_evidence(metrics: pd.DataFrame, operator: str) -> pd.DataFrame:
    selected = metrics[metrics["operator"].eq(operator)].copy()
    if selected.duplicated(["cell", "tf"]).any():
        raise ValueError(f"duplicate {operator} test result")
    keep = [
        "cell",
        "tf",
        "status",
        "n_sites",
        "n_positive",
        "n_negative",
        "auroc",
        "auprc",
        "dwm_auroc",
        "dwm_auprc",
        "raw_auroc",
        "raw_auprc",
        "auroc_gain_over_dwm",
        "relative_auprc_gain_over_dwm",
        "auroc_gain_over_raw",
        "relative_auprc_gain_over_raw",
        "score_standardized_separation",
    ]
    return selected[[column for column in keep if column in selected]].rename(
        columns={
            column: f"test_{column}"
            for column in selected
            if column not in {"cell", "tf", "operator"}
        }
    )


def bootstrap_evidence(bootstrap: pd.DataFrame, operator: str) -> pd.DataFrame:
    selected = bootstrap[bootstrap["operator"].eq(operator)].copy()
    pieces = []
    for baseline, prefix in (
        ("raw", "test_raw_bootstrap"),
        ("DWM_conventional_deviance", "test_dwm_bootstrap"),
    ):
        current = selected[selected["baseline"].eq(baseline)].copy()
        if current.duplicated(["cell", "tf"]).any():
            raise ValueError(f"duplicate {operator}/{baseline} bootstrap result")
        current = current.drop(columns=["operator", "baseline", "motif_family"])
        current = current.rename(
            columns={
                column: f"{prefix}_{column}"
                for column in current
                if column not in {"cell", "tf"}
            }
        )
        pieces.append(current)
    return pieces[0].merge(pieces[1], on=["cell", "tf"], how="outer", validate="one_to_one")


def replicate_evidence(summary: pd.DataFrame) -> pd.DataFrame:
    selected = candidate_rows(summary)
    if selected.duplicated(["cell", "tf"]).any():
        raise ValueError("duplicate frozen replicate summary")
    keep = [
        "cell",
        "tf",
        "candidate_id",
        "method",
        "observations",
        "samples",
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


def _depth_rank(value: str) -> float:
    text = str(value).lower()
    if text == "full":
        return np.inf
    if text.endswith("m"):
        return float(text[:-1])
    raise ValueError(f"unsupported depth: {value}")


def depth_evidence(metrics: pd.DataFrame) -> pd.DataFrame:
    selected = candidate_rows(metrics)
    rows = []
    for (cell, tf), group in selected.groupby(["cell", "tf"], sort=True):
        raw_both = (group["auroc_gain_over_raw"] > 0) & (
            group["relative_auprc_gain_over_raw"] > 0
        )
        dwm_both = (group["auroc_gain_over_dwm"] > 0) & (
            group["relative_auprc_gain_over_dwm"] > 0
        )
        ranks = group["depth"].map(_depth_rank)
        high = ranks >= 25.0
        low = ranks < 25.0
        rows.append(
            {
                "cell": cell,
                "tf": tf,
                "depth_levels": ",".join(
                    sorted(group["depth"].astype(str).unique(), key=_depth_rank)
                ),
                "depth_observations": int(len(group)),
                "depth_seeds": int(group["seed"].nunique()),
                "depth_both_gain_over_raw_fraction": float(raw_both.mean()),
                "depth_both_gain_over_dwm_fraction": float(dwm_both.mean()),
                "depth_high_both_gain_over_raw_fraction": (
                    float(raw_both[high].mean()) if high.any() else np.nan
                ),
                "depth_low_both_gain_over_raw_fraction": (
                    float(raw_both[low].mean()) if low.any() else np.nan
                ),
                "depth_min_auroc_gain_over_raw": float(
                    group["auroc_gain_over_raw"].min()
                ),
                "depth_min_relative_auprc_gain_over_raw": float(
                    group["relative_auprc_gain_over_raw"].min()
                ),
                "depth_mean_auroc_gain_over_raw": float(
                    group["auroc_gain_over_raw"].mean()
                ),
                "depth_mean_relative_auprc_gain_over_raw": float(
                    group["relative_auprc_gain_over_raw"].mean()
                ),
            }
        )
    return pd.DataFrame(rows)


def naked_evidence(rates: pd.DataFrame, operator: str) -> pd.DataFrame:
    selected = rates[
        rates["operator"].eq(operator) & rates["selected_for_primary"].astype(bool)
    ].copy()
    if selected.duplicated(["cell", "tf"]).any():
        raise ValueError(f"duplicate {operator} naked-DNA result")
    keep = [
        "cell",
        "tf",
        "valid_sites",
        "informative_sites",
        "calls",
        "false_positive_rate",
        "wilson_upper_95",
        "candidate_minus_dwm",
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


def classify_task(row: pd.Series) -> tuple[str, str]:
    if str(row.get("test_status", "underpowered")) != "eligible":
        return "underpowered", "collect_more_sites_or_labels"
    if pd.isna(row.get("naked_passes_safety")):
        return "safety_not_evaluable", "do_not_promote"
    if not bool(row.get("naked_passes_safety")):
        return "safety_limited", "reject_detector"

    test_raw = (
        float(row.get("test_auroc_gain_over_raw", -np.inf)) > 0
        and float(row.get("test_relative_auprc_gain_over_raw", -np.inf)) > 0
    )
    raw_ci = (
        float(row.get("test_raw_bootstrap_auroc_gain_lower_95", -np.inf)) > 0
        and float(
            row.get("test_raw_bootstrap_relative_auprc_gain_lower_95", -np.inf)
        )
        > 0
    )
    replicate_fraction = min(
        float(row.get("replicate_auroc_gain_over_raw_positive_fraction", 0.0)),
        float(row.get("replicate_auprc_gain_over_raw_positive_fraction", 0.0)),
    )
    depth_all = float(row.get("depth_both_gain_over_raw_fraction", 0.0))
    depth_high = float(row.get("depth_high_both_gain_over_raw_fraction", 0.0))
    if test_raw and raw_ci and replicate_fraction == 1.0 and depth_all == 1.0:
        return "robust_tf_specific_gain", "tf_specific_count_model_research_only"
    if test_raw and raw_ci and replicate_fraction == 1.0 and depth_high == 1.0:
        return (
            "depth_dependent_tf_specific_gain",
            "tf_specific_count_model_research_only",
        )
    if test_raw and raw_ci:
        return "support_or_depth_sensitive_gain", "retain_raw_guardrail"

    test_dwm = (
        float(row.get("test_auroc_gain_over_dwm", -np.inf)) > 0
        and float(row.get("test_relative_auprc_gain_over_dwm", -np.inf)) > 0
    )
    replicate_dwm = min(
        float(row.get("replicate_auroc_gain_positive_fraction", 0.0)),
        float(row.get("replicate_auprc_gain_positive_fraction", 0.0)),
    )
    depth_dwm = float(row.get("depth_both_gain_over_dwm_fraction", 0.0))
    if test_dwm and replicate_dwm == 1.0 and depth_dwm == 1.0:
        return "dwm_baseline_recovery_only", "retain_raw_guardrail"
    return "shape_or_assay_limited", "retain_raw_or_dwm_baseline"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--test-manifest", type=Path, required=True)
    parser.add_argument("--replicate-manifest", type=Path, required=True)
    parser.add_argument("--depth-manifest", type=Path, required=True)
    parser.add_argument("--naked-manifest", type=Path, required=True)
    parser.add_argument("--operator", default="count_only")
    parser.add_argument("--outdir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    tasks, base_study = development_tasks(args.study)
    test_doc, test, test_path = load_output(
        args.test_manifest, expected_schema=TEST_SCHEMA, output="metrics"
    )
    _test_doc, bootstrap, bootstrap_path = load_output(
        args.test_manifest, expected_schema=TEST_SCHEMA, output="bootstrap"
    )
    replicate_doc, replicates, replicate_path = load_output(
        args.replicate_manifest,
        expected_schema=DEPTH_SCHEMAS,
        output="replicate_summary",
    )
    depth_doc, depth, depth_path = load_output(
        args.depth_manifest, expected_schema=DEPTH_SCHEMAS, output="metrics"
    )
    naked_doc, naked, naked_path = load_output(
        args.naked_manifest, expected_schema=NAKED_SCHEMA, output="rates"
    )

    count_test_doc, count_test_path = checked_reference(
        test_doc["count_test_manifest"], args.test_manifest
    )
    count_naked_doc, count_naked_path = checked_reference(
        naked_doc["count_naked_manifest"], args.naked_manifest
    )
    policy_ids = {
        str(count_test_doc["policy_id"]),
        str(count_naked_doc["policy_id"]),
        str(replicate_doc["policy_id"]),
        str(depth_doc["policy_id"]),
    }
    if len(policy_ids) != 1:
        raise ValueError("detector evidence uses different frozen policies")
    consensus_ids = {str(test_doc["consensus_id"]), str(naked_doc["consensus_id"])}
    if len(consensus_ids) != 1:
        raise ValueError("test and naked-DNA evidence use different consensus freezes")
    if test_doc["selected_operator"] != args.operator:
        raise ValueError("requested operator was not selected on validation")
    if naked_doc["selected_operator"] != args.operator:
        raise ValueError("naked-DNA evidence uses a different selected operator")

    report = tasks.copy()
    frames = (
        test_evidence(test, args.operator),
        bootstrap_evidence(bootstrap, args.operator),
        replicate_evidence(replicates),
        depth_evidence(depth),
        naked_evidence(naked, args.operator),
    )
    for frame in frames:
        report = report.merge(frame, on=["cell", "tf"], how="left", validate="one_to_one")
    classifications = report.apply(classify_task, axis=1, result_type="expand")
    report["detector_classification"] = classifications[0]
    report["recommended_configuration"] = classifications[1]
    report["operator"] = args.operator
    report["policy_id"] = next(iter(policy_ids))
    report["research_only"] = True

    summary = (
        report.groupby(
            ["detector_classification", "recommended_configuration"],
            dropna=False,
            sort=True,
        )
        .size()
        .rename("tasks")
        .reset_index()
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    report_path = args.outdir / "frozen_detector_per_tf_evidence.tsv"
    summary_path = args.outdir / "frozen_detector_evidence_summary.tsv"
    report.to_csv(report_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    inputs = [
        args.study,
        base_study,
        args.test_manifest,
        test_path,
        bootstrap_path,
        count_test_path,
        args.replicate_manifest,
        replicate_path,
        args.depth_manifest,
        depth_path,
        args.naked_manifest,
        naked_path,
        count_naked_path,
    ]
    manifest = {
        "schema": SCHEMA,
        "operator": args.operator,
        "policy_id": next(iter(policy_ids)),
        "consensus_id": next(iter(consensus_ids)),
        "development_tasks": int(len(report)),
        "eligible_tasks": int(report["test_auroc"].notna().sum()),
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
    immutable_write_json(args.outdir / "frozen_detector_evidence_manifest.json", manifest)
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
