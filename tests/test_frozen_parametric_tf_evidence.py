from __future__ import annotations

import json
from pathlib import Path
import sys

import pandas as pd
import pytest


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import summarize_frozen_parametric_tf_evidence as evidence  # noqa: E402


def _row(**updates) -> pd.Series:
    values = {
        "candidate_status": "eligible",
        "candidate_auroc": 0.70,
        "candidate_auprc": 0.65,
        "candidate_auroc_gain_over_raw": 0.04,
        "candidate_relative_auprc_gain_over_raw": 0.12,
        "naked_passes_safety": True,
        "replicate_auroc_gain_over_raw_positive_fraction": 1.0,
        "replicate_auprc_gain_over_raw_positive_fraction": 1.0,
        "depth_both_gain_over_raw_fraction": 1.0,
        "raw_auroc_gain_lower_95": 0.01,
        "raw_relative_auprc_gain_lower_95": 0.02,
        "raw_auroc": 0.66,
        "dwm_auroc": 0.60,
    }
    values.update(updates)
    return pd.Series(values)


def test_classify_task_obeys_eligibility_and_safety_first() -> None:
    assert evidence.classify_task(_row(candidate_status="underpowered")) == (
        "underpowered",
        "collect_more_sites_or_labels",
    )
    assert evidence.classify_task(_row(naked_passes_safety=False)) == (
        "safety_limited",
        "reject_candidate",
    )


def test_classify_task_separates_robust_and_partial_raw_gains() -> None:
    assert evidence.classify_task(_row()) == (
        "robust_gain_over_raw",
        "frozen_tf_specific_shrinkage",
    )
    partial = _row(raw_relative_auprc_gain_lower_95=-0.001)
    assert evidence.classify_task(partial) == (
        "replicate_stable_partial_correction_gain",
        "frozen_tf_specific_shrinkage_research_only",
    )


def test_classify_task_identifies_dwm_overcorrection() -> None:
    row = _row(
        candidate_auroc_gain_over_raw=-0.002,
        candidate_relative_auprc_gain_over_raw=0.003,
        replicate_auroc_gain_over_raw_positive_fraction=0.0,
        depth_both_gain_over_raw_fraction=0.0,
    )
    assert evidence.classify_task(row) == (
        "dwm_overcorrection",
        "raw_geometry_research_baseline",
    )


def test_depth_evidence_requires_both_metrics_for_stability() -> None:
    metrics = pd.DataFrame(
        {
            "cell": ["CellA", "CellA"],
            "tf": ["TF1", "TF1"],
            "method": ["candidate", "candidate"],
            "depth": ["10m", "50m"],
            "seed": [2026, 2026],
            "auroc_gain_over_raw": [0.01, -0.01],
            "relative_auprc_gain_over_raw": [0.02, 0.02],
            "auroc_gain_over_dwm": [0.03, 0.04],
            "relative_auprc_gain_over_dwm": [0.05, 0.06],
        }
    )
    result = evidence.depth_evidence(metrics, "candidate").iloc[0]
    assert result["depth_both_gain_over_raw_fraction"] == 0.5
    assert not result["depth_all_gain_over_raw"]
    assert result["depth_high_endpoint"] == "50m"


def test_load_output_rejects_checksum_mismatch(tmp_path: Path) -> None:
    table = tmp_path / "metrics.tsv"
    table.write_text("cell\ttf\nCellA\tTF1\n", encoding="utf-8")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "schema-v1",
                "outputs": {
                    "metrics": {"path": str(table), "sha256": "0" * 64}
                },
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="checksum mismatch"):
        evidence.load_output(
            manifest,
            expected_schema="schema-v1",
            output="metrics",
        )
