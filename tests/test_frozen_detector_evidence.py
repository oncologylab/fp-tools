from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import pytest


SCRIPTS = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import summarize_frozen_detector_evidence as evidence  # noqa: E402


def eligible_row(**updates) -> pd.Series:
    values = {
        "test_status": "eligible",
        "naked_passes_safety": True,
        "shape_has_central_protection": True,
        "test_auroc_gain_over_raw": 0.04,
        "test_relative_auprc_gain_over_raw": 0.12,
        "test_raw_bootstrap_auroc_gain_lower_95": 0.01,
        "test_raw_bootstrap_relative_auprc_gain_lower_95": 0.02,
        "replicate_auroc_gain_over_raw_positive_fraction": 1.0,
        "replicate_auprc_gain_over_raw_positive_fraction": 1.0,
        "depth_both_gain_over_raw_fraction": 1.0,
        "depth_high_both_gain_over_raw_fraction": 1.0,
        "test_auroc_gain_over_dwm": 0.06,
        "test_relative_auprc_gain_over_dwm": 0.20,
        "replicate_auroc_gain_positive_fraction": 1.0,
        "replicate_auprc_gain_positive_fraction": 1.0,
        "depth_both_gain_over_dwm_fraction": 1.0,
    }
    values.update(updates)
    return pd.Series(values)


def test_classifies_robust_and_depth_dependent_gain():
    assert evidence.classify_task(eligible_row())[0] == "robust_tf_specific_gain"
    depth_limited = eligible_row(depth_both_gain_over_raw_fraction=0.8)
    assert evidence.classify_task(depth_limited)[0] == (
        "depth_dependent_tf_specific_gain"
    )


def test_raw_inconsistency_is_not_promoted_as_dwm_gain():
    row = eligible_row(
        replicate_auroc_gain_over_raw_positive_fraction=0.0,
        depth_both_gain_over_raw_fraction=0.0,
        depth_high_both_gain_over_raw_fraction=0.0,
    )
    assert evidence.classify_task(row) == (
        "support_or_depth_sensitive_gain",
        "retain_raw_guardrail",
    )


def test_accessibility_gain_without_depletion_is_not_called_a_footprint_gain():
    row = eligible_row(
        shape_has_central_protection=False,
        depth_both_gain_over_raw_fraction=0.8,
    )
    assert evidence.classify_task(row) == (
        "occupancy_signal_gain_without_footprint_protection",
        "use_as_occupancy_diagnostic_not_footprint",
    )


def test_depth_summary_separates_low_and_high_depth():
    rows = []
    for depth, raw_gain in (("10m", -0.01), ("25m", 0.02), ("50m", 0.03)):
        for seed in (2026, 2027):
            rows.append(
                {
                    "cell": "HepG2",
                    "sample": "rep1",
                    "tf": "FOXA1",
                    "depth": depth,
                    "seed": seed,
                    "method": "frozen_count",
                    "auroc_gain_over_raw": raw_gain,
                    "relative_auprc_gain_over_raw": raw_gain,
                    "auroc_gain_over_dwm": 0.1,
                    "relative_auprc_gain_over_dwm": 0.1,
                }
            )
    result = evidence.depth_evidence(pd.DataFrame(rows)).iloc[0]
    assert result["depth_both_gain_over_raw_fraction"] == pytest.approx(2 / 3)
    assert result["depth_high_both_gain_over_raw_fraction"] == 1.0
    assert result["depth_low_both_gain_over_raw_fraction"] == 0.0


def test_duplicate_frozen_candidate_is_rejected():
    frame = pd.DataFrame(
        [
            {"cell": "A", "tf": "TF", "method": "frozen_one"},
            {"cell": "A", "tf": "TF", "method": "frozen_two"},
        ]
    )
    with pytest.raises(ValueError, match="multiple frozen candidates"):
        evidence.candidate_rows(frame)


def test_shape_evidence_requires_positive_depletion_for_protection():
    frame = pd.DataFrame(
        [
            {
                "cell": "A",
                "tf": "TF",
                "method": "frozen_count",
                "position": position,
                "center": 0.2,
                "shoulders": 0.1,
                "depletion": -0.1,
                "width": 4.0,
                "shoulder_distance": 20.0,
                "asymmetry": 0.0,
                "periodicity": 0.2,
            }
            for position in (-1, 0, 1)
        ]
    )
    result = evidence.shape_evidence(frame).iloc[0]
    assert not result["shape_has_central_protection"]
