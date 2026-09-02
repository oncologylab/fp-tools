from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import pytest


SCRIPTS = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPTS))

import evaluate_shape_guarded_functional_screen as screen  # noqa: E402


def test_shape_guard_distinguishes_protection_from_center_enrichment():
    rows = []
    for tf, difference, depletion in (
        ("PROTECTED", -0.1, 0.3),
        ("ENRICHED", 0.1, -0.2),
    ):
        for position in range(-6, 7):
            rows.append(
                {
                    "cell": "A",
                    "tf": tf,
                    "method": "frozen_model",
                    "position": position,
                    "positive_minus_negative": difference,
                    "upper_95": difference + 0.02,
                    "depletion": depletion,
                    "width": 6.0,
                    "shoulder_distance": 25.0,
                    "asymmetry": 0.0,
                    "periodicity": 0.2,
                }
            )
    result = screen.task_shapes(pd.DataFrame(rows)).set_index("tf")
    assert result.loc["PROTECTED", "shape_has_strong_canonical_protection"]
    assert not result.loc["ENRICHED", "shape_has_central_protection"]


def test_common_support_validates_labels_and_coverage():
    candidate = pd.DataFrame(
        [
            {
                "cell": "A",
                "tf": "TF",
                "TFBS_chr": "chr19",
                "TFBS_start": 10,
                "TFBS_end": 12,
                "TFBS_strand": "+",
                "label": 1,
                "candidate_probability": 0.8,
                "dwm_score": 0.2,
            }
        ]
    )
    raw = candidate.rename(columns={"label": "chip_label"}).copy()
    raw["raw_score"] = 0.3
    merged, coverage = screen.common_support(
        candidate, raw, minimum_coverage=1.0
    )
    assert len(merged) == 1
    assert coverage.iloc[0]["common_support_fraction"] == 1.0
    raw["chip_label"] = 0
    with pytest.raises(ValueError, match="labels differ"):
        screen.common_support(candidate, raw, minimum_coverage=1.0)


def test_common_support_uses_authoritative_raw_guardrail_on_name_collision():
    candidate = pd.DataFrame(
        [
            {
                "cell": "A",
                "tf": "TF",
                "TFBS_chr": "chr19",
                "TFBS_start": 10,
                "TFBS_end": 12,
                "TFBS_strand": "+",
                "label": 1,
                "candidate_probability": 0.8,
                "dwm_score": 0.2,
                "raw_score": -99.0,
            }
        ]
    )
    raw = candidate.rename(columns={"label": "chip_label"}).copy()
    raw["raw_score"] = 0.3
    merged, _coverage = screen.common_support(
        candidate, raw, minimum_coverage=1.0
    )
    assert merged.iloc[0]["raw_score"] == 0.3
    assert "raw_score_candidate" not in merged


def test_classification_requires_shape_and_safety():
    base = {
        "status": "eligible",
        "naked_passes_safety": True,
        "auroc_gain_over_raw": 0.05,
        "relative_auprc_gain_over_raw": 0.1,
        "raw_auroc_gain_lower_95": 0.01,
        "raw_relative_auprc_gain_lower_95": 0.02,
        "auroc_gain_over_dwm": 0.1,
        "relative_auprc_gain_over_dwm": 0.2,
        "shape_has_central_protection": True,
        "shape_has_strong_canonical_protection": True,
    }
    assert screen.classify(base) == "strong_shape_guarded_gain"
    base["shape_has_central_protection"] = False
    base["shape_has_strong_canonical_protection"] = False
    assert screen.classify(base) == "occupancy_gain_without_protection"
    base["naked_passes_safety"] = False
    assert screen.classify(base) == "naked_dna_safety_failed_or_missing"
