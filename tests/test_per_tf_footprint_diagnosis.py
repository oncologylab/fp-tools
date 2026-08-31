from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from diagnose_per_tf_footprint_failures import (  # noqa: E402
    bias_control_table,
    classify_diagnosis,
    development_models,
    summarize_tf_consistency,
)


def diagnosis_row(**updates: float) -> pd.Series:
    values = {
        "positive_sites": 200,
        "current_auroc": 0.50,
        "current_auprc": 0.20,
        "prevalence": 0.20,
        "development_best_auroc": 0.50,
        "unlabeled_auroc": 0.50,
        "supervised_functional_relative_auprc_gain": 0.20,
        "same_bias_residual_auroc_gain": 0.0,
        "factorial_bias_model_contribution_auroc": 0.0,
    }
    values.update(updates)
    return pd.Series(values)


def test_failure_classification_distinguishes_power_and_current_detection() -> None:
    power = classify_diagnosis(diagnosis_row(positive_sites=10))[0]
    detected = classify_diagnosis(
        diagnosis_row(current_auroc=0.72, current_auprc=0.50)
    )[0]
    assert power == "power_limited"
    assert detected == "detectable_current"


def test_failure_classification_exposes_selector_gap() -> None:
    selected = classify_diagnosis(
        diagnosis_row(
            development_best_auroc=0.68,
            unlabeled_auroc=0.64,
        )
    )[0]
    failed_selector = classify_diagnosis(
        diagnosis_row(
            development_best_auroc=0.68,
            unlabeled_auroc=0.45,
        )
    )[0]
    assert selected == "detectable_with_label_free_route"
    assert failed_selector == "shape_model_selection_limited"


def test_development_model_comparison_applies_gp_gate_per_tf() -> None:
    rows = []
    for smoother, auprc, seconds in (("spline", 0.20, 1.0), ("gp", 0.24, 1.4)):
        rows.append(
            {
                "cell": "K562",
                "tf": "TF1",
                "motif_family": "F1",
                "bias_configuration": "MT_SELMA10_4m4",
                "candidate_id": smoother,
                "family": "count",
                "smoother": smoother,
                "background": "none",
                "window": 30,
                "channel": "combined_residual",
                "training_pool": "tf",
                "status": "ok",
                "converged": True,
                "selection_score": auprc,
                "auroc": 0.65,
                "auprc": auprc,
                "brier": 0.20,
                "functional_separation": 0.2,
                "depletion": 0.3,
                "width": 20,
                "shoulder_distance": 30,
                "asymmetry": 0.1,
                "periodicity": 0.01,
                "fit_seconds": seconds,
            }
        )
    _overall, comparison = development_models(pd.DataFrame(rows))
    assert comparison.loc[0, "gp_relative_auprc_gain_over_spline"] > 0.05
    assert bool(comparison.loc[0, "gp_passes_relative_auprc_gate"])
    assert comparison.loc[0, "gp_runtime_ratio_over_spline"] == 1.4


def test_bias_control_labels_encode_model_and_shift() -> None:
    selection = pd.DataFrame(
        {
            "model": ["selma10", "loglinear81"],
            "shift_reverse": [-4, -5],
            "mean_conditional_nll": [4.0, 3.9],
            "mean_nll_gain": [0.6, 0.7],
            "passed_control_likelihood": [True, True],
            "model_size_mb": [0.01, 0.03],
            "runtime_seconds": [0.1, 0.2],
        }
    )
    control = bias_control_table(selection)
    assert set(control["bias_configuration"]) == {
        "MT_SELMA10_4m4",
        "MT_LOG81_4m5",
    }


def test_tf_summary_requires_rescue_in_every_context() -> None:
    diagnosis = pd.DataFrame(
        {
            "cell": ["A", "B"],
            "tf": ["TF1", "TF1"],
            "motif_family": ["F1", "F1"],
            "current_auroc": [0.45, 0.50],
            "development_best_auroc": [0.65, 0.55],
            "development_best_auroc_gain_vs_dwm": [0.20, 0.05],
            "unlabeled_auroc": [0.62, 0.54],
            "unlabeled_auroc_gain_vs_dwm": [0.17, 0.04],
            "failure_classification": ["x", "y"],
            "development_best_bias_configuration": ["M1", "M1"],
            "development_best_candidate_id": ["C1", "C1"],
        }
    )
    summary = summarize_tf_consistency(diagnosis)
    assert bool(summary.loc[0, "consistently_hard_current"])
    assert not bool(summary.loc[0, "development_rescue_in_all_contexts"])
