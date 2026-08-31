from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from render_label_free_per_tf_aggregate_audit import (  # noqa: E402
    conservative_aggregate_difference,
    curve_summary,
    development_best,
    normal_difference_band,
    normal_mean_band,
)


def test_normal_bands_have_expected_shapes_and_contain_means() -> None:
    values = np.arange(40, dtype=float).reshape(10, 4)
    mean, lower, upper = normal_mean_band(values)
    assert mean.shape == lower.shape == upper.shape == (4,)
    assert np.all(lower < mean)
    assert np.all(mean < upper)
    difference, difference_lower, difference_upper = normal_difference_band(
        values + 2.0, values
    )
    assert np.allclose(difference, 2.0)
    assert np.all(difference_lower < difference)
    assert np.all(difference < difference_upper)


def test_conservative_stored_difference_uses_opposite_interval_edges() -> None:
    positive = pd.DataFrame(
        {
            "position": [-1, 0, 1],
            "mean": [2.0, 2.0, 2.0],
            "lower_95": [1.5, 1.5, 1.5],
            "upper_95": [2.5, 2.5, 2.5],
        }
    )
    negative = pd.DataFrame(
        {
            "position": [-1, 0, 1],
            "mean": [1.0, 1.0, 1.0],
            "lower_95": [0.75, 0.75, 0.75],
            "upper_95": [1.25, 1.25, 1.25],
        }
    )
    difference, lower, upper = conservative_aggregate_difference(positive, negative)
    assert np.allclose(difference, 1.0)
    assert np.allclose(lower, 0.25)
    assert np.allclose(upper, 1.75)


def test_curve_summary_detects_central_protection_shape() -> None:
    positions = np.arange(-100, 101)
    difference = np.zeros_like(positions, dtype=float)
    difference[np.abs(positions) <= 10] = -2.0
    lower = difference - 0.1
    upper = difference + 0.1
    summary = curve_summary(positions, difference, lower, upper)
    assert summary["center_mean_difference"] == -2.0
    assert summary["shoulder_mean_difference"] == 0.0
    assert summary["protection_shape_contrast"] == 2.0
    assert summary["visual_rms_difference"] > 0


def test_development_best_is_explicitly_marked_label_selected() -> None:
    rows = []
    for candidate, score in (("weak", 0.5), ("strong", 0.8)):
        rows.append(
            {
                "cell": "K562",
                "tf": "TF1",
                "motif_family": "F1",
                "bias_configuration": "MT_SELMA10_4m4",
                "candidate_id": candidate,
                "family": "count",
                "channel": "combined_residual",
                "status": "ok",
                "converged": True,
                "selection_score": score,
                "auroc": score,
                "auprc": score - 0.1,
            }
        )
    selected = development_best(pd.DataFrame(rows))
    assert selected.loc[0, "candidate_id"] == "strong"
    assert bool(selected.loc[0, "selection_labels_used"])
    assert "diagnostic" in selected.loc[0, "display_method"].lower()
