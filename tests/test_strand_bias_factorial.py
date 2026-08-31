from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from analyze_strand_bias_factorial import (  # noqa: E402
    ARMS,
    factorial_effects,
    select_common_settings,
)


def _metrics() -> pd.DataFrame:
    values = {
        "MT_SELMA10_4m4": (0.80, 0.70, 0.15),
        "MT_SELMA10_4m5": (0.70, 0.60, 0.20),
        "MT_LOG81_4m4": (0.65, 0.55, 0.23),
        "MT_LOG81_4m5": (0.60, 0.50, 0.25),
    }
    rows = []
    for setting, selection in (("combined", 1.0), ("shared", 0.5)):
        for arm in ARMS:
            auroc, auprc, brier = values[arm]
            rows.append(
                {
                    "cell": "K562",
                    "tf": "TF1",
                    "motif_family": "F1",
                    "training_scope": "cross_cell_tf",
                    "channel_set": setting,
                    "smoother": "spline",
                    "window_limit": 50.0,
                    "bias_configuration": arm,
                    "status": "ok",
                    "selection_score": selection,
                    "auroc": auroc,
                    "auprc": auprc,
                    "brier": brier,
                }
            )
    return pd.DataFrame(rows)


def test_common_setting_selection_is_balanced_across_arms() -> None:
    selected = select_common_settings(_metrics())
    assert len(selected) == 1
    assert selected.iloc[0]["channel_set"] == "combined"


def test_factorial_shapley_contributions_sum_to_diagonal_gain() -> None:
    metrics = _metrics()
    effects = factorial_effects(metrics, select_common_settings(metrics))
    auroc = effects[effects["metric"] == "auroc"].iloc[0]
    assert np.isclose(auroc["shift_contribution"], 0.075)
    assert np.isclose(auroc["bias_model_contribution"], 0.125)
    assert np.isclose(auroc["total_improvement"], 0.20)
    assert np.isclose(
        auroc["shift_contribution"] + auroc["bias_model_contribution"],
        auroc["total_improvement"],
    )
    brier = effects[effects["metric"] == "brier"].iloc[0]
    assert brier["total_improvement"] > 0
    assert brier["candidate_minus_reference_raw"] < 0
