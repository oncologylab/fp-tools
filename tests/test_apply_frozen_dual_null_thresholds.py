import importlib.util
from pathlib import Path
import sys

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "apply_frozen_dual_null_thresholds.py"
spec = importlib.util.spec_from_file_location(
    "apply_frozen_dual_null_thresholds", SCRIPT
)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def fixtures():
    route = {
        "cell": "K562",
        "tf": "MEF2A",
        "motif_family": "MEF2",
        "method": "frozen_policy_candidate",
        "candidate_id": "count_gp.bg_gp-long.window_30",
        "bias_configuration": "MT_SELMA10_4m4",
        "replicate": "rep2",
    }
    scores = pd.DataFrame(
        [
            {**route, "binding_probability": 0.90, "total_signal": 5, "valid": True},
            {**route, "binding_probability": 0.79, "total_signal": 4, "valid": True},
            {**route, "binding_probability": 0.99, "total_signal": 0, "valid": True},
            {**route, "binding_probability": 0.95, "total_signal": 6, "valid": False},
        ]
    )
    calibration = pd.DataFrame(
        [
            {
                key: value
                for key, value in route.items()
                if key in module.JOIN_COLUMNS
            }
            | {"dual_null_threshold": 0.80}
        ]
    )
    return scores, calibration


def test_applies_frozen_threshold_and_requires_signal():
    scores, calibration = fixtures()
    calls, rates = module.apply_thresholds(scores, calibration)
    assert calls["dual_null_call"].tolist() == [True, False, False, False]
    assert rates.loc[0, "valid"] == 3
    assert rates.loc[0, "informative"] == 2
    assert rates.loc[0, "calls"] == 1
    assert rates.loc[0, "all_site_rate"] == pytest.approx(1 / 3)
    assert rates.loc[0, "informative_rate"] == pytest.approx(1 / 2)
    assert (
        rates.loc[0, "all_site_rate_lower_95"]
        < rates.loc[0, "all_site_rate"]
        < rates.loc[0, "all_site_rate_upper_95"]
    )
    assert (
        rates.loc[0, "informative_rate_lower_95"]
        < rates.loc[0, "informative_rate"]
        < rates.loc[0, "informative_rate_upper_95"]
    )


def test_rejects_missing_or_duplicate_thresholds():
    scores, calibration = fixtures()
    with pytest.raises(ValueError, match="no frozen threshold"):
        module.apply_thresholds(
            scores.assign(tf="MEF2D"), calibration
        )
    with pytest.raises(ValueError, match="duplicate frozen thresholds"):
        module.apply_thresholds(
            scores, pd.concat([calibration, calibration], ignore_index=True)
        )
