from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "calibrate_naked_dna_posteriors.py"
spec = importlib.util.spec_from_file_location("calibrate_naked_dna_posteriors", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def _scores(probabilities: np.ndarray, chromosome: str) -> pd.DataFrame:
    count = len(probabilities)
    return pd.DataFrame(
        {
            "cell": ["K562"] * count,
            "tf": ["CTCF"] * count,
            "motif_family": ["CTCF"] * count,
            "method": ["candidate"] * count,
            "candidate_id": ["model"] * count,
            "bias_configuration": ["LOG81"] * count,
            "site_hash": np.arange(count),
            "TFBS_chr": [chromosome] * count,
            "TFBS_start": np.arange(count) * 10,
            "binding_probability": probabilities,
            "valid": [True] * count,
            "informative": [True] * count,
        }
    )


def test_conservative_threshold_controls_calibration_and_transfers() -> None:
    calibration = _scores(np.linspace(0.0, 0.99, 100), "chr1")
    validation = _scores(np.linspace(0.0, 0.98, 100), "chr16")
    summary, calls = module.calibrate(calibration, validation, alpha=0.05)
    assert summary.loc[0, "calibration_calls"] <= 5
    assert summary.loc[0, "calibration_informative_fpr"] <= 0.05
    assert summary.loc[0, "validation_informative_fpr"] <= 0.05
    assert calls["null_calibrated_call"].sum() == summary.loc[0, "validation_calls"]


def test_threshold_handles_ties_and_empty_information() -> None:
    tied = np.repeat([0.1, 0.9], [90, 10])
    threshold = module.conservative_threshold(tied, 0.05)
    assert np.sum(tied >= threshold) <= 5
    assert np.isinf(module.conservative_threshold(np.array([]), 0.05))
    with pytest.raises(ValueError, match="alpha"):
        module.conservative_threshold(tied, 0.0)


def test_calibration_rejects_mismatched_groups() -> None:
    calibration = _scores(np.linspace(0.0, 1.0, 20), "chr1")
    validation = _scores(np.linspace(0.0, 1.0, 20), "chr16")
    validation["tf"] = "MAX"
    with pytest.raises(ValueError, match="groups differ"):
        module.calibrate(calibration, validation, alpha=0.05)
