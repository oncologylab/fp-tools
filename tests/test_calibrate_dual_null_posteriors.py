from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "calibrate_dual_null_posteriors.py"
spec = importlib.util.spec_from_file_location("calibrate_dual_null_posteriors", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def _scores(values: np.ndarray, chromosome: str) -> pd.DataFrame:
    count = len(values)
    return pd.DataFrame(
        {
            "cell": "K562",
            "tf": "CTCF",
            "motif_family": "CTCF",
            "method": "candidate",
            "candidate_id": "model",
            "bias_configuration": "LOG81",
            "site_hash": np.arange(count),
            "TFBS_chr": chromosome,
            "TFBS_start": np.arange(count),
            "binding_probability": values,
            "valid": True,
            "informative": True,
        }
    )


def test_dual_threshold_uses_stricter_null() -> None:
    primary = _scores(np.linspace(0.0, 0.8, 100), "chr1")
    secondary = _scores(np.linspace(0.0, 0.99, 100), "chr1")
    validation = _scores(np.linspace(0.0, 0.98, 100), "chr16")
    summary, calls = module.calibrate_dual_null(
        primary,
        secondary,
        validation,
        primary_alpha=0.05,
        secondary_alpha=0.05,
    )
    assert summary.loc[0, "threshold_source"] == "secondary_shifted_atac"
    assert summary.loc[0, "dual_null_threshold"] == pytest.approx(
        summary.loc[0, "secondary_threshold"]
    )
    assert calls["dual_null_call"].sum() == summary.loc[0, "validation_calls"]


def test_dual_threshold_rejects_group_mismatch() -> None:
    primary = _scores(np.linspace(0.0, 1.0, 20), "chr1")
    secondary = primary.copy()
    validation = primary.copy()
    secondary["tf"] = "MAX"
    with pytest.raises(ValueError, match="groups differ"):
        module.calibrate_dual_null(
            primary,
            secondary,
            validation,
            primary_alpha=0.05,
            secondary_alpha=0.05,
        )
