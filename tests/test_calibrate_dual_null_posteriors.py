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


def test_empirical_pvalues_and_bh_are_finite_sample_safe() -> None:
    pvalues = module.empirical_upper_tail_pvalues(
        np.array([0.1, 0.2, 0.3, 0.4]),
        np.array([0.05, 0.25, 0.5, np.nan]),
    )
    assert pvalues[:3].tolist() == pytest.approx([1.0, 0.6, 0.2])
    assert np.isnan(pvalues[3])
    adjusted = module.benjamini_hochberg(np.array([0.01, 0.04, 0.03, np.nan]))
    assert adjusted[:3].tolist() == pytest.approx([0.03, 0.04, 0.04])
    assert np.isnan(adjusted[3])


def test_dual_empirical_pvalue_requires_both_nulls() -> None:
    primary = _scores(np.linspace(0.0, 0.8, 100), "chr1")
    secondary = _scores(np.linspace(0.0, 0.99, 100), "chr1")
    validation = _scores(np.array([0.95, 0.995]), "chr16")
    summary, calls = module.calibrate_dual_null(
        primary,
        secondary,
        validation,
        primary_alpha=0.05,
        secondary_alpha=0.05,
        fdr=0.05,
    )
    assert calls.loc[0, "naked_dna_pvalue"] < calls.loc[0, "shifted_atac_pvalue"]
    assert calls.loc[0, "dual_null_pvalue"] == calls.loc[0, "shifted_atac_pvalue"]
    assert summary.loc[0, "validation_fdr_calls"] == int(
        calls["dual_null_fdr_call"].sum()
    )
