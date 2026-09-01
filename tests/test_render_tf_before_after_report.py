import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "render_tf_before_after_report.py"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("render_tf_before_after_report", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_crossfit_covariate_residuals_are_deterministic_and_label_free():
    rng = np.random.default_rng(13)
    groups = np.repeat(["chr19", "chr20", "chr21", "chr22", "chrX"], 40)
    covariates = rng.normal(size=(len(groups), 3))
    values = 2.5 * covariates[:, 0] - 0.7 * covariates[:, 1]
    values += rng.normal(scale=0.05, size=len(groups))
    first = module.crossfit_covariate_residuals(
        values, covariates, groups, ridge_alpha=1e-8
    )
    second = module.crossfit_covariate_residuals(
        values, covariates, groups, ridge_alpha=1e-8
    )
    assert np.array_equal(first, second)
    assert first.shape == values.shape
    assert np.isfinite(first).all()
    assert abs(np.corrcoef(first, covariates[:, 0])[0, 1]) < 0.1


def test_crossfit_covariate_residuals_require_multiple_groups():
    with pytest.raises(ValueError, match="at least two groups"):
        module.crossfit_covariate_residuals(
            np.arange(4.0), np.ones((4, 2)), np.repeat("chr19", 4)
        )


def test_naked_dna_evidence_reports_calls_and_upper_bound():
    summary = pd.DataFrame(
        {
            "cell": ["K562"],
            "tf": ["CTCF"],
            "candidate_valid": [123],
            "candidate_calls": [1],
            "candidate_false_positive_rate": [1 / 123],
            "candidate_false_positive_rate_upper_95": [0.044617],
        }
    )
    text, metrics = module._naked_dna_evidence(summary, "K562", "CTCF")
    assert text == "1/123 calls (0.8%); upper 95% bound 4.5%"
    assert metrics["naked_dna_candidate_valid"] == 123
    assert metrics["naked_dna_candidate_calls"] == 1
