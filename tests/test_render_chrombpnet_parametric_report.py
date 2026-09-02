from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import render_chrombpnet_parametric_report as report  # noqa: E402


def _metrics() -> pd.DataFrame:
    rows = []
    for index, method in enumerate(report.METHODS):
        rows.append(
            {
                "cell": "K562",
                "tf": "CTCF",
                "method": method,
                "status": "eligible",
                "auroc": 0.7 + index * 0.02,
                "auprc": 0.6 + index * 0.03,
                "functional_separation": 0.1 + index * 0.04,
                "n_positive": 250,
                "n_negative": 250,
                "auroc_gain_over_dwm": index * 0.02,
                "relative_auprc_gain_over_dwm": index * 0.05,
            }
        )
    return pd.DataFrame(rows)


def _bootstrap() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "cell": "K562",
                "tf": "CTCF",
                "method": report.PARAMETRIC,
                "baseline": baseline,
                "auroc_gain_lower_95": lower,
                "auroc_gain_upper_95": upper,
                "relative_auprc_gain_lower_95": lower * 2,
                "relative_auprc_gain_upper_95": upper * 2,
            }
            for baseline, lower, upper in (
                (report.DWM, 0.03, 0.08),
                (report.DEEP_BIAS, -0.02, 0.01),
            )
        ]
    )


def _safety() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "cell": "K562",
                "tf": "CTCF",
                "residual": "deviance",
                "naked_sites": 200,
                "finite_support": 200,
                "false_positive_calls": 0,
                "false_positive_rate": 0.0,
                "false_positive_rate_upper_95": 0.0188,
            }
        ]
    )


def test_report_rows_require_complete_paired_comparison() -> None:
    selected, paired, safety = report.select_report_rows(
        _metrics(), _bootstrap(), _safety(), cell="K562", tf="CTCF"
    )
    assert set(selected["method"]) == set(report.METHODS)
    assert set(paired["baseline"]) == {report.DWM, report.DEEP_BIAS}
    assert int(safety["false_positive_calls"]) == 0

    table = report.build_metrics_table(selected, paired, safety)
    candidate = table[table["method"].eq(report.PARAMETRIC)].iloc[0]
    assert candidate["auroc_gain_lower_95_vs_dwm"] == pytest.approx(0.03)
    assert candidate["auroc_gain_lower_95_vs_chrombpnet_bias"] == pytest.approx(
        -0.02
    )
    np.testing.assert_allclose(table["false_positive_rate_upper_95"], 0.0188)

    incomplete = _bootstrap()[lambda frame: frame["baseline"].ne(report.DEEP_BIAS)]
    with pytest.raises(ValueError, match="lacks DWM or ChromBPNet"):
        report.select_report_rows(
            _metrics(), incomplete, _safety(), cell="K562", tf="CTCF"
        )
