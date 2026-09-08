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
                (report.RAW, 0.02, 0.06),
                (report.DEEP_BIAS, -0.02, 0.01),
            )
        ]
    )


def _safety() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                **{key: 0 for key in report.SAFETY_COLUMNS},
                "cell": "K562",
                "tf": "CTCF",
                "support": 200,
                "informative_support": 200,
                "paired_support": 200,
                "false_positive_calls": 0,
                "false_positive_rate": 0.0,
                "false_positive_rate_upper_95": 0.0188,
                "passes_candidate_screen": True,
                "dispersion": 1.6144318781609228,
            }
        ]
    )


def test_report_rows_require_complete_paired_comparison() -> None:
    selected, paired, safety = report.select_report_rows(
        _metrics(), _bootstrap(), _safety(), cell="K562", tf="CTCF"
    )
    assert set(selected["method"]) == set(report.METHODS)
    assert set(paired["baseline"]) == {report.DWM, report.RAW, report.DEEP_BIAS}
    assert int(safety["false_positive_calls"]) == 0

    table = report.build_metrics_table(selected, paired, safety)
    candidate = table[table["method"].eq(report.PARAMETRIC)].iloc[0]
    assert candidate["auroc_gain_lower_95_vs_dwm"] == pytest.approx(0.03)
    assert candidate["auroc_gain_lower_95_vs_chrombpnet_bias"] == pytest.approx(
        -0.02
    )
    assert candidate.false_positive_rate_upper_95 == 0.0188
    assert table.loc[table.method.ne(report.PARAMETRIC), "false_positive_rate_upper_95"].isna().all()

    incomplete = _bootstrap()[lambda frame: frame["baseline"].ne(report.DEEP_BIAS)]
    with pytest.raises(ValueError, match="lacks DWM, raw or ChromBPNet"):
        report.select_report_rows(
            _metrics(), incomplete, _safety(), cell="K562", tf="CTCF"
        )


def test_legacy_factorization_safety_table_is_rejected():
    legacy = pd.DataFrame([{"cell": "K562", "tf": "CTCF", "residual": "deviance",
                            "finite_support": 200, "false_positive_calls": 0}])
    with pytest.raises(ValueError, match="coverage-qualified direct-operator"):
        report.select_report_rows(_metrics(), _bootstrap(), legacy, cell="K562", tf="CTCF")


@pytest.mark.parametrize("baseline", [report.DWM, report.RAW])
def test_report_rejects_nonsignificant_gain_and_failed_safety(baseline):
    bootstrap = _bootstrap()
    bootstrap.loc[bootstrap.baseline.eq(baseline), "auroc_gain_lower_95"] = 0
    with pytest.raises(ValueError, match="significant AUROC and AUPRC"):
        report.select_report_rows(_metrics(), bootstrap, _safety(), cell="K562", tf="CTCF")
    safety = _safety(); safety["passes_candidate_screen"] = False
    with pytest.raises(ValueError, match="passing exact-operator"):
        report.select_report_rows(_metrics(), _bootstrap(), safety, cell="K562", tf="CTCF")


def test_safety_provenance_must_name_the_exact_comparison():
    config = {"cells": {"K562": {"comparison_result_id": "comparison-A"}}}
    config["configuration_id"] = report.canonical_id(config)
    source = {"schema": report.SAFETY_SCHEMA, "operator": report.OPERATOR, "replicate": "NakedDNA_rep2",
              "configuration_id": config["configuration_id"], "threshold_id": "threshold-A"}
    source["safety_id"] = report.canonical_id(source)
    coverage = {"schema": report.SCHEMA_COVERAGE, "operator": report.OPERATOR,
                "replicate": "NakedDNA_rep2", "new_holdout_labels_read": False,
                "final_untouched_safety": False, "configuration_id": config["configuration_id"],
                "threshold_id": "threshold-A"}
    coverage["coverage_id"] = report.canonical_id(coverage)
    report.validate_safety_link(coverage, source, config, {"result_id": "comparison-A"}, cell="K562")
    with pytest.raises(ValueError, match="does not match"):
        report.validate_safety_link(coverage, source, config, {"result_id": "comparison-B"}, cell="K562")
    coverage["operator"] = "factorization_expected_unbound"
    with pytest.raises(ValueError, match="direct-operator"):
        report.validate_safety_link(coverage, source, config, {"result_id": "comparison-A"}, cell="K562")
