from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))


def evidence():
    import render_direct_bias_stability_report as m
    metrics, shapes, bootstrap, coverage, safety = [], [], [], [], []
    for cell in ("HepG2", "K562"):
        runs = [(f"{cell}_rep1", depth, seed) for depth in ("10m", "25m", "50m") for seed in range(2026, 2031)]
        runs += [(f"{cell}_rep{rep}", "full", 2026) for rep in (1, 2, 3)]
        for sample, depth, seed in runs:
            key = dict(cell=cell, sample=sample, depth=depth, seed=seed, tf="CTCF", cohort=m.COHORT)
            for method in m.METHODS:
                row = {**key, "method": method, "status": "eligible"}
                for baseline in ("fixed_dwm", "sample_dwm", "raw"):
                    row[f"auroc_gain_over_{baseline}"] = .05
                    row[f"relative_auprc_gain_over_{baseline}"] = .1
                metrics.append(row)
            shapes.append({**key, "method": m.PARAMETRIC, "canonical_protection": True})
            coverage.append({**key, "observed_total_smd": .02})
            if depth == "full" or (depth == "10m" and seed == 2026):
                for baseline in (m.DWM, m.SAMPLE_DWM, m.RAW):
                    bootstrap.append({**key, "baseline": baseline, "auroc_gain_lower_95": .01,
                                      "auroc_gain_upper_95": .1, "relative_auprc_gain_lower_95": .02,
                                      "relative_auprc_gain_upper_95": .2})
        safety.append({"cell": cell, "tf": "CTCF", "passes_candidate_screen": True,
                       "false_positive_calls": 0, "informative_support": 100,
                       "informative_false_positive_rate_upper_95": .037,
                       "paired_false_positive_calls": 0, "paired_support": 100,
                       "paired_false_positive_rate_upper_95": .037,
                       "paired_false_positive_rate_increase_over_dwm": 0})
    return m, [pd.DataFrame(data) for data in (metrics, shapes, bootstrap, coverage, safety)]


def test_complete_report_preserves_balance_drift_and_missing_nonprincipal_intervals():
    m, tables = evidence()
    tables[3].loc[0, "observed_total_smd"] = .12
    result = m.select_evidence(*tables)
    assert len(result) == 36 and result.current_cut_balance_exceeds_0p1.sum() == 1
    assert result.auroc_gain_lower_95_vs_raw.notna().sum() == 8
    assert result.safety_scope.str.contains("not_final_untouched").all()


@pytest.mark.parametrize("failure", ["replicate", "shape", "raw_interval", "direction", "safety"])
def test_report_rejects_incomplete_or_failed_evidence(failure):
    m, tables = evidence()
    if failure == "replicate":
        tables[0] = tables[0].loc[~(tables[0]["sample"].eq("K562_rep3") & tables[0].depth.eq("full"))]
    elif failure == "shape":
        tables[1].loc[0, "canonical_protection"] = False
    elif failure == "raw_interval":
        selected = tables[2].depth.eq("full") & tables[2].baseline.eq(m.RAW)
        tables[2].loc[selected, "auroc_gain_lower_95"] = -.001
    elif failure == "direction":
        tables[0].loc[tables[0].method.eq(m.PARAMETRIC), "relative_auprc_gain_over_raw"] = -.01
    else:
        tables[4].loc[0, "paired_false_positive_rate_upper_95"] = .051
    with pytest.raises(ValueError):
        m.select_evidence(*tables)
