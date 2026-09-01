from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks/scripts/render_bias_shrinkage_before_after.py"
SPEC = importlib.util.spec_from_file_location("render_bias_shrinkage_report", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(module)


def _metrics() -> pd.DataFrame:
    rows = []
    for cell in ("CellA", "CellB"):
        for tf in ("SAFE", "ONE_CELL", "UNSAFE"):
            rows.append(
                {
                    "cell": cell,
                    "tf": tf,
                    "motif_family": tf,
                    "method": module.CANDIDATE_METHOD,
                    "status": "eligible",
                    "positive_sites": 250,
                    "negative_sites": 250,
                    "auroc": 0.65,
                    "auprc": 0.48,
                    "raw_auroc": 0.63,
                    "raw_auprc": 0.46,
                    "dwm_auroc": 0.55,
                    "dwm_auprc": 0.40,
                    "auroc_gain_over_dwm": 0.10,
                    "relative_auprc_gain_over_dwm": 0.20,
                }
            )
    return pd.DataFrame(rows)


def test_qualification_requires_two_significant_safe_contexts() -> None:
    metrics = _metrics()
    bootstrap_rows = []
    safety_rows = []
    for row in metrics.to_dict("records"):
        lower = 0.02
        if row["tf"] == "ONE_CELL" and row["cell"] == "CellB":
            lower = -0.01
        bootstrap_rows.append(
            {
                "cell": row["cell"],
                "tf": row["tf"],
                "method": module.CANDIDATE_METHOD,
                "baseline": module.BASELINE_METHOD,
                "auroc_gain_lower_95": lower,
                "auroc_gain_upper_95": 0.18,
                "relative_auprc_gain_lower_95": lower,
                "relative_auprc_gain_upper_95": 0.35,
            }
        )
        safety_rows.append(
            {
                "cell": row["cell"],
                "tf": row["tf"],
                "method": module.CANDIDATE_METHOD,
                "finite_sites": 200,
                "calls": 0,
                "false_positive_rate": 0.0,
                "false_positive_rate_upper_95": 0.019,
                "false_positive_rate_increase_over_dwm": 0.0,
                "passes_safety": row["tf"] != "UNSAFE",
            }
        )
    result = module.qualified_candidates(
        metrics,
        pd.DataFrame(bootstrap_rows),
        pd.DataFrame(safety_rows),
        minimum_contexts=2,
    )
    reported = set(result.loc[result["report_qualified"], "tf"])
    assert reported == {"SAFE"}


def test_renderer_writes_one_page_pdf(tmp_path: Path) -> None:
    rows = _metrics().query("tf == 'SAFE'").copy()
    rows["auroc_gain_lower_95"] = 0.02
    rows["auroc_gain_upper_95"] = 0.18
    rows["relative_auprc_gain_lower_95"] = 0.05
    rows["relative_auprc_gain_upper_95"] = 0.35
    rows["false_positive_rate"] = 0.0
    rows["false_positive_rate_upper_95"] = 0.019
    profile_rows = []
    for cell in ("CellA", "CellB"):
        for method in (module.BASELINE_METHOD, module.CANDIDATE_METHOD):
            for position in range(-100, 101):
                effect = -1.0 if abs(position) <= 10 else 0.1 * (position % 5)
                profile_rows.append(
                    {
                        "cell": cell,
                        "tf": "SAFE",
                        "method": method,
                        "position": position,
                        "positive_minus_negative": effect,
                        "lower_95": effect - 0.2,
                        "upper_95": effect + 0.2,
                    }
                )
    output = tmp_path / "report.pdf"
    module.render_report(
        rows,
        pd.DataFrame(profile_rows),
        alpha=0.8,
        source="parametric_lambda",
        output=output,
    )
    assert output.read_bytes().startswith(b"%PDF-")
