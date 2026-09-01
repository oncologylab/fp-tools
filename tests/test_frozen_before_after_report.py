from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks/scripts/render_frozen_functional_before_after.py"
SPEC = importlib.util.spec_from_file_location("render_frozen_before_after", SCRIPT)
module = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(module)


def _metrics() -> pd.DataFrame:
    rows = []
    for tf, status in (("SAFE", "eligible"), ("UNSAFE", "eligible"), ("SMALL", "underpowered")):
        rows.append(
            {
                "cell": "CellA",
                "tf": tf,
                "motif_family": tf,
                "candidate_id": "candidate",
                "method": "frozen_candidate",
                "status": status,
                "n_positive": 200,
                "n_negative": 200,
                "dwm_auroc": 0.55,
                "dwm_auprc": 0.40,
                "auroc": 0.65,
                "auprc": 0.48,
                "auroc_gain_over_dwm": 0.10,
                "relative_auprc_gain_over_dwm": 0.20,
                "functional_separation_relative_change_over_dwm": 0.30,
            }
        )
    return pd.DataFrame(rows)


def test_report_qualification_requires_significance_safety_and_power() -> None:
    bootstrap = pd.DataFrame(
        {
            "cell": ["CellA"] * 3,
            "tf": ["SAFE", "UNSAFE", "SMALL"],
            "method": ["frozen_candidate"] * 3,
            "auroc_gain_lower_95": [0.02, 0.02, 0.02],
            "auroc_gain_upper_95": [0.18, 0.18, 0.18],
            "relative_auprc_gain_lower_95": [0.05, 0.05, 0.05],
            "relative_auprc_gain_upper_95": [0.35, 0.35, 0.35],
        }
    )
    safety = pd.DataFrame(
        {
            "cell": ["CellA"] * 3,
            "tf": ["SAFE", "UNSAFE", "SMALL"],
            "method": ["paired_safety"] * 3,
            "candidate_false_positive_rate": [0.01, 0.02, 0.01],
            "candidate_wilson_upper_95": [0.04, 0.06, 0.04],
            "candidate_minus_dwm": [0.0, 0.02, 0.0],
            "passes_safety": [True, False, True],
        }
    )
    result = module.qualified_candidates(_metrics(), bootstrap, safety)
    observed = dict(zip(result["tf"], result["report_qualified"], strict=True))
    assert observed == {"SAFE": True, "SMALL": False, "UNSAFE": False}


def test_raw_guardrail_emits_only_robust_or_depth_dependent_tasks() -> None:
    summary = pd.DataFrame(
        {
            "cell": ["CellA", "CellA", "CellA"],
            "tf": ["ROBUST", "DEPTH", "SENSITIVE"],
            "report_qualified": [True, True, True],
        }
    )
    evidence = pd.DataFrame(
        {
            "cell": ["CellA", "CellA", "CellA"],
            "tf": ["ROBUST", "DEPTH", "SENSITIVE"],
            "detector_classification": [
                "robust_tf_specific_gain",
                "depth_dependent_tf_specific_gain",
                "support_or_depth_sensitive_gain",
            ],
        }
    )
    result = module.apply_raw_guardrail(summary, evidence)
    observed = dict(zip(result["tf"], result["report_qualified"], strict=True))
    assert observed == {"DEPTH": True, "ROBUST": True, "SENSITIVE": False}


def test_report_renderer_writes_a_one_page_pdf(tmp_path: Path) -> None:
    row = _metrics().iloc[0].copy()
    row["auroc_gain_lower_95"] = 0.02
    row["auroc_gain_upper_95"] = 0.18
    row["relative_auprc_gain_lower_95"] = 0.05
    row["relative_auprc_gain_upper_95"] = 0.35
    row["candidate_false_positive_rate"] = 0.01
    row["candidate_wilson_upper_95"] = 0.04
    rows = []
    for method in (module.BASELINE_METHOD, "frozen_candidate"):
        for position in range(-2, 3):
            rows.append(
                {
                    "cell": "CellA",
                    "tf": "SAFE",
                    "method": method,
                    "position": position,
                    "positive_mean": -0.2 if position == 0 else 0.1,
                    "negative_mean": 0.0,
                    "positive_minus_negative": -0.2 if position == 0 else 0.1,
                    "lower_95": -0.3 if position == 0 else 0.0,
                    "upper_95": -0.1 if position == 0 else 0.2,
                }
            )
    output = tmp_path / "report.pdf"
    module.render_report(row, pd.DataFrame(rows), output)
    assert output.read_bytes().startswith(b"%PDF-")
