from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks" / "scripts"))

import freeze_safe_parametric_configuration as freeze  # noqa: E402


def _selection() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "residual": "deviance",
                "mean_relative_auprc_gain": 0.20,
                "standard_error": 0.02,
                "passes_ctcf_gate": True,
            },
            {
                "residual": "pearson",
                "mean_relative_auprc_gain": 0.19,
                "standard_error": 0.01,
                "passes_ctcf_gate": True,
            },
            {
                "residual": "difference",
                "mean_relative_auprc_gain": 0.18,
                "standard_error": 0.01,
                "passes_ctcf_gate": True,
            },
            {
                "residual": "log-ratio",
                "mean_relative_auprc_gain": 0.17,
                "standard_error": 0.01,
                "passes_ctcf_gate": False,
            },
        ]
    )


def test_safe_selection_excludes_naked_dna_failure_before_one_se_rule() -> None:
    safety = pd.DataFrame(
        {
            "residual": ["deviance", "pearson", "difference", "log-ratio"],
            "passes_naked_dna_safety": [False, True, True, True],
        }
    )
    selected, audit = freeze.select_safe_residual(_selection(), safety)
    # Deviance has the best development result but cannot enter the freeze.
    # Pearson and difference are within one SE; fixed simplicity selects Pearson.
    assert selected == "pearson"
    indexed = audit.set_index("residual")
    assert not bool(indexed.loc["deviance", "eligible_for_freeze"])
    assert bool(indexed.loc["pearson", "selected_after_safety"])


def test_safe_selection_fails_closed_when_every_residual_is_ineligible() -> None:
    safety = pd.DataFrame(
        {
            "residual": _selection()["residual"],
            "passes_naked_dna_safety": False,
        }
    )
    with pytest.raises(RuntimeError, match="no residual passed"):
        freeze.select_safe_residual(_selection(), safety)


def test_safe_selection_rejects_incomplete_safety_rows() -> None:
    safety = pd.DataFrame(
        {"residual": ["deviance"], "passes_naked_dna_safety": [True]}
    )
    with pytest.raises(ValueError, match="incomplete"):
        freeze.select_safe_residual(_selection(), safety)
