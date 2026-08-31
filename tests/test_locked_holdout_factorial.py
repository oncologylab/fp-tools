from __future__ import annotations

import numpy as np
import pandas as pd

from benchmarks.scripts.diagnose_locked_holdout_factorial import (
    orientation_diagnostic,
    select_posthoc_winners,
)


def test_orientation_diagnostic_exposes_but_does_not_promote_inversion() -> None:
    labels = np.asarray([0, 0, 1, 1])
    probabilities = np.asarray([0.9, 0.8, 0.2, 0.1])
    result = orientation_diagnostic(labels, probabilities)
    assert result["orientation_inverted_posthoc"] is True
    assert result["oracle_orientation_auroc"] == 1.0


def test_posthoc_winners_are_explicitly_ineligible() -> None:
    metrics = pd.DataFrame(
        {
            "cell": ["A", "A", "A"],
            "tf": ["TF1", "TF1", "TF1"],
            "correction": ["new", "new", "new"],
            "candidate_id": ["b", "a", "failed"],
            "status": ["ok", "ok", "error"],
            "converged": [True, True, False],
            "auroc_gain": [0.1, 0.1, 1.0],
            "relative_auprc_gain": [0.2, 0.2, 1.0],
        }
    )
    winner = select_posthoc_winners(metrics)
    assert winner.loc[0, "candidate_id"] == "a"
    assert winner.loc[0, "selection_status"] == "posthoc_diagnostic_only"
    assert not bool(winner.loc[0, "eligible_for_promotion"])
