from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import compare_corrected_count_reassessment as m


def tables():
    common = {"cell": "HepG2", "tf": "CTCF", "sample": "HepG2_rep1", "depth": "full", "seed": 2026, "cohort": m.COHORT}
    counts = {"n_positive": 200, "n_negative": 200, "status": "eligible"}
    metric = {**common, **counts, "method": m.CANDIDATE, "auroc": .7, "auprc": .7,
              "functional_separation": .3, "auroc_gain_over_sample_dwm": .1, "relative_auprc_gain_over_sample_dwm": .15}
    return {"metrics": pd.DataFrame([metric]),
            "coverage": pd.DataFrame([{**common, **counts}, {**common, "tf": "HMGA1", "n_positive": 1, "n_negative": 1, "status": "underpowered"}]),
            "shapes": pd.DataFrame([{**common, "method": m.CANDIDATE, "canonical_protection": True}]),
            "bootstrap": pd.DataFrame([{**common, "method": m.CANDIDATE, "baseline": "DWM", "auroc_gain_lower_95": .01}])}


def test_preserves_unscored_task_and_reports_signed_change_and_shape_loss():
    old, new = tables(), tables()
    new["metrics"].loc[0, "auroc"] = .68
    new["shapes"].loc[0, "canonical_protection"] = False
    pairs, summary = m.comparisons(old, new)
    assert pairs["metrics"].auroc_change.iloc[0] == pytest.approx(-.02)
    assert summary.scored_runs.tolist() == [1, 0]
    assert summary.protection_lost_runs.tolist() == [1, 0]
    assert not summary.promotion_eligible.any()


@pytest.mark.parametrize("failure", ["task", "duplicate", "counts", "eligibility", "shape"])
def test_changed_denominators_or_ambiguous_shapes_rejected(failure):
    old, new = tables(), tables()
    if failure == "task":
        new["coverage"] = new["coverage"].iloc[:1]
    elif failure == "duplicate":
        new["metrics"] = pd.concat([new["metrics"], new["metrics"]])
    elif failure == "counts":
        new["metrics"].loc[0, "n_positive"] = 199
    elif failure == "eligibility":
        new["coverage"].loc[0, "status"] = "underpowered"
    else:
        new["shapes"]["canonical_protection"] = "unknown"
    with pytest.raises(ValueError):
        m.comparisons(old, new)
