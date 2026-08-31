from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from analyze_label_free_factorial_selection import (  # noqa: E402
    ARMS,
    factorial_effects,
    select_development_oracle_common,
    select_unlabeled_common,
    selection_rule_choices,
)


def evaluation_fixture() -> pd.DataFrame:
    values = {
        "MT_SELMA10_4m4": 0.58,
        "MT_SELMA10_4m5": 0.55,
        "MT_LOG81_4m4": 0.51,
        "MT_LOG81_4m5": 0.50,
    }
    rows = []
    for candidate, bonus in (("candidate_a", 0.0), ("candidate_b", -0.10)):
        for arm in ARMS:
            auroc = values[arm] + bonus
            rows.append(
                {
                    "cell": "K562",
                    "tf": "TF1",
                    "motif_family": "F1",
                    "bias_configuration": arm,
                    "candidate_id": candidate,
                    "family": "count",
                    "status": "ok",
                    "converged": True,
                    "selection_score": auroc,
                    "auroc": auroc,
                    "auprc": auroc - 0.1,
                    "brier": 1.0 - auroc,
                    "prevalence": 0.2,
                }
            )
    return pd.DataFrame(rows)


def cv_fixture() -> pd.DataFrame:
    rows = []
    for candidate, gain in (("candidate_a", 0.001), ("candidate_b", 0.005)):
        for index, arm in enumerate(ARMS):
            rows.append(
                {
                    "cell": "K562",
                    "tf": "TF1",
                    "motif_family": "F1",
                    "bias_configuration": arm,
                    "candidate_id": candidate,
                    "status": "ok",
                    "converged": True,
                    "profile_plausible": True,
                    "heldout_gain_per_site": gain * 100,
                    "heldout_gain_per_site_position": gain + index * 1e-5,
                    "heldout_posterior_entropy": 0.5,
                    "profile_depletion": 0.5,
                }
            )
    return pd.DataFrame(rows)


def test_balanced_factorial_contributions_sum_to_diagonal_contrast() -> None:
    metrics = evaluation_fixture()
    settings = select_development_oracle_common(metrics)
    assert settings.loc[0, "candidate_id"] == "candidate_a"
    effects = factorial_effects(metrics, settings)
    auroc = effects[effects["metric"].eq("auroc")].iloc[0]
    assert np.isclose(auroc["shift_contribution"], 0.02)
    assert np.isclose(auroc["bias_model_contribution"], 0.06)
    assert np.isclose(auroc["total_improvement"], 0.08)
    assert np.isclose(
        auroc["shift_contribution"] + auroc["bias_model_contribution"],
        auroc["total_improvement"],
    )


def test_unlabeled_common_selection_uses_heldout_likelihood() -> None:
    selected = select_unlabeled_common(cv_fixture())
    assert selected.loc[0, "candidate_id"] == "candidate_b"
    assert not bool(selected.loc[0, "selection_labels_used"])


def test_label_free_selection_rules_are_invariant_to_evaluation_labels() -> None:
    cv = cv_fixture()
    evaluation = evaluation_fixture()
    merged = cv.merge(
        evaluation,
        on=[
            "cell",
            "tf",
            "motif_family",
            "bias_configuration",
            "candidate_id",
        ],
        suffixes=("_cv", ""),
    )
    first = selection_rule_choices(merged)
    perturbed = merged.copy()
    perturbed["auroc"] = perturbed["auroc"].iloc[::-1].to_numpy()
    perturbed["auprc"] = 1.0 - perturbed["auprc"]
    perturbed["selection_score"] = -perturbed["selection_score"]
    second = selection_rule_choices(perturbed)
    keys = ["selection_rule", "bias_configuration", "candidate_id"]
    assert first[keys].sort_values(keys).reset_index(drop=True).equals(
        second[keys].sort_values(keys).reset_index(drop=True)
    )
    assert not first["selection_labels_used"].astype(bool).any()
