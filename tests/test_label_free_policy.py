from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from freeze_tf_dependent_label_free_policy import freeze_policy, map_holdout_routes  # noqa: E402


def _matrix() -> pd.DataFrame:
    rows = []
    for family, tasks, gain, ap_gain in (
        ("PASS", [("A", "T1"), ("B", "T2")], 0.06, 0.08),
        ("FAIL", [("A", "T3"), ("B", "T4")], -0.01, -0.02),
    ):
        for cell, tf in tasks:
            for bias, candidate, auc, ap in (
                ("DWM", "dwm", 0.60, 0.40),
                ("NEW", "new", 0.60 + gain, 0.40 + ap_gain),
            ):
                rows.append(
                    {
                        "cell": cell,
                        "tf": tf,
                        "motif_family": family,
                        "bias_configuration": bias,
                        "candidate_id": candidate,
                        "selection_score": auc + ap,
                        "auroc": auc,
                        "auprc": ap,
                    }
                )
    return pd.DataFrame(rows)


def test_policy_promotes_only_stable_family_gain() -> None:
    policy, contexts, global_dwm = freeze_policy(
        _matrix(),
        minimum_contexts=2,
        minimum_mean_auroc_gain=0.03,
        minimum_relative_auprc_gain=0.10,
        maximum_context_auroc_loss=0.02,
    )
    indexed = policy.set_index("motif_family")
    assert bool(indexed.loc["PASS", "new_route_passes_development_gates"])
    assert indexed.loc["PASS", "recommended_bias_configuration"] == "NEW"
    assert not bool(indexed.loc["FAIL", "new_route_passes_development_gates"])
    assert indexed.loc["FAIL", "recommended_bias_configuration"] == "DWM"
    assert len(contexts) == 4
    assert global_dwm["bias_configuration"] == "DWM"


def test_unseen_holdout_family_fails_closed_to_dwm() -> None:
    policy, _contexts, global_dwm = freeze_policy(
        _matrix(),
        minimum_contexts=2,
        minimum_mean_auroc_gain=0.03,
        minimum_relative_auprc_gain=0.10,
        maximum_context_auroc_loss=0.02,
    )
    study = {
        "tasks": [
            {
                "cell": "H",
                "tf": "KNOWN",
                "motif_id": "M1",
                "motif_family": "PASS",
                "role": "difficult",
                "split": "locked_holdout",
            },
            {
                "cell": "H",
                "tf": "UNKNOWN",
                "motif_id": "M2",
                "motif_family": "UNSEEN",
                "role": "difficult",
                "split": "locked_holdout",
            },
        ]
    }
    routes = map_holdout_routes(study, policy, global_dwm).set_index("tf")
    assert routes.loc["KNOWN", "bias_configuration"] == "NEW"
    assert routes.loc["UNKNOWN", "bias_configuration"] == "DWM"
    assert routes.loc["UNKNOWN", "route_source"] == "unseen_family_dwm_fallback"
