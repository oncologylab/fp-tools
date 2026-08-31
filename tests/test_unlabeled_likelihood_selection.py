from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from select_count_models_by_unlabeled_likelihood import (  # noqa: E402
    chromosome_masks,
    log_marginal_gain,
    select_candidates,
)


def test_log_marginal_gain_has_unbound_reference() -> None:
    neutral = log_marginal_gain(np.zeros(3), np.zeros(3))
    assert np.allclose(neutral, 0.0)
    favorable = log_marginal_gain(np.full(3, 2.0), np.zeros(3))
    unfavorable = log_marginal_gain(np.full(3, -2.0), np.zeros(3))
    assert np.all(favorable > 0)
    assert np.all(unfavorable < 0)


def test_chromosome_tuning_split_is_disjoint() -> None:
    study = {"chromosome_split": {"train": ["chr1", "chr2", "chr3", "chr4"]}}
    sites = pd.DataFrame({"TFBS_chr": ["chr1", "chr2", "chr3", "chr4"]})
    fit, tune, chromosomes = chromosome_masks(sites, study, 2)
    assert chromosomes == ["chr3", "chr4"]
    assert np.array_equal(fit, [True, True, False, False])
    assert np.array_equal(tune, [False, False, True, True])
    assert not np.any(fit & tune)


def test_candidate_selection_uses_no_label_metric() -> None:
    metrics = pd.DataFrame(
        {
            "cell": ["K562", "K562"],
            "tf": ["TF1", "TF1"],
            "status": ["ok", "ok"],
            "converged": [True, True],
            "profile_plausible": [True, True],
            "heldout_gain_per_site_position": [0.01, 0.02],
            "heldout_gain_per_site": [1.0, 2.0],
            "fit_seconds": [1.0, 2.0],
            "bias_configuration": ["A", "B"],
            "candidate_id": ["one", "two"],
            "auroc": [0.99, 0.01],
        }
    )
    selected = select_candidates(metrics)
    assert len(selected) == 1
    assert selected.iloc[0]["candidate_id"] == "two"
