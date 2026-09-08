from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import audit_direct_bias_safety_coverage as audit


def test_tf_grouped_score_output_rejoins_original_coordinate_order():
    sites = pd.DataFrame({"cell": ["K562"] * 3, "tf": ["CTCF", "MYC", "CTCF"],
                          "TFBS_chr": ["chr19"] * 3, "TFBS_start": [10, 20, 30],
                          "TFBS_end": [15, 25, 35], "TFBS_strand": ["+"] * 3})
    scores = sites.iloc[[0, 2, 1]].reset_index(drop=True)
    np.testing.assert_array_equal(audit.profile_indexes(sites, np.ones(3, dtype=bool), scores), [0, 2, 1])
    scores.loc[0, "TFBS_start"] = 99
    with pytest.raises(ValueError, match="different motif loci"):
        audit.profile_indexes(sites, np.ones(3, dtype=bool), scores)


def test_normalized_dwm_counts_are_comparable_but_missing_or_shifted_tracks_are_not():
    y = np.zeros((4, 201)); y[:, 100] = 2
    dy = y * 133.944
    de = np.ones_like(y)
    dy[1] = 0; de[1] = 0
    dy[2, 100] = 0; dy[2, 101] = 267.888
    y[3] = 0; dy[3] = 0
    common, disagreement = audit.comparable_profiles(y, dy, de)
    np.testing.assert_array_equal(common, [True, False, False, False])
    np.testing.assert_array_equal(disagreement, [False, False, True, False])


def test_zero_filled_baseline_cannot_certify_candidate_safety():
    row = audit.qualify_row(np.zeros(200), np.zeros(200), np.ones(200),
                            np.zeros(200, dtype=bool), np.zeros(200, dtype=bool),
                            candidate_threshold=1, baseline_threshold=1,
                            threshold_eligible=True, motif_residual_flag=False)
    assert row["absolute_fpr_screen_pass"]
    assert row["paired_support"] == 0
    assert not row["passes_candidate_screen"]
    assert "false_positive_rate_increase_over_dwm" not in row


def test_absolute_and_paired_rates_keep_their_own_denominators():
    y = np.zeros(200); y[0] = 2
    coverage = np.arange(200) < 114
    row = audit.qualify_row(y, np.zeros(200), np.ones(200), coverage, np.zeros(200, dtype=bool),
                            candidate_threshold=1, baseline_threshold=1,
                            threshold_eligible=True, motif_residual_flag=False)
    assert row["false_positive_rate"] == 1 / 200
    assert row["paired_false_positive_rate"] == 1 / 114
    assert row["dwm_unavailable_informative_sites"] == 86
    assert row["passes_candidate_screen"]
    row = audit.qualify_row(y, np.zeros(200), np.ones(200), coverage, np.zeros(200, dtype=bool),
                            candidate_threshold=1, baseline_threshold=1,
                            threshold_eligible=True, motif_residual_flag=True)
    assert not row["passes_candidate_screen"]
