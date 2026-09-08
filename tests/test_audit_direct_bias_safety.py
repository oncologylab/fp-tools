from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"))


def audit_module():
    import audit_direct_bias_safety
    return audit_direct_bias_safety


@pytest.mark.parametrize("dispersion", [0.0, 1.6144318781609228])
def test_direct_safety_scores_equal_frozen_comparison_operator(dispersion):
    import evaluate_chrombpnet_reference as reference
    audit = audit_module()
    positions = np.arange(-100, 101)
    counts = np.ones((3, 201)) * np.asarray([1.0, 5.0, 10.0])[:, None]
    counts[:, np.abs(positions) <= 7] *= np.asarray([0.2, 1.0, 3.0])[:, None]
    expected = np.ones_like(counts) * 5
    frozen = reference.score_methods(counts, expected, expected, expected,
                                     counts, counts, positions, dispersion=dispersion)
    scores, scaled, residuals = audit.direct_scores(counts, expected, dispersion)
    np.testing.assert_allclose(scores, frozen[reference_method()][0])
    np.testing.assert_allclose(scaled.sum(axis=1), counts.sum(axis=1))
    assert scores[0] > 0
    assert scores[1] == pytest.approx(0)
    assert scores[2] < 0


def reference_method():
    return "frozen_parametric_bias_conventional_geometry"


def test_threshold_calibration_never_uses_test_or_positive_labels():
    audit = audit_module()
    sites = pd.DataFrame({"tf": ["CTCF"] * 6, "TFBS_chr": ["chr16"] * 4 + ["chr19", "chr16"],
                          "chip_label": [0, 0, 0, 0, 0, 1]})
    scores = np.asarray([1.0, 2.0, 3.0, 4.0, 999.0, 999.0])
    row = audit.calibrate(sites, scores, np.ones(6, dtype=bool), np.ones(6),
                          target=0.4, minimum_negatives=4).iloc[0]
    assert row.validation_negatives == 4
    assert row.threshold == 3.0
    assert row.threshold_eligible


def test_unsupported_threshold_never_becomes_safety_pass():
    audit = audit_module()
    row = audit.safety_row(np.zeros(100), np.zeros(100), np.ones(100),
                           candidate_threshold=np.nan, baseline_threshold=1.0,
                           threshold_eligible=False)
    assert not row["passes_candidate_screen"]


def test_safety_checks_informative_denominator_wilson_and_dwm_increase():
    audit = audit_module()
    # The only nonzero site is called. Pooling 199 zero-cut sites would hide it.
    score = np.zeros(200); score[0] = 2
    total = np.zeros(200); total[0] = 2
    row = audit.safety_row(score, np.zeros(200), total, candidate_threshold=1,
                           baseline_threshold=1, threshold_eligible=True)
    assert row["false_positive_rate"] == 0.005
    assert row["informative_false_positive_rate"] == 1.0
    assert not row["passes_candidate_screen"]
    # Four calls among 200 are below 5%, but more than 1 percentage point above DWM.
    score = np.zeros(200); score[:4] = 2
    row = audit.safety_row(score, np.zeros(200), np.ones(200), candidate_threshold=1,
                           baseline_threshold=1, threshold_eligible=True)
    assert row["false_positive_rate_increase_over_dwm"] == 0.02
    assert not row["passes_candidate_screen"]


def test_zero_calls_with_sufficient_support_pass_only_candidate_screen():
    audit = audit_module()
    row = audit.safety_row(np.zeros(200), np.zeros(200), np.ones(200),
                           candidate_threshold=0, baseline_threshold=0,
                           threshold_eligible=True)
    assert row["false_positive_calls"] == 0  # Calls are strictly greater than cutoff.
    assert row["false_positive_rate_upper_95"] == pytest.approx(0.0188453264)
    assert row["passes_candidate_screen"]
    assert not row["final_untouched_safety"]


@pytest.mark.parametrize("change", [{"naked_replicate": "NakedDNA_rep3"},
                                    {"comparisons": [{"cell": "SK-N-SH"}]},
                                    {"new_holdout_labels_allowed": True}])
def test_plan_rejects_unopened_controls_before_reading_any_input(tmp_path, change):
    audit = audit_module()
    plan = json.loads((Path(audit.__file__).parents[1] / "manifests/chrombpnet_direct_bias_safety_v1.json").read_text())
    plan.update(change)
    with pytest.raises(ValueError, match="plan"):
        audit.resolve_plan(plan, tmp_path)


def test_control_rep3_is_rejected_before_site_table_access(tmp_path):
    audit = audit_module()
    manifest = tmp_path / "control.json"
    manifest.write_text(json.dumps({"metadata": {"labels_used": False, "flank": 100,
                                                 "bam": "NakedDNA_rep3.bam"}}))
    with pytest.raises(ValueError, match="replicate 2"):
        audit.load_artifact(audit.record(manifest), tmp_path, "K562", control=True)


def test_freezes_and_threshold_checks_fail_closed(tmp_path):
    audit = audit_module()
    table_path = tmp_path / "thresholds.tsv"
    table_path.write_text("cell\ttf\tmethod\tthreshold\nK562\tCTCF\tparametric_direct\t1.0\n")
    frozen = {"configuration_id": "frozen", "thresholds": audit.record(table_path),
              "calibration_chromosomes": ["chr16", "chr17", "chr18"], "naked_profiles_scored": False}
    frozen["threshold_id"] = audit.canonical_id(frozen)
    path = tmp_path / "thresholds.freeze.json"
    audit.write_immutable(path, frozen)
    audit.write_immutable(path, frozen)
    with pytest.raises(ValueError, match="immutable"):
        audit.write_immutable(path, {**frozen, "configuration_id": "changed"})
    _, table = audit.read_thresholds(path, "frozen", tmp_path)
    assert table.threshold.iloc[0] == 1.0
    with pytest.raises(ValueError, match="configuration mismatch"):
        audit.read_thresholds(path, "changed", tmp_path)
    table_path.write_text(table_path.read_text().replace("1.0", "2.0"))
    with pytest.raises(ValueError, match="checksum"):
        audit.read_thresholds(path, "frozen", tmp_path)
