from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))


def module():
    import audit_direct_bias_stability
    return audit_direct_bias_stability


def sites():
    return pd.DataFrame({"cell": ["K562"] * 4, "tf": ["CTCF"] * 4,
                         "motif_family": ["CTCF"] * 4, "TFBS_chr": ["chr19", "chr19", "chr20", "chr20"],
                         "TFBS_start": [10, 30, 50, 70], "TFBS_end": [20, 40, 60, 80],
                         "TFBS_strand": ["+", "-", "+", "-"], "chip_label": [0, 1, 0, 1],
                         "motif_score": [1., 1., 2., 2.], "accessibility": [1., 1., 2., 2.],
                         "site_id": ["same_motif_name"] * 4})


def matrix():
    m = module()
    rows = [(cell, f"{cell}_rep1", depth, seed) for cell in ("HepG2", "K562")
            for depth in ("10m", "25m", "50m") for seed in range(2026, 2031)]
    rows += [(cell, f"{cell}_rep{rep}", "full", 2026) for cell in ("HepG2", "K562") for rep in (1, 2, 3)]
    return pd.DataFrame(rows, columns=m.RUN_KEYS)


def test_matrix_requires_all_declared_runs_without_pseudoreplicated_full_seeds():
    m = module()
    frame = matrix()
    m.validate_matrix(frame)
    for bad in (frame.iloc[:-1], pd.concat([frame, frame.iloc[[0]]]),
                pd.concat([frame, frame.iloc[[-1]].assign(seed=2027)])):
        with pytest.raises(ValueError, match="36 prespecified"):
            m.validate_matrix(bad)


@pytest.mark.parametrize("field,value", [("cells", ["SK-N-SH"]), ("candidate_retuning_allowed", True),
                                        ("evaluation_chromosomes", ["chr16"]),
                                        ("new_holdout_labels_allowed", True)])
def test_scope_rejected_before_any_artifact_access(field, value):
    m = module()
    plan = json.loads((Path(m.__file__).parents[1] / "manifests/direct_bias_stability_v1.json").read_text())
    m.validate_plan(plan)
    plan[field] = value
    with pytest.raises(ValueError, match="development scope"):
        m.validate_plan(plan)


def test_coordinate_join_ignores_nonunique_motif_name_but_enforces_labels():
    m = module()
    frame = sites()
    selected = frame.iloc[[3, 1]].copy()
    np.testing.assert_array_equal(m.align_sites(frame, selected), [3, 1])
    selected.loc[selected.index[0], "chip_label"] = 0
    with pytest.raises(ValueError, match="loci or labels"):
        m.align_sites(frame, selected)
    with pytest.raises(ValueError, match="unique"):
        m.align_sites(pd.concat([frame, frame.iloc[[0]]]), frame)


def test_cached_expectation_preserves_each_strand_total_and_orientation():
    m = module()
    plus = np.asarray([[9., 0., 0.], [2., 0., 0.]])
    minus = np.asarray([[0., 1., 0.], [0., 0., 8.]])
    p_log = np.log([[3., 2., 1.], [3., 2., 1.]])
    m_log = np.log([[1., 2., 3.], [1., 2., 3.]])
    expected = m.expected_for_observed(plus, minus, p_log, m_log)
    np.testing.assert_allclose(expected[0], (9 * np.asarray([3, 2, 1]) + np.asarray([1, 2, 3])) / 6)
    np.testing.assert_allclose(expected.sum(axis=1), [10, 10])
    p_rev, m_rev, _ = m.orient_strand_log_bias(plus, minus, ["+", "-"])
    p_back, m_back, _ = m.orient_strand_log_bias(p_rev, m_rev, ["+", "-"])
    np.testing.assert_array_equal(p_rev[1], minus[1, ::-1])
    np.testing.assert_array_equal(p_back, plus)
    np.testing.assert_array_equal(m_back, minus)


@pytest.mark.parametrize("dispersion", [0., 1.6144318781609228])
def test_exact_operator_matches_recorded_reference_and_retains_zero_cut_sites(dispersion):
    from evaluate_chrombpnet_reference import score_methods
    m = module()
    rng = np.random.default_rng(29)
    observed = rng.poisson(.8, (4, 201)).astype(float)
    observed[0] = 0
    dwm, param, sample = [rng.uniform(.2, 2., observed.shape) for _ in range(3)]
    actual = m.score_panel(observed, param, dwm, sample, dispersion)
    expected = score_methods(observed, dwm, param, dwm, dwm, dwm, m.POSITIONS, dispersion=dispersion)
    for method in (m.RAW, m.DWM, m.PARAMETRIC):
        np.testing.assert_allclose(actual[method][0], expected[method][0])
        np.testing.assert_allclose(actual[method][1], expected[method][1])
        assert len(actual[method][0]) == 4
        assert np.isfinite(actual[method][1]).all()
    sample_score = m.direct_scores(observed, sample, dispersion)[0]
    np.testing.assert_allclose(actual[m.SAMPLE_DWM][0], sample_score)


def test_missing_dwm_expectation_is_not_treated_as_valid_bias_evidence():
    m = module()
    observed = np.ones((2, 201))
    prediction = observed.copy()
    prediction[0] = 0
    with pytest.raises(ValueError, match="no expectation"):
        m.score_panel(observed, observed, prediction, observed, 0)


def test_expected_coverage_keeps_complete_pairs_without_rematching():
    m = module()
    frame = sites()
    strict = frame.copy(); strict["pair_id"] = ["first", "first", "second", "second"]
    retained, cohorts = m.qualify_support(frame, strict, np.asarray([True, False, True, True]))
    np.testing.assert_array_equal(retained, [0, 2, 3])
    np.testing.assert_array_equal(cohorts["original_support_with_expected_coverage"], [0, 1, 2])
    np.testing.assert_array_equal(cohorts["motif_accessibility_observed_total"], [1, 2])


def test_frozen_expected_coverage_resumes_without_extracting_tracks(tmp_path, monkeypatch):
    m = module()
    frame = matrix().query("cell == 'K562'").iloc[[0]].assign(dwm_expected="existing.bw")
    source = sites()
    matched = source.assign(pair_id=["first", "first", "second", "second"], policy="motif_accessibility_observed_total")
    values = np.ones((4, 201)); values[1] = 0
    monkeypatch.setattr(m, "extract_prediction_profiles", lambda *a, **kw: (values, np.ones(4, dtype=bool)))
    first = m.freeze_support("K562", frame, source, {"counts": values}, np.ones_like(values), matched,
                             tmp_path, "frozen", tmp_path)
    def unexpected(*args, **kwargs):
        raise AssertionError("completed expectation extraction was rerun")
    monkeypatch.setattr(m, "extract_prediction_profiles", unexpected)
    again = m.freeze_support("K562", frame, source, {"counts": values}, np.ones_like(values), matched,
                             tmp_path, "frozen", tmp_path)
    pd.testing.assert_frame_equal(first[0], again[0])
    np.testing.assert_array_equal(first[3]["motif_accessibility_observed_total"], [1, 2])


def test_stage_reuse_is_bound_to_configuration_and_leaf_checksum(tmp_path):
    m = module()
    data = tmp_path / "profiles.npz"
    np.savez(data, values=np.arange(5))
    done = {"configuration_id": "frozen", "outputs": {"profiles": m.record(data)}}
    done["stage_id"] = m.canonical_id(done)
    path = tmp_path / "extraction.json"
    m.write_immutable(path, done)
    assert m.verify_stage(path, "frozen", tmp_path) == done
    with pytest.raises(ValueError, match="different frozen configuration"):
        m.verify_stage(path, "changed", tmp_path)
    with pytest.raises(ValueError, match="different sample, depth or seed"):
        m.verify_stage(path, "frozen", tmp_path, expected_run={"sample": "different"})
    np.savez(data, values=np.arange(6))
    with pytest.raises(ValueError, match="checksum"):
        m.verify_stage(path, "frozen", tmp_path)


def test_small_run_keeps_zero_sites_and_both_cohorts_then_reuses_outputs(tmp_path):
    from types import SimpleNamespace
    m = module()
    frame = sites()
    observed = np.ones((4, 201)); observed[0] = 0
    methods = m.score_panel(observed, np.ones_like(observed), np.ones_like(observed), np.ones_like(observed), 0)
    cohorts = {"original_common_support": np.arange(4), "motif_accessibility_observed_total": np.arange(4)}
    plan = {"minimum_sites_per_class": 200, "aggregate_bootstrap_iterations": 10, "seed": 2026}
    row = SimpleNamespace(cell="K562", sample="K562_rep1", depth="25m", seed=2026)
    done = m.evaluate_run(row, frame, cohorts, methods, observed, plan, tmp_path, "frozen", tmp_path)
    metrics = pd.read_csv(tmp_path / "metrics.tsv.gz", sep="\t")
    coverage = pd.read_csv(tmp_path / "coverage.tsv.gz", sep="\t")
    assert len(metrics) == 8 and metrics.n_sites.eq(4).all()
    assert metrics.status.eq("underpowered").all() and coverage.zero_cut_sites.eq(1).all()
    assert m.evaluate_run(row, None, None, None, None, None, tmp_path, "frozen", tmp_path) == done
