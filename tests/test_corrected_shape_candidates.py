from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import corrected_shape_candidates as m


def fixture(enrichment=False):
    rng = np.random.default_rng(2026)
    labels = np.tile([0]*20 + [1]*20, 2)
    sites = pd.DataFrame(dict(cell=["HepG2"]*80, tf=["A"]*40+["B"]*40,
        motif_family=["family"]*80, TFBS_chr=["chr1"]*80, chip_label=labels))
    values = rng.normal(0, 0.1, (80, 201))
    values[np.ix_(labels == 1, np.abs(m.POSITIONS) <= 7)] += 2 if enrichment else -2
    return sites, values


def test_training_only_and_finite_inputs():
    sites, values = fixture()
    sites.loc[0, "TFBS_chr"] = "chr16"
    with pytest.raises(ValueError, match="training chromosomes"):
        m.fit_candidates(sites, values, minimum_per_class=10)
    sites.loc[0, "TFBS_chr"] = "chr19"
    with pytest.raises(ValueError, match="training chromosomes"):
        m.fit_candidates(sites, values, minimum_per_class=10)
    sites.loc[0, "TFBS_chr"] = "chr1"
    values[0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        m.fit_candidates(sites, values, minimum_per_class=10)


def test_enrichment_is_not_a_footprint_candidate():
    sites, values = fixture(enrichment=True)
    models, rows = m.fit_candidates(sites, values, minimum_per_class=10)
    assert not models
    assert rows.status.eq("no_training_central_protection").all()


def test_family_prior_excludes_own_task_and_scores_are_frozen():
    sites, values = fixture()
    models, rows = m.fit_candidates(sites, values, minimum_per_class=10)
    assert len(models) == 2
    assert rows.prior_scope.eq("family").all()
    assert rows.iloc[0].prior_tasks == [["HepG2", "B"]]
    before = models[("HepG2", "A")]["spline"].discriminant_.copy()
    scores = m.score_candidates(models[("HepG2", "A")], values[:5])
    assert set(scores) == {"conventional", "trained_geometry", "hierarchical_spline"}
    m.score_candidates(models[("HepG2", "A")], values[20:30])
    np.testing.assert_array_equal(before, models[("HepG2", "A")]["spline"].discriminant_)


def test_underpowered_tasks_remain_in_denominator():
    sites, values = fixture()
    models, rows = m.fit_candidates(sites, values, minimum_per_class=200)
    assert not models and len(rows) == 2
    assert rows.status.eq("underpowered").all()
    assert set(m.score_candidates(None, values)) == {"conventional"}


def test_conventional_geometry_has_original_score_definition():
    _, values = fixture()
    expected = values[:, (np.abs(m.POSITIONS) >= 12) & (np.abs(m.POSITIONS) <= 40)].mean(axis=1) - values[:, np.abs(m.POSITIONS) <= 7].mean(axis=1)
    np.testing.assert_allclose(m.score_candidates(None, values)["conventional"], expected)
