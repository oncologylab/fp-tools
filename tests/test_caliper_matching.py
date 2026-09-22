import itertools
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.optimize import linear_sum_assignment

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
from caliper_matching import caliper_assignment


def test_rematching_retains_support_splits_and_power_status():
    import json
    import pandas as pd
    import rematch_corrected_shape_cohorts as m
    p = json.loads(m.base.PROTOCOL.read_text())
    p["minimum_sites_per_class"] = 2
    sites = pd.DataFrame(dict(cell=["HepG2"]*8, tf=["A"]*8, TFBS_chr=["chr1"]*4+["chr16"]*4,
        chip_label=[1,1,0,0]*2, motif_score=[1,2,1,2]*2))
    support = np.ones(8,dtype=bool)
    support[7] = False
    pairs, rows = m.rematch(sites,np.ones(8),support,p)
    assert rows.status.tolist() == ["eligible","underpowered"]
    assert 7 not in pairs.site_index.values
    assert (pairs.match_distance <= .5).all()
    for _, group in pairs.groupby(["split","pair_id"]):
        assert sorted(group.chip_label)==[0,1]
    sites.loc[0,"TFBS_chr"]="chr19"
    with pytest.raises(ValueError):m.rematch(sites,np.ones(8),support,p)


def test_posthoc_caliper_can_discard_every_usable_pair():
    d = np.array([[0, .36], [.36, 1.]])
    r, c = linear_sum_assignment(d)
    assert not (d[r, c] <= .25).any()
    r, c = caliper_assignment(d, .5)
    assert r.tolist() == [0] and c.tolist() == [0]


@pytest.mark.parametrize("shape", [(2, 3), (3, 2), (3, 3)])
def test_cardinality_and_cost_match_exhaustive_oracle(shape):
    rng = np.random.default_rng(2026)
    for _ in range(20):
        d = rng.random(shape)
        answers = []
        for assignments in itertools.product(range(-1, shape[1]), repeat=shape[0]):
            cols = [c for c in assignments if c >= 0]
            if len(cols) != len(set(cols)):
                continue
            edges = [(r, c) for r, c in enumerate(assignments) if c >= 0]
            if any(d[r, c] > .25 for r, c in edges):
                continue
            answers.append((-len(edges), sum(d[r, c] for r, c in edges)))
        expected = min(answers)
        r, c = caliper_assignment(d, .5)
        assert -len(r) == expected[0]
        np.testing.assert_allclose(d[r, c].sum(), expected[1])


def test_empty_and_no_allowed_edges():
    assert len(caliper_assignment(np.empty((0, 2)), .5)[0]) == 0
    assert len(caliper_assignment(np.ones((2, 3)), .5)[0]) == 0


@pytest.mark.parametrize("distance", [np.array([[np.nan]]), np.array([[-1.]]), np.ones(3)])
def test_invalid_distances_rejected(distance):
    with pytest.raises(ValueError):
        caliper_assignment(distance, .5)
