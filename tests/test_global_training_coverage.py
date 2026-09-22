from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import summarize_global_training_coverage as coverage
from global_parametric_expectation import expected_track


def test_prefix_sums_match_explicit_neighborhood_and_annulus():
    x = np.random.default_rng(2026).integers(0, 5, (3, 401), dtype=np.uint32)
    full, annulus = coverage.neighborhood_sum(x), coverage.annular_sum(x)
    for i in range(201):
        np.testing.assert_array_equal(full[:, i], x[:, i:i+201].sum(axis=1))
        np.testing.assert_array_equal(annulus[:, i], x[:, i:i+51].sum(axis=1)+x[:, i+150:i+201].sum(axis=1))


def fixture():
    sites = pd.DataFrame(dict(cell=['HepG2']*3, tf=['A', 'A', 'B']))
    plus, minus = np.ones((3, 401)), np.ones((3, 401))
    bias = np.zeros((3, 401)); bias[0, 390] = np.nan
    plus[1] = 0
    valid = np.array([expected_track(a, b)[1][100:301] & expected_track(c, b)[1][100:301]
                      for a, b, c in zip(plus, bias, minus)])
    return sites, dict(plus_observed=plus, minus_observed=minus, plus_log_bias=bias, minus_log_bias=bias,
        position_valid=valid, valid=valid.all(axis=1), bam_valid=np.ones(3, bool), observed=(plus+minus)[:, 100:301])


def test_missing_outer_positions_do_not_hide_supported_fixed_geometry():
    sites, arrays = fixture()
    rows, positions = coverage.summarize(sites, arrays)
    a = rows.set_index('tf').loc['A']
    assert a.sites == a.positive_observed_sites == 2
    assert a.full_profile_sites == 0 and a.conventional_geometry_sites == 1
    assert a.central_empty_plus_background == 1 and a.central_empty_minus_background == 0
    assert len(positions) == 402


def test_stored_validity_and_label_columns_are_checked():
    sites, arrays = fixture()
    arrays['position_valid'][0, 200] = True
    with pytest.raises(AssertionError):
        coverage.summarize(sites, arrays)
    sites, arrays = fixture(); sites['chip_label'] = 0
    with pytest.raises(ValueError, match='labels'):
        coverage.summarize(sites, arrays)
