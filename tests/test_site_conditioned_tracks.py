from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
from audit_site_conditioned_tracks import adjacent_overlaps
from build_strand_functional_profiles import _expected_from_scores


def fixture():
    sites = pd.DataFrame(dict(TFBS_chr=['chr1']*2, TFBS_start=[2, 4], TFBS_end=[3, 5],
        TFBS_strand=['+', '+'], tf=['A', 'B']))
    observed = np.array([[10, 0, 0, 2, 0], [0, 2, 0, 0, 0]], dtype=float)
    expected = np.array([_expected_from_scores(row, np.zeros(5), np.ones(5, bool)) for row in observed])
    return sites, observed, expected


def test_local_total_normalization_cannot_always_be_a_single_track():
    sites, observed, expected = fixture()
    row = adjacent_overlaps(sites, observed, expected, [True, True], flank=2).iloc[0]
    assert row.observed_disagreements == 0
    assert row.overlap_positions == row.expected_disagreements == 3
    assert row.maximum_expected_difference == pytest.approx(2.)


def test_reverse_motif_orientation_is_restored():
    sites, observed, expected = fixture()
    before = adjacent_overlaps(sites, observed, expected, [True, True], flank=2)
    sites.loc[1, 'TFBS_strand'] = '-'
    observed[1], expected[1] = observed[1, ::-1], expected[1, ::-1]
    after = adjacent_overlaps(sites, observed, expected, [True, True], flank=2)
    pd.testing.assert_frame_equal(before, after)


def test_invalid_rows_and_nonoverlapping_windows_are_not_compared():
    sites, observed, expected = fixture()
    assert adjacent_overlaps(sites, observed, expected, [True, False], flank=2).empty
    sites.loc[1, ['TFBS_start', 'TFBS_end']] = [20, 21]
    assert adjacent_overlaps(sites, observed, expected, [True, True], flank=2).empty


def test_a_consistent_expected_track_has_no_conflict():
    sites, observed, expected = fixture()
    expected[:] = 1
    row = adjacent_overlaps(sites, observed, expected, [True, True], flank=2).iloc[0]
    assert row.expected_disagreements == 0


def test_bad_shape_or_nonfinite_input_rejected():
    sites, observed, expected = fixture()
    with pytest.raises(ValueError, match='shape'):
        adjacent_overlaps(sites, observed[:, :3], expected, [True, True], flank=2)
    expected[0, 0] = np.nan
    with pytest.raises(ValueError, match='nonfinite'):
        adjacent_overlaps(sites, observed, expected, [True, True], flank=2)
