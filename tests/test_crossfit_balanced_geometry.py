from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
from crossfit_balanced_geometry import fitting_mask, fit_geometry, geometry_contrast, GEOMETRIES


def test_fitting_excludes_chromosome_and_cross_chromosome_pair_mates():
    sites = pd.DataFrame(dict(TFBS_chr=['chr1', 'chr2', 'chr2', 'chr3'], pair_id=[0, 0, 1, 1]))
    np.testing.assert_array_equal(fitting_mask(sites, 'chr1'), [False, False, True, True])


def test_geometry_requires_protection_and_retains_underpowered_fallback():
    y = np.array([0, 0, 1, 1])
    values = np.zeros((4, 201))
    values[2:, 97:104] = 1
    assert fit_geometry(values, y, np.zeros(4, bool), minimum=2)[1] == 'no_protection_fallback'
    values[2:, 97:104] = -1
    g, reason = fit_geometry(values, y, np.zeros(4, bool), minimum=2)
    assert reason == 'selected' and g == GEOMETRIES[0]  # exact AUROC ties keep first
    assert geometry_contrast(values[2], g)[3]
    assert fit_geometry(values, y, np.zeros(4, bool))[1] == 'underpowered_fallback'


def test_held_labels_cannot_change_fitting_geometry():
    sites = pd.DataFrame(dict(TFBS_chr=['chr1', 'chr1'] + ['chr2'] * 4, pair_id=[0, 0, 1, 1, 2, 2]))
    y = np.array([0, 1, 0, 1, 0, 1]); values = np.zeros((6, 201))
    values[y == 1, 95:106] = -1
    fit = fitting_mask(sites, 'chr1').to_numpy()
    before = fit_geometry(values[fit], y[fit], np.zeros(fit.sum(), bool), minimum=2)
    y[:2] = y[:2][::-1]; values[:2] = 100
    assert fit_geometry(values[fit], y[fit], np.zeros(fit.sum(), bool), minimum=2) == before
