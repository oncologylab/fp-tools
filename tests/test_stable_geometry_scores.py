from pathlib import Path
import sys

import numpy as np
import pytest
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
from stable_geometry_scores import stable_scores, rank_auc
from crossfit_stable_geometry import fit_geometry, geometry_contrast
from corrected_shape_candidates import GEOMETRIES, POSITIONS, geometry_weights
from fixed_shape_normalization import normalize_fixed


def test_sparse_flat_profiles_solver_and_batch_parity():
    rng = np.random.default_rng(2026)
    values = rng.poisson(.05, size=(200, 201)).astype(float)
    values[:4] = np.array([np.zeros(201), np.ones(201), 2 + .01 * POSITIONS, np.zeros(201)])
    outer = abs(POSITIONS) >= 60
    design = np.column_stack([np.ones(201), POSITIONS])
    fit = np.linalg.lstsq(design[outer], values[:, outer].T, rcond=None)[0]
    residual = values - (design @ fit).T
    scale = np.sqrt(np.mean(residual[:, outer] ** 2, axis=1))
    alternative = np.clip(residual / np.where(scale > 1e-8, scale, 1)[:, None], -12, 12)
    normal = normalize_fixed(values, POSITIONS, 1.)
    zero = values.sum(1) == 0
    for g in GEOMETRIES:
        w = geometry_weights(g)
        expected = stable_scores(normal, w, zero)
        np.testing.assert_array_equal(expected, stable_scores(alternative, w, zero))
        individual = [stable_scores(normal[i:i+1], w, zero[i:i+1])[0] for i in range(len(values))]
        np.testing.assert_array_equal(expected, individual)
        np.testing.assert_array_equal(expected[:4], 0)


def test_precision_preserves_resolvable_scores_and_exact_tied_auc():
    values = np.array([[0.], [1e-14], [-1e-14], [1e-8]])
    scores = stable_scores(values, np.ones(1), np.zeros(4, bool))
    np.testing.assert_array_equal(scores, [0, 0, 0, 1e-8])
    y = np.array([0, 1, 1, 0])
    assert rank_auc(y, scores) == roc_auc_score(y, scores)
    with pytest.raises(ValueError, match='both classes'):
        rank_auc(np.zeros(4), scores)


def test_tiny_contrast_cannot_pass_protection_and_true_signal_selects():
    values = np.zeros((4, 201)); y = np.array([0, 0, 1, 1])
    values[2:, 97:104] = -1e-14
    assert not geometry_contrast(values[2], GEOMETRIES[0])[3]
    assert fit_geometry(values, y, np.zeros(4, bool), minimum=2)[1] == 'no_protection_fallback'
    values[2:, 97:104] = -1
    assert fit_geometry(values, y, np.zeros(4, bool), minimum=2) == (GEOMETRIES[0], 'selected')
