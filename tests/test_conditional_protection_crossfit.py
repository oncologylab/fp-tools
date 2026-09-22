from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'benchmarks/scripts'))
from crossfit_conditional_protection import select_geometry, prior


def test_selection_preserves_underpowered_and_nonprotective_fallbacks():
    width = len(prior.GEOMETRIES)
    assert select_geometry(np.zeros((2, 201)), np.zeros((2, width)), np.array([0, 1])) == (0, 'underpowered_fallback')
    labels = np.repeat([0, 1], 200)
    assert select_geometry(np.zeros((400, 201)), np.zeros((400, width)), labels) == (0, 'no_protection_fallback')


def test_selection_uses_count_scores_only_after_protection_check():
    labels = np.repeat([0, 1], 200)
    profiles = np.zeros((400, 201))
    profiles[200:, abs(prior.POSITIONS) <= 3] = -1
    scores = np.zeros((400, len(prior.GEOMETRIES)))
    scores[:, 2] = labels
    assert select_geometry(profiles, scores, labels) == (2, 'selected')
