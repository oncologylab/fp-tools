from __future__ import annotations

from pathlib import Path
import sys

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_bias_motif_leakage import (  # noqa: E402
    score_sequence_profiles,
    summarize_response,
)
from fp_tools.tools.parametric_bias import (  # noqa: E402
    BiasFeatureSpec,
    ConditionalSequenceBiasModel,
)


def test_vectorized_sequence_profile_scoring() -> None:
    model = ConditionalSequenceBiasModel(BiasFeatureSpec.selma10())
    model.main[4, 0] = 2.0
    model.main[4] -= model.main[4].mean()
    sequences = np.asarray(["ACGT" * 40, "TGCA" * 40])
    profiles = score_sequence_profiles(model, sequences, width=51, margin=41, batch_size=1)
    assert profiles.shape == (2, 51)
    assert np.isfinite(profiles).all()
    assert np.std(profiles[0]) > 0


def test_response_summary_flags_broad_reproducible_effect() -> None:
    rng = np.random.default_rng(9)
    positions = np.arange(-100, 101)
    motif = rng.normal(scale=0.05, size=(200, len(positions)))
    control = rng.normal(scale=0.05, size=(200, len(positions)))
    motif[:, np.abs(positions) <= 15] -= 0.8
    summary, curves = summarize_response(
        motif,
        control,
        positions,
        bootstraps=100,
        seed=4,
        review_threshold=0.25,
    )
    assert summary["potential_motif_response_requires_review"]
    assert summary["center_flank_log_bias_effect"] > 0.5
    assert len(curves) == len(positions)
    assert curves.loc[curves["position"] == 0, "response"].iloc[0] < -0.5
