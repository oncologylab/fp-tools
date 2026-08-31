import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "evaluate_common_support_sensitivity.py"
spec = importlib.util.spec_from_file_location("evaluate_common_support_sensitivity", SCRIPT)
common = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = common
spec.loader.exec_module(common)


def fixture_scores():
    rows = []
    for label in (0, 1):
        for index in range(40):
            rows.append(
                {
                    "cell": "Fixture",
                    "tf": "MAX",
                    "label": label,
                    "TFBS_chr": f"chr{19 + index % 2}",
                    "TFBS_start": index * 100 + label,
                    "TFBS_end": index * 100 + label + 6,
                    "TFBS_strand": "+" if index % 2 else "-",
                    "motif_score": 8.0 + (index % 3) * 0.01,
                    "log_accessibility": 2.0 + (index % 8) * 0.2 + label * 0.01,
                    "gc_fraction": 0.4 + (index % 5) * 0.04 + label * 0.001,
                    "peak_position_signed": -0.8 + (index % 10) * 0.16,
                    "peak_position_abs": abs(-0.8 + (index % 10) * 0.16),
                    "artifact_index": index + label * 40,
                    "candidate_probability": index / 40 + label * 0.05,
                    "reference_probability": index / 40,
                }
            )
    return pd.DataFrame(rows)


def test_common_support_selection_does_not_depend_on_model_scores():
    frame = fixture_scores()
    first_grid, first, first_summary = common.evaluate_grid(
        frame,
        accessibility_bins=[2, 3],
        gc_bins=[2, 3],
        motif_bins=2,
        peak_position_bins=2,
        minimum_pairs=10,
        maximum_smd=0.25,
        seed=7,
    )
    changed = frame.copy()
    changed["candidate_probability"] = np.random.default_rng(3).random(len(changed))
    changed["reference_probability"] = np.random.default_rng(4).random(len(changed))
    second_grid, second, second_summary = common.evaluate_grid(
        changed,
        accessibility_bins=[2, 3],
        gc_bins=[2, 3],
        motif_bins=2,
        peak_position_bins=2,
        minimum_pairs=10,
        maximum_smd=0.25,
        seed=7,
    )
    assert first_summary["accessibility_bins"] == second_summary["accessibility_bins"]
    assert first_summary["gc_bins"] == second_summary["gc_bins"]
    assert first["artifact_index"].tolist() == second["artifact_index"].tolist()
    assert len(first_grid) == len(second_grid) == 4


def test_common_support_is_balanced_and_deterministic():
    frame = fixture_scores()
    first = common.coarsened_exact_match(
        frame, accessibility_bins=3, gc_bins=3, peak_position_bins=2, seed=11
    )
    second = common.coarsened_exact_match(
        frame, accessibility_bins=3, gc_bins=3, peak_position_bins=2, seed=11
    )
    assert first["label"].value_counts().to_dict()[0] == first["label"].value_counts().to_dict()[1]
    pd.testing.assert_frame_equal(first, second)
