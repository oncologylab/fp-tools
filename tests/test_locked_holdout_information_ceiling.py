from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "diagnose_locked_holdout_information_ceiling.py"
spec = importlib.util.spec_from_file_location(
    "diagnose_locked_holdout_information_ceiling", SCRIPT
)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_cross_chromosome_predictions_are_oof_and_deterministic():
    rng = np.random.default_rng(4)
    chromosomes = np.repeat(["chr19", "chr20", "chr21", "chr22"], 40)
    feature = rng.normal(size=len(chromosomes))
    labels = (feature + rng.normal(scale=0.2, size=len(feature)) > 0).astype(int)
    features = np.column_stack((feature, rng.normal(size=len(feature))))
    first, folds = module.cross_chromosome_predictions(
        features, labels, chromosomes, seed=7
    )
    second, second_folds = module.cross_chromosome_predictions(
        features, labels, chromosomes, seed=7
    )
    assert folds == second_folds == 4
    assert np.isfinite(first).all()
    assert np.array_equal(first, second)
    assert module.scored_metrics(labels, first)["auroc"] > 0.9


def test_artifact_channels_exposes_strand_channels_only_when_available():
    combined = type("Artifact", (), {"profiles": {"combined_residual": np.ones((2, 3))}})()
    strand = type(
        "Artifact",
        (),
        {
            "profiles": {
                "combined_residual": np.ones((2, 3)),
                "shared_strand_residual": np.ones((2, 3)),
                "antisymmetric_strand_residual": np.ones((2, 3)),
            }
        },
    )()
    assert module.artifact_channels(combined) == ["combined_residual"]
    assert module.artifact_channels(strand) == [
        "combined_residual",
        "shared_strand_residual",
        "antisymmetric_strand_residual",
    ]
