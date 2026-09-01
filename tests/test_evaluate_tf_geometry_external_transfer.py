import importlib.util
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT / "benchmarks" / "scripts" / "evaluate_tf_geometry_external_transfer.py"
)
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location(
    "evaluate_tf_geometry_external_transfer", SCRIPT
)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_conventional_profile_scores_returns_center_score():
    profiles = np.zeros((2, 201), dtype=float)
    profiles[0, :50] = 2.0
    profiles[0, 151:] = 2.0
    scores = module.conventional_profile_scores(profiles, 100)
    assert scores.shape == (2,)
    assert scores[0] > scores[1]


def test_chromosome_block_bootstrap_is_deterministic():
    labels = np.tile([0, 1], 50)
    conventional = labels * 0.2 + np.linspace(0, 0.1, len(labels))
    candidate = labels * 0.8 + np.linspace(0, 0.1, len(labels))
    chromosomes = np.repeat(["chr19", "chr20", "chr21", "chr22", "chrX"], 20)
    first = module.chromosome_block_bootstrap(
        labels,
        conventional,
        candidate,
        chromosomes,
        iterations=50,
        seed=7,
    )
    second = module.chromosome_block_bootstrap(
        labels,
        conventional,
        candidate,
        chromosomes,
        iterations=50,
        seed=7,
    )
    assert first == second
    assert first["bootstrap_successful"] == 50
