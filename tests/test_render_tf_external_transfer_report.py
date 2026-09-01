import importlib.util
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "render_tf_external_transfer_report.py"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location(
    "render_tf_external_transfer_report", SCRIPT
)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_conventional_score_profiles_agrees_at_center():
    profiles = np.zeros((3, 201), dtype=float)
    profiles[:, :60] = 1.0
    profiles[:, 141:] = 1.0
    scored = module.conventional_score_profiles(profiles)
    centered = module.conventional_profile_scores(profiles, 100)
    assert scored.shape == profiles.shape
    assert np.allclose(scored[:, 100], centered)
