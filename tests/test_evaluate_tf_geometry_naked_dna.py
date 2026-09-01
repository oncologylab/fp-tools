import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "evaluate_tf_geometry_naked_dna.py"
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("evaluate_tf_geometry_naked_dna", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_frozen_upper_tail_threshold_respects_call_budget_with_ties():
    values = np.array([9.0, 9.0, 8.0, 7.0, 6.0] + list(range(15)))
    threshold = module.frozen_upper_tail_threshold(values, 0.10)
    assert np.sum(values >= threshold) <= 2


def test_frozen_upper_tail_threshold_rejects_empty_and_bad_alpha():
    with pytest.raises(ValueError, match="zero finite"):
        module.frozen_upper_tail_threshold(np.array([np.nan]), 0.05)
    with pytest.raises(ValueError, match="between zero and one"):
        module.frozen_upper_tail_threshold(np.array([1.0]), 1.0)
