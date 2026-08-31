import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "build_combined_functional_profiles.py"
spec = importlib.util.spec_from_file_location("build_combined_functional_profiles", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def sites() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell": ["K562", "K562"],
            "tf": ["CTCF", "CTCF"],
            "TFBS_chr": ["chr1", "chr1"],
            "TFBS_start": [100, 200],
            "TFBS_end": [110, 210],
            "TFBS_strand": ["+", "-"],
        }
    )


def test_validate_unlabeled_sites_rejects_label_columns():
    frame = sites()
    frame["chip_label"] = 0
    with pytest.raises(ValueError, match="refuses columns"):
        module.validate_unlabeled_sites(frame, "sites.tsv")


def test_write_artifact_preserves_orientation_ready_arrays(tmp_path):
    frame = sites()
    observed = np.arange(18, dtype=float).reshape(2, 9)
    expected = np.full((2, 9), 2.0)
    valid = np.array([True, False])
    npz_path, json_path, sites_path = module.write_artifact(
        tmp_path / "fixture",
        frame,
        observed,
        expected,
        valid,
        dispersion=0.1,
        metadata={"labels_used": False},
    )
    assert json_path.is_file() and sites_path.is_file()
    with np.load(npz_path, allow_pickle=False) as arrays:
        np.testing.assert_allclose(arrays["observed"], observed)
        np.testing.assert_allclose(arrays["expected"], expected)
        assert arrays["combined_residual"].shape == observed.shape
        assert arrays["valid"].tolist() == [True, False]
