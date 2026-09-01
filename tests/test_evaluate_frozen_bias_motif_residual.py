from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "evaluate_frozen_bias_motif_residual.py"
spec = importlib.util.spec_from_file_location(
    "evaluate_frozen_bias_motif_residual", SCRIPT
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_residual_summary_flags_unexplained_center_shape() -> None:
    positions = np.arange(-100, 101)
    expected = np.full((100, len(positions)), 5.0)
    observed = expected.copy()
    observed[:, np.abs(positions) <= 15] = 2.0
    arrays = {
        "plus_observed": observed / 2,
        "minus_observed": observed / 2,
        "plus_expected": expected / 2,
        "minus_expected": expected / 2,
        "valid": np.ones(100, dtype=bool),
    }
    sites = pd.DataFrame(
        {
            "tf": np.repeat("TF1", 100),
            "motif": np.repeat("MA0001.1", 100),
            "motif_family": np.repeat("family", 100),
        }
    )
    summary, curves = module.summarize_artifact(
        "candidate",
        arrays,
        sites,
        threshold=0.25,
        bootstraps=100,
        seed=2026,
    )
    assert bool(summary.loc[0, "motif_residual_flag"])
    assert len(curves) == len(positions)


def test_profile_loader_rejects_chipped_labels(tmp_path) -> None:
    prefix = tmp_path / "profiles"
    arrays = {
        "plus_observed": np.ones((2, 201)),
        "minus_observed": np.ones((2, 201)),
        "plus_expected": np.ones((2, 201)),
        "minus_expected": np.ones((2, 201)),
        "valid": np.ones(2, dtype=bool),
    }
    np.savez_compressed(prefix.with_suffix(".npz"), **arrays)
    sites = pd.DataFrame({"tf": ["TF", "TF"], "chip_label": [0, 1]})
    sites.to_csv(Path(str(prefix) + ".sites.tsv.gz"), sep="\t", index=False)
    document = {
        "profiles_sha256": module.file_sha256(prefix.with_suffix(".npz")),
        "sites_sha256": module.file_sha256(Path(str(prefix) + ".sites.tsv.gz")),
    }
    prefix.with_suffix(".json").write_text(__import__("json").dumps(document))
    try:
        module.load_profile_artifact(prefix)
    except ValueError as error:
        assert "cannot contain ChIP labels" in str(error)
    else:
        raise AssertionError("ChIP labels were accepted")
