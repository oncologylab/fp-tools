import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "benchmarks" / "scripts"


def load(name):
    path = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


combined_builder = load("build_combined_functional_profiles")
strand_builder = load("build_strand_functional_profiles")
pooler = load("pool_functional_profile_artifacts")


def sites():
    return pd.DataFrame(
        {
            "cell": "Fixture",
            "tf": ["A", "A"],
            "motif_id": "M",
            "motif_family": "F",
            "TFBS_chr": "chr1",
            "TFBS_start": [100, 200],
            "TFBS_end": [110, 210],
            "TFBS_strand": ["+", "-"],
            "motif_score": 1.0,
            "peak_start": [50, 150],
            "peak_end": [150, 250],
            "chromosome_split": "train",
        }
    )


def combined_artifact(tmp_path, name, observed):
    prefix = tmp_path / name
    expected = np.ones_like(observed)
    _npz, document, _sites = combined_builder.write_artifact(
        prefix,
        sites(),
        observed,
        expected,
        np.array([True, True]),
        dispersion=0.0,
        metadata={"labels_used": False},
    )
    return document


def test_combined_mean_pool_recomputes_residual(tmp_path):
    first = combined_artifact(tmp_path, "r1", np.ones((2, 9)))
    second = combined_artifact(tmp_path, "r2", np.full((2, 9), 3.0))
    npz, document, _sites = pooler.pool_artifacts(
        [first, second], tmp_path / "pooled", mode="mean", dispersion=0.0
    )
    with np.load(npz, allow_pickle=False) as payload:
        np.testing.assert_allclose(payload["observed"], 2.0)
        assert payload["combined_residual"].shape == (2, 9)
    metadata = json.loads(document.read_text(encoding="utf-8"))["metadata"]
    assert metadata["replicate_weights"] == [0.5, 0.5]
    assert metadata["labels_used"] is False


def test_library_equalized_pool_gives_replicates_equal_total_weight(tmp_path):
    first = combined_artifact(tmp_path, "r1", np.ones((2, 9)))
    second = combined_artifact(tmp_path, "r2", np.full((2, 9), 4.0))
    npz, document, _sites = pooler.pool_artifacts(
        [first, second],
        tmp_path / "pooled",
        mode="library-equalized-mean",
        dispersion=0.0,
    )
    metadata = json.loads(document.read_text(encoding="utf-8"))["metadata"]
    weights = np.asarray(metadata["replicate_weights"])
    totals = np.asarray([18.0, 72.0])
    np.testing.assert_allclose(weights * totals, [22.5, 22.5])
    with np.load(npz, allow_pickle=False) as payload:
        assert np.isclose(payload["observed"].sum(), 45.0)


def test_pool_rejects_changed_site_order(tmp_path):
    first = combined_artifact(tmp_path, "r1", np.ones((2, 9)))
    second = combined_artifact(tmp_path, "r2", np.ones((2, 9)))
    document = json.loads(second.read_text(encoding="utf-8"))
    npz = Path(document["profiles_npz"])
    with np.load(npz, allow_pickle=False) as payload:
        arrays = {key: payload[key] for key in payload.files}
    arrays["site_hash"] = arrays["site_hash"][::-1]
    np.savez_compressed(npz, **arrays)
    document["profiles_sha256"] = pooler.file_sha256(npz)
    second.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(ValueError, match="site hashes or order differ"):
        pooler.pool_artifacts(
            [first, second], tmp_path / "pooled", mode="sum", dispersion=0.0
        )
