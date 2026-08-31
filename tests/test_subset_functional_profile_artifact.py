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
subsetter = load("subset_functional_profile_artifact")


def fixture_sites():
    return pd.DataFrame(
        {
            "cell": "Fixture",
            "tf": ["A", "A", "A"],
            "motif_id": "M",
            "motif_family": "F",
            "TFBS_chr": ["chr1", "chr16", "chr19"],
            "TFBS_start": [100, 200, 300],
            "TFBS_end": [110, 210, 310],
            "TFBS_strand": ["+", "-", "+"],
            "chromosome_split": ["train", "validation", "test"],
        }
    )


def fixture_artifact(tmp_path):
    observed = np.arange(27, dtype=float).reshape(3, 9)
    _npz, document, _sites = combined_builder.write_artifact(
        tmp_path / "parent",
        fixture_sites(),
        observed,
        np.ones_like(observed),
        np.array([True, False, True]),
        dispersion=0.0,
        metadata={"labels_used": False},
    )
    return document, observed


def test_subset_keeps_only_requested_rows_and_recomputes_hashes(tmp_path):
    parent, observed = fixture_artifact(tmp_path)
    npz, document, sites_path = subsetter.subset_artifact(
        parent,
        tmp_path / "train",
        [("chromosome_split", "train")],
    )
    sites = pd.read_csv(sites_path, sep="\t")
    assert sites["chromosome_split"].tolist() == ["train"]
    with np.load(npz, allow_pickle=False) as payload:
        np.testing.assert_array_equal(payload["observed"], observed[[0]])
        np.testing.assert_array_equal(payload["site_hash"], subsetter.site_hashes(sites))
    manifest = json.loads(document.read_text(encoding="utf-8"))
    assert manifest["sites_total"] == 1
    assert manifest["sites_valid"] == 1
    assert manifest["metadata"]["labels_used"] is False
    assert manifest["metadata"]["row_filters"] == [
        {"column": "chromosome_split", "equals": "train"}
    ]


def test_subset_rejects_empty_and_label_filters(tmp_path):
    parent, _observed = fixture_artifact(tmp_path)
    with pytest.raises(ValueError, match="selected no sites"):
        subsetter.subset_artifact(
            parent, tmp_path / "empty", [("chromosome_split", "missing")]
        )
    with pytest.raises(ValueError, match="label-derived"):
        subsetter.subset_artifact(parent, tmp_path / "bad", [("chip_label", "1")])
