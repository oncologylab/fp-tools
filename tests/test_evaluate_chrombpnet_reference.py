from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import evaluate_chrombpnet_reference as reference  # noqa: E402
from fp_tools.utils import bigwig as pyBigWig  # noqa: E402


def test_scale_to_observed_total_is_shape_and_scale_invariant() -> None:
    observed = np.asarray([[2.0, 3.0, 5.0], [4.0, 6.0, 10.0]])
    prediction = np.asarray([[1.0, 2.0, 1.0], [2.0, 1.0, 1.0]])
    scaled = reference.scale_to_observed_total(observed, prediction)
    np.testing.assert_allclose(scaled.sum(axis=1), observed.sum(axis=1))
    np.testing.assert_allclose(
        scaled,
        reference.scale_to_observed_total(observed, prediction * 37.0),
    )


def _write_bigwig(path: Path) -> None:
    handle = pyBigWig.open(str(path), "w")
    handle.addHeader([("chr1", 1000)])
    starts = list(range(98, 103))
    handle.addEntries(
        ["chr1"] * len(starts),
        starts,
        ends=[value + 1 for value in starts],
        values=[1.0, 2.0, 3.0, 4.0, 5.0],
    )
    handle.close()


def test_dense_prediction_extraction_rejects_unpredicted_sites_and_orients(
    tmp_path: Path,
) -> None:
    signal = tmp_path / "prediction.bw"
    _write_bigwig(signal)
    sites = pd.DataFrame(
        {
            "TFBS_chr": ["chr1", "chr1", "chr1"],
            "TFBS_start": [100, 100, 200],
            "TFBS_end": [101, 101, 201],
            "TFBS_strand": ["+", "-", "+"],
        }
    )
    profiles, valid = reference.extract_prediction_profiles(
        sites, signal, 2, require_dense=True
    )
    assert valid.tolist() == [True, True, False]
    np.testing.assert_allclose(profiles[0], [1, 2, 3, 4, 5])
    np.testing.assert_allclose(profiles[1], [5, 4, 3, 2, 1])
    sparse, sparse_valid = reference.extract_prediction_profiles(
        sites, signal, 2, require_dense=False
    )
    assert sparse_valid.tolist() == [True, True, True]
    np.testing.assert_allclose(sparse[2], 0.0)
    assert reference.validate_bigwig(signal)["covered_bases"] == 5


def test_bias_prediction_equal_to_observed_has_zero_residual_geometry() -> None:
    positions = np.arange(-100, 101, dtype=float)
    bias_shape = 1.0 + 2.0 * np.exp(-0.5 * np.square(positions / 7.0))
    observed = np.vstack([bias_shape * 10.0, bias_shape * 20.0])
    uniform = np.ones_like(observed)
    methods = reference.score_methods(
        observed,
        uniform,
        observed * 2.0,
        observed * 4.0,
        observed,
        observed,
        positions,
        dispersion=0.0,
    )
    score, profile = methods["ChromBPNet_bias_conventional_geometry"]
    # Total matching and the deviance square root can leave machine-scale
    # roundoff even when the predicted and observed shapes are identical.
    np.testing.assert_allclose(score, 0.0, atol=5e-8)
    np.testing.assert_allclose(profile, 0.0, atol=2e-7)
    parametric_score, parametric_profile = methods[
        "frozen_parametric_bias_conventional_geometry"
    ]
    np.testing.assert_allclose(parametric_score, 0.0, atol=5e-8)
    np.testing.assert_allclose(parametric_profile, 0.0, atol=2e-7)
    assert np.any(np.abs(methods["DWM_conventional_geometry"][0]) > 0)


def test_reference_manifest_requires_pin_completion_and_checksums(
    tmp_path: Path,
) -> None:
    output = tmp_path / "prediction.bw"
    output.write_bytes(b"fixture")
    manifest = tmp_path / "manifest.json"
    document = {
        "schema": reference.REFERENCE_SCHEMA,
        "source_commit": reference.PINNED_COMMIT,
        "completed": True,
        "stage": "predict",
        "outputs": [
            {
                "path": str(output),
                "sha256": reference.file_sha256(output),
            }
        ],
    }
    manifest.write_text(json.dumps(document), encoding="utf-8")
    assert reference.validate_reference_manifest(manifest)["stage"] == "predict"
    output.write_bytes(b"changed")
    with pytest.raises(ValueError, match="checksum mismatch"):
        reference.validate_reference_manifest(manifest)
