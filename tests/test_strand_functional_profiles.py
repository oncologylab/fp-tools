from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from build_strand_functional_profiles import (  # noqa: E402
    predict_strand_expected_profiles,
    site_hashes,
    write_profiles,
)
from fp_tools.tools.functional_footprints import (  # noqa: E402
    construct_strand_functional_profiles,
)
from fp_tools.tools.parametric_bias import (  # noqa: E402
    BiasFeatureSpec,
    ConditionalSequenceBiasModel,
)


def _sites() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "TFBS_chr": ["chr1", "chr1"],
            "TFBS_start": [480, 680],
            "TFBS_end": [490, 690],
            "TFBS_strand": ["+", "-"],
            "tf": ["A", "A"],
        }
    )


def test_strand_expected_profiles_preserve_each_strand_total(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "ACGT" * 300 + "\n", encoding="utf-8")
    pysam.faidx(str(genome))
    model = ConditionalSequenceBiasModel(BiasFeatureSpec.selma10())
    model.main[4, 0] = 1.5
    model.main[4] -= model.main[4].mean()
    plus = np.ones((2, 101), dtype=float)
    minus = np.full((2, 101), 2.0)
    plus_expected, minus_expected, valid = predict_strand_expected_profiles(
        _sites(), plus, minus, model, genome, flank=50
    )
    assert valid.all()
    assert np.allclose(plus_expected.sum(axis=1), plus.sum(axis=1))
    assert np.allclose(minus_expected.sum(axis=1), minus.sum(axis=1))
    assert np.std(plus_expected[0]) > 0


def test_strand_profile_artifact_is_safe_and_hashed(tmp_path: Path) -> None:
    sites = _sites()
    plus = np.ones((2, 21), dtype=float)
    minus = np.full((2, 21), 2.0)
    profiles = construct_strand_functional_profiles(
        plus, minus, plus, minus, sites["TFBS_strand"]
    )
    npz, metadata, site_table = write_profiles(
        tmp_path / "profiles", sites, profiles, np.asarray([True, True]), {"test": True}
    )
    assert npz.is_file() and metadata.is_file() and site_table.is_file()
    with np.load(npz, allow_pickle=False) as arrays:
        assert np.array_equal(arrays["site_hash"], site_hashes(sites))
        assert "antisymmetric_strand_residual" in arrays
