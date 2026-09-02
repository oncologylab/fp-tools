from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from build_strand_functional_profiles import (  # noqa: E402
    orient_strand_log_bias,
    predict_strand_expected_profiles,
    site_hashes,
    write_profiles,
)
from evaluate_strand_functional_templates import (  # noqa: E402
    load_artifact,
    parse_artifact,
    stack_channels,
)
from evaluate_strand_label_free_models import (  # noqa: E402
    candidate_grid as label_free_candidate_grid,
    validate_unlabeled_training_sites,
)
from render_functional_aggregate_comparison import (  # noqa: E402
    combined_strand_shape,
    difference_band,
    mean_band,
    smooth_profiles,
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


def test_strand_log_bias_orientation_swaps_reverse_motifs() -> None:
    plus = np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    minus = np.asarray([[7.0, 8.0, 9.0], [10.0, 11.0, 12.0]])
    oriented_plus, oriented_minus, combined = orient_strand_log_bias(
        plus,
        minus,
        ["+", "-"],
    )
    assert oriented_plus[0].tolist() == [1.0, 2.0, 3.0]
    assert oriented_minus[0].tolist() == [7.0, 8.0, 9.0]
    assert oriented_plus[1].tolist() == [12.0, 11.0, 10.0]
    assert oriented_minus[1].tolist() == [6.0, 5.0, 4.0]
    assert np.allclose(combined, np.logaddexp(oriented_plus, oriented_minus) - np.log(2.0))


def test_strand_log_bias_orientation_silences_nonfinite_sentinels() -> None:
    plus = np.asarray([[0.0, np.nan, -np.inf]])
    minus = np.asarray([[1.0, 0.0, -np.inf]])

    with np.errstate(invalid="raise"):
        oriented_plus, oriented_minus, combined = orient_strand_log_bias(
            plus, minus, ["+"]
        )

    assert np.isnan(combined[0, 1])
    assert np.isneginf(combined[0, 2])
    assert np.array_equal(oriented_plus, plus, equal_nan=True)
    assert np.array_equal(oriented_minus, minus, equal_nan=True)


def test_strand_profile_artifact_is_safe_and_hashed(tmp_path: Path) -> None:
    sites = _sites()
    plus = np.ones((2, 21), dtype=float)
    minus = np.full((2, 21), 2.0)
    profiles = construct_strand_functional_profiles(
        plus, minus, plus, minus, sites["TFBS_strand"]
    )
    log_bias = (
        np.zeros_like(plus),
        np.ones_like(plus),
        np.full_like(plus, np.logaddexp(0.0, 1.0) - np.log(2.0)),
    )
    npz, metadata, site_table = write_profiles(
        tmp_path / "profiles",
        sites,
        profiles,
        np.asarray([True, True]),
        {"test": True},
        log_bias=log_bias,
    )
    assert npz.is_file() and metadata.is_file() and site_table.is_file()
    with np.load(npz, allow_pickle=False) as arrays:
        assert np.array_equal(arrays["site_hash"], site_hashes(sites))
        assert "antisymmetric_strand_residual" in arrays
        assert "combined_log_bias" in arrays

    dotted_npz, dotted_metadata, dotted_sites = write_profiles(
        tmp_path / "K562.selma10.shift_4_-4",
        sites,
        profiles,
        np.asarray([True, True]),
        {"test": True},
    )
    assert dotted_npz.name == "K562.selma10.shift_4_-4.npz"
    assert dotted_metadata.name == "K562.selma10.shift_4_-4.json"
    assert dotted_sites.name == "K562.selma10.shift_4_-4.sites.tsv.gz"


def test_loaded_strand_artifact_includes_raw_counts(tmp_path: Path) -> None:
    sites = _sites().assign(cell="K562")
    plus = np.ones((2, 21), dtype=float)
    minus = np.full((2, 21), 2.0)
    profiles = construct_strand_functional_profiles(
        plus, minus, plus, minus, sites["TFBS_strand"]
    )
    _npz, metadata, _site_table = write_profiles(
        tmp_path / "profiles",
        sites,
        profiles,
        np.asarray([True, True]),
        {"labels_used": False},
    )
    study = {
        "chromosome_split": {
            "train": ["chr1"],
            "validation": [],
            "internal_test": [],
        }
    }
    loaded_sites, loaded_profiles, _document = load_artifact(metadata, "K562", study)
    assert len(loaded_sites) == 2
    for name in ("plus_observed", "minus_observed", "plus_expected", "minus_expected"):
        assert loaded_profiles[name].shape == (2, 21)


def test_label_free_grid_and_label_firewall() -> None:
    candidates = label_free_candidate_grid()
    assert len(candidates) == 102
    assert {candidate.family for candidate in candidates} == {
        "count",
        "conditional",
        "fda",
        "hybrid",
        "anchored-fda",
        "residualized-fda",
    }
    protected = [
        candidate
        for candidate in candidates
        if candidate.candidate_id.startswith("conditional-protected_")
    ]
    assert len(protected) == 18
    assert {candidate.family for candidate in protected} == {"conditional"}
    anchored = [candidate for candidate in candidates if candidate.family == "anchored-fda"]
    assert len(anchored) == 18
    assert {candidate.anchor_strength for candidate in anchored} == {0.5, 1.0, 2.0}
    residualized = [
        candidate for candidate in candidates if candidate.family == "residualized-fda"
    ]
    assert len(residualized) == 18
    assert {candidate.covariate_ridge for candidate in residualized} == {
        1.0,
        10.0,
        100.0,
    }
    validate_unlabeled_training_sites(
        pd.DataFrame({"tf": ["A"], "chromosome_split": ["train"]}),
        "safe.tsv",
    )
    with pytest.raises(ValueError, match="chip_label"):
        validate_unlabeled_training_sites(
            pd.DataFrame(
                {"tf": ["A"], "chip_label": [1], "chromosome_split": ["train"]}
            ),
            "leaky.tsv",
        )
    with pytest.raises(ValueError, match="only train chromosomes"):
        validate_unlabeled_training_sites(
            pd.DataFrame({"tf": ["A"], "chromosome_split": ["validation"]}),
            "leaky-split.tsv",
        )


def test_strand_channel_sets_preserve_site_position_axes() -> None:
    values = {
        "combined_residual": np.zeros((5, 21)),
        "shared_strand_residual": np.ones((5, 21)),
        "antisymmetric_strand_residual": np.full((5, 21), 2.0),
    }
    combined = stack_channels(values, "combined")
    all_channels = stack_channels(values, "all")
    assert combined.shape == (5, 1, 21)
    assert all_channels.shape == (5, 3, 21)
    assert np.all(all_channels[:, 2, :] == 2.0)
    assert parse_artifact("SELMA,K562,profiles.json") == (
        "SELMA",
        "K562",
        Path("profiles.json"),
    )


def test_blinded_aggregate_combines_channels_and_bootstraps() -> None:
    rng = np.random.default_rng(3)
    values = {
        "combined_residual": rng.normal(size=(20, 21)),
        "shared_strand_residual": rng.normal(size=(20, 21)),
        "antisymmetric_strand_residual": rng.normal(size=(20, 21)),
    }
    winner = pd.Series(
        {
            "channel_set": "shared_antisymmetric",
            "channel_weights": "2,-1",
        }
    )
    combined = combined_strand_shape(values, winner, np.arange(-10, 11))
    assert combined.shape == (20, 21)
    mean, lower, upper = mean_band(combined, bootstraps=25, seed=8)
    assert mean.shape == lower.shape == upper.shape == (21,)
    assert np.all(lower <= upper)
    smoothed = smooth_profiles(combined, 2.0)
    assert smoothed.shape == combined.shape
    labels = np.asarray([0, 1] * 10)
    difference, lower, upper = difference_band(
        smoothed,
        labels,
        1,
        0,
        bootstraps=25,
        seed=9,
    )
    assert difference.shape == lower.shape == upper.shape == (21,)
    assert np.all(lower <= upper)
