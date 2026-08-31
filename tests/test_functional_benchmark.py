from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_footprints import (  # noqa: E402
    build_unlabeled_training_sites,
    chromosome_split,
    classify_failures,
    cluster_functional_phenotypes,
    derive_parametric_expected_profiles,
    fit_supervised_ceiling,
    residual_score,
    site_hashes,
    stable_seed,
    summarize_aggregate_profiles,
    validate_sites,
)
from search_functional_model_grid import candidate_grid  # noqa: E402
from fp_tools.tools.parametric_bias import (  # noqa: E402
    BiasFeatureSpec,
    ConditionalSequenceBiasModel,
)


SPEC = Path(__file__).resolve().parents[1] / "benchmarks" / "manifests" / "footprint_functional_v1.spec.json"


def test_functional_spec_is_locked_and_complete() -> None:
    study = json.loads(SPEC.read_text(encoding="utf-8"))
    assert study["status"] == "development_locked_holdout_unscored"
    assert study["cut_shift_arms"] == {"legacy": [4, -5], "aligned": [4, -4]}
    assert study["locked_holdout_cells"] == ["GM12878", "IMR-90"]
    assert {task["cell"] for task in study["tasks"] if task["split"] == "locked_holdout"} == {
        "GM12878",
        "IMR-90",
    }
    assert study["promotion_gates"]["minimum_gp_relative_auprc_gain_over_spline"] == 0.05
    assert chromosome_split("chr16", study) == "validation"
    assert chromosome_split("chr19", study) == "test"


def test_functional_hyperparameter_grid_covers_prespecified_ablations() -> None:
    compact = candidate_grid("compact")
    full = candidate_grid("full")
    assert len(compact) == 30
    assert len(full) == 77
    assert len({candidate.candidate_id for candidate in full}) == len(full)
    assert {candidate.background for candidate in compact} == {
        "none",
        "linear",
        "quadratic",
        "gp-long",
    }
    assert {candidate.prior_constraint for candidate in compact} == {
        "none",
        "motif-accessibility",
    }
    assert {candidate.likelihood_limit for candidate in compact if candidate.likelihood_limit} == {
        30.0,
        50.0,
        80.0,
    }
    gp_scales = {
        (candidate.long_length_scale, candidate.short_length_scale)
        for candidate in full
        if candidate.family == "gp"
    }
    assert {(long, short) for long in (30.0, 50.0, 80.0) for short in (3.0, 6.0, 10.0, 15.0)}.issubset(
        gp_scales
    )


def test_unlabeled_training_sites_never_read_labels(tmp_path: Path) -> None:
    study = json.loads(SPEC.read_text(encoding="utf-8"))
    tasks = pd.DataFrame(
        [
            {
                "cell": "K562",
                "tf": "CTCF",
                "motif_id": "MA0139.2",
                "motif_family": "CTCF",
            }
        ]
    )
    source = pd.DataFrame(
        {
            "motif": ["CTCF_MA0139.2"] * 12,
            "TFBS_chr": ["chr1"] * 8 + ["chr19"] * 4,
            "TFBS_start": np.arange(12) * 100 + 1000,
            "TFBS_end": np.arange(12) * 100 + 1015,
            "TFBS_strand": ["+", "-"] * 6,
            "TFBS_score": np.linspace(1, 2, 12),
        }
    )
    source_path = tmp_path / "motifs.tsv"
    source.to_csv(source_path, sep="\t", index=False)
    sampled = build_unlabeled_training_sites(
        source_path,
        "K562",
        tasks,
        study,
        maximum_per_tf=5,
        seed=2026,
    )
    assert len(sampled) == 5
    assert set(sampled["chromosome_split"]) == {"train"}
    assert "chip_label" not in sampled

    source["chip_label"] = 1
    source.to_csv(source_path, sep="\t", index=False)
    with pytest.raises(ValueError, match="must not contain"):
        build_unlabeled_training_sites(
            source_path,
            "K562",
            tasks,
            study,
            maximum_per_tf=5,
            seed=2026,
        )


def test_site_hash_and_seed_are_deterministic() -> None:
    frame = pd.DataFrame(
        {
            "TFBS_chr": ["chr1", "chr2"],
            "TFBS_start": [10, 20],
            "TFBS_end": [20, 30],
            "TFBS_strand": ["+", "-"],
            "tf": ["A", "B"],
        }
    )
    assert np.array_equal(site_hashes(frame), site_hashes(frame.copy()))
    assert stable_seed("K562", "CTCF") == stable_seed("K562", "CTCF")
    assert stable_seed("K562", "CTCF") != stable_seed("HepG2", "CTCF")


def test_legacy_site_table_gets_coverage_placeholder() -> None:
    frame = pd.DataFrame(
        {
            "cell": ["K562", "K562"],
            "tf": ["CTCF", "CTCF"],
            "motif_family": ["CTCF", "CTCF"],
            "TFBS_chr": ["chr1", "chr1"],
            "TFBS_start": [10, 30],
            "TFBS_end": [25, 45],
            "TFBS_strand": ["+", "-"],
            "motif_score": [1.0, 2.0],
            "chip_label": [0, 1],
            "chromosome_split": ["train", "train"],
        }
    )
    validated = validate_sites(frame, "memory")
    assert validated["accessibility"].tolist() == [0.0, 0.0]


def test_residual_factorial_scores_protection() -> None:
    observed = np.full((20, 101), 4.0)
    expected = np.full_like(observed, 4.0)
    observed[:10, 45:56] = 1.0
    for mode in ("difference", "pearson", "deviance", "log-ratio", "nb-likelihood"):
        residual, score = residual_score(observed, expected, mode, dispersion=0.05)
        assert residual.shape == observed.shape
        assert np.mean(score[:10]) > np.mean(score[10:])


def test_supervised_ceiling_uses_functional_shape() -> None:
    rng = np.random.default_rng(8)
    width = 101
    x = np.arange(width) - width // 2
    shape = -np.exp(-0.5 * np.square(x / 6.0))

    def make_frame(n: int) -> tuple[np.ndarray, pd.DataFrame]:
        labels = rng.integers(0, 2, size=n)
        profiles = labels[:, None] * shape + rng.normal(scale=0.4, size=(n, width))
        sites = pd.DataFrame(
            {
                "motif_score": rng.normal(size=n),
                "accessibility": rng.uniform(1, 20, size=n),
                "chip_label": labels,
            }
        )
        return profiles, sites

    train_profiles, train_sites = make_frame(300)
    validation_profiles, validation_sites = make_frame(180)
    probabilities, selected, fpca = fit_supervised_ceiling(
        train_profiles,
        validation_profiles,
        train_sites,
        validation_sites,
        seed=4,
    )
    assert selected["auroc"] > 0.8
    assert selected["converged"] is True
    assert selected["converged_candidates"] > 0
    assert selected["candidate_count"] == 15
    assert selected["iterations"] < 20000
    assert probabilities.shape == (180,)
    assert fpca.components_.shape[0] <= 20


def test_failure_classification_separates_assay_and_shape_limits() -> None:
    rows = []
    for tf, baseline, ceiling, label_free in (
        ("NO_INFO", 0.50, 0.52, 0.51),
        ("DETECTABLE", 0.50, 0.80, 0.70),
        ("MODEL_LIMIT", 0.50, 0.80, 0.52),
    ):
        for method, auprc in (
            ("supervised_baseline", baseline),
            ("supervised_fpca", ceiling),
            ("gp", label_free),
        ):
            rows.append(
                {
                    "split": "validation",
                    "cell": "K562",
                    "tf": tf,
                    "motif_family": tf,
                    "correction": "DWM",
                    "method": method,
                    "positive_sites": 500,
                    "prevalence": 0.5,
                    "auroc": auprc,
                    "auprc": auprc,
                }
            )
    classified = classify_failures(pd.DataFrame(rows)).set_index("tf")
    assert classified.loc["NO_INFO", "classification"] == "assay_limited_or_motif_ambiguous"
    assert classified.loc["DETECTABLE", "classification"] == "detectable"
    assert classified.loc["MODEL_LIMIT", "classification"] == "shape_model_limited"


def test_parametric_expected_profiles_preserve_site_totals(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "ACGT" * 300 + "\n", encoding="utf-8")
    pysam.faidx(str(genome))
    model = ConditionalSequenceBiasModel(BiasFeatureSpec.selma10())
    model.main[4, 0] = 1.5
    model.main[4] -= model.main[4].mean()
    model_path, _ = model.save(tmp_path / "bias")
    sites = pd.DataFrame(
        {
            "TFBS_chr": ["chr1", "chr1"],
            "TFBS_start": [480, 680],
            "TFBS_end": [490, 690],
            "TFBS_strand": ["+", "-"],
            "tf": ["A", "A"],
        }
    )
    observed = np.zeros((2, 101), dtype=np.float32)
    observed[0, 20:80] = 2.0
    observed[1, 10:90] = 1.0
    expected = derive_parametric_expected_profiles(
        sites,
        observed,
        model_path,
        genome,
        tmp_path / "expected.npz",
        flank=50,
    )
    assert np.allclose(expected.sum(axis=1), observed.sum(axis=1))
    assert np.std(expected[0]) > 0
    cached = derive_parametric_expected_profiles(
        sites,
        observed,
        model_path,
        genome,
        tmp_path / "expected.npz",
        flank=50,
    )
    assert np.array_equal(cached, expected)


def test_aggregate_profile_summary_and_unsupervised_clustering() -> None:
    rng = np.random.default_rng(22)
    positions = np.arange(-50, 51)
    expected = np.full((120, len(positions)), 4.0)
    observed = rng.poisson(expected).astype(float)
    labels = np.repeat([0, 1], 60)
    observed[labels == 1, 45:56] *= 0.2
    aggregate_frames = []
    descriptor_frames = []
    for index, tf in enumerate(("A", "B", "C", "D")):
        aggregate, descriptors = summarize_aggregate_profiles(
            observed,
            expected,
            labels,
            split="validation",
            cell="K562",
            tf=tf,
            motif_family=f"family_{index}",
            correction="DWM",
            dispersion=0.05,
            positions=positions,
            bootstraps=20,
            seed=5,
        )
        aggregate_frames.append(aggregate)
        descriptor_frames.append(descriptors)
    aggregates = pd.concat(aggregate_frames, ignore_index=True)
    descriptors = pd.concat(descriptor_frames, ignore_index=True)
    positive = descriptors[descriptors["group"] == "chip_positive"]
    negative = descriptors[descriptors["group"] == "matched_negative"]
    assert positive["depletion"].mean() > negative["depletion"].mean()
    assert {"mean", "lower_95", "upper_95"}.issubset(aggregates.columns)
    clusters = cluster_functional_phenotypes(aggregates, descriptors, maximum_clusters=3)
    assert len(clusters) == 4
    assert clusters["functional_cluster"].notna().all()
    assert clusters["phenotype"].notna().all()
