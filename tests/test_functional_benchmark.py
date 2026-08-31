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
    fit_supervised_ceiling,
    residual_score,
    site_hashes,
    stable_seed,
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
                    "auroc": auprc,
                    "auprc": auprc,
                }
            )
    classified = classify_failures(pd.DataFrame(rows)).set_index("tf")
    assert classified.loc["NO_INFO", "classification"] == "assay_limited_or_motif_ambiguous"
    assert classified.loc["DETECTABLE", "classification"] == "detectable"
    assert classified.loc["MODEL_LIMIT", "classification"] == "shape_model_limited"
