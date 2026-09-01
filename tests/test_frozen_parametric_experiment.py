from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from types import SimpleNamespace
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import freeze_parametric_holdouts  # noqa: E402
import freeze_label_free_functional_models  # noqa: E402
import freeze_functional_call_thresholds  # noqa: E402
import evaluate_frozen_functional_naked_dna  # noqa: E402
import evaluate_frozen_functional_depth_matrix  # noqa: E402
import evaluate_frozen_functional_information_ceiling  # noqa: E402
import evaluate_frozen_bias_shrinkage  # noqa: E402
import evaluate_frozen_functional_policy  # noqa: E402
import evaluate_parametric_factorization  # noqa: E402
import evaluate_strand_label_free_models  # noqa: E402
import derive_label_free_family_selection  # noqa: E402
import sample_label_free_motif_sites  # noqa: E402
import run_frozen_parametric_experiment  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402


def _record(
    accession: str,
    *,
    step: str,
    preferred: bool = False,
    major: int = 1,
    minor: int = 0,
) -> dict:
    return {
        "accession": accession,
        "status": "released",
        "assembly": "GRCh38",
        "file_format": "bed",
        "output_type": "conservative IDR thresholded peaks",
        "biological_replicates": [1, 2],
        "preferred_default": preferred,
        "analysis_step_version": {
            "minor_version": minor,
            "name": f"{step}-v-{major}-{minor}",
            "analysis_step": {
                "step_label": step,
                "major_version": major,
                "pipelines": [],
            },
        },
    }


def test_replicate_aware_selector_excludes_pseudoreplicates() -> None:
    replicated = _record("REPLICATED", step="tf-chip-seq-replicated-idr-step")
    pseudo = _record("PSEUDO", step="tf-chip-seq-pooled-pseudoreplicated-idr-step")
    assert freeze_parametric_holdouts.is_replicate_aware_idr(replicated)
    assert not freeze_parametric_holdouts.is_replicate_aware_idr(pseudo)
    replicated["biological_replicates"] = [1]
    assert not freeze_parametric_holdouts.is_replicate_aware_idr(replicated)


def test_locked_candidate_tie_break_is_preferred_version_accession() -> None:
    preferred = _record("A", step="tf-idr-step", preferred=True, major=1)
    newer = _record("B", step="tf-idr-step", preferred=False, major=2)
    assert (
        max([preferred, newer], key=freeze_parametric_holdouts.candidate_rank)
        is preferred
    )
    newer_preferred = _record("B", step="tf-idr-step", preferred=True, major=2)
    assert (
        max([preferred, newer_preferred], key=freeze_parametric_holdouts.candidate_rank)
        is newer_preferred
    )
    same = _record("Z", step="tf-idr-step", preferred=True, major=2)
    assert (
        max([newer_preferred, same], key=freeze_parametric_holdouts.candidate_rank)
        is same
    )


def test_committed_holdout_freeze_locks_metadata_before_labels() -> None:
    freeze_path = (
        ROOT
        / "benchmarks"
        / "manifests"
        / "frozen_parametric_factorization_v1.freeze.json"
    )
    document = json.loads(freeze_path.read_text(encoding="utf-8"))
    assert document["chipped_peak_contents_read"] is False
    assert document["holdout_labels_scored"] is False
    assert document["manifest"]["rows"] == 17
    for key in ("study", "selector", "manifest"):
        path = Path(document[key]["path"])
        if not path.is_absolute():
            path = ROOT / path
        digest = sha256(path.read_bytes()).hexdigest()
        assert digest == document[key]["sha256"]
    decisions = {
        (row["cell"], row["tf"]): row["status"] for row in document["chip_decisions"]
    }
    assert decisions[("GM23338", "REST")].startswith("ineligible")
    assert decisions[("GM23338", "NANOG")] == "ineligible_no_jaspar_motif"
    assert decisions[("SK-N-SH", "MEF2A")] == "selected_pending_power_check"


def test_factorization_geometry_score_detects_center_protection() -> None:
    positions = np.arange(-100, 101, dtype=float)
    profiles = np.zeros((2, len(positions)), dtype=float)
    profiles[0, np.abs(positions) <= 7] = -2.0
    profiles[1, np.abs(positions) <= 7] = 2.0
    scores = evaluate_parametric_factorization.geometry_score(profiles, positions)
    assert scores[0] > 0
    assert scores[1] < 0


def test_factorization_metrics_mark_small_tasks_underpowered() -> None:
    positions = np.arange(-100, 101, dtype=float)
    sites = pd.DataFrame(
        {
            "cell": ["K562"] * 4,
            "tf": ["TF"] * 4,
            "motif": ["M1"] * 4,
            "motif_family": ["family"] * 4,
            "role": ["weak_shape"] * 4,
            "chip_label": [0, 0, 1, 1],
        }
    )
    row = evaluate_parametric_factorization.metric_row(
        sites,
        np.asarray([0.1, 0.2, 0.8, 0.9]),
        np.zeros((4, len(positions))),
        cell="K562",
        tf="TF",
        method="candidate",
        split="test",
        positions=positions,
        minimum_sites_per_class=200,
    )
    assert row is not None
    assert row["status"] == "underpowered"
    assert row["n_positive"] == row["n_negative"] == 2


def test_frozen_functional_aggregate_bootstrap_is_deterministic() -> None:
    profiles = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 3.0, 0.0],
        ]
    )
    labels = np.asarray([0, 0, 1, 1])
    chromosomes = np.asarray(["chr19", "chr20", "chr19", "chr20"])
    first = evaluate_frozen_functional_policy.aggregate_curve(
        profiles,
        labels,
        chromosomes,
        iterations=100,
        seed=17,
    )
    second = evaluate_frozen_functional_policy.aggregate_curve(
        profiles,
        labels,
        chromosomes,
        iterations=100,
        seed=17,
    )
    assert np.allclose(first["difference"], [0.0, 2.0, 0.0])
    assert np.array_equal(first["lower_95"], second["lower_95"])
    assert np.array_equal(first["upper_95"], second["upper_95"])


def test_frozen_functional_model_names_cannot_create_suffix_collisions() -> None:
    first = freeze_label_free_functional_models.safe_token("fda.shared.pool_tf")
    second = freeze_label_free_functional_models.safe_token("fda.antisymmetric.pool_tf")
    assert first == "fda_shared_pool_tf"
    assert second == "fda_antisymmetric_pool_tf"
    assert first != second


def test_cached_bias_depth_profiles_preserve_each_site_total() -> None:
    observed = np.asarray([[2.0, 3.0, 5.0], [0.0, 7.0, 1.0]])
    log_bias = np.log(np.asarray([[1.0, 2.0, 1.0], [3.0, 1.0, 2.0]]))
    expected, valid = (
        evaluate_frozen_functional_depth_matrix.expected_from_cached_log_bias(
            observed,
            log_bias,
        )
    )
    assert valid.tolist() == [True, True]
    assert np.allclose(expected.sum(axis=1), observed.sum(axis=1))
    assert not np.isclose(expected[0].sum(), observed.sum())


def test_depth_dwm_scaling_matches_local_totals() -> None:
    observed = np.asarray([[2.0, 3.0, 5.0], [0.0, 7.0, 1.0]])
    expected = np.asarray([[10.0, 20.0, 10.0], [3.0, 1.0, 2.0]])
    scaled = evaluate_frozen_functional_depth_matrix.scale_expected_to_observed(
        observed,
        expected,
    )
    assert np.allclose(scaled.sum(axis=1), observed.sum(axis=1))
    assert np.allclose(scaled[0] / scaled[0].sum(), expected[0] / expected[0].sum())


def test_depth_matrix_applies_frozen_global_and_tf_shrinkage() -> None:
    positions = np.arange(-100, 101, dtype=float)
    observed = np.ones((1, len(positions)), dtype=float)
    observed[0, 95:106] = 5.0
    direct = np.ones_like(observed) * observed.sum() / observed.shape[1]
    methods = evaluate_frozen_functional_depth_matrix.frozen_shrinkage_methods(
        observed,
        direct,
        np.zeros_like(observed),
        positions,
        bias_strength=0.75,
        global_choice={"source": "parametric_direct", "alpha": 0.25},
        tf_choice={"source": "raw", "alpha": 0.0},
    )
    values = {method: profiles for method, _score, profiles, _seconds, _model in methods}
    assert np.allclose(
        values["frozen_global_shrinkage"], observed - 0.25 * direct
    )
    assert np.array_equal(values["frozen_tf_specific_shrinkage"], observed)


def test_depth_metrics_keep_both_dwm_and_raw_guardrails() -> None:
    rows = [
        {
            "method": "DWM_conventional_geometry",
            "auroc": 0.60,
            "auprc": 0.25,
            "functional_separation": 2.0,
        },
        {
            "method": "raw_geometry",
            "auroc": 0.65,
            "auprc": 0.30,
            "functional_separation": 2.5,
        },
        {
            "method": "frozen_candidate",
            "auroc": 0.70,
            "auprc": 0.36,
            "functional_separation": 3.0,
        },
    ]
    result = evaluate_frozen_functional_depth_matrix.add_depth_baseline_deltas(rows)
    candidate = result[-1]
    assert candidate["auroc_gain_over_dwm"] == pytest.approx(0.10)
    assert candidate["relative_auprc_gain_over_dwm"] == pytest.approx(0.44)
    assert candidate["auroc_gain_over_raw"] == pytest.approx(0.05)
    assert candidate["relative_auprc_gain_over_raw"] == pytest.approx(0.20)
    assert candidate["functional_separation_relative_change_over_dwm"] == pytest.approx(0.5)
    assert candidate["functional_separation_relative_change_over_raw"] == pytest.approx(0.2)


def test_depth_discovery_and_artifacts_preserve_biological_samples(
    tmp_path: Path,
) -> None:
    samples = [("CellA", "CellA_rep1"), ("CellA", "CellA_rep2")]
    for _cell, sample in samples:
        directory = tmp_path / sample / "full" / "seed_2026"
        (directory / "fp_tools_dwm").mkdir(parents=True)
        (directory / f"{sample}.full.s2026.bam").write_bytes(b"bam")
        (directory / f"{sample}.full.s2026.bam.bai").write_bytes(b"bai")
        (directory / "fp_tools_dwm" / f"{sample}.full.s2026_expected.bw").write_bytes(
            b"bigwig"
        )
    discovered = evaluate_frozen_functional_depth_matrix.discover_signals(
        tmp_path,
        samples,
        ["full"],
        [2026],
        allow_incomplete=False,
    )
    assert discovered[["cell", "sample"]].to_records(index=False).tolist() == samples
    prefixes = {
        evaluate_frozen_functional_depth_matrix.artifact_prefix(
            tmp_path,
            SimpleNamespace(cell=cell, sample=sample, depth="full", seed=2026),
        )
        for cell, sample in samples
    }
    assert len(prefixes) == 2


def test_depth_replicate_summary_reports_direction_stability() -> None:
    rows = []
    for sample, gain in (("CellA_rep1", 0.1), ("CellA_rep2", -0.05)):
        rows.append(
            {
                "cell": "CellA",
                "sample": sample,
                "tf": "TF1",
                "motif_family": "FAMILY1",
                "candidate_id": "candidate",
                "method": "frozen_candidate",
                "depth": "full",
                "seed": 2026,
                "auroc": 0.6 + gain,
                "auprc": 0.4 + gain,
                "auroc_gain_over_dwm": gain,
                "relative_auprc_gain_over_dwm": gain,
                "auroc_gain_over_raw": gain - 0.02,
                "relative_auprc_gain_over_raw": gain - 0.02,
                "functional_separation": 1.0 + gain,
                "functional_separation_relative_change_over_dwm": gain,
                "functional_separation_relative_change_over_raw": gain - 0.02,
                "brier": 0.2,
                "calibration_error": 0.1,
                "prediction_seconds": 0.01,
            }
        )
    metrics = pd.DataFrame(rows)
    per_sample = evaluate_frozen_functional_depth_matrix.summarize_metrics(metrics)
    across = evaluate_frozen_functional_depth_matrix.summarize_replicates(metrics)
    assert len(per_sample) == 2
    assert across.loc[0, "samples"] == 2
    assert across.loc[0, "observations"] == 2
    assert across.loc[0, "auroc_gain_positive_fraction"] == 0.5
    assert across.loc[0, "auroc_gain_over_raw_positive_fraction"] == 0.5


def test_depth_classification_prefers_full_endpoint() -> None:
    summary = pd.DataFrame(
        {
            "cell": ["CellA"] * 3,
            "tf": ["TF1"] * 3,
            "motif_family": ["FAMILY1"] * 3,
            "candidate_id": ["candidate"] * 3,
            "method": ["frozen_candidate"] * 3,
            "depth": ["10m", "50m", "full"],
            "auroc_mean": [0.55, 0.60, 0.70],
            "auroc_gain_over_dwm_mean": [0.01, 0.02, 0.04],
            "relative_auprc_gain_over_dwm_mean": [0.01, 0.02, 0.12],
            "auroc_gain_over_raw_mean": [0.01, 0.02, 0.04],
            "relative_auprc_gain_over_raw_mean": [0.01, 0.02, 0.12],
            "auroc_gain_positive_fraction": [0.6, 0.8, 1.0],
            "auroc_gain_over_raw_positive_fraction": [0.6, 0.8, 1.0],
            "auprc_gain_over_raw_positive_fraction": [0.6, 0.8, 1.0],
        }
    )
    result = evaluate_frozen_functional_depth_matrix.classify_depth(summary)
    assert result.loc[0, "high_depth"] == "full"
    assert result.loc[0, "high_auroc"] == 0.70
    assert result.loc[0, "classification"] == "detectable_above_raw_at_high_depth"


def test_depth_classification_keeps_frozen_methods_separate() -> None:
    rows = []
    for method, high_auroc in (
        ("frozen_global_shrinkage", 0.62),
        ("frozen_tf_specific_shrinkage", 0.71),
    ):
        for depth, auroc in (("10m", 0.55), ("50m", high_auroc)):
            rows.append(
                {
                    "cell": "CellA",
                    "tf": "TF1",
                    "motif_family": "FAMILY1",
                    "candidate_id": "shrinkage",
                    "method": method,
                    "depth": depth,
                    "auroc_mean": auroc,
                    "auroc_gain_over_dwm_mean": 0.04,
                    "relative_auprc_gain_over_dwm_mean": 0.12,
                    "auroc_gain_over_raw_mean": 0.03,
                    "relative_auprc_gain_over_raw_mean": 0.08,
                    "auroc_gain_positive_fraction": 1.0,
                    "auroc_gain_over_raw_positive_fraction": 1.0,
                    "auprc_gain_over_raw_positive_fraction": 1.0,
                }
            )

    result = evaluate_frozen_functional_depth_matrix.classify_depth(
        pd.DataFrame(rows)
    ).set_index("method")

    assert set(result.index) == {
        "frozen_global_shrinkage",
        "frozen_tf_specific_shrinkage",
    }
    assert result.loc["frozen_tf_specific_shrinkage", "high_auroc"] == 0.71


def test_depth_classification_identifies_dwm_overcorrection() -> None:
    summary = pd.DataFrame(
        {
            "cell": ["CellA", "CellA"],
            "tf": ["TF1", "TF1"],
            "motif_family": ["FAMILY1", "FAMILY1"],
            "candidate_id": ["candidate", "candidate"],
            "method": ["frozen_candidate", "frozen_candidate"],
            "depth": ["10m", "50m"],
            "auroc_mean": [0.52, 0.60],
            "auroc_gain_over_dwm_mean": [0.05, 0.10],
            "relative_auprc_gain_over_dwm_mean": [0.10, 0.20],
            "auroc_gain_over_raw_mean": [-0.01, -0.01],
            "relative_auprc_gain_over_raw_mean": [0.01, 0.01],
            "auroc_gain_positive_fraction": [1.0, 1.0],
            "auroc_gain_over_raw_positive_fraction": [0.0, 0.0],
            "auprc_gain_over_raw_positive_fraction": [1.0, 1.0],
        }
    )

    result = evaluate_frozen_functional_depth_matrix.classify_depth(summary)

    assert result.loc[0, "classification"] == (
        "dwm_overcorrection_recovered_toward_raw"
    )
    assert result.loc[0, "high_raw_guard_positive_fraction"] == 0.0


def test_frozen_site_score_frame_preserves_artifact_indexes() -> None:
    sites = pd.DataFrame(
        {
            "TFBS_chr": ["chr19", "chr20", "chr21"],
            "TFBS_start": [10, 20, 30],
            "TFBS_end": [15, 25, 35],
            "TFBS_strand": ["+", "-", "+"],
            "motif_score": [7.0, 8.0, 9.0],
            "accessibility": [0.0, 9.0, 3.0],
            "chip_label": [0, 1, 0],
            "motif": ["M1", "M1", "M1"],
        }
    )
    candidate = evaluate_strand_label_free_models.Candidate(
        candidate_id="count_spline.bg_none.window_30",
        family="count",
        smoother="spline",
        background="none",
        window=30.0,
        channel="combined_residual",
        training_pool="tf",
    )
    result = evaluate_frozen_functional_policy.site_score_frame(
        record={
            "cell": "CellA",
            "tf": "TF1",
            "motif_family": "F1",
            "bias_configuration": "LOG21",
        },
        candidate=candidate,
        sites=sites,
        indexes=np.asarray([2, 0]),
        candidate_score=np.asarray([0.8, 0.2]),
        dwm_score=np.asarray([1.5, -0.5]),
        raw_score=np.asarray([2.0, 0.1]),
        direct_score=np.asarray([1.2, -0.2]),
    )
    assert result["artifact_index"].tolist() == [2, 0]
    assert result["TFBS_start"].tolist() == [30, 10]
    assert result["label"].tolist() == [0, 0]
    assert result["candidate_probability"].tolist() == [0.8, 0.2]
    assert result["raw_score"].tolist() == [2.0, 0.1]
    assert result["log_accessibility"].tolist() == pytest.approx(
        [np.log1p(3.0), 0.0]
    )


@pytest.mark.parametrize(
    ("supervised", "auroc", "auprc", "expected"),
    [
        (0.05, 0.10, 0.20, "assay_limited"),
        (0.20, 0.01, 0.05, "shape_model_limited"),
        (0.20, 0.03, 0.10, "detectable"),
        (np.nan, 0.10, 0.20, "insufficient_supervised_folds"),
    ],
)
def test_frozen_information_ceiling_failure_classification(
    supervised: float,
    auroc: float,
    auprc: float,
    expected: str,
) -> None:
    assert (
        evaluate_frozen_functional_information_ceiling.failure_classification(
            supervised_relative_auprc_gain=supervised,
            label_free_auroc_gain=auroc,
            label_free_relative_auprc_gain=auprc,
        )
        == expected
    )


def test_frozen_information_ceiling_marks_unstable_classifier() -> None:
    assert (
        evaluate_frozen_functional_information_ceiling.failure_classification(
            supervised_relative_auprc_gain=0.25,
            label_free_auroc_gain=0.05,
            label_free_relative_auprc_gain=0.20,
            supervised_converged=False,
        )
        == "supervised_fit_unstable"
    )


@pytest.mark.parametrize(
    (
        "supervised",
        "signal_panel",
        "functional",
        "label_free_auroc",
        "label_free_auprc",
        "expected",
    ),
    [
        (0.05, 0.02, 0.01, 0.01, 0.02, "assay_limited_relative_to_raw"),
        (0.20, 0.02, 0.15, 0.01, 0.02, "shape_model_limited"),
        (0.20, 0.15, 0.02, 0.01, 0.02, "signal_combination_limited"),
        (0.20, 0.05, 0.05, 0.01, 0.02, "covariate_or_shape_model_limited"),
        (0.20, 0.05, 0.05, 0.03, 0.10, "detectable_above_raw"),
        (np.nan, 0.05, 0.05, 0.03, 0.10, "insufficient_supervised_folds"),
    ],
)
def test_raw_guarded_information_ceiling_classification(
    supervised: float,
    signal_panel: float,
    functional: float,
    label_free_auroc: float,
    label_free_auprc: float,
    expected: str,
) -> None:
    assert (
        evaluate_frozen_functional_information_ceiling.raw_guarded_failure_classification(
            supervised_relative_auprc_gain_over_raw=supervised,
            signal_panel_relative_auprc_gain_over_raw=signal_panel,
            functional_relative_auprc_gain_over_signal_panel=functional,
            label_free_auroc_gain_over_raw=label_free_auroc,
            label_free_relative_auprc_gain_over_raw=label_free_auprc,
        )
        == expected
    )


def test_raw_guarded_information_ceiling_marks_unstable_classifier() -> None:
    assert (
        evaluate_frozen_functional_information_ceiling.raw_guarded_failure_classification(
            supervised_relative_auprc_gain_over_raw=0.20,
            signal_panel_relative_auprc_gain_over_raw=0.15,
            functional_relative_auprc_gain_over_signal_panel=0.05,
            label_free_auroc_gain_over_raw=0.05,
            label_free_relative_auprc_gain_over_raw=0.20,
            supervised_converged=False,
        )
        == "supervised_fit_unstable"
    )


def test_bias_shrinkage_preserves_raw_at_zero_and_partially_subtracts() -> None:
    positions = np.arange(-40, 41, dtype=float)
    observed = np.full((1, len(positions)), 2.0)
    observed[:, np.abs(positions) <= 7] = 0.5
    expected = np.full((1, len(positions)), 0.5)
    expected[:, np.abs(positions) <= 7] = 1.0
    panel = {
        "observed": observed,
        "expected": {"parametric_lambda": expected},
        "positions": positions,
    }
    raw_score, raw_profile = evaluate_frozen_bias_shrinkage.choice_profile(
        panel,
        np.asarray([0]),
        {"source": "raw", "alpha": 0.0},
    )
    partial_score, partial_profile = evaluate_frozen_bias_shrinkage.choice_profile(
        panel,
        np.asarray([0]),
        {"source": "parametric_lambda", "alpha": 0.5},
    )
    assert np.array_equal(raw_profile, observed)
    assert np.allclose(partial_profile, observed - 0.5 * expected)
    assert raw_score.shape == partial_score.shape == (1,)


def test_bias_shrinkage_policy_selects_global_and_tf_abstention() -> None:
    rows = []
    task_values = (
        ("CellA", "CTCF", "CTCF", "positive_control"),
        ("CellB", "CTCF", "CTCF", "positive_control"),
        ("CellA", "MEF2A", "MEF2", "difficult"),
        ("CellB", "MEF2A", "MEF2", "difficult"),
        ("CellA", "MYC", "MYC_MAX", "difficult"),
        ("CellB", "MYC", "MYC_MAX", "difficult"),
    )
    for cell, tf, family, role in task_values:
        common = {
            "cell": cell,
            "tf": tf,
            "motif_family": family,
            "role": role,
            "status": "eligible",
        }
        rows.append(
            {
                **common,
                "method": "raw",
                "source": "raw",
                "alpha": 0.0,
                "auroc_gain_over_raw": 0.0,
                "relative_auprc_gain_over_raw": 0.0,
            }
        )
        for source, alpha, mef2_gain, myc_gain in (
            ("parametric_direct", 0.5, 0.04, -0.03),
            ("parametric_lambda", 0.8, 0.08, -0.01),
        ):
            gain = 0.01 if tf == "CTCF" else (mef2_gain if tf == "MEF2A" else myc_gain)
            rows.append(
                {
                    **common,
                    "method": "partial_bias_subtraction",
                    "source": source,
                    "alpha": alpha,
                    "auroc_gain_over_raw": gain,
                    "relative_auprc_gain_over_raw": gain,
                }
            )
    _global_rows, _tf_rows, global_choice, tf_choices = (
        evaluate_frozen_bias_shrinkage.select_policy_rows(pd.DataFrame(rows))
    )
    assert global_choice["source"] == "parametric_lambda"
    assert global_choice["alpha"] == 0.8
    assert tf_choices["MEF2A"]["source"] == "parametric_lambda"
    assert tf_choices["MYC"] == {
        "source": "raw",
        "alpha": 0.0,
        "selection_scope": "tf_across_cells",
        "validation_cells": 2,
        "mean_auroc_gain_over_raw": 0.0,
        "mean_relative_auprc_gain_over_raw": 0.0,
    }


def test_bias_shrinkage_naked_dna_uses_frozen_thresholds_and_zero_cut_support() -> None:
    positions = np.arange(-40, 41, dtype=float)
    observed = np.zeros((100, len(positions)), dtype=float)
    expected = np.zeros_like(observed)
    sites = pd.DataFrame(
        {
            "tf": ["TF1"] * 100,
            "TFBS_chr": ["chr1"] * 100,
            "TFBS_start": np.arange(100),
            "TFBS_end": np.arange(100) + 1,
            "TFBS_strand": ["+"] * 100,
        }
    )
    panel = {
        "sites": sites,
        "observed": observed,
        "expected": {
            "parametric_direct": expected,
            "parametric_lambda": expected,
        },
        "dwm_expected": expected,
        "positions": positions,
    }
    methods = (
        "raw",
        "DWM_conventional_deviance",
        "frozen_global_shrinkage",
        "frozen_tf_specific_shrinkage",
    )
    thresholds = pd.DataFrame(
        {
            "cell": ["CellA"] * len(methods),
            "tf": ["TF1"] * len(methods),
            "method": methods,
            "threshold": [1.0] * len(methods),
        }
    )
    policy = {
        "global_choice": {"source": "raw", "alpha": 0.0},
        "per_tf_choices": {"TF1": {"source": "raw", "alpha": 0.0}},
        "dispersion": 1.0,
    }
    rates, scores = evaluate_frozen_bias_shrinkage.evaluate_naked_dna(
        {"CellA": panel},
        policy,
        thresholds,
        replicate="NakedDNA_rep2",
    )
    assert len(scores) == 100
    assert not any("label" in column.lower() for column in scores)
    assert set(rates["finite_sites"]) == {100}
    assert set(rates["informative_sites"]) == {0}
    assert set(rates["calls"]) == {0}
    assert rates["passes_safety"].all()


def test_derived_family_selection_filters_candidates_and_profiles() -> None:
    metrics = pd.DataFrame(
        {
            "cell": ["CellA"] * 3,
            "tf": ["TF1"] * 3,
            "bias_configuration": ["LOG21"] * 3,
            "candidate_id": ["count_a", "count_b", "fda_a"],
            "family": ["count", "count", "fda"],
            "status": ["ok", "ok", "ok"],
            "converged": [True, True, True],
            "selection_score": [0.8, 0.9, 1.0],
            "auprc": [0.6, 0.7, 0.8],
            "auroc": [0.6, 0.7, 0.8],
            "validation_positive_sites": [200, 200, 200],
            "validation_negative_sites": [200, 200, 200],
        }
    )
    profiles = pd.DataFrame(
        {
            "candidate_id": ["count_a", "count_b", "fda_a"],
            "position": [0, 0, 0],
        }
    )
    selected, selected_profiles, winners = (
        derive_label_free_family_selection.select_families(
            metrics,
            profiles,
            {"count"},
            minimum_sites_per_class=200,
        )
    )
    assert set(selected["candidate_id"]) == {"count_a", "count_b"}
    assert set(selected_profiles["candidate_id"]) == {"count_a", "count_b"}
    assert winners["candidate_id"].tolist() == ["count_b"]


def test_frozen_call_threshold_respects_ties_and_target_rate() -> None:
    scores = np.asarray([0.1, 0.2, 0.3, 0.8, 0.8, 0.9])
    threshold, calls = freeze_functional_call_thresholds.upper_tail_threshold(
        scores,
        0.34,
    )
    assert threshold == 0.9
    assert calls == 1
    assert np.sum(scores >= threshold) / len(scores) <= 0.34


def test_validation_site_score_frame_preserves_common_site_hashes() -> None:
    sites = pd.DataFrame(
        {
            "TFBS_chr": ["chr16", "chr17", "chr18"],
            "chip_label": [0, 1, 0],
        }
    )
    result = freeze_functional_call_thresholds.validation_site_score_frame(
        record={
            "cell": "CellA",
            "tf": "TF1",
            "motif_family": "F1",
            "bias_configuration": "LOG21",
        },
        candidate_id="candidate",
        sites=sites,
        indexes=np.asarray([2, 0]),
        site_hash=np.asarray([91, 17], dtype=np.uint64),
        candidate_score=np.asarray([0.9, 0.1]),
        dwm_score=np.asarray([1.2, -0.4]),
    )
    assert result["site_hash"].tolist() == [91, 17]
    assert result["label"].tolist() == [0, 0]
    assert result["candidate_score"].tolist() == [0.9, 0.1]


def test_naked_dna_rate_keeps_zero_cut_sites_in_finite_denominator() -> None:
    score = np.asarray([0.9, 0.1, 0.0, 0.0])
    valid = np.ones(4, dtype=bool)
    informative = np.asarray([True, True, False, False])
    record, calls = evaluate_frozen_functional_naked_dna.rate_record(
        score,
        valid,
        informative,
        threshold=0.8,
    )
    assert record["valid_sites"] == 4
    assert record["informative_sites"] == 2
    assert record["calls"] == 1
    assert record["false_positive_rate"] == 0.25
    assert record["informative_false_positive_rate"] == 0.5
    assert calls.tolist() == [True, False, False, False]


def test_label_free_motif_pool_is_deterministic_and_checksum_locked(
    tmp_path: Path,
) -> None:
    study = tmp_path / "study.json"
    study.write_text(
        json.dumps(
            {
                "tasks": [
                    {
                        "cell": "CellA",
                        "tf": "TF1",
                        "motif_id": "MA0001.1",
                        "motif_family": "FAMILY1",
                        "split": "development",
                    }
                ],
                "chromosome_split": {
                    "train": ["chr1"],
                    "validation": ["chr2"],
                    "test": ["chr3"],
                },
            }
        ),
        encoding="utf-8",
    )
    source = tmp_path / "motifs.tsv"
    pd.DataFrame(
        {
            "motif": ["TF1_MA0001.1"] * 120,
            "TFBS_chr": ["chr1"] * 120,
            "TFBS_start": np.arange(120) * 10,
            "TFBS_end": np.arange(120) * 10 + 6,
            "TFBS_strand": ["+"] * 120,
            "TFBS_score": np.linspace(1.0, 2.0, 120),
        }
    ).to_csv(source, sep="\t", index=False)
    outdir = tmp_path / "pools"
    arguments = [
        "--study",
        str(study),
        "--source",
        f"CellA={source}",
        "--maximum-per-tf",
        "100",
        "--outdir",
        str(outdir),
    ]
    assert sample_label_free_motif_sites.main(arguments) == 0
    output = outdir / "CellA.unlabeled_training_sites.tsv.gz"
    first_hash = sample_label_free_motif_sites.file_sha256(output)
    assert len(pd.read_csv(output, sep="\t")) == 100
    assert sample_label_free_motif_sites.main(arguments) == 0
    assert sample_label_free_motif_sites.file_sha256(output) == first_hash
    output.write_bytes(output.read_bytes() + b"tampered")
    with pytest.raises(ValueError, match="checksum mismatch"):
        sample_label_free_motif_sites.main(arguments)


def test_factorization_dwm_loader_accepts_verified_cache_and_rejects_parametric_label(
    tmp_path: Path,
) -> None:
    cache = tmp_path / "DWM.expected.npz"
    np.savez_compressed(
        cache,
        profiles=np.ones((2, 5)),
        valid=np.ones(2, dtype=bool),
        site_hash=np.asarray([1, 2], dtype=np.uint64),
        signal_identity=np.asarray("/run/fp_tools_dwm/sample_expected.bw:1:2"),
    )
    baseline, inputs = evaluate_parametric_factorization.load_dwm_baseline(cache)
    assert baseline["expected"].shape == (2, 5)
    assert not bool(baseline["orientation_aligned"])
    assert inputs == [cache]

    invalid = tmp_path / "LOG81.expected.npz"
    np.savez_compressed(
        invalid,
        profiles=np.ones((2, 5)),
        valid=np.ones(2, dtype=bool),
        site_hash=np.asarray([1, 2], dtype=np.uint64),
        signal_identity=np.asarray("/run/loglinear81/sample_expected.bw:1:2"),
    )
    try:
        evaluate_parametric_factorization.load_dwm_baseline(invalid)
    except ValueError as error:
        assert "conventional DWM" in str(error)
    else:
        raise AssertionError("parametric baseline was mislabeled as DWM")


def test_factorization_pwm_loader_requires_pwm_identity(tmp_path: Path) -> None:
    cache = tmp_path / "PWM.expected.npz"
    np.savez_compressed(
        cache,
        profiles=np.ones((2, 5)),
        valid=np.ones(2, dtype=bool),
        site_hash=np.asarray([1, 2], dtype=np.uint64),
        signal_identity=np.asarray("/run/fp_tools_pwm/sample_expected.bw:1:2"),
    )
    baseline, inputs = evaluate_parametric_factorization.load_pwm_baseline(cache)
    assert baseline["expected"].shape == (2, 5)
    assert not bool(baseline["orientation_aligned"])
    assert inputs == [cache]


def test_direct_baseline_is_oriented_after_hash_alignment(tmp_path: Path) -> None:
    sites = pd.DataFrame(
        {
            "cell": ["K562", "K562"],
            "tf": ["TF", "TF"],
            "motif_family": ["family", "family"],
            "motif": ["M1", "M1"],
            "role": ["difficult", "difficult"],
            "TFBS_chr": ["chr1", "chr1"],
            "TFBS_start": [100, 200],
            "TFBS_end": [101, 201],
            "TFBS_strand": ["+", "-"],
            "motif_score": [1.0, 1.0],
            "chip_label": [0, 1],
            "chromosome_split": ["train", "train"],
        }
    )
    hashes = np.asarray([1, 2], dtype=np.uint64)
    cache = tmp_path / "K562.DWM.npz"
    np.savez_compressed(
        cache,
        profiles=np.asarray([[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]], dtype=float),
        valid=np.ones(2, dtype=bool),
        site_hash=hashes,
        signal_identity=np.asarray("/run/fp_tools_dwm/sample_expected.bw:1:2"),
    )
    baseline, _inputs = evaluate_parametric_factorization.load_dwm_baseline(cache)
    expected, _valid = evaluate_parametric_factorization.align_baseline(
        {"site_hash": hashes}, baseline
    )
    oriented = evaluate_parametric_factorization.orient_aligned_baseline(
        expected, baseline, sites
    )
    assert oriented.tolist() == [
        [1.0, 2.0, 3.0, 4.0, 5.0],
        [10.0, 9.0, 8.0, 7.0, 6.0],
    ]


def test_residual_selector_uses_difficult_tasks_and_ctcf_gate() -> None:
    rows = []
    for residual, difficult_ap, ctcf_auc in (
        ("deviance", 0.70, 0.79),
        ("pearson", 0.71, 0.75),
        ("difference", 0.69, 0.80),
        ("log-ratio", 0.68, 0.80),
        ("nb-center-flank", 0.67, 0.80),
    ):
        rows.extend(
            [
                {
                    "cell": "K562",
                    "tf": "MEF2A",
                    "role": "weak_shape",
                    "method": f"factorized_residual_{residual}",
                    "auroc": 0.7,
                    "auprc": difficult_ap,
                },
                {
                    "cell": "K562",
                    "tf": "CTCF",
                    "role": "positive_control",
                    "method": f"factorized_residual_{residual}",
                    "auroc": ctcf_auc,
                    "auprc": 0.8,
                },
            ]
        )
    rows.extend(
        [
            {
                "cell": "K562",
                "tf": "MEF2A",
                "role": "weak_shape",
                "method": "DWM",
                "auroc": 0.6,
                "auprc": 0.5,
            },
            {
                "cell": "K562",
                "tf": "CTCF",
                "role": "positive_control",
                "method": "DWM",
                "auroc": 0.8,
                "auprc": 0.8,
            },
        ]
    )
    selected, summary = evaluate_parametric_factorization.select_residual(
        pd.DataFrame(rows)
    )
    # Pearson has the best difficult-task AP but violates the -0.02 CTCF gate.
    assert selected == "deviance"
    assert (
        summary.loc[summary["residual"] == "pearson", "passes_ctcf_gate"].item()
        is False
    )


def test_block_bootstrap_caches_equivalent_chromosome_multisets() -> None:
    sites = pd.DataFrame(
        {
            "TFBS_chr": np.repeat(["chr16", "chr17", "chr18"], 20),
            "chip_label": np.tile([0, 1], 30),
        }
    )
    baseline = np.linspace(0.0, 1.0, len(sites))
    candidate = baseline + np.tile([-0.05, 0.05], 30)
    result = evaluate_parametric_factorization.block_bootstrap_delta(
        sites,
        candidate,
        baseline,
        iterations=1000,
        seed=2026,
    )
    assert result["bootstrap_successful"] == 1000
    # Three chromosomes sampled three times have only C(5, 3)=10 multisets.
    assert result["bootstrap_unique_resamples"] <= 10


def test_test_input_freeze_is_immutable_and_records_candidate_signature(
    tmp_path: Path,
) -> None:
    candidate = tmp_path / "K562.json"
    candidate.write_text(
        json.dumps(
            {
                "schema": "fp-tools-strand-functional-profiles-v1",
                "metadata": {
                    "bias_model_sha256": "bias",
                    "genome_sha256": "genome",
                    "read_shift": [4, -5],
                    "flank": 100,
                },
            }
        )
    )
    signature = evaluate_parametric_factorization.profile_model_signature(
        tmp_path / "K562"
    )
    assert signature["read_shift"] == [4, -5]
    test_input = tmp_path / "test-input.npz"
    np.savez_compressed(test_input, values=np.ones(2))
    configuration = tmp_path / "safe.json"
    configuration.write_text("{}")
    output = tmp_path / "test.freeze.json"
    document = evaluate_parametric_factorization.write_test_input_freeze(
        output,
        configuration_path=configuration,
        configuration={"configuration_id": "locked"},
        inputs=[test_input],
        candidate_signatures={"K562": signature},
    )
    assert document["refitted"] is False
    assert document["thresholds_changed"] is False
    evaluate_parametric_factorization.write_test_input_freeze(
        output,
        configuration_path=configuration,
        configuration={"configuration_id": "locked"},
        inputs=[test_input],
        candidate_signatures={"K562": signature},
    )
    np.savez_compressed(test_input, values=np.zeros(2))
    try:
        evaluate_parametric_factorization.write_test_input_freeze(
            output,
            configuration_path=configuration,
            configuration={"configuration_id": "locked"},
            inputs=[test_input],
            candidate_signatures={"K562": signature},
        )
    except ValueError as error:
        assert "immutable test-input freeze differs" in str(error)
    else:  # pragma: no cover
        raise AssertionError("changed test input unexpectedly reused its freeze")


def test_bigwig_integrity_rejects_empty_and_accepts_covered(tmp_path: Path) -> None:
    pybigwig = __import__("pyBigWig")
    empty = tmp_path / "empty.bw"
    handle = pybigwig.open(str(empty), "w")
    handle.addHeader([("chr1", 100)])
    handle.close()
    try:
        run_frozen_parametric_experiment.validate_bigwig(empty)
    except ValueError as error:
        assert "no covered bases" in str(error)
    else:  # pragma: no cover - documents the fail-closed contract
        raise AssertionError("empty bigWig unexpectedly passed validation")

    covered = tmp_path / "covered.bw"
    handle = pybigwig.open(str(covered), "w")
    handle.addHeader([("chr1", 100)])
    handle.addEntries(["chr1"], [10], ends=[20], values=[1.5])
    handle.close()
    result = run_frozen_parametric_experiment.validate_bigwig(covered)
    assert result["covered_bases"] == 10


def test_stage_runner_completes_and_resumes_synthetic(tmp_path: Path) -> None:
    runner = run_frozen_parametric_experiment.StageRunner(
        study_path=ROOT
        / "benchmarks/manifests/frozen_parametric_factorization_v1.spec.json",
        holdout_freeze_path=ROOT
        / "benchmarks/manifests/frozen_parametric_factorization_v1.freeze.json",
        registry_path=ROOT
        / "benchmarks/manifests/frozen_parametric_local_inputs_v1.json",
        outdir=tmp_path / "run",
    )
    assert runner.run_stage("synthetic") == "completed"
    assert runner.run_stage("synthetic") == "resumed_verified"
    completion = json.loads(
        (tmp_path / "run/synthetic/stage.complete.json").read_text()
    )
    assert (
        completion["schema"] == run_frozen_parametric_experiment.STAGE_COMPLETE_SCHEMA
    )
    assert all(Path(row["path"]).is_file() for row in completion["outputs"])
