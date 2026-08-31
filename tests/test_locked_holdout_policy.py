import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pysam
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "evaluate_locked_holdout_policy.py"
spec = importlib.util.spec_from_file_location("evaluate_locked_holdout_policy", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_frozen_routes_resolve_to_implemented_candidate_families():
    compact = ROOT / "benchmarks" / "manifests" / "compact"
    routes = pd.read_csv(compact / "functional_holdout_routes_v1.tsv", sep="\t")
    policy = pd.read_csv(compact / "functional_detector_policy_v1.tsv", sep="\t")
    strand = module.strand_candidates()
    combined = module.combined_candidates()
    for route in routes.itertuples(index=False):
        if route.bias_configuration == "DWM":
            assert route.candidate_id in combined
        else:
            assert route.candidate_id in strand
        assert module.reference_candidate(policy, route) in combined


def test_validate_tables_requires_exact_preregistered_accessions():
    study = json.loads(
        (ROOT / "benchmarks" / "manifests" / "footprint_functional_v1.spec.json")
        .read_text(encoding="utf-8")
    )
    compact = ROOT / "benchmarks" / "manifests" / "compact"
    routes = pd.read_csv(compact / "functional_holdout_routes_v1.tsv", sep="\t")
    policy = pd.read_csv(compact / "functional_detector_policy_v1.tsv", sep="\t")
    chip = pd.read_csv(compact / "functional_holdout_chip_peaks.tsv", sep="\t")
    tasks = module.validate_tables(study, routes, policy, chip)
    assert len(tasks) == len(routes) == len(chip) == 15
    tampered = chip.copy()
    tampered.loc[0, "file_accession"] = "ENCFF000AAA"
    with pytest.raises(ValueError, match="preregistered study"):
        module.validate_tables(study, routes, policy, tampered)


def _synthetic_sites() -> pd.DataFrame:
    centers = [500, 800, 1100, 1400, 2500, 2800, 3100, 3400, 3700, 4000]
    return pd.DataFrame(
        {
            "cell": "Fixture",
            "tf": "TFX",
            "motif_id": "MA0000.1",
            "motif_family": "FIX",
            "TFBS_chr": "chr1",
            "TFBS_start": [center - 5 for center in centers],
            "TFBS_end": [center + 5 for center in centers],
            "TFBS_strand": ["+", "-", "+", "-", "+", "-", "+", "-", "+", "-"],
            "motif_score": np.linspace(5.0, 6.0, len(centers)),
            "peak_start": [center - 100 for center in centers],
            "peak_end": [center + 100 for center in centers],
            "chromosome_split": "test",
        }
    )


def test_build_matched_test_sites_uses_summits_and_all_frozen_covariates(tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 1300 + "\n", encoding="utf-8")
    pysam.faidx(str(fasta))
    peaks = tmp_path / "chip.narrowPeak"
    lines = []
    for center in (500, 800, 1100, 1400):
        start = center - 50
        lines.append(
            f"chr1\t{start}\t{center + 50}\tpeak\t1000\t.\t10\t10\t10\t50\n"
        )
    peaks.write_text("".join(lines), encoding="utf-8")
    sites = _synthetic_sites()
    observed = np.ones((len(sites), 201), dtype=float)
    observed *= np.linspace(1.0, 1.1, len(sites))[:, None]
    matched, diagnostics = module.build_matched_test_sites(
        sites,
        np.ones(len(sites), dtype=bool),
        observed,
        peaks,
        fasta,
        positive_summit_distance=25,
        negative_peak_distance=500,
        seed=7,
    )
    assert diagnostics["positive_summit_supported"] == 4
    assert diagnostics["matched_per_class"] == 4
    assert matched["label"].value_counts().to_dict() == {0: 4, 1: 4}
    assert set(module.MATCH_FEATURES).issubset(matched.columns)
    assert np.isfinite(matched[list(module.MATCH_FEATURES)].to_numpy()).all()


def test_bootstrap_delta_is_paired_by_chromosome():
    frame = pd.DataFrame(
        {
            "TFBS_chr": np.repeat(["chr19", "chr20", "chr21"], 4),
            "label": [0, 0, 1, 1] * 3,
        }
    )
    candidate = np.tile([0.1, 0.2, 0.8, 0.9], 3)
    reference = np.tile([0.2, 0.8, 0.3, 0.7], 3)
    summary = module.bootstrap_delta(
        frame, candidate, reference, iterations=50, seed=4
    )
    assert summary["bootstrap_successful"] == 50
    assert summary["auroc_gain_bootstrap_probability_positive"] == 1.0
    assert summary["relative_auprc_gain_bootstrap_probability_positive"] == 1.0


def test_freeze_identifier_detects_tampering():
    document = {
        "schema": module.FREEZE_SCHEMA,
        "locked_holdout_labels_read": False,
        "options": {"seed": 2026},
    }
    identifier = module.canonical_hash(document)
    assert identifier == module.canonical_hash(document)
    document["options"]["seed"] = 2027
    assert identifier != module.canonical_hash(document)


def test_candidate_artifact_may_be_an_exact_dwm_site_subset():
    frame = _synthetic_sites()
    hashes = module.site_hashes(frame)
    arrays = {"valid": np.ones(len(frame), dtype=bool), "site_hash": hashes}
    reference = module.Artifact(
        "DWM",
        "Fixture",
        Path("dwm.json"),
        {"schema": "fp-tools-combined-functional-profiles-v1"},
        frame,
        arrays,
    )
    subset_rows = np.array([0, 2, 5, 8])
    candidate = module.Artifact(
        "NEW",
        "Fixture",
        Path("new.json"),
        {"schema": "fp-tools-strand-functional-profiles-v1"},
        frame.iloc[subset_rows].reset_index(drop=True),
        {"valid": np.ones(len(subset_rows), dtype=bool), "site_hash": hashes[subset_rows]},
    )
    module.validate_exact_site_alignment(
        {("DWM", "Fixture"): reference, ("NEW", "Fixture"): candidate}
    )
    mapping = module.reference_to_candidate_indexes(reference, candidate)
    assert mapping.tolist() == [0, -1, 1, -1, -1, 2, -1, -1, 3, -1]
    projected = module.valid_on_reference(reference, candidate)
    assert projected.tolist() == [True, False, True, False, False, True, False, False, True, False]


def test_replicate_groups_require_two_complete_site_aligned_replicates():
    frame = _synthetic_sites()
    hashes = module.site_hashes(frame)
    arrays = {"valid": np.ones(len(frame), dtype=bool), "site_hash": hashes}
    pooled = module.Artifact(
        "DWM",
        "Fixture",
        Path("pooled.json"),
        {"schema": "fp-tools-combined-functional-profiles-v1"},
        frame,
        arrays,
    )
    rep1 = module.Artifact(
        "DWM",
        "Fixture",
        Path("rep1.json"),
        {"schema": "fp-tools-combined-functional-profiles-v1"},
        frame.copy(),
        {"valid": np.ones(len(frame), dtype=bool), "site_hash": hashes.copy()},
    )
    rep2 = module.Artifact(
        "DWM",
        "Fixture",
        Path("rep2.json"),
        {"schema": "fp-tools-combined-functional-profiles-v1"},
        frame.copy(),
        {"valid": np.ones(len(frame), dtype=bool), "site_hash": hashes.copy()},
    )
    routes = pd.DataFrame(
        {
            "cell": ["Fixture"],
            "tf": ["TFX"],
            "bias_configuration": ["DWM"],
        }
    )
    module.validate_replicate_artifacts(
        {("DWM", "Fixture"): pooled},
        {
            ("rep1", "DWM", "Fixture"): rep1,
            ("rep2", "DWM", "Fixture"): rep2,
        },
        routes,
    )
    mapped = module.map_indexes_by_hash(pooled, rep2, np.array([1, 4, 7]))
    assert mapped.tolist() == [1, 4, 7]
    with pytest.raises(ValueError, match="at least two"):
        module.validate_replicate_artifacts(
            {("DWM", "Fixture"): pooled},
            {("rep1", "DWM", "Fixture"): rep1},
            routes,
        )


def test_frozen_combined_count_and_strand_routes_fit_without_labels():
    rng = np.random.default_rng(12)
    rows, width = 180, 201
    positions = np.arange(-100, 101, dtype=float)
    sites = pd.DataFrame(
        {
            "cell": "Fixture",
            "tf": "TFX",
            "motif_id": "MA0000.1",
            "motif_family": "FIX",
            "TFBS_chr": "chr1",
            "TFBS_start": np.arange(rows) * 10 + 1000,
            "TFBS_end": np.arange(rows) * 10 + 1005,
            "TFBS_strand": "+",
            "motif_score": rng.normal(size=rows),
            "peak_start": np.arange(rows) * 10 + 900,
            "peak_end": np.arange(rows) * 10 + 1100,
            "chromosome_split": ["train"] * 150 + ["test"] * 30,
        }
    )
    expected = np.full((rows, width), 2.0)
    observed = rng.poisson(expected).astype(float)
    observed[:75, 97:104] = rng.poisson(0.5, size=(75, 7))
    residual = module.deviance_profiles(observed, expected, 0.0)
    common = {
        "combined_residual": residual,
        "valid": np.ones(rows, dtype=bool),
        "site_hash": np.arange(rows, dtype=np.uint64),
    }
    combined = module.Artifact(
        "DWM",
        "Fixture",
        Path("fixture.json"),
        {"schema": "fp-tools-combined-functional-profiles-v1"},
        sites,
        {"observed": observed, "expected": expected, **common},
    )
    strand = module.Artifact(
        "MT_SELMA10_4m4",
        "Fixture",
        Path("fixture.json"),
        {"schema": "fp-tools-strand-functional-profiles-v1"},
        sites,
        {
            "plus_observed": observed / 2,
            "minus_observed": observed / 2,
            "plus_expected": expected / 2,
            "minus_expected": expected / 2,
            "shared_strand_residual": residual,
            "antisymmetric_strand_residual": residual,
            **common,
        },
    )
    evaluation = np.arange(150, 180)
    combined_result = module.fit_combined_route(
        combined,
        "spline.bg_gp-long.prior_free.pen_10.shrink_50",
        "TFX",
        "FIX",
        evaluation,
        positions,
        10_000,
        2026,
        {},
    )
    strand_result = module.fit_strand_route(
        strand,
        "count_gp.bg_gp-long.window_30",
        "TFX",
        "FIX",
        evaluation,
        positions,
        10_000,
        2026,
        {},
    )
    assert combined_result.probabilities.shape == (30,)
    assert strand_result.probabilities.shape == (30,)
    assert np.isfinite(combined_result.probabilities).all()
    assert np.isfinite(strand_result.probabilities).all()
    assert combined_result.converged and strand_result.converged


def test_promotion_tables_preserve_paired_metrics_and_fail_closed_stability():
    metrics = pd.DataFrame(
        [
            {
                "cell": "Fixture",
                "tf": "TFX",
                "motif_id": "MA0000.1",
                "motif_family": "FIX",
                "role": "difficult",
                "status": "ok",
                "bias_configuration": "NEW",
                "candidate_id": "hybrid.fixture",
                "reference_candidate_id": "spline.fixture",
                "candidate_auroc": 0.72,
                "candidate_auprc": 0.41,
                "candidate_brier": 0.20,
                "candidate_ece": 0.06,
                "reference_auroc": 0.66,
                "reference_auprc": 0.34,
                "reference_brier": 0.23,
                "reference_ece": 0.09,
                "candidate_positive_depletion": 1.4,
                "candidate_negative_depletion": 0.4,
                "reference_positive_depletion": 1.0,
                "reference_negative_depletion": 0.5,
                "replicate_direction_stable": False,
                "biological_replicates": 2,
            }
        ]
    )
    long_metrics = module.promotion_metric_table(metrics)
    assert long_metrics["method"].tolist() == [
        module.PROMOTION_CANDIDATE,
        module.PROMOTION_REFERENCE,
    ]
    assert long_metrics["auroc"].tolist() == [0.72, 0.66]
    descriptors = module.promotion_descriptor_table(metrics)
    assert len(descriptors) == 4
    stability = module.promotion_stability_table(metrics)
    assert stability["candidate_id"].tolist() == [module.PROMOTION_CANDIDATE]
    assert not bool(stability.loc[0, "direction_consistent"])
