from __future__ import annotations

from pathlib import Path
import json
import sys

import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from evaluate_functional_promotion import (  # noqa: E402
    evaluate_promotion,
    prepare_pairs,
)


SPEC = Path(__file__).resolve().parents[1] / "benchmarks" / "manifests" / "footprint_functional_v1.spec.json"


def _study_and_pairs() -> tuple[dict, pd.DataFrame]:
    study = json.loads(SPEC.read_text(encoding="utf-8"))
    tasks = pd.DataFrame(study["tasks"])
    tasks = (
        tasks[tasks["split"] == "development"]
        .groupby("cell", sort=True, group_keys=False)
        .head(4)
        .reset_index(drop=True)
    )
    tasks.loc[:, "role"] = ["positive_control", "difficult"] * 4
    tasks.loc[:, "motif_family"] = [f"family_{index}" for index in range(8)]
    metrics = []
    for row in tasks.itertuples(index=False):
        for candidate, auroc, auprc in (
            ("DWM:spline", 0.70, 0.30),
            ("LOG:gp", 0.76, 0.38),
        ):
            metrics.append(
                {
                    "cell": row.cell,
                    "tf": row.tf,
                    "candidate_id": candidate,
                    "auroc": auroc,
                    "auprc": auprc,
                    "split": "validation",
                }
            )
    custom_study = dict(study)
    custom_study["tasks"] = tasks.to_dict("records")
    return custom_study, prepare_pairs(
        pd.DataFrame(metrics), custom_study, "LOG:gp", "DWM:spline", "development"
    )


def _evidence(pairs: pd.DataFrame):
    descriptors = []
    for row in pairs.itertuples(index=False):
        for correction, positive, negative in (("DWM", 1.0, 0.5), ("LOG", 1.8, 0.5)):
            for group, depletion in (("chip_positive", positive), ("matched_negative", negative)):
                descriptors.append(
                    {
                        "cell": row.cell,
                        "tf": row.tf,
                        "motif_family": row.motif_family,
                        "correction": correction,
                        "group": group,
                        "depletion": depletion,
                    }
                )
    negative = pd.DataFrame(
        {"candidate_id": ["DWM:spline", "LOG:gp"], "false_positive_rate": [0.03, 0.035]}
    )
    resources = pd.DataFrame(
        {
            "candidate_id": ["DWM:spline", "LOG:gp"],
            "runtime_seconds": [100.0, 150.0],
            "peak_memory_mb": [1000.0, 1200.0],
            "model_size_mb": [0.0, 2.0],
        }
    )
    uncertainty = pd.DataFrame(
        {"candidate_id": ["LOG:gp"], "empirical_coverage": [0.90]}
    )
    stability = pd.DataFrame(
        {"candidate_id": ["LOG:gp", "LOG:gp"], "direction_consistent": [True, True]}
    )
    leakage = pd.DataFrame(
        {
            "candidate_id": ["LOG:gp", "LOG:gp"],
            "potential_motif_response_requires_review": [False, False],
        }
    )
    complexity = pd.DataFrame(
        {
            "candidate_id": ["LOG:gp", "LOG:gp"],
            "relative_auprc_gain_over_spline": [0.08, 0.09],
            "uncertainty_calibration_improved": [False, False],
        }
    )
    return (
        pd.DataFrame(descriptors),
        negative,
        resources,
        uncertainty,
        stability,
        leakage,
        complexity,
    )


def test_functional_promotion_passes_only_with_complete_evidence() -> None:
    study, pairs = _study_and_pairs()
    evidence = _evidence(pairs)
    _pairs, separation, summary = evaluate_promotion(
        pairs,
        study,
        candidate="LOG:gp",
        baseline="DWM:spline",
        descriptors=evidence[0],
        negative_controls=evidence[1],
        resources=evidence[2],
        uncertainty=evidence[3],
        stability=evidence[4],
        leakage=evidence[5],
        complexity=evidence[6],
        locked_holdout_scored=True,
        bootstrap=200,
        seed=4,
    )
    assert len(separation)
    assert summary["passed"]
    assert all(summary["checks"].values())


def test_functional_promotion_fails_closed_when_evidence_is_missing() -> None:
    study, pairs = _study_and_pairs()
    _pairs, _separation, summary = evaluate_promotion(
        pairs,
        study,
        candidate="LOG:gp",
        baseline="DWM:spline",
        bootstrap=50,
        seed=4,
    )
    assert not summary["passed"]
    assert not summary["checks"]["naked_dna_false_positive_control"]
    assert not summary["checks"]["functional_uncertainty_coverage"]


def test_prepare_pairs_rejects_incomplete_frozen_task_coverage() -> None:
    study, _pairs = _study_and_pairs()
    tasks = pd.DataFrame(study["tasks"])
    metrics = []
    for row in tasks.iloc[:-1].itertuples(index=False):
        for candidate, auroc, auprc in (
            ("DWM:spline", 0.70, 0.30),
            ("LOG:gp", 0.76, 0.38),
        ):
            metrics.append(
                {
                    "cell": row.cell,
                    "tf": row.tf,
                    "candidate_id": candidate,
                    "auroc": auroc,
                    "auprc": auprc,
                    "split": "validation",
                }
            )
    with pytest.raises(ValueError, match="exactly cover"):
        prepare_pairs(
            pd.DataFrame(metrics),
            study,
            "LOG:gp",
            "DWM:spline",
            "development",
        )
