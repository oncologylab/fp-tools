from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))


def audit_module():
    import audit_chrombpnet_comparison

    return audit_chrombpnet_comparison


def test_checksum_validation_uses_declared_data_root_and_rejects_mutation(tmp_path):
    audit = audit_module()
    path = tmp_path / "data.tsv"
    path.write_bytes(b"frozen input\n")
    record = {
        "path": "data.tsv",
        "sha256": hashlib.sha256(b"frozen input\n").hexdigest(),
    }
    assert audit.checked_path(record, tmp_path) == path
    path.write_bytes(b"changed input\n")
    with pytest.raises(ValueError, match="checksum"):
        audit.checked_path(record, tmp_path)


def test_audit_rejects_unopened_cells_before_reading_label_outputs(tmp_path):
    audit = audit_module()
    freeze = {
        "schema": "fp-tools-chrombpnet-comparison-inputs-v1",
        "cells": ["SK-N-SH"],
        "test_labels_previously_opened_for_other_models": False,
    }
    freeze["input_id"] = audit.canonical_id(freeze)
    (tmp_path / "chrombpnet_comparison_inputs.freeze.json").write_text(
        json.dumps(freeze)
    )
    path = tmp_path / "chrombpnet_comparison_manifest.json"
    document = {
        "schema": "fp-tools-chrombpnet-frozen-comparison-v1",
        "input_id": freeze["input_id"],
        "outputs": {"site_scores": {"path": "DO_NOT_OPEN_LABELS"}},
    }
    document["result_id"] = audit.canonical_id(document)
    path.write_text(json.dumps(document))
    with pytest.raises(ValueError, match="previously opened development"):
        audit.load_comparison(path, tmp_path)


def profile_rows(center, shoulders):
    positions = np.arange(-100, 101)
    values = np.zeros(len(positions))
    values[np.abs(positions) <= 5] = center
    values[(np.abs(positions) >= 15) & (np.abs(positions) <= 40)] = shoulders
    return pd.DataFrame(
        {
            "cell": "HepG2",
            "tf": "CTCF",
            "method": "frozen_parametric_bias_conventional_geometry",
            "position": positions,
            "positive_mean": values,
            "negative_mean": 0.0,
            "positive_minus_negative": values,
            "lower_95": values - 0.05,
            "upper_95": values + 0.05,
        }
    )


def test_center_enrichment_cannot_pass_shape_guard_despite_positive_depletion():
    audit = audit_module()
    row = audit.shape_audit(profile_rows(center=0.2, shoulders=0.4)).iloc[0]
    assert row["depletion"] == pytest.approx(0.2)
    assert not row["canonical_protection"]
    assert row["significant_center_points"] == 0


def test_shape_guard_preserves_signed_center_and_both_shoulders():
    audit = audit_module()
    row = audit.shape_audit(profile_rows(center=-0.2, shoulders=0.1)).iloc[0]
    assert row["canonical_protection"]
    assert row["center_contrast"] == pytest.approx(-0.2)
    assert row["depletion"] == pytest.approx(0.3)
    assert row["left_shoulder"] == pytest.approx(0.1)
    assert row["right_shoulder"] == pytest.approx(0.1)
    assert row["significant_center_points"] == 11


def test_balance_joins_sites_by_identity_and_reports_label_dependent_attrition():
    audit = audit_module()
    source = pd.DataFrame(
        {
            "cell": ["HepG2"] * 4,
            "tf": ["CTCF"] * 4,
            # In actual benchmark tables site_id is the motif name, not a
            # unique genomic locus. Identity includes coordinates and strand.
            "site_id": ["CTCF_MA0139.2"] * 4,
            "TFBS_chr": ["chr19"] * 4,
            "TFBS_start": [100, 200, 300, 400],
            "TFBS_end": [119, 219, 319, 419],
            "TFBS_strand": ["+", "-", "+", "-"],
            "chip_label": [1, 1, 0, 0],
            "motif_score": [2.0, 4.0, 2.0, 4.0],
            "accessibility": [10.0, 20.0, 10.0, 20.0],
        }
    )
    retained = source.iloc[[2, 1, 0]].drop(columns=["motif_score", "accessibility"])
    totals = np.asarray([1.0, 2.0, 10.0, 20.0])
    table = audit.balance_audit(source, retained, totals)
    row = table.iloc[0]
    assert row["positive_retained_fraction"] == 1.0
    assert row["negative_retained_fraction"] == 0.5
    assert row["retention_difference"] == 0.5
    assert row["motif_score_smd_before"] == 0.0
    assert row["raw_total_positive_mean"] == 1.5
    assert row["raw_total_negative_mean"] == 10.0
    duplicate = pd.concat([retained, retained.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate"):
        audit.balance_audit(source, duplicate, totals)
    inconsistent = retained.copy()
    inconsistent.iloc[0, inconsistent.columns.get_loc("chip_label")] = 1
    with pytest.raises(ValueError, match="label"):
        audit.balance_audit(source, inconsistent, totals)


def test_zero_variance_covariate_difference_is_not_reported_as_balanced():
    audit = audit_module()
    assert np.isinf(audit.standardized_difference([2.0, 2.0], [1.0, 1.0]))
    assert audit.standardized_difference([1.0, 1.0], [1.0, 1.0]) == 0.0


def test_large_nondirectional_separation_without_protection_cannot_pass_review():
    audit = audit_module()
    metrics = pd.DataFrame(
        [{
            "cell": "HepG2", "tf": "FOXA1", "status": "eligible",
            "method": "frozen_parametric_bias_conventional_geometry",
            "auroc_gain_over_dwm": 0.05, "relative_auprc_gain_over_dwm": 0.20,
            "auroc_gain_over_raw": 0.04, "relative_auprc_gain_over_raw": 0.15,
            "functional_separation_relative_change_over_dwm": 1.5,
        }]
    )
    shapes = pd.DataFrame(
        [{"cell": "HepG2", "tf": "FOXA1", "method": metrics.iloc[0]["method"],
          "canonical_protection": False}]
    )
    bootstrap = pd.DataFrame(
        [dict(cell="HepG2", tf="FOXA1", method=metrics.iloc[0]["method"],
              baseline=baseline, auroc_gain_lower_95=0.02,
              relative_auprc_gain_lower_95=0.05)
         for baseline in ["DWM_conventional_geometry", "raw_geometry"]]
    )
    decisions = audit.task_decisions(metrics, shapes, bootstrap)
    assert decisions.iloc[0]["interpretation"] == "no_canonical_protection"
    assert not decisions.iloc[0]["promotion_eligible"]
    assert decisions.iloc[0]["exact_operator_safety"] == "not_established"
    with pytest.raises(ValueError, match="missing.*shape"):
        audit.task_decisions(metrics, shapes.iloc[:0], bootstrap)


def test_missing_aggregate_values_cannot_be_silently_averaged():
    audit = audit_module()
    profiles = profile_rows(center=-0.2, shoulders=0.1)
    profiles.loc[profiles.position.eq(0), "positive_minus_negative"] = np.nan
    with pytest.raises(ValueError, match="finite"):
        audit.shape_audit(profiles)
