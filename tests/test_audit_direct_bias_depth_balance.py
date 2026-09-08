from copy import deepcopy
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import audit_direct_bias_depth_balance as m


def sites(cell="K562"):
    return pd.DataFrame({"cell": [cell] * 8, "tf": ["CTCF"] * 8, "motif_family": ["CTCF"] * 8,
                         "TFBS_chr": ["chr19"] * 4 + ["chr20"] * 4,
                         "TFBS_start": np.arange(8) * 30 + 100, "TFBS_end": np.arange(8) * 30 + 110,
                         "TFBS_strand": ["+", "-"] * 4, "chip_label": [0, 1] * 4,
                         "motif_score": np.repeat(np.arange(4) + 1., 2), "accessibility": np.repeat(np.arange(4) + 1., 2),
                         "site_id": ["repeated_motif_name"] * 8})


def pairs(frame):
    return frame.assign(pair_id=np.repeat(["a", "b", "c", "d"], 2), policy="motif_accessibility_observed_total")


def observed():
    result = np.zeros((8, 201))
    result[:, 100] = [0, 0, 10, 11, 10, 100, 20, 20]
    return result


def scales(sd=1.):
    return {"CTCF": {"log_observed_total": {"sd": sd}}}


def matrix():
    rows = [(cell, f"{cell}_rep1", depth, seed) for cell in ("HepG2", "K562")
            for depth in ("10m", "25m", "50m") for seed in range(2026, 2031)]
    rows += [(cell, f"{cell}_rep{rep}", "full", 2026) for cell in ("HepG2", "K562") for rep in (1, 2, 3)]
    return pd.DataFrame(rows, columns=m.base.RUN_KEYS)


def document(path, payload, field="stage_id"):
    payload = deepcopy(payload)
    payload[field] = m.canonical_id(payload)
    m.write_immutable(path, payload)
    return m.record(path)


@pytest.mark.parametrize("field,value", [("cells", ["GM23338"]), ("caliper_standard_deviations", .3),
                                        ("candidate_retuning_allowed", True), ("rematching_allowed", True),
                                        ("primary_results_recomputed", True), ("new_holdout_labels_allowed", True)])
def test_scope_and_single_caliper_are_fixed(field, value):
    plan = json.loads((m.REPOSITORY / "benchmarks/manifests/direct_bias_depth_balance_v1.json").read_text())
    m.validate_plan(plan)
    plan[field] = value
    with pytest.raises(ValueError, match="fixed post-hoc"):
        m.validate_plan(plan)


def test_symmetric_current_count_caliper_keeps_original_pairs_and_zero_ties():
    frame = sites()
    source_pairs = m.qualify_pairs(frame, pairs(frame))
    selected = m.select_pairs(frame, source_pairs, observed(), scales())
    assert selected.selected.tolist() == [True, True, False, True]
    np.testing.assert_array_equal(m.member_indexes(source_pairs, selected), [0, 1, 2, 3, 6, 7])
    assert selected.loc[selected.pair_id.eq("a"), "positive_log_total"].iloc[0] == 0
    # A candidate's score cannot alter an existing pair or its count mask.
    decorated = frame.assign(candidate_score=np.arange(8) * 1000.)
    decorated_pairs = source_pairs.assign(candidate_score=-np.arange(8) * 1000.)
    pd.testing.assert_frame_equal(selected, m.select_pairs(decorated, decorated_pairs, observed(), scales()))
    flipped = frame.assign(chip_label=1 - frame.chip_label)
    flipped_selection = m.select_pairs(flipped, m.qualify_pairs(flipped, pairs(flipped)), observed(), scales())
    np.testing.assert_array_equal(selected.selected, flipped_selection.selected)
    np.testing.assert_allclose(selected.standardized_gap, flipped_selection.standardized_gap)


def test_expectation_support_removes_whole_pairs_without_repair():
    frame = sites()
    reduced = frame.drop(index=5).reset_index(drop=True)
    retained = m.qualify_pairs(reduced, pairs(frame))
    assert set(retained.pair_id) == {"a", "b", "d"} and len(retained) == 6
    with pytest.raises(ValueError, match="complete pairs"):
        m.qualify_pairs(frame, pairs(frame).iloc[:-1])
    wrong = pairs(frame).copy()
    wrong.loc[1, "TFBS_chr"] = "chr20"
    changed_sites = frame.copy()
    changed_sites.loc[1, "TFBS_chr"] = "chr20"
    with pytest.raises(ValueError, match="one TF and chromosome"):
        m.qualify_pairs(changed_sites, wrong)


def test_zero_training_scale_only_accepts_equal_current_pair_totals():
    frame = sites()
    selected = m.select_pairs(frame, m.qualify_pairs(frame, pairs(frame)), observed(), scales(0.))
    assert selected.selected.tolist() == [True, False, False, True]
    assert np.isinf(selected.standardized_gap.iloc[1])
    for sd in (-1., float("nan")):
        with pytest.raises(ValueError, match="scale is invalid"):
            m.select_pairs(frame, m.qualify_pairs(frame, pairs(frame)), observed(), scales(sd))


def test_bad_counts_and_adaptive_caliper_are_rejected():
    frame = sites()
    source_pairs = m.qualify_pairs(frame, pairs(frame))
    for value in (-1., float("nan")):
        bad = observed()
        bad[0, 0] = value
        with pytest.raises(ValueError, match="fixed caliper"):
            m.select_pairs(frame, source_pairs, bad, scales())
    with pytest.raises(ValueError, match="fixed caliper"):
        m.select_pairs(frame, source_pairs, observed(), scales(), .3)


def test_retention_reports_tfs_without_pairs_and_unassessable_balance():
    frame = pd.concat([sites(), sites().iloc[:2].assign(tf="ZNF384")], ignore_index=True)
    counts = np.concatenate([observed(), np.zeros((2, 201))])
    source_pairs = m.qualify_pairs(frame, pairs(sites()))
    selected = m.select_pairs(frame, source_pairs, counts, scales())
    result = m.retention_rows(frame, source_pairs, selected, counts, 200).set_index("tf")
    assert result.loc["CTCF", "original_pairs"] == 4 and result.loc["CTCF", "retained_pairs"] == 3
    assert result.loc["CTCF", "zero_cut_sites_retained"] == 2
    assert result.loc["ZNF384", "retained_pairs"] == 0
    assert result.loc["ZNF384", "status"] == "underpowered"
    assert not result.loc["ZNF384", "all_measured_covariates_balanced"]
    assert np.isnan(result.loc["ZNF384", "observed_total_smd_retained"])


def test_all_36_selections_are_frozen_before_scoring_and_resume_is_immutable(tmp_path, monkeypatch):
    frame_by_cell = {cell: sites(cell) for cell in ("HepG2", "K562")}
    matched = pd.concat([pairs(frame) for frame in frame_by_cell.values()], ignore_index=True)
    matched_path = tmp_path / "matched.tsv"
    matched.to_csv(matched_path, sep="\t", index=False)
    config = {"configuration_id": "new", "settings": {"cells": ["HepG2", "K562"]}}
    source = {"configuration_id": "old", "settings": {"minimum_sites_per_class": 200}}
    result = {"stages": []}
    for cell, frame in frame_by_cell.items():
        path = tmp_path / (cell + ".npz")
        # Selection must not inspect model expectations, even if unusable.
        np.savez_compressed(path, observed=observed(), expected=np.full((8, 201), np.nan),
                            site_hash=m.base.site_hashes(frame))
        for row in matrix().query("cell == @cell").itertuples(index=False):
            directory = tmp_path / row.sample / row.depth / str(row.seed)
            directory.mkdir(parents=True)
            item = document(directory / "extraction.json", {"configuration_id": "old", "run": row._asdict(),
                                                              "outputs": {"profiles": m.record(path)}})
            result["stages"].append({"extraction": item})
    def load(cell, *args):
        frame = frame_by_cell[cell]
        return frame, m.qualify_pairs(frame, pairs(frame)), np.ones((8, 201))
    monkeypatch.setattr(m, "load_cell", load)
    def forbidden(*args, **kwargs):
        raise AssertionError("a model was scored during cohort selection")
    monkeypatch.setattr(m.base, "score_panel", forbidden)
    selection = {"outputs": {"cohorts": m.record(matched_path)}}
    scale_document = {"scales": {cell: scales() for cell in frame_by_cell}}
    directory = tmp_path / "new"
    directory.mkdir()
    _, done = m.freeze_cohorts(config, source, result, selection, scale_document, directory, tmp_path)
    assert len(done["runs"]) == 36 and not done["new_evaluation_started"]
    retained = pd.read_csv(directory / "retention.tsv.gz", sep="\t")
    assert len(retained) == 36 and retained.retained_pairs.eq(3).all()
    monkeypatch.setattr(m, "select_pairs", forbidden)
    assert m.freeze_cohorts(config, source, result, selection, scale_document, directory, tmp_path)[1] == done
    (directory / "pairs.tsv.gz").write_bytes(b"changed selection")
    with pytest.raises(ValueError, match="checksum"):
        m.freeze_cohorts(config, source, result, selection, scale_document, directory, tmp_path)


def test_zero_pair_evaluation_retains_coverage_without_inventing_metrics(tmp_path):
    result = {"outputs": {}}
    for name in ("metrics", "bootstrap", "profiles", "shapes", "coverage"):
        path = tmp_path / ("source_" + name + ".tsv")
        pd.DataFrame(columns=["cell", "tf", "cohort", "method"]).to_csv(path, sep="\t", index=False)
        result["outputs"][name] = m.record(path)
    frame = sites()
    source_pairs = m.qualify_pairs(frame, pairs(frame))
    selected = m.select_pairs(frame, source_pairs, observed(), scales())
    selected["selected"] = False
    retention = m.retention_rows(frame, source_pairs, selected, observed(), 200)
    run = {"cell": "K562", "sample": "K562_rep1", "depth": "full", "seed": 2026}
    output = tmp_path / "empty"
    output.mkdir()
    done = m.empty_evaluation(run, result, retention, output, "new", tmp_path)
    assert done["status"] == "no_task_with_two_retained_pairs"
    assert pd.read_csv(output / "metrics.tsv.gz", sep="\t").empty
    coverage = pd.read_csv(output / "coverage.tsv.gz", sep="\t")
    assert len(coverage) == 1 and coverage.n_positive.iloc[0] == coverage.n_negative.iloc[0] == 0
    assert coverage.status.iloc[0] == "underpowered"
    m.base.verify_stage(output / "evaluation.json", "new", tmp_path, expected_run=run)


def test_secondary_evaluation_preserves_primary_values_and_uses_only_archived_profiles(tmp_path, monkeypatch):
    frame = sites()
    counts = observed()
    expected = np.repeat(counts.sum(axis=1)[:, None] / 201, 201, axis=1)
    dwm = np.ones((8, 201))
    source_pairs = m.qualify_pairs(frame, pairs(frame))
    selected = m.select_pairs(frame, source_pairs, counts, scales())
    run = {"cell": "K562", "sample": "K562_rep1", "depth": "full", "seed": 2026}
    row = SimpleNamespace(**run)
    settings = {"minimum_sites_per_class": 200, "aggregate_bootstrap_iterations": 2, "seed": 2026,
                "dispersion": {"K562": 1.6144318781609228}}
    primary = tmp_path / "primary"
    primary.mkdir()
    methods = m.base.score_panel(counts, expected, dwm, dwm, settings["dispersion"]["K562"])
    original = m.base.evaluate_run(row, frame, {"motif_accessibility_observed_total": np.arange(8)}, methods,
                                   counts, settings, primary, "old", tmp_path)
    original_records = {name: (entry, Path(entry["path"]).stat().st_mtime_ns) for name, entry in original["outputs"].items()}
    archive_path = primary / "profiles.npz"
    np.savez_compressed(archive_path, observed=counts, expected=expected, sample_dwm=dwm, site_hash=m.base.site_hashes(frame))
    source_extraction = document(primary / "extraction.json", {"configuration_id": "old", "run": run,
                                                               "outputs": {"profiles": m.record(archive_path)}})
    result = {"stages": [{"extraction": source_extraction}], "outputs": original["outputs"]}
    config = {"configuration_id": "new", "settings": {"source_cohort": "motif_accessibility_observed_total", "evidence": "fixture"}}
    output = tmp_path / "secondary"
    output.mkdir()
    records = {}
    for name, table in (("pairs", selected.assign(**run)),
                        ("retention", m.retention_rows(frame, source_pairs, selected, counts, 200).assign(**run))):
        path = output / (name + ".tsv.gz")
        table.to_csv(path, sep="\t", index=False)
        records[name] = m.record(path)
    cohort_record = document(output / "cohorts.freeze.json", {"configuration_id": "new", "outputs": records})
    cohorts = m.read_document(Path(cohort_record["path"]), "stage_id")
    def forbidden(*args, **kwargs):
        raise AssertionError("a BAM or sequence model was re-extracted")
    monkeypatch.setattr(m.base, "extract_strand_cut_profiles", forbidden)
    done = m.evaluate(config, {"settings": settings}, result, {"K562": (frame, source_pairs, dwm)}, cohorts, output, tmp_path)
    assert done["primary_results_recomputed"] is False and done["promotion_eligible"] is False
    combined = pd.read_csv(output / "balance_metrics.tsv.gz", sep="\t")
    preserved = combined.loc[combined.evidence_role.eq("preserved_primary")].drop(columns="evidence_role").reset_index(drop=True)
    pd.testing.assert_frame_equal(preserved, pd.read_csv(primary / "metrics.tsv.gz", sep="\t"), check_exact=False, rtol=1e-14)
    current = combined.loc[combined.cohort.eq(m.SECONDARY)]
    assert len(current) == 4 and current.n_sites.eq(6).all() and current.status.eq("underpowered").all()
    for entry, mtime in original_records.values():
        assert m.record(Path(entry["path"])) == entry
        assert Path(entry["path"]).stat().st_mtime_ns == mtime
    monkeypatch.setattr(m.base, "score_panel", forbidden)
    assert m.evaluate(config, {"settings": settings}, result, {"K562": (frame, source_pairs, dwm)}, cohorts, output, tmp_path) == done
