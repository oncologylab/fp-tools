from copy import deepcopy
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import extend_direct_bias_stability as m


def matrix():
    rows = [(cell, f"{cell}_rep{rep}", depth, seed) for cell in ("HepG2", "K562") for rep in (1, 2, 3)
            for depth in ("10m", "25m", "50m") for seed in range(2026, 2031)]
    rows += [(cell, f"{cell}_rep{rep}", "full", 2026) for cell in ("HepG2", "K562") for rep in (1, 2, 3)]
    return pd.DataFrame(rows, columns=m.base.RUN_KEYS).sort_values(m.base.RUN_KEYS).reset_index(drop=True)


def split_matrix():
    frame = matrix()
    old = frame[frame["sample"].str.endswith("rep1") | frame.depth.eq("full")]
    return old.copy(), frame.drop(old.index).copy()


def sites(cell="K562"):
    return pd.DataFrame({"cell": [cell] * 4, "tf": ["CTCF"] * 4, "motif_family": ["CTCF"] * 4,
                         "TFBS_chr": ["chr19", "chr19", "chr20", "chr20"], "TFBS_start": [10, 30, 50, 70],
                         "TFBS_end": [20, 40, 60, 80], "TFBS_strand": ["+", "-", "+", "-"],
                         "chip_label": [0, 1, 0, 1], "motif_score": [1., 1., 2., 2.],
                         "accessibility": [1., 1., 2., 2.], "site_id": ["repeated_motif_name"] * 4})


def document(path, payload, field="stage_id"):
    payload = deepcopy(payload)
    payload[field] = m.canonical_id(payload)
    m.write_immutable(path, payload)
    return m.record(path)


def test_all_replicates_depths_and_seeds_without_full_depth_pseudoreplicates():
    old, new = split_matrix()
    actual = m.validate_merge(old, new)
    pd.testing.assert_frame_equal(actual, matrix())
    for bad in (actual.iloc[:-1], pd.concat([actual, actual.iloc[[0]]]),
                pd.concat([actual, actual.query("depth == 'full'").iloc[[0]].assign(seed=2027)])):
        with pytest.raises(ValueError, match="96 development"):
            m.validate_matrix(bad)
    with pytest.raises(ValueError, match="missing 60"):
        m.validate_merge(old, new.iloc[:-1])
    with pytest.raises(ValueError, match="96 development"):
        m.validate_merge(old, pd.concat([new.iloc[:-1], old.iloc[[0]]]))


@pytest.mark.parametrize("field,value", [("cells", ["SK-N-SH"]), ("candidate_retuning_allowed", True),
                                        ("new_holdout_labels_allowed", True), ("final_untouched_control_allowed", True),
                                        ("evaluation_chromosomes", ["chr16"]), ("replicates", [1])])
def test_scope_cannot_change(field, value):
    plan = json.loads((m.REPOSITORY / "benchmarks/manifests/direct_bias_stability_96_v1.json").read_text())
    m.validate_plan(plan)
    plan[field] = value
    with pytest.raises(ValueError, match="fixed 96-run"):
        m.validate_plan(plan)


def test_missing_prerequisite_stops_before_reading_archived_labels(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("artifact access before the prerequisite completed")
    monkeypatch.setattr(m, "checked_path", forbidden)
    with pytest.raises(ValueError, match="preparation is incomplete"):
        m.freeze_inputs({"prepared_inputs": {"path": str(tmp_path / "missing.json")}}, tmp_path)


def test_frozen_plan_pins_the_current_evaluator():
    path = m.REPOSITORY / "benchmarks/manifests/direct_bias_stability_96_v1.json"
    plan = json.loads(path.read_text())
    assert m.checked_path(plan["evaluator_code"], m.REPOSITORY).resolve() == Path(m.__file__).resolve()
    queue = json.loads((m.REPOSITORY / "benchmarks/manifests/direct_bias_stability_96_queue_v1.json").read_text())
    command = queue["stages"][0]["command"]
    assert command[command.index("--plan-sha256") + 1] == m.file_sha256(path)
    assert queue["prerequisite"]["manifest"] == plan["prepared_inputs"]["path"]


def test_changed_queue_plan_stops_before_input_access(tmp_path):
    path = tmp_path / "changed.json"
    path.write_text("{}")
    with pytest.raises(ValueError, match="queued plan checksum changed"):
        m.main(["--plan", str(path), "--plan-sha256", "0" * 64,
                "--data-root", str(tmp_path), "--outdir", str(tmp_path / "output")])
    assert not (tmp_path / "output").exists()


def test_output_lock_prevents_two_evaluators(tmp_path):
    with m.exclusive_run(tmp_path):
        with pytest.raises(ValueError, match="another 96-run evaluator"):
            with m.exclusive_run(tmp_path):
                raise AssertionError("second evaluator acquired the output directory")
    with m.exclusive_run(tmp_path):
        pass


def test_prepared_table_is_bound_to_all_60_correction_inputs():
    _, new = split_matrix()
    new = new.assign(bam_sha256="bam", dwm_expected_sha256="dwm")
    stages = [{"sample": row.sample, "depth": int(row.depth[:-1]) * 1_000_000, "seed": row.seed,
               "input_bam": {"sha256": "bam"}, "outputs": {"sample_expected.bw": {"sha256": "dwm"}}}
              for row in new.itertuples(index=False)]
    m.bind_prepared_rows(new, stages)
    for field in ("bam_sha256", "dwm_expected_sha256"):
        changed = new.copy()
        changed.loc[changed.index[0], field] = "substituted"
        with pytest.raises(ValueError, match="substituted"):
            m.bind_prepared_rows(changed, stages)
    with pytest.raises(ValueError, match="duplicate"):
        m.bind_prepared_rows(new, stages[:-1] + stages[:1])
    changed = deepcopy(stages)
    changed[0]["depth"] += 1
    with pytest.raises(ValueError, match="undeclared depth"):
        m.bind_prepared_rows(new, changed)


def test_correction_command_must_preserve_the_frozen_dwm_template(tmp_path):
    template = ["python", "-m", "fp_tools.cli", "--bams", "{bam}", "--prefix", "{prefix}",
                "--outdir", "{outdir}", "--score_mat", "DWM", "--cores", "4"]
    stage = {"sample": "HepG2_rep2", "depth": 10_000_000, "seed": 2026,
             "input_bam": {"path": str(tmp_path / "sample.bam")}}
    replacements = {"{bam}": stage["input_bam"]["path"], "{prefix}": "HepG2_rep2.10m.s2026", "{outdir}": str(tmp_path)}
    stage["command"] = [replacements.get(item, item) for item in template]
    m.bind_correction_command(stage, {"correction_template": template}, tmp_path / "correction.json")
    stage["command"] += ["--bias-pkl", "full_depth_model.pickle"]
    with pytest.raises(ValueError, match="changed the frozen DWM command"):
        m.bind_correction_command(stage, {"correction_template": template}, tmp_path / "correction.json")


@pytest.mark.parametrize("field,value", [("source_bam_sha256", "wrong"), ("output_bai_sha256", "wrong"),
                                        ("available_fragments", 60_000_001), ("target_fragments", 25_000_000),
                                        ("seed", 2027), ("nested_with_same_seed", False)])
def test_nested_metadata_cannot_substitute_source_count_or_sampling(tmp_path, field, value):
    row = SimpleNamespace(depth="10m", seed=2026, bam_sha256="bam", bai_sha256="bai")
    metadata = {"schema": "fp-tools-nested-depth-matrix-v1", "source_bam_sha256": "source",
                "output_bam_sha256": "bam", "output_bai_sha256": "bai", "available_fragments": 60_000_000,
                "target_fragments": 10_000_000, "seed": 2026, "sampling_fraction": 1 / 6,
                "selection_rule": "blake2b(seed,query_name)<sampling_fraction", "nested_with_same_seed": True}
    sampling = {"outputs": {"10000000:2026:bam": {"sha256": "bam"}, "10000000:2026:bai": {"sha256": "bai"},
                            "10000000:2026:metadata": document(tmp_path / "good.json", metadata, "metadata_id")}}
    check = lambda entry: m.checked_path(entry, tmp_path)
    args = (sampling, {"available_fragments": 60_000_000}, {"bam": {"sha256": "source"}}, check)
    m.bind_sampling_row(row, *args)
    metadata[field] = value
    sampling["outputs"]["10000000:2026:metadata"] = document(tmp_path / "changed.json", metadata, "metadata_id")
    with pytest.raises(ValueError, match="metadata differs"):
        m.bind_sampling_row(row, *args)


def archive(tmp_path, frame):
    observed = np.arange(len(frame) * 201).reshape(len(frame), 201).astype(float)
    observed[0] = 0
    expected, dwm = observed + 1, observed + 2
    path = tmp_path / "old.npz"
    np.savez_compressed(path, observed=observed, expected=expected, sample_dwm=dwm, site_hash=m.base.site_hashes(frame))
    return {"extraction": document(tmp_path / "old.json", {"configuration_id": "old",
                                                          "outputs": {"profiles": m.record(path)}})}, observed, expected, dwm


def test_archived_profiles_reuse_coordinate_order_and_keep_zero_cut_sites(tmp_path):
    frame = sites()
    archived, observed, expected, dwm = archive(tmp_path, frame)
    selected = frame.iloc[[3, 0]].copy()
    actual = m.subset_archived_profiles(archived, frame, selected, dwm[[3, 0]], tmp_path)
    for values, source in zip(actual, (observed, expected, dwm)):
        np.testing.assert_array_equal(values, source[[3, 0]])
    assert actual[0][1].sum() == 0
    with pytest.raises(ValueError, match="current track"):
        m.subset_archived_profiles(archived, frame, selected, dwm[[3, 0]] + 1, tmp_path)
    selected.loc[selected.index[0], "chip_label"] = 0
    with pytest.raises(ValueError, match="loci or labels"):
        m.subset_archived_profiles(archived, frame, selected, dwm[[3, 0]], tmp_path)
    with pytest.raises(ValueError, match="old site order"):
        m.subset_archived_profiles(archived, frame.iloc[::-1], frame, dwm, tmp_path)
    (tmp_path / "old.npz").write_bytes(b"changed leaf")
    with pytest.raises(ValueError, match="checksum"):
        m.subset_archived_profiles(archived, frame, frame, dwm, tmp_path)


def test_metric_reuse_requires_identical_sites_and_both_cohort_memberships(tmp_path):
    frame = sites()
    cohorts = {"original_support_with_expected_coverage": np.arange(4),
               "motif_accessibility_observed_total": np.arange(4)}
    identity = m.cohort_identity(frame, cohorts)
    assert identity != m.cohort_identity(frame.iloc[::-1], cohorts)
    assert identity != m.cohort_identity(frame, {**cohorts, "motif_accessibility_observed_total": np.arange(2)})
    row = SimpleNamespace(cell="K562", sample="K562_rep1", depth="10m", seed=2026)
    run = vars(row)
    data = tmp_path / "metrics.tsv"
    data.write_text("auroc\n0.5\n")
    archived = {"evaluation": document(tmp_path / "old.json", {"configuration_id": "old", "run": run,
                                                               "outputs": {"metrics": m.record(data)}, "bootstrap_run": True})}
    membership = {"archived_metrics_reusable": True, "archived_identity": identity, "site_and_cohort_identity": identity}
    result = m.reuse_evaluation(row, archived, membership, "old", tmp_path, "new", tmp_path)
    assert result["outputs"]["metrics"] == m.record(data)
    before = data.stat().st_mtime_ns
    assert m.reuse_evaluation(row, archived, membership, "old", tmp_path, "new", tmp_path) == result
    assert data.stat().st_mtime_ns == before
    with pytest.raises(ValueError, match="site or cohort change"):
        m.reuse_evaluation(row, archived, {**membership, "site_and_cohort_identity": "changed"}, "old", tmp_path, "new", tmp_path)
    data.write_text("changed\n")
    with pytest.raises(ValueError, match="checksum"):
        m.reuse_evaluation(row, archived, membership, "old", tmp_path, "new", tmp_path)


def test_joint_freeze_retains_pairs_in_both_cells_before_scoring(tmp_path, monkeypatch):
    old, _ = split_matrix()
    full = matrix().assign(dwm_expected="sample.bw")
    matched = pd.concat([sites(cell).assign(pair_id=["a", "a", "b", "b"], policy="motif_accessibility_observed_total")
                         for cell in ("HepG2", "K562")], ignore_index=True)
    values = np.ones((4, 201))
    monkeypatch.setattr(m.base, "extract_prediction_profiles", lambda *a, **kw: (values.copy(), np.ones(4, dtype=bool)))
    supports = {}
    for cell in ("HepG2", "K562"):
        directory = tmp_path / ("old_" + cell)
        directory.mkdir()
        m.base.freeze_support(cell, old.query("cell == @cell").assign(dwm_expected="sample.bw"), sites(cell),
                              {"counts": values}, values, matched, directory, "old", tmp_path)
        supports[cell] = m.record(directory / "support.json")
    visited = []
    def load_reference(cell, *args):
        visited.append(cell)
        return sites(cell), {"counts": values.copy()}, values.copy(), {}
    monkeypatch.setattr(m.base, "load_reference", load_reference)
    missing = values.copy()
    missing[1] = 0
    monkeypatch.setattr(m.base, "extract_prediction_profiles", lambda *a, **kw: (missing.copy(), np.ones(4, dtype=bool)))
    directory = tmp_path / "new"
    directory.mkdir()
    contexts, joint = m.freeze_supports(full, matched, {"configuration_id": "old", "references": {"HepG2": {}, "K562": {}}},
                                        {"supports": supports}, directory, "new", tmp_path)
    assert visited == ["HepG2", "K562"]
    assert set(joint["supports"]) == set(contexts) == {"HepG2", "K562"}
    assert (directory / "cohorts.freeze.json").is_file()
    for cell, context in contexts.items():
        assert len(context[1]) == 3
        np.testing.assert_array_equal(context[4]["motif_accessibility_observed_total"], [1, 2])
        assert joint["memberships"][cell]["archived_metrics_reusable"] is False


def test_completed_stage_matrix_rejects_missing_or_duplicate_runs(tmp_path, monkeypatch):
    full = matrix()
    old, _ = split_matrix()
    archives = [{"run": row._asdict(), "extraction": {"path": f"archive_{index}"}}
                for index, row in enumerate(old.itertuples(index=False, name="Run"))]
    archived_by_run = {tuple(item["run"].values()): item for item in archives}
    stages, documents = [], {}
    for index, row in enumerate(full.itertuples(index=False, name="Run")):
        run = row._asdict()
        extraction = {"run": run}
        if tuple(run.values()) in archived_by_run:
            extraction["reused_from"] = archived_by_run[tuple(run.values())]["extraction"]
        documents[f"extract_{index}"] = extraction
        documents[f"evaluate_{index}"] = {"run": run}
        stages.append({"extraction": {"path": f"extract_{index}"}, "evaluation": {"path": f"evaluate_{index}"}})
    monkeypatch.setattr(m, "checked_path", lambda entry, root: Path(entry["path"]))
    monkeypatch.setattr(m.base, "verify_stage", lambda path, *args, **kw: documents[str(path)])
    config = {"configuration_id": "new", "archives": archives}
    assert m.verify_completed_stages(stages, config, {}, tmp_path) == 0
    with pytest.raises(ValueError, match="96-run matrix"):
        m.verify_completed_stages(stages[:-1], config, {}, tmp_path)
    with pytest.raises(ValueError, match="96 development"):
        m.verify_completed_stages(stages[:-1] + stages[:1], config, {}, tmp_path)


def test_complete_small_96_run_execution_reuses_archives_and_rechecks_changed_support(tmp_path, monkeypatch):
    for name, value in m.THREADS.items():
        monkeypatch.setenv(name, value)
    old, new = split_matrix()
    full = m.validate_merge(old, new).assign(bam="fixture.bam", bai="fixture.bai", dwm_expected="fixture.bw")
    full["bam_sha256"] = np.where(full["sample"].str.endswith("rep1") | full.depth.eq("full"), "reference", "new")
    old = full[full.bam_sha256.eq("reference")]
    frames = {cell: pd.concat([sites(cell), sites(cell).assign(TFBS_start=lambda f: f.TFBS_start + 1000,
                                                              TFBS_end=lambda f: f.TFBS_end + 1000)], ignore_index=True)
              for cell in ("HepG2", "K562")}
    matched = pd.concat([frame.assign(pair_id=np.repeat(np.arange(4), 2), policy="motif_accessibility_observed_total")
                         for frame in frames.values()], ignore_index=True)
    values = np.ones((8, 201))
    plus = values.copy()
    plus[0] = 0
    minus = plus.copy()
    arrays = {"plus_observed": plus, "minus_observed": minus,
              "plus_log_bias": np.zeros_like(plus), "minus_log_bias": np.zeros_like(minus)}
    reference = {"profile_manifest": {"metadata": {"bam_sha256": "reference"}}}
    settings = {"minimum_sites_per_class": 200, "aggregate_bootstrap_iterations": 2, "seed": 2026,
                "dispersion": {"HepG2": 0., "K562": 1.6144318781609228}}
    source_config = {"settings": settings, "runtime": m.runtime(), "references": {cell: reference for cell in frames}, "code": []}
    source_config["configuration_id"] = m.canonical_id(source_config)
    source_config_record = document(tmp_path / "source_configuration.json", {k: v for k, v in source_config.items()
                                                                            if k != "configuration_id"}, "configuration_id")
    monkeypatch.setattr(m.base, "extract_prediction_profiles", lambda *a, **kw: (values.copy(), np.ones(8, dtype=bool)))
    supports, archives = {}, []
    for cell, frame in frames.items():
        directory = tmp_path / ("old_" + cell)
        directory.mkdir()
        _, _, _, cohorts = m.base.freeze_support(cell, old.query("cell == @cell"), frame, arrays, values, matched,
                                                 directory, source_config["configuration_id"], tmp_path)
        supports[cell] = m.record(directory / "support.json")
        for row in old.query("cell == @cell").itertuples(index=False):
            run_dir = directory / row.sample / row.depth / str(row.seed)
            run_dir.mkdir(parents=True)
            observed, expected, dwm = m.base.extract_run(row, frame, arrays, reference, values, run_dir,
                                                        source_config["configuration_id"], tmp_path)
            methods = m.base.score_panel(observed, expected, values, dwm, settings["dispersion"][cell])
            m.base.evaluate_run(row, frame, cohorts, methods, observed, settings, run_dir, source_config["configuration_id"], tmp_path)
            archives.append({"run": {key: getattr(row, key) for key in m.base.RUN_KEYS},
                             **{name: m.record(run_dir / (name + ".json")) for name in ("extraction", "evaluation")}})
    source_result = {"supports": supports}
    preparation_record = document(tmp_path / "preparation.json", {"code": {}}, "configuration_id")
    plan = json.loads((m.REPOSITORY / "benchmarks/manifests/direct_bias_stability_96_v1.json").read_text())
    plan.update(source_configuration=source_config_record, preparation_configuration=preparation_record,
                evaluator_code=m.record(Path(m.__file__)))
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan))
    monkeypatch.setattr(m, "freeze_inputs", lambda *args: (full, matched, source_config, source_result, archives,
                                                           [source_config_record, preparation_record]))
    monkeypatch.setattr(m.base, "load_reference", lambda cell, *args: (frames[cell], arrays, values.copy(), {}))
    def prediction(frame, *args, **kwargs):
        result = values.copy()
        if frame.cell.iloc[0] == "HepG2":
            result[1] = 0
        return result, np.ones(8, dtype=bool)
    monkeypatch.setattr(m.base, "extract_prediction_profiles", prediction)
    output = tmp_path / "extension"
    output.mkdir()
    extraction_calls, evaluation_calls = [], []
    def extract(frame, *args, **kwargs):
        joint = m.read_document(output / "cohorts.freeze.json", "stage_id")
        assert set(joint["supports"]) == {"HepG2", "K562"}
        extraction_calls.append(len(frame))
        order = m.base.align_sites(frames[frame.cell.iloc[0]], frame)
        return plus[order], minus[order], np.ones(len(frame), dtype=bool)
    monkeypatch.setattr(m.base, "extract_strand_cut_profiles", extract)
    evaluate = m.base.evaluate_run
    def evaluate_checked(row, *args, **kwargs):
        assert (output / "cohorts.freeze.json").is_file()
        evaluation_calls.append(m.run_key(row))
        return evaluate(row, *args, **kwargs)
    monkeypatch.setattr(m.base, "evaluate_run", evaluate_checked)
    result = m.execute(plan, plan_path, tmp_path, output)
    assert result["archived_extractions_reused"] == 36
    assert result["archived_evaluations_reused"] == 18
    assert len(extraction_calls) == 60 and len(evaluation_calls) == 78
    assert result["promotion_eligible"] is False and result["new_holdout_labels_read"] is False
    assert m.verify_result(output, tmp_path) == result
    before = len(extraction_calls), len(evaluation_calls)
    assert m.main(["--plan", str(plan_path), "--data-root", str(tmp_path), "--outdir", str(output)]) == 0
    assert before == (len(extraction_calls), len(evaluation_calls))
    coverage = pd.read_csv(output / "stability_coverage.tsv.gz", sep="\t")
    assert len(coverage) == 192 and coverage.zero_cut_sites.max() == 1
    # A completed manifest cannot hide a changed nested output.
    leaf = Path(json.loads(Path(result["stages"][0]["extraction"]["path"]).read_text())["outputs"]["profiles"]["path"])
    leaf.write_bytes(b"changed nested extraction")
    with pytest.raises(ValueError, match="checksum"):
        m.verify_result(output, tmp_path)
