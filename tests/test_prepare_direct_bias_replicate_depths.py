from pathlib import Path
import json
import shlex
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import prepare_direct_bias_replicate_depths as module


def plan():
    return json.loads((module.REPOSITORY / "benchmarks/manifests/direct_bias_replicate_depths_v1.json").read_text())


@pytest.mark.parametrize("key,value", [("samples", ["SK-N-SH_rep1"]), ("candidate_retuning_allowed", True),
                                      ("new_holdout_labels_allowed", True), ("final_untouched_control_allowed", True),
                                      ("depths", [10_000_000]), ("dwm_comparator", "reuse_full_model")])
def test_plan_rejects_changed_scope_sampling_or_comparator(key, value):
    document = plan()
    module.validate_plan(document)
    document[key] = value
    with pytest.raises(ValueError, match="fixed missing development"):
        module.validate_plan(document)


def test_source_registry_requires_exact_full_replicate_identity():
    frame = pd.DataFrame([{"cell": sample.split("_")[0], "sample": sample, "depth": "full", "seed": 2026}
                          for sample in module.SAMPLES])
    assert len(module.source_rows(frame)) == 4
    for bad in (frame.iloc[:-1], pd.concat([frame, frame.iloc[[0]]]), frame.assign(cell="wrong"),
                frame.assign(depth="50m"), frame.assign(seed=2027)):
        with pytest.raises(ValueError, match="four full development"):
            module.source_rows(bad)


def test_native_count_matches_original_rule_for_every_sam_flag_and_keeps_mapq_zero(tmp_path):
    pysam = pytest.importorskip("pysam")
    path = tmp_path / "flags.bam"
    expected = 0
    with pysam.AlignmentFile(path, "wb", header={"SQ": [{"SN": "chr1", "LN": 10000}]}) as bam:
        for flag in range(4096):
            read = pysam.AlignedSegment()
            read.query_name = f"flag_{flag}"
            read.query_sequence = "A" * 20
            read.flag = flag
            read.reference_id = 0
            read.reference_start = flag
            read.mapping_quality = 0
            read.cigarstring = "20M"
            expected += int(module.fragment_representative(read))
            bam.write(read)
    assert expected == 96
    assert module.native_fragment_count(path) == expected


def test_correction_template_preserves_original_options_and_never_reuses_full_bias(tmp_path):
    common = {key: {"path": str((tmp_path / key).resolve())} for key in ("genome", "peaks", "blacklist")}
    rows = []
    for cell in ("HepG2", "K562"):
        args = [sys.executable, "-m", "fp_tools.cli", "--bams", f"{cell}.bam",
                "--genome", common["genome"]["path"], "--peaks", common["peaks"]["path"],
                "--blacklist", common["blacklist"]["path"], "--outdir", f"{cell}/dwm", "--prefix", cell,
                "--score_mat", "DWM", "--write-tracks", "all", "--cores", "4"]
        rows.append({"job_id": f"correct:{cell}_rep1.10m.s2026:fp_tools_dwm", "command": shlex.join(args)})
    frame = pd.DataFrame(rows)
    result = module.correction_template(frame, tmp_path, plan(), common)
    assert result[result.index("--bams")+1] == "{bam}"
    assert result[result.index("--write-tracks")+1] == "all"
    changed = frame.copy()
    changed.loc[0, "command"] += " --bias-pkl full.pickle"
    with pytest.raises(ValueError, match="refitted conventional"):
        module.correction_template(changed, tmp_path, plan(), common)
    changed = frame.copy()
    changed.loc[0, "command"] = changed.loc[0, "command"].replace("--cores 4", "--cores 8")
    with pytest.raises(ValueError, match="CPU settings differ"):
        module.correction_template(changed, tmp_path, plan(), common)


def test_sampling_resume_verifies_count_manifest_transitively(tmp_path):
    config = {"configuration_id": "frozen"}
    count_path = tmp_path / "count.json"
    document = {"configuration_id": "frozen", "sample": "HepG2_rep2", "stage": "fragment_count", "available_fragments": 100}
    module.write_json(count_path, document, "stage_id")
    sampling = tmp_path / "sampling.json"
    module.write_json(sampling, {"configuration_id": "frozen", "sample": "HepG2_rep2", "stage": "sampling",
                                 "outputs": {}, "count_manifest": module.record(count_path)}, "stage_id")
    module.verify_stage(sampling, config, "HepG2_rep2", "sampling")
    changed = json.loads(count_path.read_text())
    changed["available_fragments"] = 200
    changed["stage_id"] = module.identity(changed, "stage_id")
    count_path.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="checksum"):
        module.verify_stage(sampling, config, "HepG2_rep2", "sampling")


def test_correction_resume_rejects_another_depth_or_seed(tmp_path):
    path = tmp_path / "correction.json"
    config = {"configuration_id": "frozen"}
    module.write_json(path, {"configuration_id": "frozen", "sample": "K562_rep3", "stage": "correction",
                            "depth": 10_000_000, "seed": 2026, "outputs": {}}, "stage_id")
    with pytest.raises(ValueError, match="different frozen request"):
        module.verify_stage(path, config, "K562_rep3", "correction", depth=25_000_000, seed=2026)
    with pytest.raises(ValueError, match="different frozen request"):
        module.verify_stage(path, config, "K562_rep3", "correction", depth=10_000_000, seed=2027)


@pytest.mark.skipif(sys.platform != "linux", reason="preparation requires a Linux file lock")
def test_lock_uses_kernel_state_and_releases_after_failure(tmp_path):
    with module.exclusive_run(tmp_path):
        with pytest.raises(RuntimeError, match="another input preparation"):
            with module.exclusive_run(tmp_path):
                pass
    with pytest.raises(ValueError):
        with module.exclusive_run(tmp_path):
            raise ValueError("interrupted stage")
    with module.exclusive_run(tmp_path):
        pass
