import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import benchmark_direct_bias_inference as module


def plan():
    return json.loads((module.REPOSITORY / "benchmarks/manifests/direct_bias_inference_benchmark_v1.json").read_text())


@pytest.mark.parametrize("key,value", [("promotion_evidence", True), ("new_holdout_labels_allowed", True),
                                      ("model_training_allowed", True), ("chromosomes", ["chr16"]),
                                      ("methods", ["python_dwm_reference", module.METHODS[1]])])
def test_scope_cannot_be_promoted_or_swapped_to_a_slow_reference(key, value):
    document = plan()
    module.validate_plan(document)
    document[key] = value
    with pytest.raises(ValueError, match="component diagnostic"):
        module.validate_plan(document)


def test_site_selection_ignores_labels_order_and_numeric_scores():
    frame = pd.DataFrame({"cell": ["HepG2"] * 10, "tf": ["CTCF"] * 10,
                          "TFBS_chr": ["chr19"] * 10, "TFBS_start": np.arange(10) * 10,
                          "TFBS_end": np.arange(10) * 10 + 5, "TFBS_strand": ["+", "-"] * 5,
                          "all_run_dwm_present": [True] * 9 + [False],
                          "chip_label": np.arange(10) % 2, "motif_score": np.arange(10)})
    selected = module.select_sites(frame, 5, 2026)
    changed = frame.sample(frac=1, random_state=7).assign(chip_label=99, motif_score=-1000)
    pd.testing.assert_frame_equal(selected, module.select_sites(changed, 5, 2026))
    assert list(selected.columns) == module.KEYS + ["selection_hash"]
    assert 90 not in selected.TFBS_start.tolist()
    with pytest.raises(ValueError, match="insufficient"):
        module.select_sites(frame, 10, 2026)
    with pytest.raises(ValueError, match="duplicate"):
        module.select_sites(pd.concat([frame, frame.iloc[[0]]]), 5, 2026)


def test_compiled_dwm_uses_native_encoding_reverse_complement_and_exact_positions():
    model = module.AtacBias(L=5, stype="DWM")
    rng = np.random.default_rng(17)
    for strand in ("forward", "reverse"):
        for _ in range(40):
            codes = rng.integers(0, 4, size=5, dtype=np.int64)
            model.bias[strand].add_sequence(codes, float(rng.integers(1, 8)))
            model.bias[strand].add_background(codes, float(rng.integers(1, 8)))
        model.bias[strand].prepare_mat()
    sequence = "ATCGGCTAATCGATGCATCGG"
    width, margin = 7, 7
    observed = module.sequence_scores(module.METHODS[0], model, [sequence], width, margin)[0]
    complement = str.maketrans("ATCGN", "TAGCN")
    expected = np.empty_like(observed)
    for i in range(width):
        context = sequence[margin+i-2:margin+i+3]
        expected[0, i] = model.bias["forward"].score_sequence(module.native_sequences.nuc_to_num(context))[2]
        expected[1, i] = model.bias["reverse"].score_sequence(module.native_sequences.nuc_to_num(context.translate(complement)[::-1]))[2]
    np.testing.assert_allclose(observed, expected, atol=1e-12, rtol=1e-12)
    assert not np.allclose(observed[0], observed[1])


def test_parametric_scoring_preserves_recorded_operator_and_unknown_bases():
    from fp_tools.tools.parametric_bias import BiasFeatureSpec
    model = module.ConditionalSequenceBiasModel(BiasFeatureSpec.loglinear21())
    rng = np.random.default_rng(8)
    model.main[:] = rng.normal(size=model.main.shape)
    sequences = ["".join(rng.choice(list("ACGT"), size=283)), "N" * 283]
    result = module.sequence_scores(module.METHODS[1], model, sequences, 201, 41)
    plus, minus, _ = module.strand_log_bias(model, sequences[0], 41 + np.arange(201))
    np.testing.assert_array_equal(result[0], np.stack([plus, minus]))
    assert np.isnan(result[1]).all()


def rows():
    return [{"cell": "HepG2", "method": method, "repetition": repetition,
             "inference_wall_seconds": 10 if method == module.METHODS[0] else 1,
             "inference_cpu_seconds": 9 if method == module.METHODS[0] else 0.9,
             "worker_peak_rss_mib": 100, "prediction_id": method}
            for repetition in range(6) for method in module.METHODS]


def test_component_ratios_never_certify_whole_pipeline_and_require_complete_repeats():
    data = rows()
    result = module.summarize(data)[0]
    assert result["component_mean_wall_ratio"] == pytest.approx(0.1)
    assert result["component_max_rss_ratio"] == 1
    assert result["promotion_evidence"] is False
    assert result["end_to_end_runtime_measured"] is False
    assert result["end_to_end_memory_measured"] is False
    assert "runtime_seconds" not in result and "peak_memory_mb" not in result
    with pytest.raises(ValueError, match="six paired"):
        module.summarize(data[:-1])
    with pytest.raises(ValueError, match="six paired"):
        module.summarize(data[:-1] + [data[0]])
    data[0]["prediction_id"] = "different"
    with pytest.raises(ValueError, match="predictions changed"):
        module.summarize(data)


def test_worker_resume_binds_repetition_before_opening_predictions(tmp_path):
    path = tmp_path / "worker.json"
    module.write_json(path, {"configuration_id": "config", "cell": "HepG2", "method": module.METHODS[0],
                            "sites": 512, "promotion_evidence": False, "repetition": 0}, "worker_id")
    with pytest.raises(ValueError, match="frozen benchmark"):
        module.verify_worker(path, {"configuration_id": "config", "plan": {"sites_per_cell": 512}},
                             "HepG2", module.METHODS[0], 1)


def test_checksums_reject_changed_nested_output(tmp_path):
    path = tmp_path / "scores.npz"
    path.write_bytes(b"before")
    entry = module.record(path)
    path.write_bytes(b"after")
    with pytest.raises(ValueError, match="checksum mismatch"):
        module.checked(entry, tmp_path)


@pytest.mark.skipif(sys.platform != "linux", reason="declared Linux VmHWM measurement")
def test_fresh_worker_peak_does_not_inherit_a_previous_large_child():
    def peak(megabytes):
        command = ("from pathlib import Path\n"
                   f"allocation = bytearray({megabytes} * 1024 * 1024)\n"
                   "print(next(int(line.split()[1]) for line in Path('/proc/self/status').read_text().splitlines() "
                   "if line.startswith('VmHWM:')))\n")
        return int(subprocess.check_output([sys.executable, "-c", command], text=True))
    high = peak(80)
    low = peak(1)
    assert high - low > 50 * 1024


def test_json_identity_and_array_digest_detect_mutation(tmp_path):
    path = tmp_path / "config.json"
    document = module.write_json(path, {"scope": "component"}, "configuration_id")
    assert module.read_json(path, "configuration_id") == document
    document["scope"] = "promotion"
    path.write_text(json.dumps(document))
    with pytest.raises(ValueError, match="invalid configuration_id"):
        module.read_json(path, "configuration_id")
    a = np.asarray([1., np.nan])
    b = a.copy()
    assert module.array_identity(a) == module.array_identity(b)
    b[0] = 2
    assert module.array_identity(a) != module.array_identity(b)
