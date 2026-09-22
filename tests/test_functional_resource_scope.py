from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
from evaluate_functional_promotion import resource_evidence


GATES = {"maximum_runtime_ratio_to_dwm": 2.,
         "maximum_memory_ratio_to_dwm": 1.5, "maximum_model_size_mb": 25.}


def evidence():
    return pd.DataFrame([
        {"method": method, "runtime_seconds": wall, "cpu_seconds": cpu,
         "peak_memory_mb": memory, "model_size_mb": 2.,
         "benchmark_id": "a" * 64, "workload_id": "b" * 64,
         "repetition": repetition, "cores": 1, "cache_policy": "fresh_process_warm_os_cache",
         "pipeline_mode": "fit_and_apply", "resource_scope": "complete_pipeline",
         "cpu_scope": "process_tree", "memory_scope": "absolute_process_tree_peak"}
        for repetition in range(6)
        for method, wall, cpu, memory in (("DWM", 100., 80., 100.), ("candidate", 90., 100., 120.))])


def test_legacy_component_table_cannot_pass_end_to_end_gate():
    frame = evidence()[["method", "runtime_seconds", "peak_memory_mb", "model_size_mb"]]
    assert not resource_evidence(frame, "candidate", "DWM", GATES)[0]


@pytest.mark.parametrize("field,value", [
    ("resource_scope", "sequence_inference_component_only"),
    ("memory_scope", "cumulative_children_high_water_mark"),
    ("cpu_scope", "parent_only"), ("pipeline_mode", "apply_only"),
    ("workload_id", "c" * 64), ("cores", 4), ("cpu_seconds", -1),
    ("peak_memory_mb", float("nan")), ("model_size_mb", -1),
])
def test_incomplete_or_unpaired_measurement_fails(field, value):
    frame = evidence()
    frame.loc[1, field] = value
    assert not resource_evidence(frame, "candidate", "DWM", GATES)[0]


def test_cpu_budget_cannot_be_replaced_by_wall_time():
    frame = evidence()
    frame.loc[frame.method.eq("candidate"), "cpu_seconds"] = 200.
    assert not resource_evidence(frame, "candidate", "DWM", GATES)[0]


def test_valid_paired_cpu_evidence_and_duplicates():
    frame = evidence()
    passed, result = resource_evidence(frame, "candidate", "DWM", GATES)
    assert passed and result["runtime_ratio"] == 1.25
    assert not resource_evidence(pd.concat([frame, frame.iloc[[1]]]), "candidate", "DWM", GATES)[0]


def test_fast_case_cannot_hide_slow_case():
    first = evidence()
    second = evidence()
    first.loc[first.method.eq("DWM"), "cpu_seconds"] = 10000
    second["workload_id"] = "c" * 64
    second.loc[second.method.eq("candidate"), "cpu_seconds"] = 200
    assert not resource_evidence(pd.concat([first, second]), "candidate", "DWM", GATES)[0]


def test_missing_repetition_in_both_methods_fails():
    frame = evidence().query("repetition != 5")
    assert not resource_evidence(frame, "candidate", "DWM", GATES)[0]
