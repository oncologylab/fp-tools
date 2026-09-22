from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import probe_process_tree_resources as probe


def measurements():
    base = {"measurement_complete": True, "descendants_remaining": False, "cpu_seconds": .5}
    return {"serial": {**base, "peak_memory_bytes": 42_000_000},
            "parallel": {**base, "peak_memory_bytes": 80_000_000}}


def test_kernel_counters_require_integer_values(tmp_path):
    path = tmp_path / "counter"
    path.write_text("1024\n")
    assert probe.counter(path) == 1024
    for value in ("max", "-1", "nan", ""):
        path.write_text(value)
        with pytest.raises(ValueError):
            probe.counter(path)


def test_simultaneous_memory_and_descendant_cpu_qualified():
    probe.validate_probes(measurements())


@pytest.mark.parametrize("mode,field,value", [
    ("serial", "cpu_seconds", .1),
    ("parallel", "measurement_complete", False),
    ("parallel", "descendants_remaining", True),
    ("serial", "peak_memory_bytes", 80_000_000),
    ("parallel", "peak_memory_bytes", 43_000_000),
])
def test_inadequate_measurement_rejected(mode, field, value):
    records = measurements()
    records[mode][field] = value
    with pytest.raises(ValueError):
        probe.validate_probes(records)


def test_foreign_cgroup_rejected_before_writes(tmp_path, monkeypatch):
    monkeypatch.setattr(probe, "current_cgroup", lambda: tmp_path / "actual")
    with pytest.raises(RuntimeError, match="refusing"):
        probe.measure_tree([sys.executable], tmp_path, {}, tmp_path / "foreign", tmp_path)
    assert not (tmp_path / "foreign").exists()


def test_probe_commands_are_direct_argument_arrays_and_failure_is_real():
    for mode in ("serial", "parallel", "failure", "orphan"):
        args = probe.probe_command(mode)
        assert args[:2] == [sys.executable, "-c"]
    assert subprocess.run(probe.probe_command("failure")).returncode == 7
    with pytest.raises(ValueError):
        probe.probe_command("scientific_workload")
