from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "scripts" / "run_chrombpnet_reference.py"
spec = importlib.util.spec_from_file_location("run_chrombpnet_reference", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_container_paths_are_confined_to_repository(tmp_path) -> None:
    path = ROOT / "test_data/genome.fa.gz"
    assert module.container_path(path).startswith("/work/")
    try:
        module.container_path(tmp_path / "outside")
    except ValueError as error:
        assert "must stay" in str(error)
    else:
        raise AssertionError("outside path was accepted")


def test_smoke_command_uses_pinned_source_and_manual_gpu_mounts(monkeypatch) -> None:
    monkeypatch.setattr(module, "driver_library", lambda name: Path("/driver") / name)
    command = module.wrapped_container_command(module.smoke_arguments())
    joined = " ".join(str(value) for value in command)
    assert module.PINNED_IMAGE in command
    assert "--device" in command
    assert "/dev/nvidia0" in command
    assert "host_libnvidia-ptxjitcompiler.so" in joined
    assert "host_libnvidia-nvvm.so" in joined
    assert (
        "PYTHONPATH=/work/benchmarks/results/footprint_external_references/chrombpnet"
        in command
    )
    assert "tensorflow" in joined
    assert "tf.random.uniform" in joined


def test_regulatory_stage_keeps_models_external(tmp_path) -> None:
    parser = module.build_parser()
    args = parser.parse_args(
        [
            "regulatory",
            "--genome",
            "test_data/genome.fa.gz",
            "--chrom-sizes",
            "test_data/chrom_sizes.txt",
            "--peaks",
            "test_data/merged_peaks.bed",
            "--fold",
            "benchmarks/manifests/frozen_parametric_factorization_v1.spec.json",
            "--bam",
            "test_data/Bcell.bam",
            "--nonpeaks",
            "test_data/blacklist.bed",
            "--bias-model",
            "test_data/fake_bias.h5",
            "--output-dir",
            "benchmarks/results/footprint_external_references/fixture",
        ]
    )
    arguments = module.stage_arguments(args)
    assert arguments[:2] == ["chrombpnet", "pipeline"]
    assert "-b" in arguments
    assert all(
        str(value).startswith("/work/")
        for value in arguments
        if str(value).startswith("/work")
    )


def test_core_only_stages_skip_optional_interpretation() -> None:
    parser = module.build_parser()
    common = [
        "--genome",
        "test_data/genome.fa.gz",
        "--chrom-sizes",
        "test_data/chrom_sizes.txt",
        "--peaks",
        "test_data/merged_peaks.bed",
        "--fold",
        "benchmarks/manifests/frozen_parametric_factorization_v1.spec.json",
        "--bam",
        "test_data/Bcell.bam",
        "--nonpeaks",
        "test_data/blacklist.bed",
        "--output-dir",
        "benchmarks/results/footprint_external_references/fixture",
        "--core-only",
    ]
    bias = parser.parse_args(["bias", *common])
    regulatory = parser.parse_args(
        ["regulatory", *common, "--bias-model", "test_data/fake_bias.h5"]
    )
    assert module.stage_arguments(bias)[:3] == ["chrombpnet", "bias", "train"]
    assert module.stage_arguments(regulatory)[:2] == ["chrombpnet", "train"]


def test_prep_output_discovery_ignores_auxiliary_directories(
    tmp_path, monkeypatch
) -> None:
    prefix = tmp_path / "K562"
    expected = tmp_path / "K562_filtered.nonpeaks.bed"
    expected.write_text("chr1\t100\t200\n", encoding="utf-8")
    (tmp_path / "K562_auxiliary").mkdir()
    monkeypatch.setattr(
        module,
        "checked_repository_path",
        lambda path, *, must_exist: Path(path),
    )

    outputs = module.discover_outputs(
        SimpleNamespace(stage="prep-nonpeaks", output_prefix=prefix)
    )

    assert outputs == [expected]
