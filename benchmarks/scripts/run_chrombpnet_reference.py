#!/usr/bin/env python3
"""Run the pinned official ChromBPNet GPU reference in an isolated container.

The wrapper verifies both the source commit and container digest, places the
pinned checkout first on ``PYTHONPATH``, and records checksums for every input
and output.  It works around hosts where the NVIDIA Docker prestart hook lacks
NVML permission by mounting only the GPU devices and driver libraries needed
by the TensorFlow 2.8 image.  No ChromBPNet dependency enters fp-tools.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shlex
import subprocess
from typing import Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
PINNED_SOURCE = (
    REPOSITORY / "benchmarks/results/footprint_external_references/chrombpnet"
)
PINNED_COMMIT = "09938fdb4397ec0006510e5251e48920a505d4de"
PINNED_IMAGE = (
    "kundajelab/chrombpnet@"
    "sha256:6f41e0f59fc025285645e2cbfd1bb6347431b4e0ab6088804d11843cd6aed169"
)
SCHEMA = "fp-tools-chrombpnet-reference-run-v1"


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_repository_path(path: str | Path, *, must_exist: bool) -> Path:
    value = Path(path)
    if not value.is_absolute():
        value = REPOSITORY / value
    value = value.resolve()
    try:
        value.relative_to(REPOSITORY.resolve())
    except ValueError as exc:
        raise ValueError(
            f"ChromBPNet reference paths must stay in {REPOSITORY}"
        ) from exc
    if must_exist and not value.exists():
        raise FileNotFoundError(value)
    return value


def container_path(path: str | Path) -> str:
    value = checked_repository_path(path, must_exist=False)
    return str(Path("/work") / value.relative_to(REPOSITORY.resolve()))


def verify_source(source: Path = PINNED_SOURCE) -> str:
    commit = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if commit != PINNED_COMMIT:
        raise ValueError(
            f"ChromBPNet source commit changed: {commit} != {PINNED_COMMIT}"
        )
    status = subprocess.run(
        ["git", "-C", str(source), "status", "--porcelain"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if status:
        raise ValueError("ChromBPNet source checkout is dirty")
    return commit


def docker_image_record(image: str = PINNED_IMAGE) -> dict[str, object]:
    result = subprocess.run(
        ["docker", "image", "inspect", image],
        check=True,
        capture_output=True,
        text=True,
    )
    document = json.loads(result.stdout)[0]
    repo_digests = document.get("RepoDigests", [])
    expected_digest = image.split("@", 1)[1]
    if not any(str(value).endswith("@" + expected_digest) for value in repo_digests):
        raise ValueError("local ChromBPNet image does not match the pinned digest")
    return {
        "image": image,
        "image_id": document["Id"],
        "repo_digests": repo_digests,
        "size_bytes": int(document["Size"]),
    }


def driver_library(name: str) -> Path:
    result = subprocess.run(
        ["bash", "-lc", f"readlink -f /lib/x86_64-linux-gnu/{shlex.quote(name)}"],
        check=True,
        capture_output=True,
        text=True,
    )
    path = Path(result.stdout.strip())
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def docker_prefix(image: str = PINNED_IMAGE) -> list[str]:
    cuda = driver_library("libcuda.so.1")
    nvml = driver_library("libnvidia-ml.so.1")
    source_in_container = container_path(PINNED_SOURCE)
    return [
        "docker",
        "run",
        "--rm",
        "--device",
        "/dev/nvidia0",
        "--device",
        "/dev/nvidiactl",
        "--device",
        "/dev/nvidia-uvm",
        "--device",
        "/dev/nvidia-uvm-tools",
        "-v",
        f"{cuda}:/usr/local/nvidia/lib64/host_libcuda.so:ro",
        "-v",
        f"{nvml}:/usr/local/nvidia/lib64/host_libnvidia-ml.so:ro",
        "-v",
        f"{REPOSITORY.resolve()}:/work",
        "-w",
        "/work",
        "-e",
        f"PYTHONPATH={source_in_container}",
        "-e",
        "TF_FORCE_GPU_ALLOW_GROWTH=true",
        image,
    ]


def wrapped_container_command(arguments: Sequence[str]) -> list[str]:
    preamble = (
        "ln -sf /usr/local/nvidia/lib64/host_libcuda.so "
        "/usr/local/nvidia/lib64/libcuda.so.1 && "
        "ln -sf /usr/local/nvidia/lib64/host_libnvidia-ml.so "
        "/usr/local/nvidia/lib64/libnvidia-ml.so.1 && exec "
    )
    return docker_prefix() + ["bash", "-lc", preamble + shlex.join(arguments)]


def smoke_arguments() -> list[str]:
    code = (
        "import json,tensorflow as tf,chrombpnet.CHROMBPNET as c;"
        "g=tf.config.list_physical_devices('GPU');"
        "assert tf.__version__=='2.8.0' and g;"
        "print(json.dumps({'tensorflow':tf.__version__,'gpu':str(g[0]),"
        "'entrypoint':c.__file__},sort_keys=True))"
    )
    return ["python", "-c", code]


def stage_arguments(args: argparse.Namespace) -> list[str]:
    if args.stage == "smoke":
        return smoke_arguments()
    common = [
        "-g",
        container_path(args.genome),
        "-c",
        container_path(args.chrom_sizes),
        "-p",
        container_path(args.peaks),
        "-fl",
        container_path(args.fold),
    ]
    if args.stage == "prep-nonpeaks":
        return [
            "chrombpnet",
            "prep",
            "nonpeaks",
            *common,
            "-o",
            container_path(args.output_prefix),
            "-br",
            container_path(args.blacklist),
            "-s",
            str(args.seed),
        ]
    training = [
        "-ibam",
        container_path(args.bam),
        "-d",
        "ATAC",
        *common,
        "-n",
        container_path(args.nonpeaks),
        "-o",
        container_path(args.output_dir),
        "-s",
        str(args.seed),
        "-e",
        str(args.epochs),
        "-bs",
        str(args.batch_size),
    ]
    if args.stage == "bias":
        return ["chrombpnet", "bias", "pipeline", *training, "-b", "0.5"]
    if args.stage == "regulatory":
        return [
            "chrombpnet",
            "pipeline",
            *training,
            "-b",
            container_path(args.bias_model),
        ]
    if args.stage == "predict":
        return [
            "chrombpnet",
            "pred_bw",
            "-bm",
            container_path(args.bias_model),
            "-cm",
            container_path(args.chrombpnet_model),
            "-cmb",
            container_path(args.chrombpnet_nobias_model),
            "-r",
            container_path(args.peaks),
            "-g",
            container_path(args.genome),
            "-c",
            container_path(args.chrom_sizes),
            "-op",
            container_path(args.output_prefix),
            "-bs",
            str(args.batch_size),
        ]
    raise ValueError(f"unknown stage: {args.stage}")


def declared_inputs(args: argparse.Namespace) -> list[Path]:
    names = {
        "prep-nonpeaks": ("genome", "chrom_sizes", "peaks", "fold", "blacklist"),
        "bias": ("bam", "genome", "chrom_sizes", "peaks", "nonpeaks", "fold"),
        "regulatory": (
            "bam",
            "genome",
            "chrom_sizes",
            "peaks",
            "nonpeaks",
            "fold",
            "bias_model",
        ),
        "predict": (
            "genome",
            "chrom_sizes",
            "peaks",
            "bias_model",
            "chrombpnet_model",
            "chrombpnet_nobias_model",
        ),
    }.get(args.stage, ())
    return [
        checked_repository_path(getattr(args, name), must_exist=True) for name in names
    ]


def output_manifest_path(args: argparse.Namespace) -> Path:
    if args.stage in {"prep-nonpeaks", "predict"}:
        prefix = checked_repository_path(args.output_prefix, must_exist=False)
        return prefix.parent / f"{prefix.name}.{args.stage}.manifest.json"
    if args.stage in {"bias", "regulatory"}:
        directory = checked_repository_path(args.output_dir, must_exist=False)
        return directory.parent / f"{directory.name}.{args.stage}.manifest.json"
    return checked_repository_path(
        "benchmarks/results/footprint_external_references/chrombpnet_smoke.json",
        must_exist=False,
    )


def discover_outputs(args: argparse.Namespace) -> list[Path]:
    if args.stage == "smoke":
        return []
    if args.stage == "prep-nonpeaks":
        prefix = checked_repository_path(args.output_prefix, must_exist=False)
        return sorted(
            path for path in prefix.parent.glob(prefix.name + "*") if path.is_file()
        )
    if args.stage in {"bias", "regulatory"}:
        directory = checked_repository_path(args.output_dir, must_exist=True)
        return sorted(path for path in directory.rglob("*") if path.is_file())
    prefix = checked_repository_path(args.output_prefix, must_exist=False)
    return sorted(
        path for path in prefix.parent.glob(prefix.name + "*") if path.is_file()
    )


def add_common_reference_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path, required=True)
    parser.add_argument("--peaks", type=Path, required=True)
    parser.add_argument("--fold", type=Path, required=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="stage", required=True)
    subparsers.add_parser("smoke")
    prep = subparsers.add_parser("prep-nonpeaks")
    add_common_reference_arguments(prep)
    prep.add_argument("--blacklist", type=Path, required=True)
    prep.add_argument("--output-prefix", type=Path, required=True)
    prep.add_argument("--seed", type=int, default=2026)
    for name in ("bias", "regulatory"):
        stage = subparsers.add_parser(name)
        add_common_reference_arguments(stage)
        stage.add_argument("--bam", type=Path, required=True)
        stage.add_argument("--nonpeaks", type=Path, required=True)
        stage.add_argument("--output-dir", type=Path, required=True)
        stage.add_argument("--seed", type=int, default=2026)
        stage.add_argument("--epochs", type=int, default=50)
        stage.add_argument("--batch-size", type=int, default=64)
        if name == "regulatory":
            stage.add_argument("--bias-model", type=Path, required=True)
    predict = subparsers.add_parser("predict")
    predict.add_argument("--genome", type=Path, required=True)
    predict.add_argument("--chrom-sizes", type=Path, required=True)
    predict.add_argument("--peaks", type=Path, required=True)
    predict.add_argument("--bias-model", type=Path, required=True)
    predict.add_argument("--chrombpnet-model", type=Path, required=True)
    predict.add_argument("--chrombpnet-nobias-model", type=Path, required=True)
    predict.add_argument("--output-prefix", type=Path, required=True)
    predict.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    source_commit = verify_source()
    image = docker_image_record()
    inputs = declared_inputs(args)
    command = wrapped_container_command(stage_arguments(args))
    manifest_path = output_manifest_path(args)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    planned = {
        "schema": SCHEMA,
        "stage": args.stage,
        "source": str(PINNED_SOURCE),
        "source_commit": source_commit,
        "container": image,
        "command": command,
        "inputs": [
            {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": file_sha256(path),
            }
            for path in inputs
        ],
        "dry_run": bool(args.dry_run),
        "completed": False,
    }
    if args.dry_run:
        manifest_path.write_text(
            json.dumps(planned, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(shlex.join(command))
        return 0
    result = subprocess.run(
        command,
        cwd=REPOSITORY,
        check=True,
        text=True,
        capture_output=args.stage == "smoke",
    )
    if args.stage == "smoke":
        smoke = json.loads(result.stdout.strip().splitlines()[-1])
        planned["smoke"] = smoke
    outputs = discover_outputs(args)
    planned["outputs"] = [
        {"path": str(path), "bytes": path.stat().st_size, "sha256": file_sha256(path)}
        for path in outputs
    ]
    planned["completed"] = True
    manifest_path.write_text(
        json.dumps(planned, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
