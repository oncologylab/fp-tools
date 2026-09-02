from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "benchmarks" / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import run_background_experiment_queue as queue  # noqa: E402


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _write_manifest(root: Path, output: Path, *, checksum: str | None = None) -> Path:
    manifest = root / "stage.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "test-stage-v1",
                "stage": "prediction",
                "completed": True,
                "source_commit": "abc123",
                "outputs": [
                    {
                        "path": str(output),
                        "sha256": checksum or _sha(output),
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    return manifest


def test_stage_manifest_requires_matching_outputs(tmp_path: Path) -> None:
    output = tmp_path / "output.bw"
    output.write_bytes(b"nonempty")
    manifest = _write_manifest(tmp_path, output)
    specification = {
        "schema": "test-stage-v1",
        "stage": "prediction",
        "source_commit": "abc123",
        "require_completed": True,
    }
    document = queue.validate_stage_manifest(
        manifest, specification, repository=tmp_path
    )
    assert document["completed"] is True

    output.write_bytes(b"changed")
    with pytest.raises(ValueError, match="checksum mismatch"):
        queue.validate_stage_manifest(manifest, specification, repository=tmp_path)


def test_stage_manifest_rejects_paths_outside_repository(tmp_path: Path) -> None:
    repository = tmp_path / "repository"
    repository.mkdir()
    outside = tmp_path / "outside.tsv"
    outside.write_text("value\n1\n", encoding="utf-8")
    manifest = _write_manifest(repository, outside)
    with pytest.raises(ValueError, match="leaves repository"):
        queue.validate_stage_manifest(
            manifest,
            {"schema": "test-stage-v1"},
            repository=repository,
        )


def test_command_expansion_uses_current_interpreter(tmp_path: Path) -> None:
    command = queue.expand_command(
        ["{python}", "tool.py", "--root", "{root}"], repository=tmp_path
    )
    assert command[0] == sys.executable
    assert command[-1] == str(tmp_path.resolve())
