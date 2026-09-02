#!/usr/bin/env python3
"""Run a checksum-validated sequence of research stages in the background.

The queue is intentionally small and generic: a JSON plan declares one
prerequisite manifest and an ordered set of commands.  Existing stages are
resumed only when every artifact declared by their manifest still matches its
checksum.  This is benchmark infrastructure, not a public fp-tools command.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY = SCRIPT_DIR.parent.parent
PLAN_SCHEMA = "fp-tools-background-experiment-plan-v1"
CONTROLLER_SCHEMA = "fp-tools-background-experiment-controller-v1"


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def repository_path(path: str | Path, repository: Path = REPOSITORY) -> Path:
    value = Path(path)
    if not value.is_absolute():
        value = repository / value
    value = value.resolve()
    try:
        value.relative_to(repository.resolve())
    except ValueError as exc:
        raise ValueError(f"experiment path leaves repository: {value}") from exc
    return value


def output_records(document: dict[str, Any]) -> list[dict[str, Any]]:
    outputs = document.get("outputs")
    if isinstance(outputs, list):
        records = outputs
    elif isinstance(outputs, dict):
        records = list(outputs.values())
    else:
        raise ValueError("stage manifest does not declare outputs")
    if not records or not all(isinstance(record, dict) for record in records):
        raise ValueError("stage manifest has no valid output records")
    return records


def validate_stage_manifest(
    path: str | Path,
    specification: dict[str, Any],
    *,
    repository: Path = REPOSITORY,
) -> dict[str, Any]:
    manifest_path = repository_path(path, repository)
    document = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_schema = specification.get("schema")
    if expected_schema and document.get("schema") != expected_schema:
        raise ValueError(f"unexpected schema in {manifest_path}")
    expected_stage = specification.get("stage")
    if expected_stage and document.get("stage") != expected_stage:
        raise ValueError(f"unexpected stage in {manifest_path}")
    expected_commit = specification.get("source_commit")
    if expected_commit and document.get("source_commit") != expected_commit:
        raise ValueError(f"unexpected source commit in {manifest_path}")
    if specification.get("require_completed", False) and document.get("completed") is not True:
        raise ValueError(f"stage is not marked complete in {manifest_path}")
    for record in output_records(document):
        output = repository_path(record.get("path", ""), repository)
        if not output.is_file() or output.stat().st_size <= 0:
            raise ValueError(f"missing or empty stage output: {output}")
        expected = record.get("sha256")
        if not expected or file_sha256(output) != expected:
            raise ValueError(f"stage output checksum mismatch: {output}")
    return document


def process_is_alive(pid: int | None) -> bool:
    if pid is None:
        return True
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def wait_for_prerequisite(
    specification: dict[str, Any],
    *,
    poll_seconds: float,
    timeout_seconds: float,
    watch_pid: int | None,
    repository: Path = REPOSITORY,
) -> dict[str, Any]:
    manifest = repository_path(specification["manifest"], repository)
    started = time.monotonic()
    while not manifest.is_file():
        elapsed = time.monotonic() - started
        if timeout_seconds > 0 and elapsed >= timeout_seconds:
            raise TimeoutError(f"timed out waiting for {manifest}")
        if not process_is_alive(watch_pid):
            time.sleep(min(2.0, poll_seconds))
            if not manifest.is_file():
                raise RuntimeError(
                    f"watched process {watch_pid} exited without creating {manifest}"
                )
        print(f"waiting for prerequisite: {manifest} ({elapsed:.0f}s)", flush=True)
        time.sleep(poll_seconds)
    return validate_stage_manifest(manifest, specification, repository=repository)


def expand_command(
    command: Sequence[str], *, repository: Path = REPOSITORY
) -> list[str]:
    replacements = {
        "{python}": sys.executable,
        "{root}": str(repository.resolve()),
    }
    expanded: list[str] = []
    for item in command:
        value = str(item)
        for token, replacement in replacements.items():
            value = value.replace(token, replacement)
        expanded.append(value)
    return expanded


def atomic_json(path: Path, document: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def run_plan(
    plan_path: Path,
    *,
    poll_seconds: float,
    timeout_seconds: float,
    watch_pid: int | None,
    repository: Path = REPOSITORY,
) -> dict[str, Any]:
    plan_path = repository_path(plan_path, repository)
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("schema") != PLAN_SCHEMA:
        raise ValueError("unsupported background experiment plan")
    if not plan.get("stages"):
        raise ValueError("background experiment plan has no stages")

    plan_digest = file_sha256(plan_path)
    wait_for_prerequisite(
        plan["prerequisite"],
        poll_seconds=poll_seconds,
        timeout_seconds=timeout_seconds,
        watch_pid=watch_pid,
        repository=repository,
    )
    stage_records = []
    for stage in plan["stages"]:
        manifest = repository_path(stage["manifest"], repository)
        if manifest.is_file():
            validate_stage_manifest(manifest, stage, repository=repository)
            disposition = "resumed"
        else:
            command = expand_command(stage["command"], repository=repository)
            print(f"starting stage {stage['name']}: {' '.join(command)}", flush=True)
            subprocess.run(command, cwd=repository, check=True)
            validate_stage_manifest(manifest, stage, repository=repository)
            disposition = "executed"
        stage_records.append(
            {
                "name": stage["name"],
                "disposition": disposition,
                "manifest": str(manifest),
                "manifest_sha256": file_sha256(manifest),
            }
        )

    controller = {
        "schema": CONTROLLER_SCHEMA,
        "name": plan.get("name", plan_path.stem),
        "plan": str(plan_path),
        "plan_sha256": plan_digest,
        "completed": True,
        "stages": stage_records,
    }
    controller_path = repository_path(plan["controller_manifest"], repository)
    atomic_json(controller_path, controller)
    print(controller_path)
    return controller


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--poll-seconds", type=float, default=30.0)
    parser.add_argument("--timeout-seconds", type=float, default=43_200.0)
    parser.add_argument("--watch-pid", type=int)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.poll_seconds <= 0 or args.timeout_seconds < 0:
        raise ValueError("poll and timeout intervals are invalid")
    run_plan(
        args.plan,
        poll_seconds=args.poll_seconds,
        timeout_seconds=args.timeout_seconds,
        watch_pid=args.watch_pid,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
