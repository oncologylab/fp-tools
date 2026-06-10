"""Optional YAML-driven batch/config runner for fp-tools.

Direct CLI usage remains primary. This module adds an extra path for:
- replaying GUI-saved configs
- running batch sample lists
- running batch BINDetect comparison lists
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from fp_tools.gui_config import JobSpec, canonical_tool_name, dump_yaml_config, expand_jobs, load_yaml_config, normalize_config


def main() -> None:
    parser = argparse.ArgumentParser(description="Run fp-tools jobs from a YAML config file.")
    parser.add_argument("--config", required=True, help="Path to YAML config.")
    parser.add_argument("--run-root", default=None, help="Optional directory for run metadata/logs.")
    parser.add_argument("--only", nargs="*", default=None, help="Optional tool filter, e.g. BINDetect.")
    parser.add_argument("--dry-run", action="store_true", help="Print expanded commands without running.")
    parser.add_argument("--list-jobs", action="store_true", help="List expanded jobs and exit.")
    parser.add_argument("--fail-fast", action="store_true", help="Stop at first failed job.")
    args = parser.parse_args()

    only_tools = {canonical_tool_name(tool_name) for tool_name in (args.only or [])}
    exit_code = run_config_file(
        args.config,
        run_root=args.run_root,
        only_tools=only_tools or None,
        dry_run=args.dry_run,
        list_jobs=args.list_jobs,
        fail_fast=args.fail_fast,
    )
    raise SystemExit(exit_code)


def run_config_file(
    config_path: str | os.PathLike[str],
    run_root: str | os.PathLike[str] | None = None,
    only_tools: set[str] | None = None,
    dry_run: bool = False,
    list_jobs: bool = False,
    fail_fast: bool = False,
) -> int:
    config = normalize_config(load_yaml_config(config_path))
    jobs = expand_jobs(config, only_tools=only_tools)
    if not jobs:
        print("No jobs matched the current config/filter.", file=sys.stderr)
        return 1

    if list_jobs or dry_run:
        for job in jobs:
            print(f"[{job.tool}] {job.job_id}: {' '.join(job.command)}")
        if dry_run or list_jobs:
            return 0

    root = Path(run_root or config.get("run_root") or _default_run_root()).expanduser()
    root.mkdir(parents=True, exist_ok=True)

    batch_index = root / "batch_index.tsv"
    with batch_index.open("w", encoding="utf-8") as handle:
        handle.write("job_id\ttool\tstatus\texit_code\trun_dir\n")

    exit_code = 0
    for job in jobs:
        code = run_job(job, root)
        status = "succeeded" if code == 0 else "failed"
        with batch_index.open("a", encoding="utf-8") as handle:
            handle.write(f"{job.job_id}\t{job.tool}\t{status}\t{code}\t{root / job.job_id}\n")
        if code != 0:
            exit_code = code
            if fail_fast:
                break
    return exit_code


def run_job(job: JobSpec, run_root: Path) -> int:
    run_dir = run_root / job.job_id
    run_dir.mkdir(parents=True, exist_ok=True)

    dump_yaml_config(
        {
            "version": 1,
            "run_mode": "single",
            "defaults": {},
            "samples" if job.section == "samples" else "comparisons": [
                {"job_id": job.job_id, "tool": job.tool, **job.params}
            ],
            "comparisons" if job.section == "samples" else "samples": [],
        },
        run_dir / "config.yml",
    )

    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    status_path = run_dir / "status.json"
    cache_dir = run_dir / ".cache"
    mpl_dir = run_dir / ".mplconfig"
    cache_dir.mkdir(parents=True, exist_ok=True)
    mpl_dir.mkdir(parents=True, exist_ok=True)

    status = {
        "tool": job.tool,
        "job_id": job.job_id,
        "status": "running",
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "command": job.command,
    }
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")

    command = list(job.command)
    command[0] = _resolve_executable(command[0])
    status["command"] = command
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
    (run_dir / "command.txt").write_text(" ".join(command) + "\n", encoding="utf-8")

    env = os.environ.copy()
    env.setdefault("XDG_CACHE_HOME", str(cache_dir))
    env.setdefault("MPLCONFIGDIR", str(mpl_dir))

    with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open("w", encoding="utf-8") as stderr_handle:
        process = subprocess.run(command, stdout=stdout_handle, stderr=stderr_handle, text=True, env=env)

    status["status"] = "succeeded" if process.returncode == 0 else "failed"
    status["finished_at"] = datetime.now().isoformat(timespec="seconds")
    status["exit_code"] = process.returncode
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
    return process.returncode


def _default_run_root() -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path.cwd() / f"fp-tools-batch-{stamp}"


def _resolve_executable(name: str) -> str:
    local = Path(sys.executable).parent / name
    if local.exists():
        return str(local)
    found = shutil.which(name)
    if found:
        return found
    return name


if __name__ == "__main__":
    main()
