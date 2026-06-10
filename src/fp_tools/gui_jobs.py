"""Helpers for GUI-managed run folders and config execution."""

from __future__ import annotations

import os
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

from fp_tools.gui_config import dump_yaml_config, normalize_config


def default_gui_run_dir() -> Path:
    return Path.home() / "fp-tools-gui-runs"


def materialize_run_config(config: dict[str, Any], run_root: str | os.PathLike[str] | None = None, label: str = "run") -> tuple[Path, Path]:
    root = Path(run_root or default_gui_run_dir()).expanduser()
    root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_label = "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in label).strip("_") or "run"
    run_dir = root / f"{stamp}_{safe_label}"
    run_dir.mkdir(parents=True, exist_ok=True)
    config_path = run_dir / "config.yml"
    dump_yaml_config(normalize_config(config), config_path)
    return run_dir, config_path


def run_config_sync(config_path: str | os.PathLike[str], run_root: str | os.PathLike[str] | None = None) -> subprocess.CompletedProcess[str]:
    command = [sys.executable, "-m", "fp_tools.cli_batch", "--config", str(config_path)]
    if run_root:
        command.extend(["--run-root", str(run_root)])
    return subprocess.run(command, capture_output=True, text=True)


def launch_config_async(
    config_path: str | os.PathLike[str],
    run_root: str | os.PathLike[str],
    tool_label: str,
) -> tuple[Path, int]:
    run_root = Path(run_root).expanduser()
    stdout_path = run_root / "launcher_stdout.log"
    stderr_path = run_root / "launcher_stderr.log"
    status_path = run_root / "status.json"
    command = [sys.executable, "-m", "fp_tools.cli_batch", "--config", str(config_path), "--run-root", str(run_root)]

    with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open("w", encoding="utf-8") as stderr_handle:
        process = subprocess.Popen(
            command,
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
            start_new_session=True,
            cwd=Path.cwd(),
        )

    status = {
        "tool": tool_label,
        "job_id": run_root.name,
        "status": "running",
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "command": command,
        "pid": process.pid,
    }
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
    return status_path, process.pid


def refresh_run_status(run_dir: str | os.PathLike[str]) -> dict[str, Any] | None:
    run_dir = Path(run_dir).expanduser()
    status_path = run_dir / "status.json"
    if not status_path.exists():
        return None
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("status") != "running":
        return status

    pid = status.get("pid")
    if pid and _pid_is_alive(int(pid)):
        return status

    batch_index = run_dir / "batch_index.tsv"
    if not batch_index.exists():
        status["status"] = "failed"
        status["finished_at"] = datetime.now().isoformat(timespec="seconds")
        status.setdefault("exit_code", 1)
        status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
        return status

    lines = [line.strip().split("\t") for line in batch_index.read_text(encoding="utf-8").splitlines() if line.strip()]
    rows = lines[1:] if len(lines) > 1 else []
    failed = any(len(row) >= 4 and row[2] == "failed" for row in rows)
    still_blank = any(len(row) >= 4 and row[2] == "" for row in rows)
    if still_blank:
        return status
    status["status"] = "failed" if failed else "succeeded"
    status["finished_at"] = datetime.now().isoformat(timespec="seconds")
    status["exit_code"] = 1 if failed else 0
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")
    return status


def _pid_is_alive(pid: int) -> bool:
    proc_stat = Path(f"/proc/{pid}/stat")
    if proc_stat.exists():
        try:
            state = proc_stat.read_text(encoding="utf-8").split()[2]
            if state == "Z":
                return False
        except Exception:
            pass
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False
