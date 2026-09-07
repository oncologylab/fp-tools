"""Helpers for GUI-managed run folders and config execution."""

from __future__ import annotations

import csv
import os
import json
import subprocess
import tempfile
import threading
from datetime import datetime
from pathlib import Path
from typing import Any

from fp_tools.gui_config import (
    dump_yaml_config,
    expand_jobs,
    load_yaml_config,
    normalize_config,
)
from fp_tools.utils.subprocess_commands import fp_tools_subprocess_command


_ACTIVE_PROCESSES: dict[Path, subprocess.Popen] = {}
_PROCESS_LOCK = threading.RLock()


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
    command = fp_tools_subprocess_command("run-yaml-workflow", ["--config", str(config_path)])
    if run_root:
        command.extend(["--run-root", str(run_root)])
    return subprocess.run(command, capture_output=True, text=True)


def _write_status_atomic(path: Path, status: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(status, handle, indent=2)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _read_status(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return None
    return value if isinstance(value, dict) else None


def _expected_job_count(config_path: Path) -> int:
    config = normalize_config(load_yaml_config(config_path))
    return len(expand_jobs(config))


def _batch_progress(batch_index: Path) -> dict[str, Any] | None:
    if not batch_index.is_file():
        return None
    try:
        with batch_index.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, UnicodeError, csv.Error):
        return None
    completed = [
        row for row in rows if str(row.get("status") or "") in {"succeeded", "failed"}
    ]
    failed = [row for row in completed if row.get("status") == "failed"]
    failed_codes = []
    for row in failed:
        try:
            failed_codes.append(int(str(row.get("exit_code") or "1")))
        except ValueError:
            failed_codes.append(1)
    return {
        "rows": len(rows),
        "completed": len(completed),
        "failed": len(failed),
        "exit_code": next((code for code in failed_codes if code), 1),
    }


def _finish_status(
    path: Path,
    status: dict[str, Any],
    *,
    state: str,
    exit_code: int,
    completed_jobs: int,
    detail: str | None = None,
) -> dict[str, Any]:
    status["status"] = state
    status["exit_code"] = int(exit_code)
    status["completed_jobs"] = int(completed_jobs)
    status["finished_at"] = datetime.now().isoformat(timespec="seconds")
    if detail:
        status["status_detail"] = detail
    else:
        status.pop("status_detail", None)
    _write_status_atomic(path, status)
    return status


def _reconcile_status(
    run_dir: Path,
    status: dict[str, Any],
    *,
    process_finished: bool = False,
    process_returncode: int | None = None,
) -> dict[str, Any]:
    status_path = run_dir / "status.json"
    try:
        expected = max(0, int(status.get("expected_jobs") or 0))
    except (TypeError, ValueError):
        expected = 0
    if expected == 0:
        config_path = run_dir / "config.yml"
        if config_path.is_file():
            try:
                expected = _expected_job_count(config_path)
            except (OSError, ValueError):
                expected = 0
        if expected:
            status["expected_jobs"] = expected

    progress = _batch_progress(run_dir / "batch_index.tsv")
    completed = int(progress["completed"]) if progress else 0
    status["completed_jobs"] = completed
    if expected and completed >= expected:
        failed = int(progress["failed"]) if progress else 0
        return _finish_status(
            status_path,
            status,
            state="failed" if failed else "succeeded",
            exit_code=int(progress["exit_code"]) if failed else 0,
            completed_jobs=completed,
        )

    if process_finished:
        return_code = int(process_returncode or 0)
        detail = (
            f"Runner exited after {completed} of {expected} expected jobs were recorded."
            if expected
            else "Runner exited before a complete batch record was written."
        )
        return _finish_status(
            status_path,
            status,
            state="failed",
            exit_code=return_code or 1,
            completed_jobs=completed,
            detail=detail,
        )

    _write_status_atomic(status_path, status)
    return status


def _watch_process(run_dir: Path, process: subprocess.Popen) -> None:
    return_code = int(process.wait())
    key = run_dir.resolve()
    with _PROCESS_LOCK:
        if _ACTIVE_PROCESSES.get(key) is process:
            _ACTIVE_PROCESSES.pop(key, None)
        status_path = run_dir / "status.json"
        status = _read_status(status_path)
        if status is not None and status.get("status") == "running":
            _reconcile_status(
                run_dir,
                status,
                process_finished=True,
                process_returncode=return_code,
            )


def launch_config_async(
    config_path: str | os.PathLike[str],
    run_root: str | os.PathLike[str],
    tool_label: str,
) -> tuple[Path, int]:
    run_root = Path(run_root).expanduser()
    stdout_path = run_root / "launcher_stdout.log"
    stderr_path = run_root / "launcher_stderr.log"
    status_path = run_root / "status.json"
    expected_jobs = _expected_job_count(Path(config_path).expanduser())
    if expected_jobs < 1:
        raise ValueError("The GUI configuration does not contain any runnable jobs")
    command = fp_tools_subprocess_command(
        "run-yaml-workflow",
        ["--config", str(config_path), "--run-root", str(run_root)],
    )

    with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open("w", encoding="utf-8") as stderr_handle:
        session_options = (
            {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP}
            if os.name == "nt"
            else {"start_new_session": True}
        )
        process = subprocess.Popen(
            command,
            stdout=stdout_handle,
            stderr=stderr_handle,
            text=True,
            cwd=Path.cwd(),
            **session_options,
        )

    status = {
        "tool": tool_label,
        "job_id": run_root.name,
        "status": "running",
        "started_at": datetime.now().isoformat(timespec="seconds"),
        "command": command,
        "pid": process.pid,
        "expected_jobs": expected_jobs,
        "completed_jobs": 0,
    }
    key = run_root.resolve()
    with _PROCESS_LOCK:
        _ACTIVE_PROCESSES[key] = process
        _write_status_atomic(status_path, status)
    watcher = threading.Thread(
        target=_watch_process,
        args=(run_root, process),
        name=f"fp-tools-gui-job-{process.pid}",
        daemon=True,
    )
    watcher.start()
    return status_path, process.pid


def refresh_run_status(run_dir: str | os.PathLike[str]) -> dict[str, Any] | None:
    run_dir = Path(run_dir).expanduser()
    status_path = run_dir / "status.json"
    if not status_path.exists():
        return None
    status = _read_status(status_path)
    if status is None:
        return None
    if status.get("status") != "running":
        return status

    key = run_dir.resolve()
    with _PROCESS_LOCK:
        status = _reconcile_status(run_dir, status)
        if status.get("status") != "running":
            return status
        process = _ACTIVE_PROCESSES.get(key)
        if process is not None:
            return_code = process.poll()
            if return_code is None:
                return status
            _ACTIVE_PROCESSES.pop(key, None)
            return _reconcile_status(
                run_dir,
                status,
                process_finished=True,
                process_returncode=int(return_code),
            )

        pid = status.get("pid")
        if pid:
            try:
                if _pid_is_alive(int(pid)):
                    return status
            except (TypeError, ValueError):
                pass
        return _reconcile_status(
            run_dir,
            status,
            process_finished=True,
            process_returncode=int(status.get("exit_code") or 1),
        )


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
