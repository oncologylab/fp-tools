import json
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest import mock

from fp_tools import gui_jobs
from fp_tools.gui_config import dump_yaml_config


class _ControlledProcess:
    def __init__(self, pid: int, returncode: int):
        self.pid = pid
        self.returncode = returncode
        self.finished = threading.Event()
        self.wait_called = threading.Event()

    def wait(self):
        self.wait_called.set()
        self.finished.wait(timeout=5)
        return self.returncode

    def poll(self):
        return self.returncode if self.finished.is_set() else None


def _write_batch(path: Path, rows: list[tuple[str, str, str, int]]) -> None:
    text = "job_id\ttool\tstatus\texit_code\trun_dir\n"
    for job_id, tool, status, exit_code in rows:
        text += f"{job_id}\t{tool}\t{status}\t{exit_code}\t{path / job_id}\n"
    (path / "batch_index.tsv").write_text(text, encoding="utf-8")


class GuiJobStatusTest(unittest.TestCase):
    def setUp(self):
        self.processes: list[_ControlledProcess] = []
        with gui_jobs._PROCESS_LOCK:
            gui_jobs._ACTIVE_PROCESSES.clear()

    def tearDown(self):
        for process in self.processes:
            process.finished.set()
        deadline = time.monotonic() + 1
        while time.monotonic() < deadline:
            with gui_jobs._PROCESS_LOCK:
                if not gui_jobs._ACTIVE_PROCESSES:
                    break
            time.sleep(0.01)
        with gui_jobs._PROCESS_LOCK:
            gui_jobs._ACTIVE_PROCESSES.clear()

    def _config(self, path: Path) -> Path:
        config = path / "config.yml"
        dump_yaml_config(
            {
                "version": 1,
                "run_mode": "batch",
                "defaults": {},
                "samples": [
                    {"sample_id": "one", "tool": "call-footprints"},
                    {"sample_id": "two", "tool": "call-footprints"},
                ],
                "comparisons": [],
            },
            config,
        )
        return config

    def _launch(self, root: Path, returncode: int = 0) -> _ControlledProcess:
        process = _ControlledProcess(41001 + len(self.processes), returncode)
        self.processes.append(process)
        with mock.patch.object(gui_jobs.subprocess, "Popen", return_value=process):
            gui_jobs.launch_config_async(self._config(root), root, "bulk-footprinting")
        return process

    def test_multi_job_success_is_reaped_and_parent_status_succeeds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            process = self._launch(root)
            initial = json.loads((root / "status.json").read_text(encoding="utf-8"))
            self.assertEqual(initial["expected_jobs"], 2)
            self.assertEqual(initial["completed_jobs"], 0)
            self.assertTrue(process.wait_called.wait(timeout=1))
            _write_batch(
                root,
                [
                    ("one", "call-footprints", "succeeded", 0),
                    ("two", "call-footprints", "succeeded", 0),
                ],
            )
            process.finished.set()
            status = gui_jobs.refresh_run_status(root)
            self.assertEqual(status["status"], "succeeded")
            self.assertEqual(status["completed_jobs"], 2)
            self.assertEqual(status["exit_code"], 0)
            deadline = time.monotonic() + 1
            while time.monotonic() < deadline:
                with gui_jobs._PROCESS_LOCK:
                    if root.resolve() not in gui_jobs._ACTIVE_PROCESSES:
                        break
                time.sleep(0.01)
            self.assertNotIn(root.resolve(), gui_jobs._ACTIVE_PROCESSES)
            self.assertFalse(list(root.glob(".status.json.*.tmp")))

    def test_failed_child_marks_parent_failed_with_its_exit_code(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            process = self._launch(root, returncode=9)
            _write_batch(
                root,
                [
                    ("one", "call-footprints", "succeeded", 0),
                    ("two", "call-footprints", "failed", 9),
                ],
            )
            process.finished.set()
            status = gui_jobs.refresh_run_status(root)
            self.assertEqual(status["status"], "failed")
            self.assertEqual(status["completed_jobs"], 2)
            self.assertEqual(status["exit_code"], 9)

    def test_restart_reconciles_complete_batch_before_live_pid(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "status.json").write_text(
                json.dumps(
                    {
                        "status": "running",
                        "pid": 12345,
                        "expected_jobs": 2,
                    }
                ),
                encoding="utf-8",
            )
            _write_batch(
                root,
                [
                    ("one", "call-footprints", "succeeded", 0),
                    ("two", "call-footprints", "succeeded", 0),
                ],
            )
            with mock.patch.object(gui_jobs, "_pid_is_alive") as pid_check:
                status = gui_jobs.refresh_run_status(root)
            pid_check.assert_not_called()
            self.assertEqual(status["status"], "succeeded")

    def test_restart_keeps_partial_batch_running_when_pid_is_alive(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "status.json").write_text(
                json.dumps(
                    {
                        "status": "running",
                        "pid": 12345,
                        "expected_jobs": 2,
                    }
                ),
                encoding="utf-8",
            )
            _write_batch(
                root, [("one", "call-footprints", "succeeded", 0)]
            )
            with mock.patch.object(gui_jobs, "_pid_is_alive", return_value=True):
                status = gui_jobs.refresh_run_status(root)
            self.assertEqual(status["status"], "running")
            self.assertEqual(status["completed_jobs"], 1)

    def test_stale_pid_with_incomplete_batch_fails_clearly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "status.json").write_text(
                json.dumps(
                    {
                        "status": "running",
                        "pid": 999999,
                        "expected_jobs": 2,
                    }
                ),
                encoding="utf-8",
            )
            _write_batch(
                root, [("one", "call-footprints", "succeeded", 0)]
            )
            with mock.patch.object(gui_jobs, "_pid_is_alive", return_value=False):
                status = gui_jobs.refresh_run_status(root)
            self.assertEqual(status["status"], "failed")
            self.assertEqual(status["completed_jobs"], 1)
            self.assertIn("1 of 2", status["status_detail"])

    def test_clean_runner_exit_without_complete_index_is_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            process = self._launch(root)
            process.finished.set()
            deadline = time.monotonic() + 1
            status = None
            while time.monotonic() < deadline:
                status = gui_jobs.refresh_run_status(root)
                if status and status.get("status") != "running":
                    break
                time.sleep(0.01)
            self.assertEqual(status["status"], "failed")
            self.assertEqual(status["exit_code"], 1)
            self.assertIn("0 of 2", status["status_detail"])


if __name__ == "__main__":
    unittest.main()
