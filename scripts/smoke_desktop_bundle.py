#!/usr/bin/env python3
"""Smoke-test a frozen fp-tools GUI executable."""

from __future__ import annotations

import argparse
import os
import socket
import signal
import subprocess
import tempfile
import time
import urllib.request
from pathlib import Path

def free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as handle:
        handle.bind(("127.0.0.1", 0))
        return int(handle.getsockname()[1])


def stop_process_tree(process: subprocess.Popen[str]) -> str:
    """Stop the GUI and any server child while preserving captured output."""

    if os.name == "nt":
        subprocess.run(
            ["taskkill", "/PID", str(process.pid), "/T", "/F"],
            capture_output=True,
            text=True,
            check=False,
        )
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return ""
    try:
        return process.communicate(timeout=20)[0] or ""
    except subprocess.TimeoutExpired:
        if os.name == "nt":
            process.kill()
        else:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        return process.communicate(timeout=20)[0] or ""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable")
    parser.add_argument("--timeout", type=float, default=90.0)
    args = parser.parse_args()

    executable = Path(args.executable).resolve()
    help_run = subprocess.run(
        [str(executable), "--fp-tools-internal-smoke-command-helps"],
        capture_output=True,
        text=True,
        timeout=args.timeout,
    )
    if help_run.returncode != 0 or "Validated " not in help_run.stdout:
        raise SystemExit(
            f"Internal command dispatch audit failed:\n{help_run.stdout}\n{help_run.stderr}"
        )

    raw_read_run = subprocess.run(
        [str(executable), "--fp-tools-internal-command", "prepare-atac", "--help"],
        capture_output=True,
        text=True,
        timeout=args.timeout,
    )
    if raw_read_run.returncode == 0:
        raise SystemExit("Desktop bundle unexpectedly exposes prepare-atac")

    examples_run = subprocess.run(
        [str(executable), "--fp-tools-internal-list-gui-examples"],
        capture_output=True,
        text=True,
        timeout=args.timeout,
    )
    if examples_run.returncode != 0 or "bulk_footprinting_bam.yml" not in examples_run.stdout:
        raise SystemExit(
            f"Desktop bundle is missing GUI examples:\n{examples_run.stdout}\n{examples_run.stderr}"
        )

    with tempfile.TemporaryDirectory(prefix="fp-tools-desktop-smoke-") as run_dir:
        run_path = Path(run_dir)
        fixture_root = Path(__file__).resolve().parents[1] / "test_data"
        call_output = run_path / "call_footprints"
        call_run = subprocess.run(
            [
                str(executable),
                "--fp-tools-internal-command",
                "call-footprints",
                "--signals",
                str(fixture_root / "Bcell_corrected.bw"),
                "--regions",
                str(fixture_root / "plot_regions.bed"),
                "--outdir",
                str(call_output),
                "--score",
                "footprint",
                "--cores",
                "1",
            ],
            capture_output=True,
            text=True,
            timeout=args.timeout,
            cwd=run_dir,
        )
        if call_run.returncode != 0 or not any(call_output.glob("*.bw")):
            raise SystemExit(
                f"Frozen call-footprints smoke failed:\n{call_run.stdout}\n{call_run.stderr}"
            )

        normalize_output = run_path / "normalize_bigwig"
        normalize_run = subprocess.run(
            [
                str(executable),
                "--fp-tools-internal-command",
                "normalize-bigwig",
                "--bigwigs",
                str(fixture_root / "Bcell_corrected.bw"),
                str(fixture_root / "Tcell_corrected.bw"),
                "--background",
                str(fixture_root / "plot_regions.bed"),
                "--outdir",
                str(normalize_output),
                "--method",
                "background-scale",
                "--stat",
                "q95",
                "--target",
                "median",
                "--workers",
                "2",
            ],
            capture_output=True,
            text=True,
            timeout=args.timeout,
            cwd=run_dir,
        )
        if normalize_run.returncode != 0 or len(list(normalize_output.glob("*.bw"))) != 2:
            raise SystemExit(
                "Frozen multiprocessing normalize-bigwig smoke failed:\n"
                f"{normalize_run.stdout}\n{normalize_run.stderr}"
            )

        port = free_port()
        process = subprocess.Popen(
            [str(executable), "--no-browser", "--port", str(port), "--run-dir", run_dir],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=run_dir,
            env={**os.environ, "STREAMLIT_BROWSER_GATHER_USAGE_STATS": "false"},
            start_new_session=os.name != "nt",
        )
        try:
            deadline = time.monotonic() + args.timeout
            last_error: Exception | None = None
            while time.monotonic() < deadline:
                if process.poll() is not None:
                    output = process.stdout.read() if process.stdout else ""
                    raise SystemExit(f"GUI exited before becoming ready ({process.returncode}):\n{output}")
                try:
                    with urllib.request.urlopen(f"http://127.0.0.1:{port}/_stcore/health", timeout=2) as response:
                        if response.status == 200:
                            print(f"Desktop GUI ready on port {port}")
                            return 0
                except Exception as exc:  # startup polling
                    last_error = exc
                    time.sleep(0.5)
            raise SystemExit(f"GUI did not become ready within {args.timeout:g} seconds: {last_error}")
        finally:
            output = stop_process_tree(process)
            if process.returncode not in (0, -15, 1):
                print(output)


if __name__ == "__main__":
    raise SystemExit(main())
