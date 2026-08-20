#!/usr/bin/env python3
"""Smoke-test a frozen fp-tools GUI executable."""

from __future__ import annotations

import argparse
import os
import socket
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
        process.terminate()
    try:
        return process.communicate(timeout=20)[0] or ""
    except subprocess.TimeoutExpired:
        process.kill()
        return process.communicate(timeout=20)[0] or ""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable")
    parser.add_argument("--timeout", type=float, default=90.0)
    args = parser.parse_args()

    executable = Path(args.executable).resolve()
    help_run = subprocess.run(
        [str(executable), "--fp-tools-internal-command", "summarize-motifs", "--help"],
        capture_output=True,
        text=True,
        timeout=args.timeout,
    )
    if help_run.returncode != 0 or "Summarize MEME/Tomtom" not in help_run.stdout:
        raise SystemExit(f"Internal command dispatch failed:\n{help_run.stdout}\n{help_run.stderr}")

    port = free_port()
    with tempfile.TemporaryDirectory(prefix="fp-tools-desktop-smoke-") as run_dir:
        process = subprocess.Popen(
            [str(executable), "--no-browser", "--port", str(port), "--run-dir", run_dir],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env={**os.environ, "STREAMLIT_BROWSER_GATHER_USAGE_STATS": "false"},
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
