"""Launcher for the isolated fp-tools Streamlit GUI."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import socket
import subprocess
import sys
from pathlib import Path

from fp_tools.gui_jobs import default_gui_run_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Launch the fp-tools Streamlit GUI.")
    parser.add_argument("--host", default="0.0.0.0", help="Bind address for the GUI server.")
    parser.add_argument("--port", type=int, default=None, help="Optional fixed port.")
    parser.add_argument("--run-dir", default=str(default_gui_run_dir()), help="Directory for GUI-managed runs.")
    args = parser.parse_args()

    if importlib.util.find_spec("streamlit") is None:
        raise SystemExit(
            "Streamlit is not installed in the current environment. "
            'Install the optional GUI extra to use fp-tools-gui: pip install "fp-tools-bio[gui]".'
        )

    port = args.port if args.port is not None else _find_free_port()
    _write_state(args.host, port, args.run_dir)

    display_host = _display_host(args.host)
    print(f"fp-tools GUI running at http://{display_host}:{port}")
    print(f"Run directory: {Path(args.run_dir).expanduser()}")

    app_path = Path(__file__).with_name("gui_app.py")
    env = os.environ.copy()
    env["FP_TOOLS_GUI_RUN_DIR"] = str(Path(args.run_dir).expanduser())
    command = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(app_path),
        "--server.address",
        args.host,
        "--server.port",
        str(port),
        "--server.headless",
        "true",
        "--browser.gatherUsageStats",
        "false",
    ]
    raise SystemExit(subprocess.run(command, env=env).returncode)


def _find_free_port(start: int = 8891, end: int = 8999) -> int:
    for port in range(start, end + 1):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            if sock.connect_ex(("127.0.0.1", port)) != 0:
                return port
    raise SystemExit("No free port found in the range 8891-8999.")


def _write_state(host: str, port: int, run_dir: str) -> None:
    cache_dir = _state_dir(Path(run_dir).expanduser())
    state_path = cache_dir / "gui.json"
    state_path.write_text(
        json.dumps({"host": host, "port": port, "run_dir": str(Path(run_dir).expanduser())}, indent=2),
        encoding="utf-8",
    )


def _state_dir(run_dir: Path) -> Path:
    candidates = [
        Path.home() / ".cache" / "fp-tools",
        run_dir / ".gui-state",
        Path("/tmp") / "fp-tools-gui-state",
    ]
    for path in candidates:
        try:
            path.mkdir(parents=True, exist_ok=True)
            return path
        except OSError:
            continue
    raise SystemExit("Unable to create a writable state directory for fp-tools-gui.")


def _display_host(host: str) -> str:
    if host not in {"0.0.0.0", "::"}:
        return host
    try:
        return socket.gethostbyname(socket.gethostname())
    except OSError:
        return "127.0.0.1"
