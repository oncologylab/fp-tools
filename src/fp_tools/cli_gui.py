"""Launcher for the command-backed fp-tools Streamlit GUI."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import socket
import subprocess
import sys
import tempfile
import threading
import time
import webbrowser
from pathlib import Path

from fp_tools.gui_jobs import default_gui_run_dir


STREAMLIT_LIGHT_THEME = {
    "theme.base": "light",
    "theme.primaryColor": "#2563eb",
    "theme.backgroundColor": "#f3f6fa",
    "theme.secondaryBackgroundColor": "#ffffff",
    "theme.textColor": "#111827",
    "theme.font": "sans serif",
}


def _streamlit_theme_cli_args() -> list[str]:
    """Return explicit CLI flags so user-level settings cannot enable dark mode."""

    return [part for key, value in STREAMLIT_LIGHT_THEME.items() for part in (f"--{key}", value)]


def _streamlit_theme_bootstrap_options() -> dict[str, str]:
    """Return Streamlit bootstrap keys for the frozen desktop server."""

    return {key.replace(".", "_"): value for key, value in STREAMLIT_LIGHT_THEME.items()}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Launch the fp-tools browser interface.")
    parser.add_argument("--host", default="127.0.0.1", help="Bind address (default: 127.0.0.1).")
    parser.add_argument("--port", type=int, default=None, help="Optional fixed port (default: first free port from 8891).")
    parser.add_argument("--run-dir", default=str(default_gui_run_dir()), help="Directory for GUI-managed runs.")
    parser.add_argument("--no-browser", action="store_true", help="Do not open a local browser automatically.")
    args = parser.parse_args(argv)

    if importlib.util.find_spec("streamlit") is None:
        raise SystemExit("Streamlit is missing. Reinstall with: python -m pip install --upgrade fp-tools-bio")

    port = args.port if args.port is not None else _find_free_port()
    access_urls = _access_urls(args.host, port)
    _write_state(args.host, port, args.run_dir, access_urls)
    for message in _startup_messages(args.host, port, Path(args.run_dir).expanduser()):
        print(message, flush=True)

    local_url = f"http://127.0.0.1:{port}"
    if not args.no_browser and args.host in {"127.0.0.1", "localhost", "::1"}:
        threading.Thread(target=_open_browser_when_ready, args=(port, local_url), daemon=True).start()

    app_path = Path(__file__).with_name("gui_app.py")
    env = os.environ.copy()
    env["FP_TOOLS_GUI_RUN_DIR"] = str(Path(args.run_dir).expanduser())
    try:
        if getattr(sys, "frozen", False):
            os.environ["STREAMLIT_GLOBAL_DEVELOPMENT_MODE"] = "false"
            from streamlit.web import bootstrap

            os.environ.update(env)
            flag_options = {
                "global_developmentMode": False,
                "server_address": args.host,
                "server_port": port,
                "server_headless": True,
                "browser_gatherUsageStats": False,
                **_streamlit_theme_bootstrap_options(),
            }
            bootstrap.load_config_options(flag_options)
            bootstrap.run(
                str(app_path),
                False,
                [],
                flag_options,
            )
            return_code = 0
        else:
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
                *_streamlit_theme_cli_args(),
            ]
            return_code = subprocess.run(command, env=env).returncode
    except KeyboardInterrupt:
        return_code = 130
    raise SystemExit(return_code)


def _startup_messages(host: str, port: int, run_dir: Path) -> list[str]:
    messages = [
        "fp-tools GUI",
        f"Open locally: http://127.0.0.1:{port}",
        f"Run directory: {run_dir}",
        "Stop the server with Ctrl+C.",
    ]
    if host in {"127.0.0.1", "localhost", "::1"}:
        messages.extend(
            [
                "Remote Linux server: keep this process running, then run the following on your computer:",
                f"  ssh -N -L {port}:127.0.0.1:{port} USER@SERVER",
                f"Then open http://127.0.0.1:{port} on your computer.",
            ]
        )
    elif host in {"0.0.0.0", "::"}:
        messages.extend(
            [
                "Network mode is enabled. Open http://SERVER_IP:%d from an allowed network." % port,
                "Warning: fp-tools does not add authentication; use a firewall, VPN, or SSH tunnel.",
            ]
        )
    else:
        messages.append(f"Open from a reachable client: http://{host}:{port}")
    return messages


def _open_browser_when_ready(port: int, url: str, timeout: float = 20.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=0.25):
                webbrowser.open(url)
                return
        except OSError:
            time.sleep(0.15)


def _find_free_port(start: int = 8891, end: int = 8999) -> int:
    for port in range(start, end + 1):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            if sock.connect_ex(("127.0.0.1", port)) != 0:
                return port
    raise SystemExit("No free port found in the range 8891-8999.")


def _write_state(host: str, port: int, run_dir: str, access_urls: list[str]) -> None:
    cache_dir = _state_dir(Path(run_dir).expanduser())
    (cache_dir / "gui.json").write_text(
        json.dumps({"host": host, "port": port, "run_dir": str(Path(run_dir).expanduser()), "access_urls": access_urls}, indent=2),
        encoding="utf-8",
    )


def _state_dir(run_dir: Path) -> Path:
    candidates = [
        Path(os.environ.get("LOCALAPPDATA", "")) / "fp-tools" if os.name == "nt" else Path.home() / ".cache" / "fp-tools",
        run_dir / ".gui-state",
        Path(tempfile.gettempdir()) / "fp-tools-gui-state",
    ]
    for path in candidates:
        try:
            path.mkdir(parents=True, exist_ok=True)
            return path
        except OSError:
            continue
    raise SystemExit("Unable to create a writable state directory for fp-tools-gui.")


def _access_urls(host: str, port: int) -> list[str]:
    if host in {"0.0.0.0", "::"}:
        return [f"http://127.0.0.1:{port}", f"http://SERVER_IP:{port}"]
    return [f"http://{host}:{port}"]


if __name__ == "__main__":
    main()
