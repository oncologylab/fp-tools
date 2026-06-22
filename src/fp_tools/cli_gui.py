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
from urllib.error import URLError
from urllib.request import Request, urlopen

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
    access_urls = _access_urls(args.host, port)
    _write_state(args.host, port, args.run_dir, access_urls)

    print(f"fp-tools GUI bind address: {args.host}:{port}")
    print("Access URLs:")
    for url in access_urls:
        print(f"  {url}")
    if args.host in {"0.0.0.0", "::"}:
        print(f"External access requires TCP port {port} to be open in the host firewall or cloud security group.")
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


def _write_state(host: str, port: int, run_dir: str, access_urls: list[str]) -> None:
    cache_dir = _state_dir(Path(run_dir).expanduser())
    state_path = cache_dir / "gui.json"
    state_path.write_text(
        json.dumps(
            {
                "host": host,
                "port": port,
                "run_dir": str(Path(run_dir).expanduser()),
                "access_urls": access_urls,
            },
            indent=2,
        ),
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


def _access_urls(host: str, port: int) -> list[str]:
    hosts = [host]
    if host in {"0.0.0.0", "::"}:
        hosts = ["127.0.0.1", *_candidate_public_hosts()]
    urls: list[str] = []
    for candidate in hosts:
        if not candidate or candidate in {"0.0.0.0", "::"}:
            continue
        url = f"http://{candidate}:{port}"
        if url not in urls:
            urls.append(url)
    return urls or [f"http://127.0.0.1:{port}"]


def _candidate_public_hosts() -> list[str]:
    candidates: list[str] = []
    for candidate in _hostname_ips() + [_outbound_ip()] + _metadata_public_ips():
        if not candidate:
            continue
        if candidate.startswith("127.") or candidate == "0.0.0.0":
            continue
        if candidate not in candidates:
            candidates.append(candidate)
    return candidates


def _hostname_ips() -> list[str]:
    names = {socket.gethostname()}
    try:
        names.add(socket.getfqdn())
    except OSError:
        pass
    ips: list[str] = []
    for name in names:
        try:
            for family, _, _, _, sockaddr in socket.getaddrinfo(name, None, socket.AF_INET):
                if family == socket.AF_INET and sockaddr[0] not in ips:
                    ips.append(sockaddr[0])
        except OSError:
            continue
    return ips


def _outbound_ip() -> str:
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as sock:
            sock.connect(("8.8.8.8", 80))
            return str(sock.getsockname()[0])
    except OSError:
        return ""


def _metadata_public_ips() -> list[str]:
    urls = [
        ("http://169.254.169.254/latest/meta-data/public-ipv4", 2.0),
        ("http://169.254.169.254/metadata/v1/interfaces/public/0/ipv4/address", 0.35),
        ("http://169.254.169.254/metadata/instance/network/interface/0/ipv4/ipAddress/0/publicIpAddress?api-version=2021-02-01&format=text", 0.35),
    ]
    ips: list[str] = []
    for url, timeout in urls:
        request = Request(url, headers={"Metadata": "true"})
        try:
            with urlopen(request, timeout=timeout) as response:
                value = response.read(128).decode("utf-8", errors="ignore").strip()
        except (OSError, URLError, TimeoutError):
            continue
        if _is_ipv4(value) and value not in ips:
            ips.append(value)
    return ips


def _is_ipv4(value: str) -> bool:
    try:
        socket.inet_aton(value)
    except OSError:
        return False
    return value.count(".") == 3
