"""Native window for the frozen Windows and macOS fp-tools applications."""

from __future__ import annotations

import argparse
import os
import signal
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import IO, Sequence

from fp_tools.cli_gui import _find_free_port, _state_dir
from fp_tools.gui_jobs import default_gui_run_dir


class DesktopLaunchError(RuntimeError):
    """Raised when the bundled GUI cannot be started safely."""


def launch_native_gui(
    argv: Sequence[str] | None = None,
    *,
    auto_close: bool = False,
) -> int:
    """Run the command-backed Streamlit GUI in a native Qt window."""

    args = _parse_desktop_args(argv)
    run_dir = Path(args.run_dir).expanduser().resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    port = args.port if args.port is not None else _find_free_port()
    url = f"http://127.0.0.1:{port}"
    state_dir = _state_dir(run_dir)
    log_path = state_dir / "desktop-server.log"

    with log_path.open("w", encoding="utf-8") as log_handle:
        process = _start_server(port, run_dir, log_handle)
        try:
            _wait_for_server(process, url, log_path)
            return _run_window(url, state_dir, auto_close=auto_close)
        finally:
            _stop_process_tree(process)


def show_desktop_error(message: str) -> int:
    """Display a launch failure in a native dialog when Qt is available."""

    try:
        from PySide6.QtWidgets import QApplication, QMessageBox
    except ImportError:
        print(f"fp-tools could not start: {message}", file=sys.stderr)
        return 1

    application = QApplication.instance() or QApplication(["fp-tools"])
    QMessageBox.critical(None, "fp-tools", message)
    if QApplication.instance() is application:
        application.processEvents()
    return 1


def _parse_desktop_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Launch the fp-tools desktop application."
    )
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--port",
        type=int,
        default=None,
        help="Optional local port (default: first free port from 8891).",
    )
    parser.add_argument(
        "--run-dir",
        default=str(default_gui_run_dir()),
        help="Directory for GUI-managed runs.",
    )
    parser.add_argument("--no-browser", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(list(argv or []))
    if args.host not in {"127.0.0.1", "localhost", "::1"}:
        raise DesktopLaunchError(
            "The desktop application only accepts a local connection. "
            "Use the Python fp-tools-gui command for remote or network access."
        )
    return args


def _server_command(port: int, run_dir: Path) -> list[str]:
    return [
        sys.executable,
        "--fp-tools-internal-gui-server",
        "--host",
        "127.0.0.1",
        "--port",
        str(port),
        "--run-dir",
        str(run_dir),
        "--no-browser",
    ]


def _start_server(
    port: int, run_dir: Path, log_handle: IO[str]
) -> subprocess.Popen[str]:
    kwargs: dict[str, object] = {
        "cwd": str(run_dir),
        "env": {
            **os.environ,
            "STREAMLIT_BROWSER_GATHER_USAGE_STATS": "false",
        },
        "stdout": log_handle,
        "stderr": subprocess.STDOUT,
        "text": True,
    }
    if os.name == "nt":
        kwargs["creationflags"] = getattr(
            subprocess, "CREATE_NEW_PROCESS_GROUP", 0
        ) | getattr(subprocess, "CREATE_NO_WINDOW", 0)
    else:
        kwargs["start_new_session"] = True
    return subprocess.Popen(_server_command(port, run_dir), **kwargs)


def _wait_for_server(
    process: subprocess.Popen[str],
    url: str,
    log_path: Path,
    *,
    timeout: float = 120.0,
) -> None:
    deadline = time.monotonic() + timeout
    last_error: Exception | None = None
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise DesktopLaunchError(
                "The local GUI server stopped during startup. "
                f"Details were written to {log_path}."
            )
        try:
            with urllib.request.urlopen(f"{url}/_stcore/health", timeout=2) as response:
                if response.status == 200:
                    return
        except Exception as exc:  # normal startup polling
            last_error = exc
            time.sleep(0.25)
    raise DesktopLaunchError(
        "The local GUI server did not become ready within two minutes. "
        f"Details were written to {log_path}. Last error: {last_error}"
    )


def _run_window(url: str, state_dir: Path, *, auto_close: bool) -> int:
    if auto_close:
        os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--disable-gpu")
        if sys.platform != "darwin":
            os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    try:
        from PySide6.QtCore import QTimer, QUrl
        from PySide6.QtGui import QIcon
        from PySide6.QtWebEngineCore import QWebEngineProfile
        from PySide6.QtWebEngineWidgets import QWebEngineView
        from PySide6.QtWidgets import QApplication, QMainWindow
    except ImportError as exc:
        raise DesktopLaunchError(
            "The native window component is missing from this desktop bundle."
        ) from exc

    application = QApplication.instance() or QApplication(["fp-tools"])
    application.setApplicationDisplayName("fp-tools")
    application.setApplicationName("fp-tools")
    icon_path = _bundled_icon_path()
    if icon_path is not None:
        application.setWindowIcon(QIcon(str(icon_path)))

    profile = QWebEngineProfile.defaultProfile()
    profile.setCachePath(str(state_dir / "web-cache"))
    profile.setPersistentStoragePath(str(state_dir / "web-storage"))

    window = QMainWindow()
    window.setWindowTitle("fp-tools")
    if icon_path is not None:
        window.setWindowIcon(QIcon(str(icon_path)))
    view = QWebEngineView(window)
    window.setCentralWidget(view)
    window.setMinimumSize(900, 640)
    geometry = application.primaryScreen().availableGeometry()
    window.resize(
        min(1440, int(geometry.width() * 0.92)), min(960, int(geometry.height() * 0.92))
    )

    load_failed = [False]

    def loaded(ok: bool) -> None:
        if not ok:
            load_failed[0] = True
            application.quit()
        elif auto_close:
            QTimer.singleShot(750, application.quit)

    view.loadFinished.connect(loaded)
    view.setUrl(QUrl(url))
    window.show()
    if auto_close:
        QTimer.singleShot(45_000, application.quit)
    return_code = application.exec()
    if load_failed[0]:
        raise DesktopLaunchError(
            "The fp-tools interface could not be rendered in the desktop window."
        )
    return int(return_code)


def _bundled_icon_path() -> Path | None:
    root = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parents[2]))
    candidates = [
        root / "fp_tools" / "resources" / "fp_tools_logo_icon_1024.png",
        root / "docs" / "assets" / "fp_tools_logo_icon_1024.png",
    ]
    return next((path for path in candidates if path.is_file()), None)


def _stop_process_tree(process: subprocess.Popen[str]) -> None:
    if process.poll() is not None:
        return
    if os.name == "nt":
        subprocess.run(
            ["taskkill", "/PID", str(process.pid), "/T", "/F"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
    try:
        process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        if os.name == "nt":
            process.kill()
        else:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        process.wait(timeout=15)
