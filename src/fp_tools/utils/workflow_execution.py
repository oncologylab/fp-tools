"""Stream workflow diagnostics while preserving per-stage log files."""

from __future__ import annotations

import codecs
import os
import shlex
import subprocess
import sys
import threading
import time
from collections.abc import Mapping, Sequence
from contextlib import ExitStack
from pathlib import Path
from typing import IO, Any


_OUTPUT_LOCK = threading.Lock()


def run_logged(
    command: Sequence[str],
    *,
    stdout_log: str | Path | None = None,
    stderr_log: str | Path | None = None,
    env: Mapping[str, str] | None = None,
    cwd: str | Path | None = None,
    stdout_target: IO[Any] | None = None,
    input_command: Sequence[str] | None = None,
    capture_stdout: bool = False,
    append: bool = False,
    label: str | None = None,
) -> subprocess.CompletedProcess:
    """Run a command (or two-command pipe), teeing diagnostics in bounded chunks.

    Data redirected to ``stdout_target`` or the next process is never echoed.
    ``capture_stdout`` is reserved for small machine-readable command results.
    Log files retain the original bytes; decoding applies only to console text.
    """
    command = [str(part) for part in command]
    environment = dict(os.environ if env is None else env)
    environment["PYTHONUNBUFFERED"] = "1"
    environment.setdefault("PYTHONIOENCODING", "utf-8")
    stage = label or Path(command[0]).name
    log_paths = list(dict.fromkeys(str(path) for path in (stdout_log, stderr_log) if path is not None))
    log_hint = "; logs: " + ", ".join(log_paths) if log_paths else ""
    display = shlex.join(command)
    if input_command is not None:
        display = shlex.join([str(part) for part in input_command]) + " | " + display
    print(f"[start] {stage}: {display}", file=sys.stderr, flush=True)
    started = time.monotonic()
    processes = []
    readers = []
    errors = []
    captured = bytearray()
    lock = _OUTPUT_LOCK

    with ExitStack() as stack:
        handles = {}

        def open_log(path):
            if path is None:
                return None
            path = Path(path).resolve()
            if path not in handles:
                path.parent.mkdir(parents=True, exist_ok=True)
                handles[path] = stack.enter_context(path.open("ab" if append else "wb"))
            return handles[path]

        out_file, err_file = open_log(stdout_log), open_log(stderr_log)
        if append:
            for handle in handles.values():
                handle.write(("$ " + display + "\n").encode("utf-8"))
                handle.flush()

        def drain(pipe, log, console, collect=False):
            decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")
            at_line_start = True

            def echo(text):
                nonlocal at_line_start
                if append and label and text:
                    prefix = f"[{Path(str(label)).stem}] "
                    parts = text.splitlines(keepends=True)
                    text = ""
                    for part in parts:
                        text += (prefix if at_line_start else "") + part
                        at_line_start = part.endswith(("\n", "\r"))
                try:
                    console.write(text)
                except UnicodeEncodeError:
                    encoding = getattr(console, "encoding", None) or "utf-8"
                    console.write(text.encode(encoding, errors="replace").decode(encoding))
                console.flush()
            try:
                while True:
                    block = pipe.read1(65536)
                    if not block:
                        break
                    if collect:
                        captured.extend(block)
                    try:
                        with lock:
                            if log is not None:
                                log.write(block)
                                log.flush()
                            if console is not None:
                                echo(decoder.decode(block))
                    except (OSError, ValueError, UnicodeError) as exc:
                        errors.append(exc)
                        # Keep draining so a reporting failure cannot deadlock a child.
                        log = console = None
                if console is not None:
                    with lock:
                        echo(decoder.decode(b"", final=True))
            except (OSError, ValueError) as exc:
                errors.append(exc)
            finally:
                pipe.close()

        def reader(pipe, log, console, collect=False):
            thread = threading.Thread(target=drain, args=(pipe, log, console, collect), daemon=True)
            readers.append(thread)
            thread.start()

        try:
            # The launcher owns the process group. Inherit it so terminal
            # interrupts and desktop group termination reach nested stages.
            left = None
            if input_command is not None:
                left = subprocess.Popen(
                    input_command, cwd=cwd, env=environment, stdin=subprocess.DEVNULL,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                )
                processes.append(left)
                reader(left.stderr, err_file, sys.stderr)
            process = subprocess.Popen(
                command, cwd=cwd, env=environment,
                stdin=left.stdout if left is not None else None,
                stdout=stdout_target if stdout_target is not None else subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            processes.append(process)
            if left is not None:
                left.stdout.close()
            if stdout_target is None:
                reader(process.stdout, None if capture_stdout else out_file,
                       None if capture_stdout else sys.stdout, capture_stdout)
            reader(process.stderr, err_file, sys.stderr)
            code = process.wait()
            if left is not None:
                left_code = left.wait()
                code = left_code or code
            for thread in readers:
                thread.join()
            if errors:
                raise OSError(f"Could not stream workflow output: {errors[0]}") from errors[0]
        except BaseException:
            for child in reversed(processes):
                if child.poll() is None:
                    if os.name == "nt":
                        subprocess.run(["taskkill", "/PID", str(child.pid), "/T", "/F"],
                                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
                    else:
                        child.terminate()
            for child in reversed(processes):
                try:
                    child.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    child.kill()
                    child.wait()
            if input_command is not None and processes and processes[0].stdout:
                processes[0].stdout.close()
            for thread in readers:
                thread.join(timeout=5)
            print(f"[failed] {stage}: interrupted or could not execute{log_hint}", file=sys.stderr, flush=True)
            raise

    state = "done" if code == 0 else "failed"
    print(f"[{state}] {stage}: exit {code}, {time.monotonic() - started:.1f}s{log_hint if code else ''}", file=sys.stderr, flush=True)
    return subprocess.CompletedProcess(command, code, stdout=bytes(captured) if capture_stdout else None)
