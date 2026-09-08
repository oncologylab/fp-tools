#!/usr/bin/env python3
"""Supervise one pinned Linux research command with durable exit evidence.

Launch this supervisor in a detached session for long computations. A process
exit receipt does not replace the experiment's scientific output verifier.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fcntl
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time


SCHEMA = "fp-tools-logged-experiment-job-v1"
ENVIRONMENT_KEYS = {"PYTHONPATH", "PYTHONUNBUFFERED", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                    "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "CUDA_VISIBLE_DEVICES"}


def file_record(path):
    path = Path(path).resolve()
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path), "sha256": digest.hexdigest(), "size_bytes": path.stat().st_size}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def process_identity(pid):
    fields = Path(f"/proc/{pid}/stat").read_text().rpartition(")")[2].split()
    return {"pid": pid, "start_ticks": int(fields[19]),
            "boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip()}


def write_receipt(path, document, identity_field):
    document = dict(document)
    document[identity_field] = sha256(json.dumps(document, sort_keys=True, separators=(",", ":"),
                                               allow_nan=False).encode()).hexdigest()
    temporary = path.with_name(path.name + ".partial")
    if path.exists():
        raise ValueError(f"refusing to replace an existing execution receipt: {path}")
    with temporary.open("x") as handle:
        handle.write(json.dumps(document, sort_keys=True, indent=2, allow_nan=False) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    temporary.rename(path)
    directory_fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    return document


def read_spec(path, expected_sha256):
    if not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        raise ValueError("a complete lowercase SHA-256 is required")
    payload = path.read_bytes()
    if sha256(payload).hexdigest() != expected_sha256:
        raise ValueError("job specification checksum mismatch")
    spec = json.loads(payload)
    if (set(spec) != {"schema", "name", "command", "cwd", "environment", "outdir"}
            or spec["schema"] != SCHEMA or not isinstance(spec["name"], str) or not spec["name"]
            or not isinstance(spec["command"], list) or not spec["command"]
            or any(not isinstance(value, str) or not value or "\0" in value for value in spec["command"])
            or not isinstance(spec["environment"], dict)
            or set(spec["environment"]) - ENVIRONMENT_KEYS
            or any(not isinstance(value, str) or "\0" in value for value in spec["environment"].values())):
        raise ValueError("invalid locked research job specification")
    for name in ("cwd", "outdir"):
        if not isinstance(spec[name], str) or not Path(spec[name]).is_absolute():
            raise ValueError(f"{name} must be an absolute path")
    if not Path(spec["cwd"]).is_dir() or not Path(spec["command"][0]).is_absolute():
        raise ValueError("working directory and absolute executable are required")
    return spec, payload


def supervise(spec_path, expected_sha256):
    spec, payload = read_spec(spec_path, expected_sha256)
    outdir = Path(spec["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / ".supervisor.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if any((outdir / name).exists() for name in ("launch.json", "launch.json.partial", "exit.json",
                                                    "exit.json.partial", "stdout.log", "specification.json")):
            raise ValueError("execution evidence already exists; use a new attempt directory")
        snapshot = outdir / "specification.json"
        with snapshot.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        return run_child(spec, snapshot, outdir)


def run_child(spec, spec_path, outdir):
    child = None
    received = []

    def send_to_child(signum):
        if child is not None and child.poll() is None:
            try:
                os.killpg(child.pid, signum)
            except ProcessLookupError:
                pass

    def forward(signum, _frame):
        received.append({"signal": signal.Signals(signum).name, "received_utc": utc_now()})
        send_to_child(signum)

    previous = {signum: signal.signal(signum, forward) for signum in (signal.SIGTERM, signal.SIGINT)}
    previous[signal.SIGHUP] = signal.signal(signal.SIGHUP, signal.SIG_IGN)
    started, started_utc = time.monotonic(), utc_now()
    launch = {"schema": "fp-tools-logged-experiment-launch-v1", "name": spec["name"],
              "specification": file_record(spec_path), "supervisor_code": file_record(__file__),
              "supervisor": process_identity(os.getpid()), "started_utc": started_utc,
              "command": spec["command"], "cwd": spec["cwd"], "environment_overrides": spec["environment"]}
    try:
        with (outdir / "stdout.log").open("xb") as log:
            try:
                child = subprocess.Popen(spec["command"], cwd=spec["cwd"],
                                         env={**os.environ, **spec["environment"]}, stdin=subprocess.DEVNULL,
                                         stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            except OSError as exc:
                launch.update(child=None, spawn_error=str(exc))
                returncode = 127
            else:
                launch["child"] = process_identity(child.pid)
                for item in list(received):
                    send_to_child(signal.Signals[item["signal"]])
            write_receipt(outdir / "launch.json", launch, "launch_id")
            if child is not None:
                returncode = child.wait()
            log.flush()
            os.fsync(log.fileno())
        write_receipt(outdir / "exit.json", {
            "schema": "fp-tools-logged-experiment-exit-v1", "launch": file_record(outdir / "launch.json"),
            "returncode": returncode, "finished_utc": utc_now(), "elapsed_seconds": time.monotonic() - started,
            "signals_received_by_supervisor": list(received), "stdout": file_record(outdir / "stdout.log"),
            "scientific_outputs_verified": False,
        }, "exit_id")
        return returncode if returncode >= 0 else 128 - returncode
    except BaseException:
        # Do not leave an unrecorded computation running if receipt I/O fails.
        if child is not None and child.poll() is None:
            try:
                os.killpg(child.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                child.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
        raise
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--spec-sha256", required=True)
    args = parser.parse_args(argv)
    return supervise(args.spec.resolve(), args.spec_sha256)


if __name__ == "__main__":
    raise SystemExit(main())
