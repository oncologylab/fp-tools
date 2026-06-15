#!/usr/bin/env python
"""Run a command and record engineering benchmark metadata.

This helper is intentionally tool-agnostic. It can wrap fp-tools, TOBIAS, HINT,
or any other command used in future head-to-head benchmark runs.
"""

from __future__ import annotations

import argparse
import csv
import os
import resource
import shlex
import subprocess
import time
from pathlib import Path


def run_benchmark(command: list[str], output: Path, label: str, cores: int | None = None) -> dict[str, object]:
    start = time.perf_counter()
    completed = subprocess.run(command, check=False)
    elapsed = time.perf_counter() - start
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    row = {
        "label": label,
        "command": " ".join(shlex.quote(part) for part in command),
        "exit_code": completed.returncode,
        "wall_seconds": round(elapsed, 3),
        "peak_rss_kb": int(usage.ru_maxrss),
        "cores": cores if cores is not None else "",
        "cwd": os.getcwd(),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    exists = output.exists()
    with output.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerow(row)
    return row


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True, help="Benchmark row label.")
    parser.add_argument("--out", type=Path, required=True, help="Output TSV.")
    parser.add_argument("--cores", type=int, help="Cores assigned to the wrapped command.")
    parser.add_argument("command", nargs=argparse.REMAINDER, help="Command to run after --.")
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("provide a command to benchmark after --")
    row = run_benchmark(command, args.out, args.label, cores=args.cores)
    return int(row["exit_code"])


if __name__ == "__main__":
    raise SystemExit(main())
