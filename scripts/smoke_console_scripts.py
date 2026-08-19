#!/usr/bin/env python3
"""Run ``--help`` for every console script declared in pyproject.toml."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import shutil
import sys
import sysconfig
import tomllib


ROOT = Path(__file__).resolve().parents[1]


def console_scripts() -> list[str]:
    data = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    scripts = data["project"]["scripts"]
    poetry_scripts = data["tool"]["poetry"]["scripts"]
    if scripts != poetry_scripts:
        raise RuntimeError("[project.scripts] and [tool.poetry.scripts] differ")
    return list(scripts)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bin-dir",
        type=Path,
        default=Path(sysconfig.get_path("scripts")),
        help="Directory containing installed console scripts.",
    )
    args = parser.parse_args()

    failures = []
    env = os.environ.copy()
    for command in console_scripts():
        executable = shutil.which(command, path=str(args.bin_dir))
        if executable is None:
            failures.append((command, 127, f"not found in {args.bin_dir}"))
            continue
        result = subprocess.run(
            [executable, "--help"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
            check=False,
        )
        if result.returncode:
            failures.append((command, result.returncode, result.stderr.strip()))
        else:
            print(f"OK\t{command}")
    if failures:
        for command, returncode, stderr in failures:
            print(f"FAIL\t{command}\t{returncode}\t{stderr}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
