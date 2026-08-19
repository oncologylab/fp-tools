#!/usr/bin/env python3
"""Generate the public command reference from installed console scripts."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tomllib


ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "pyproject.toml"
OUTPUT = ROOT / "docs" / "api.md"


def main() -> int:
    config = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))
    commands = config["tool"]["fp-tools"]["public-console-scripts"]
    scripts = config["project"]["scripts"]
    missing = [command for command in commands if command not in scripts]
    if missing:
        raise SystemExit(f"Public commands missing from [project.scripts]: {', '.join(missing)}")

    rows = [
        "---",
        "hide:",
        "  - navigation",
        "---",
        "",
        "# API Reference",
        "",
        "Direct CLI commands are the primary interface. The GUI and YAML runner call these same commands.",
        "",
        "| Command | Reference |",
        "| --- | --- |",
    ]
    rows.extend(f"| `{command}` | [Options](#{command}) |" for command in commands)
    rows.append("")
    bin_dir = Path(sys.executable).parent
    for command in commands:
        executable = bin_dir / command
        result = subprocess.run(
            [str(executable), "--help"],
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        if result.returncode:
            raise SystemExit(f"{command} --help failed:\n{result.stdout}")
        help_text = "\n".join(line.rstrip() for line in result.stdout.strip().splitlines())
        rows.extend(
            [
                f"## `{command}`",
                "",
                "```text",
                help_text,
                "```",
                "",
            ]
        )
    output = "\n".join(rows)
    lowered = output.lower()
    for forbidden in ("nutrient", "legacy"):
        if forbidden in lowered:
            raise SystemExit(f"Generated API reference contains forbidden term: {forbidden}")
    OUTPUT.write_text(output, encoding="utf-8")
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
