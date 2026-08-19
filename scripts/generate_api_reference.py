#!/usr/bin/env python3
"""Generate the public command reference from installed console scripts."""

from __future__ import annotations

from pathlib import Path
import os
import re
import subprocess
import sys
import tempfile
import tomllib


ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "pyproject.toml"
OUTPUT = ROOT / "docs" / "api.md"
GUIDES_DIR = ROOT / "docs" / "get-started" / "commands"


def _guide_content(command: str) -> tuple[str, str]:
    """Return the command purpose and reusable guide content for the API page."""
    path = GUIDES_DIR / f"{command}.md"
    if not path.is_file():
        raise SystemExit(f"Missing command guide: {path}")

    source = path.read_text(encoding="utf-8").strip()
    if source.startswith("---\n"):
        _, separator, source = source[4:].partition("\n---\n")
        if not separator:
            raise SystemExit(f"Command guide has unterminated front matter: {path}")
        source = source.strip()
    lines = source.splitlines()
    if not lines or not lines[0].startswith("# "):
        raise SystemExit(f"Command guide must start with an H1: {path}")

    blocks = "\n".join(lines[1:]).strip().split("\n\n")
    blocks = [block for block in blocks if "../../api.md" not in block]
    if not blocks:
        raise SystemExit(f"Command guide has no reusable content: {path}")

    purpose = " ".join(blocks[0].split())
    body = "\n\n".join(blocks)
    for heading in ("Example command", "Primary inputs", "Main outputs"):
        body = body.replace(f"## {heading}", f"**{heading}**")
    def api_link(match: re.Match[str]) -> str:
        label, target = match.groups()
        if target.startswith(("http://", "https://", "mailto:", "#")):
            return match.group(0)
        path_text, marker, fragment = target.partition("#")
        resolved = (path.parent / path_text).resolve()
        if resolved.parent == GUIDES_DIR.resolve() and resolved.suffix == ".md":
            rewritten = f"#{resolved.stem}"
        else:
            try:
                rewritten = resolved.relative_to((ROOT / "docs").resolve()).as_posix()
            except ValueError:
                return match.group(0)
            if marker:
                rewritten += f"#{fragment}"
        return f"[{label}]({rewritten})"

    body = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", api_link, body)
    return purpose, body


def main() -> int:
    config = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))
    commands = config["tool"]["fp-tools"]["public-console-scripts"]
    scripts = config["project"]["scripts"]
    missing = [command for command in commands if command not in scripts]
    if missing:
        raise SystemExit(f"Public commands missing from [project.scripts]: {', '.join(missing)}")

    guides = {command: _guide_content(command) for command in commands}

    rows = [
        "---",
        "hide:",
        "  - navigation",
        "---",
        "",
        "# API Reference",
        "",
        "Direct CLI commands are the primary interface. Each reference includes a method summary, practical example, primary inputs, outputs, and the complete command options.",
        "",
        "| Command | Purpose |",
        "| --- | --- |",
    ]
    rows.extend(
        f"| [`{command}`](#{command}) | {guides[command][0]} |"
        for command in commands
    )
    rows.append("")
    bin_dir = Path(sys.executable).parent
    with tempfile.TemporaryDirectory(prefix="fp-tools-api-reference-") as cache:
        environment = os.environ.copy()
        environment["MPLCONFIGDIR"] = str(Path(cache) / "matplotlib")
        environment["XDG_CACHE_HOME"] = str(Path(cache) / "xdg")
        for command in commands:
            _, guide = guides[command]
            executable = bin_dir / command
            result = subprocess.run(
                [str(executable), "--help"],
                cwd=ROOT,
                env=environment,
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
                    '<div class="fp-api-card" markdown="1">',
                    "",
                    guide,
                    "",
                    "**Complete options**",
                    "",
                    "```text",
                    help_text,
                    "```",
                    "",
                    "</div>",
                    "",
                ]
            )
    output = "\n".join(rows)
    lowered = output.lower()
    for forbidden in ("nutrient", "legacy", "/home/", "169.254.169.254"):
        if forbidden in lowered:
            raise SystemExit(f"Generated API reference contains forbidden term: {forbidden}")
    OUTPUT.write_text(output, encoding="utf-8")
    print(f"Wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
