"""Resolve fp-tools child commands in Python and frozen installations."""

from __future__ import annotations

import shutil
import sys
from collections.abc import Sequence
from pathlib import Path

from fp_tools.command_registry import COMMAND_TARGETS


def is_frozen() -> bool:
    return bool(getattr(sys, "frozen", False))


def fp_tools_subprocess_command(name: str, argv: Sequence[str] = ()) -> list[str]:
    """Build a child-process command that works inside a PyInstaller bundle."""

    arguments = [str(value) for value in argv]
    if is_frozen():
        from fp_tools.desktop import INTERNAL_COMMAND_FLAG

        return [sys.executable, INTERNAL_COMMAND_FLAG, name, *arguments]
    local = Path(sys.executable).parent / name
    executable = str(local) if local.exists() else (shutil.which(name) or name)
    return [executable, *arguments]


def python_script_subprocess_command(script: str | Path, argv: Sequence[str] = ()) -> list[str]:
    """Run a Python helper script with the bundled or installed interpreter."""

    arguments = [str(value) for value in argv]
    if is_frozen():
        from fp_tools.desktop import INTERNAL_PYTHON_SCRIPT_FLAG

        return [sys.executable, INTERNAL_PYTHON_SCRIPT_FLAG, str(script), *arguments]
    return [sys.executable, str(script), *arguments]


def resolve_fp_tools_subprocess(command: Sequence[str]) -> list[str]:
    """Resolve the first element when it names a packaged fp-tools command."""

    resolved = [str(value) for value in command]
    if not resolved:
        raise ValueError("Cannot resolve an empty command")
    if resolved[0] in COMMAND_TARGETS:
        return fp_tools_subprocess_command(resolved[0], resolved[1:])
    return resolved
