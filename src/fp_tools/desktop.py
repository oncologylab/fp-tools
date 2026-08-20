"""Entry point for the self-contained fp-tools desktop executable."""

from __future__ import annotations

import runpy
import sys
from collections.abc import Sequence

from fp_tools.command_registry import dispatch_command

INTERNAL_COMMAND_FLAG = "--fp-tools-internal-command"
INTERNAL_MATCH_BEDS_FLAG = "--fp-tools-internal-match-beds"
INTERNAL_PYTHON_SCRIPT_FLAG = "--fp-tools-internal-python-script"
INTERNAL_LIST_EXAMPLES_FLAG = "--fp-tools-internal-list-gui-examples"
INTERNAL_SMOKE_HELPS_FLAG = "--fp-tools-internal-smoke-command-helps"


def main(argv: Sequence[str] | None = None) -> int:
    """Launch the GUI, or dispatch a child command inside a frozen bundle."""

    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments and arguments[0] == INTERNAL_SMOKE_HELPS_FLAG:
        if len(arguments) != 1:
            raise SystemExit(f"{INTERNAL_SMOKE_HELPS_FLAG} does not accept arguments")
        from fp_tools.command_registry import COMMAND_TARGETS

        for command in sorted(COMMAND_TARGETS):
            try:
                dispatch_command(command, ["--help"])
            except SystemExit as exc:
                if exc.code not in (None, 0):
                    raise
        print(f"Validated {len(COMMAND_TARGETS)} desktop commands")
        return 0
    if arguments and arguments[0] == INTERNAL_LIST_EXAMPLES_FLAG:
        if len(arguments) != 1:
            raise SystemExit(f"{INTERNAL_LIST_EXAMPLES_FLAG} does not accept arguments")
        from importlib import resources

        root = resources.files("fp_tools.resources.gui_configs")
        for path in sorted(root.iterdir(), key=lambda value: value.name):
            if path.is_file() and path.name.endswith(".yml"):
                print(path.name)
        return 0
    if arguments and arguments[0] == INTERNAL_COMMAND_FLAG:
        if len(arguments) < 2:
            raise SystemExit(f"{INTERNAL_COMMAND_FLAG} requires a command name")
        return dispatch_command(arguments[1], arguments[2:])
    if arguments and arguments[0] == INTERNAL_MATCH_BEDS_FLAG:
        if len(arguments) != 2:
            raise SystemExit(f"{INTERNAL_MATCH_BEDS_FLAG} requires one payload path")
        from fp_tools.tools.diff_footprints import _materialize_match_motif_beds_payload

        _materialize_match_motif_beds_payload(arguments[1])
        return 0
    if arguments and arguments[0] == INTERNAL_PYTHON_SCRIPT_FLAG:
        if len(arguments) < 2:
            raise SystemExit(f"{INTERNAL_PYTHON_SCRIPT_FLAG} requires a script path")
        previous = sys.argv
        sys.argv = arguments[1:]
        try:
            runpy.run_path(arguments[1], run_name="__main__")
        finally:
            sys.argv = previous
        return 0

    from fp_tools.cli_gui import main as gui_main

    gui_main(arguments)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
