"""Lazy command dispatch used by the standalone desktop executable.

The normal console scripts remain the public interface.  A frozen executable
cannot rely on those sibling scripts being installed, so child processes route
back through the same executable and this registry.
"""

from __future__ import annotations

import importlib
import sys
from collections.abc import Callable, Sequence


COMMAND_TARGETS: dict[str, str] = {
    "bulk-footprinting": "fp_tools.tools.bulk_footprinting:main",
    "atac-correct": "fp_tools.cli:main",
    "call-footprints": "fp_tools.cli_scorebigwig:main",
    "match-motifs": "fp_tools.tools.match_motifs:main",
    "diff-footprints": "fp_tools.tools.diff_footprints:diff_footprints_cli",
    "normalize-bigwig": "fp_tools.tools.normalize_bigwig:main",
    "plot-aggregate": "fp_tools.cli_plotaggregate:main",
    "review-multi-comparisons": "fp_tools.tools.review_multi_comparisons:main",
    "run-yaml-workflow": "fp_tools.cli_batch:main",
    "fp-tools-gui": "fp_tools.cli_gui:main",
    "fp-tools-runtime": "fp_tools.cli_runtime:main",
    "discover-motifs": "fp_tools.tools.motif_discovery:motif_discovery_plan_main",
    "summarize-motifs": "fp_tools.tools.motif_discovery:motif_report_main",
    "pseudobulk-fragments": "fp_tools.tools.pseudobulk:main",
    "find-signature-fp": "fp_tools.tools.find_signature_fp:main",
    "sc-footprinting": "fp_tools.tools.pseudobulk_footprints:main",
    # Compatibility entry points shipped by the Python package.
    "plot-motif-aggregate-grid": "fp_tools.cli_compat:plot_motif_aggregate_grid_main",
    "run-workflow": "fp_tools.cli_compat:run_workflow_main",
    "motif-discovery": "fp_tools.cli_compat:motif_discovery_main",
    "motif-summary": "fp_tools.cli_compat:motif_summary_main",
    "pseudobulk-footprints": "fp_tools.cli_compat:pseudobulk_footprints_main",
}


def resolve_command(name: str) -> Callable[[], object]:
    """Return the callable for a packaged command without importing it early."""

    try:
        target = COMMAND_TARGETS[name]
    except KeyError as exc:
        choices = ", ".join(sorted(COMMAND_TARGETS))
        raise ValueError(f"Unknown fp-tools command {name!r}. Available commands: {choices}") from exc
    module_name, function_name = target.split(":", 1)
    module = importlib.import_module(module_name)
    return getattr(module, function_name)


def dispatch_command(name: str, argv: Sequence[str] = ()) -> int:
    """Run one packaged command with console-script-compatible ``sys.argv``."""

    function = resolve_command(name)
    previous = sys.argv
    sys.argv = [name, *map(str, argv)]
    try:
        result = function()
    finally:
        sys.argv = previous
    return int(result) if isinstance(result, int) else 0
