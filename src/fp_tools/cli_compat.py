"""Temporary compatibility shims for renamed fp-tools commands."""

from __future__ import annotations

import sys


def _warn(old: str, new: str) -> None:
    print(
        f"WARNING: '{old}' is deprecated; use '{new}' instead.",
        file=sys.stderr,
    )


def plot_motif_aggregate_grid_main() -> int:
    _warn("plot-motif-aggregate-grid", "plot-aggregate")
    from fp_tools.tools.motif_aggregate_grid import main

    return main()


def run_workflow_main() -> None:
    _warn("run-workflow", "run-yaml-workflow")
    from fp_tools.cli_batch import main

    main()


def motif_discovery_main() -> int:
    _warn("motif-discovery", "discover-motifs")
    from fp_tools.tools.motif_discovery import motif_discovery_plan_main

    return motif_discovery_plan_main()


def motif_summary_main() -> int:
    _warn("motif-summary", "summarize-motifs")
    from fp_tools.tools.motif_discovery import motif_report_main

    return motif_report_main()


def pseudobulk_footprints_main() -> int:
    _warn("pseudobulk-footprints", "sc-footprinting")
    from fp_tools.tools.pseudobulk_footprints import main

    return main()
