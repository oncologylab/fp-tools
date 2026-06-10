"""Shared Matplotlib style for paper figures: bold Arial/Helvetica, size 11.

Arial and Helvetica are preferred; Liberation Sans (metric-compatible with Arial)
and Nimbus Sans (a Helvetica clone) are the free fallbacks bundled on most Linux
systems. Every text element is bold. Import and call :func:`apply_style` at the top
of a figure script, optionally lowering the base size if a dense panel needs it.
"""

from __future__ import annotations

import matplotlib as mpl

SANS = ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]


def apply_style(base_size: int = 11) -> None:
    """Set a bold sans-serif (Arial/Helvetica) style at the given base font size."""

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": SANS,
            "font.size": base_size,
            "font.weight": "bold",
            "axes.titlesize": base_size,
            "axes.titleweight": "bold",
            "axes.labelsize": base_size,
            "axes.labelweight": "bold",
            "axes.linewidth": 1.2,
            "xtick.labelsize": base_size,
            "ytick.labelsize": base_size,
            "legend.fontsize": max(base_size - 1, 7),
            "figure.titlesize": base_size,
            "figure.titleweight": "bold",
            "savefig.dpi": 300,
            "pdf.fonttype": 42,  # embed TrueType so text stays editable/searchable
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def bold_all_text(ax) -> None:
    """Force every text artist on an Axes (incl. tick labels, legend) to bold."""

    items = [ax.title, ax.xaxis.label, ax.yaxis.label]
    items += ax.get_xticklabels() + ax.get_yticklabels()
    legend = ax.get_legend()
    if legend is not None:
        items += legend.get_texts()
    for item in items:
        item.set_fontweight("bold")
