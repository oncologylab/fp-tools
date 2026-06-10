"""
Shared plotting style helpers for fp-tools PDF outputs.

This module centralizes the preferred PDF font setup, bold text defaults, and
ASCII-safe tick formatting used across plotting commands in this package.
"""

from __future__ import annotations

import matplotlib as mpl
from matplotlib import ticker

PDF_FONT_SIZE = 9
PDF_SANS_FONTS = ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"]


def apply_pdf_style() -> None:
    """Configure a consistent bold sans-serif style for matplotlib PDF plots."""

    # Prefer Helvetica/Arial when available, but keep the family generic so
    # matplotlib falls back silently to the next installed sans-serif font.
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = PDF_SANS_FONTS
    mpl.rcParams["font.size"] = PDF_FONT_SIZE
    mpl.rcParams["font.weight"] = "bold"
    mpl.rcParams["axes.labelweight"] = "bold"
    mpl.rcParams["axes.titleweight"] = "bold"
    mpl.rcParams["axes.labelsize"] = PDF_FONT_SIZE
    mpl.rcParams["axes.titlesize"] = PDF_FONT_SIZE
    mpl.rcParams["axes.unicode_minus"] = False
    mpl.rcParams["xtick.labelsize"] = PDF_FONT_SIZE
    mpl.rcParams["ytick.labelsize"] = PDF_FONT_SIZE
    mpl.rcParams["legend.fontsize"] = PDF_FONT_SIZE
    mpl.rcParams["figure.titlesize"] = PDF_FONT_SIZE
    mpl.rcParams["pdf.use14corefonts"] = True
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42


def ascii_tick_formatter(decimals: int | None = None) -> ticker.FuncFormatter:
    """Return a formatter that always uses ASCII hyphen-minus."""

    def _format(value: float, _pos: int) -> str:
        if decimals is None:
            rounded = round(value)
            if abs(value - rounded) < 1e-9:
                return f"{int(rounded)}"
            return f"{value:g}".replace("−", "-")
        return f"{value:.{decimals}f}".replace("−", "-")

    return ticker.FuncFormatter(_format)


def apply_ascii_minus_to_figure(fig) -> None:
    """Force ASCII hyphen-minus on all axes in a figure."""
    for ax in getattr(fig, "axes", []):
        if hasattr(ax, "xaxis"):
            ax.xaxis.set_major_formatter(ascii_tick_formatter())
        if hasattr(ax, "yaxis"):
            ax.yaxis.set_major_formatter(ascii_tick_formatter())
