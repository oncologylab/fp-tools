#!/usr/bin/env python
"""Plot Fig5 engineering-performance and GUI-support panels."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.text as mtext
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = REPO_ROOT / "manuscript/figures/Fig5_engineering_gui_source.tsv"
DEFAULT_OUT = REPO_ROOT / "manuscript/figures/Fig5.svg"
COLORS = {"baseline": "#9CA3AF", "current": "#2563EB"}


def apply_svg_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "font.weight": "bold",
            "axes.titlesize": 9,
            "axes.titleweight": "bold",
            "axes.labelsize": 9,
            "axes.labelweight": "bold",
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "svg.fonttype": "none",
        }
    )


def force_arial_bold(fig: plt.Figure) -> None:
    for text in fig.findobj(match=mtext.Text):
        text.set_fontfamily("Arial")
        text.set_fontweight("bold")


def load_source(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    required = {
        "panel",
        "metric",
        "task",
        "baseline_label",
        "current_label",
        "baseline_value",
        "current_value",
        "unit",
        "direction",
    }
    missing = required.difference(table.columns)
    if missing:
        raise SystemExit(f"Fig5 source table is missing required columns: {', '.join(sorted(missing))}")
    return table


def plot_metric_bars(ax: plt.Axes, table: pd.DataFrame, panel: str, title: str) -> None:
    subset = table[table["panel"] == panel].copy()
    if subset.empty:
        raise SystemExit(f"No source rows found for panel {panel}")
    baseline_label = str(subset["baseline_label"].iloc[0])
    current_label = str(subset["current_label"].iloc[0])
    unit = str(subset["unit"].iloc[0])
    baseline = float(pd.to_numeric(subset["baseline_value"], errors="raise").mean())
    current = float(pd.to_numeric(subset["current_value"], errors="raise").mean())
    values = np.asarray([baseline, current], dtype=float)
    labels = [baseline_label, current_label]
    x = np.arange(len(values), dtype=float)
    ymax = max(float(values.max()), 1.0) * 1.22

    ax.bar(x, values, width=0.52, color=[COLORS["baseline"], COLORS["current"]])
    for xpos, value in zip(x, values, strict=True):
        ax.text(xpos, value + ymax * 0.025, f"{value:.1f}", ha="center", va="bottom")

    ax.set_title(f"{panel}. {title}", pad=5)
    ax.set_xticks(x, labels=labels)
    ax.set_ylabel(unit)
    ax.set_ylim(0, ymax)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E5E7EB", linewidth=0.7)
    ax.set_axisbelow(True)


def _rounded_box(
    ax: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    *,
    face: str = "#FFFFFF",
    edge: str = "#D1D5DB",
    lw: float = 0.8,
    radius: float = 0.012,
) -> mpatches.FancyBboxPatch:
    box = mpatches.FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle=f"round,pad=0.006,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
    )
    ax.add_patch(box)
    return box


def _gui_text(
    ax: plt.Axes,
    x: float,
    y: float,
    text: str,
    *,
    size: float = 6.5,
    color: str = "#111827",
    ha: str = "left",
    va: str = "top",
) -> None:
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        fontsize=size,
        color=color,
        ha=ha,
        va=va,
        fontweight="bold",
        fontfamily="Arial",
    )


def _sidebar_button(ax: plt.Axes, x: float, y: float, width: float, label: str, *, active: bool = False) -> None:
    face = "#DBEAFE" if active else "#FFFFFF"
    edge = "#93C5FD" if active else "#E5E7EB"
    color = "#1D4ED8" if active else "#111827"
    _rounded_box(ax, x, y, width, 0.039, face=face, edge=edge, lw=0.7, radius=0.008)
    _gui_text(ax, x + 0.014, y + 0.027, label, size=6.1, color=color)


def _small_status(ax: plt.Axes, x: float, y: float, label: str, value: str, color: str = "#2563EB") -> None:
    _rounded_box(ax, x, y, 0.103, 0.044, face="#FFFFFF", edge="#E5E7EB", lw=0.7, radius=0.008)
    ax.add_patch(
        mpatches.Circle((x + 0.017, y + 0.022), 0.008, transform=ax.transAxes, facecolor=color, edgecolor=color, linewidth=0.5)
    )
    _gui_text(ax, x + 0.032, y + 0.032, label, size=5.1, color="#475569")
    _gui_text(ax, x + 0.032, y + 0.017, value, size=6.0, color="#111827")


def _workflow_card(ax: plt.Axes, x: float, y: float, width: float, title: str, body: str, color: str) -> None:
    _rounded_box(ax, x, y, width, 0.118, face="#FFFFFF", edge="#D1D5DB", lw=0.8, radius=0.01)
    ax.add_patch(
        mpatches.Rectangle((x, y + 0.092), width, 0.026, transform=ax.transAxes, facecolor=color, edgecolor=color, linewidth=0)
    )
    _gui_text(ax, x + 0.012, y + 0.078, title, size=6.2, color="#111827")
    _gui_text(ax, x + 0.012, y + 0.052, body, size=5.4, color="#475569")


def plot_vector_gui_panel(ax: plt.Axes) -> None:
    """Draw a vector equivalent of the GUI home page for editable manuscript output."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    _rounded_box(ax, 0.0, 0.0, 1.0, 1.0, face="#F8FAFC", edge="#CBD5E1", lw=0.9, radius=0.0)

    sidebar_w = 0.275
    ax.add_patch(
        mpatches.Rectangle(
            (0.0, 0.0),
            sidebar_w,
            1.0,
            transform=ax.transAxes,
            facecolor="#F1F5F9",
            edgecolor="#CBD5E1",
            linewidth=0.8,
        )
    )

    _rounded_box(ax, 0.023, 0.926, sidebar_w - 0.046, 0.056, face="#0F172A", edge="#0F172A", radius=0.01)
    _gui_text(ax, 0.040, 0.962, "fp-tools", size=8.6, color="#FFFFFF")
    _gui_text(ax, 0.040, 0.939, "GUI + CLI", size=5.8, color="#CBD5E1")

    _gui_text(ax, 0.025, 0.902, "Workspace", size=6.3, color="#334155")
    _rounded_box(ax, 0.023, 0.853, sidebar_w - 0.046, 0.044, face="#FFFFFF", edge="#D1D5DB", radius=0.008)
    _gui_text(ax, 0.039, 0.881, "Example run", size=5.8, color="#111827")
    _gui_text(ax, 0.039, 0.864, "inputs ready", size=5.2, color="#16A34A")

    _gui_text(ax, 0.025, 0.823, "Commands", size=6.3, color="#334155")
    nav = [
        ("Home", True),
        ("ATACorrect", False),
        ("FootprintScores", False),
        ("BINDetect", False),
        ("PlotAggregate", False),
        ("Pseudobulk", False),
        ("Normalize bigWig", False),
        ("Variants", False),
        ("Motif discovery", False),
        ("History", False),
    ]
    y = 0.756
    for label, active in nav:
        _sidebar_button(ax, 0.023, y, sidebar_w - 0.046, label, active=active)
        y -= 0.044

    _gui_text(ax, 0.025, 0.326, "Shared options", size=6.3, color="#334155")
    options = [
        ("Genome", "genome.fa.gz"),
        ("Motifs", "JASPAR / HOCOMOCO"),
        ("Cores", "1 demo; scale up"),
        ("Output", "examples/gui_outputs"),
    ]
    y = 0.268
    for title, value in options:
        _rounded_box(ax, 0.023, y, sidebar_w - 0.046, 0.043, face="#FFFFFF", edge="#E5E7EB", lw=0.7, radius=0.007)
        _gui_text(ax, 0.037, y + 0.030, title, size=5.2, color="#475569")
        _gui_text(ax, 0.037, y + 0.014, value, size=5.0, color="#111827")
        y -= 0.048

    main_x = sidebar_w + 0.026
    main_w = 1.0 - main_x - 0.025
    _gui_text(ax, main_x, 0.976, "GUI home", size=8.4, color="#111827")
    _gui_text(
        ax,
        main_x,
        0.952,
        "Save configs, run commands, and review outputs.",
        size=5.8,
        color="#475569",
    )

    _small_status(ax, main_x, 0.899, "Inputs", "ready", "#16A34A")
    _small_status(ax, main_x + 0.108, 0.899, "Config", "saved", "#2563EB")
    _small_status(ax, main_x + 0.216, 0.899, "Logs", "live", "#7C3AED")
    _small_status(ax, main_x + 0.324, 0.899, "Plots", "SVG", "#EA580C")

    _rounded_box(ax, main_x, 0.742, main_w, 0.137, face="#EFF6FF", edge="#BFDBFE", lw=0.8, radius=0.012)
    _gui_text(ax, main_x + 0.018, 0.854, "Guided tour", size=7.2, color="#1E3A8A")
    _gui_text(
        ax,
        main_x + 0.018,
        0.826,
        "1  Select command    2  Load example    3  Preview command    4  Run",
        size=5.8,
        color="#1E40AF",
    )
    _gui_text(
        ax,
        main_x + 0.018,
        0.799,
        "Same config for GUI and CLI.",
        size=5.8,
        color="#1E40AF",
    )
    _rounded_box(ax, main_x + main_w - 0.165, 0.764, 0.136, 0.041, face="#2563EB", edge="#2563EB", radius=0.009)
    _gui_text(ax, main_x + main_w - 0.097, 0.791, "Start tour", size=6.0, color="#FFFFFF", ha="center")

    _rounded_box(ax, main_x, 0.548, main_w, 0.168, face="#FFFFFF", edge="#D1D5DB", lw=0.8, radius=0.012)
    _gui_text(ax, main_x + 0.018, 0.690, "CLI-ready config", size=7.0, color="#111827")
    command_lines = [
        "fp-tools-run examples/gui_configs/pseudobulk_footprints_dry_run.yml",
        "ATACorrect -> FootprintScores -> BINDetect",
        "paths, samples, motifs, and outputs are saved once",
    ]
    y = 0.656
    for line in command_lines:
        _rounded_box(ax, main_x + 0.018, y - 0.025, main_w - 0.036, 0.034, face="#F8FAFC", edge="#E5E7EB", lw=0.55, radius=0.006)
        _gui_text(ax, main_x + 0.032, y - 0.002, line, size=5.3, color="#0F172A")
        y -= 0.043

    _gui_text(ax, main_x, 0.517, "Core actions", size=7.1, color="#111827")
    card_gap = 0.014
    card_w = (main_w - 3 * card_gap) / 4
    cards = [
        ("Correct", "Tn5 bias", "#2563EB"),
        ("Score", "footprints", "#7C3AED"),
        ("Detect", "TF activity", "#EA580C"),
        ("Report", "plots + tables", "#16A34A"),
    ]
    for idx, (title, body, color) in enumerate(cards):
        _workflow_card(ax, main_x + idx * (card_w + card_gap), 0.380, card_w, title, body, color)

    _rounded_box(ax, main_x, 0.200, main_w, 0.148, face="#FFFFFF", edge="#D1D5DB", lw=0.8, radius=0.012)
    _gui_text(ax, main_x + 0.018, 0.322, "Templates", size=7.0, color="#111827")
    examples = [
        ("atacorrect.yml", "correction"),
        ("footprints.yml", "scores"),
        ("bindetect.yml", "TF activity"),
        ("aggregate.yml", "profiles"),
        ("pseudobulk.yml", "single-cell"),
        ("variants.yml", "variants"),
    ]
    for idx, (name, label) in enumerate(examples):
        row = idx // 3
        col = idx % 3
        box_w = (main_w - 0.054) / 3
        x = main_x + 0.018 + col * (box_w + 0.018)
        yy = 0.268 - row * 0.051
        _rounded_box(ax, x, yy, box_w, 0.041, face="#F8FAFC", edge="#E5E7EB", lw=0.6, radius=0.007)
        _gui_text(ax, x + 0.010, yy + 0.028, name, size=5.0, color="#111827")
        _gui_text(ax, x + 0.010, yy + 0.012, label, size=4.8, color="#64748B")

    _rounded_box(ax, main_x, 0.050, main_w, 0.118, face="#F8FAFC", edge="#D1D5DB", lw=0.8, radius=0.012)
    _gui_text(ax, main_x + 0.018, 0.142, "Outputs", size=7.0, color="#111827")
    statuses = [
        ("config", "#16A34A"),
        ("logs", "#2563EB"),
        ("tables", "#7C3AED"),
        ("SVG figures", "#EA580C"),
    ]
    for idx, (label, color) in enumerate(statuses):
        x = main_x + 0.020 + idx * ((main_w - 0.050) / 4)
        ax.add_patch(mpatches.Circle((x, 0.091), 0.0085, transform=ax.transAxes, facecolor=color, edgecolor=color, linewidth=0.5))
        _gui_text(ax, x + 0.014, 0.103, label, size=5.5, color="#111827")


def plot_usability_strip(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.text(0.0, 0.96, "User flow", transform=ax.transAxes, ha="left", va="top")
    items = [
        ("1", "Select", "sidebar"),
        ("2", "Guide", "tour"),
        ("3", "Configure", "inputs"),
        ("4", "Run", "results"),
    ]
    box_width = 0.235
    gap = 0.02
    y0 = 0.05
    height = 0.72
    for idx, (number, title, body) in enumerate(items):
        x0 = idx * (box_width + gap)
        rect = mpatches.FancyBboxPatch(
            (x0, y0),
            box_width,
            height,
            boxstyle="round,pad=0.012,rounding_size=0.02",
            transform=ax.transAxes,
            facecolor="#F8FAFC",
            edgecolor="#D1D5DB",
            linewidth=0.8,
        )
        ax.add_patch(rect)
        circle = mpatches.Circle(
            (x0 + 0.035, y0 + height - 0.18),
            0.035,
            transform=ax.transAxes,
            facecolor="#2563EB",
            edgecolor="#2563EB",
            linewidth=0.8,
        )
        ax.add_patch(circle)
        ax.text(x0 + 0.035, y0 + height - 0.18, number, transform=ax.transAxes, color="white", ha="center", va="center")
        ax.text(x0 + 0.085, y0 + height - 0.12, title, transform=ax.transAxes, ha="left", va="top")
        ax.text(x0 + 0.085, y0 + height - 0.36, body, transform=ax.transAxes, ha="left", va="top")


def plot_fig5(source: Path, output: Path) -> None:
    apply_svg_style()
    table = load_source(source)
    fig = plt.figure(figsize=(8.5, 11))
    outer = fig.add_gridspec(
        3,
        1,
        height_ratios=[0.18, 0.675, 0.145],
        left=0.065,
        right=0.985,
        top=0.945,
        bottom=0.022,
        hspace=0.10,
    )
    top_grid = outer[0].subgridspec(
        1,
        2,
        wspace=0.28,
    )
    axes = [fig.add_subplot(top_grid[0, 0]), fig.add_subplot(top_grid[0, 1])]
    plot_metric_bars(axes[0], table, "A", "Runtime")
    plot_metric_bars(axes[1], table, "B", "Peak memory")

    gui_ax = fig.add_subplot(outer[1])
    plot_vector_gui_panel(gui_ax)

    gui_ax.text(-0.025, 0.995, "C", transform=gui_ax.transAxes, ha="left", va="top")
    strip_ax = fig.add_subplot(outer[2])
    plot_usability_strip(strip_ax)
    fig.suptitle("Fig5. Improved performance and GUI support", y=0.982)
    force_arial_bold(fig)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, format="svg")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    plot_fig5(args.source, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
