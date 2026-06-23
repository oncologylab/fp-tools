#!/usr/bin/env python
"""Plot compact Fig5 engineering-performance panels."""

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
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
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
    y = np.arange(len(values), dtype=float)
    ymax = max(float(values.max()), 1.0) * 1.22

    ax.barh(y, values, height=0.32, color=[COLORS["baseline"], COLORS["current"]])
    for idx, (ypos, value) in enumerate(zip(y, values, strict=True)):
        label_color = "#111827" if idx == 0 else "#FFFFFF"
        ax.text(value - ymax * 0.035, ypos, f"{value:.1f}", ha="right", va="center", color=label_color)

    ax.text(
        -0.22,
        1.03,
        panel,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=18,
        fontweight="bold",
        clip_on=False,
    )
    ax.set_title(title, pad=7)
    ax.set_yticks(y, labels=labels)
    ax.set_xlabel(unit)
    ax.set_xlim(0, ymax)
    ax.set_ylim(-0.55, len(values) - 0.45)
    ax.invert_yaxis()
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.65)
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
    """Draw a vector equivalent of the current GUI home page."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    _rounded_box(ax, 0.0, 0.0, 1.0, 1.0, face="#F3F6FA", edge="#CBD5E1", lw=0.9, radius=0.0)

    sidebar_w = 0.255
    ax.add_patch(
        mpatches.Rectangle(
            (0.0, 0.0),
            sidebar_w,
            1.0,
            transform=ax.transAxes,
            facecolor="#0F172A",
            edgecolor="#0F172A",
            linewidth=0,
        )
    )
    _rounded_box(ax, 0.022, 0.900, sidebar_w - 0.044, 0.072, face="#101D33", edge="#2B3A55", radius=0.010)
    _gui_text(ax, 0.040, 0.952, "fp-tools", size=10.6, color="#FFFFFF")
    _gui_text(ax, 0.040, 0.922, "Command-ready workflows", size=6.6, color="#AEBBD0")

    def nav_item(y: float, label: str, *, active: bool = False) -> None:
        face = "#1D4ED8" if active else "#0F172A"
        edge = "#1D4ED8" if active else "#0F172A"
        _rounded_box(ax, 0.022, y, sidebar_w - 0.044, 0.040, face=face, edge=edge, lw=0.5, radius=0.008)
        if active:
            ax.add_patch(
                mpatches.Rectangle(
                    (0.026, y + 0.006),
                    0.004,
                    0.028,
                    transform=ax.transAxes,
                    facecolor="#72E0B2",
                    edgecolor="#72E0B2",
                    linewidth=0,
                )
            )
        _gui_text(ax, 0.040, y + 0.027, label, size=6.6, color="#E5EDF6")

    y = 0.848
    for group, labels in [
        ("Overview", ["Home", "Run History"]),
        ("Core Workflow", ["atac-correct", "call-footprints", "match-motifs", "diff-footprints"]),
        ("Reports", ["normalize-bigwig", "plot-aggregate", "plot-aggregate-batch"]),
        ("Single-cell ATAC", ["pseudobulk-fragments", "pseudobulk-footprints"]),
    ]:
        _gui_text(ax, 0.030, y + 0.017, group.upper(), size=5.3, color="#93A4B8")
        y -= 0.040
        for label in labels:
            nav_item(y, label, active=label == "Home")
            y -= 0.043
        y -= 0.010

    main_x = sidebar_w + 0.030
    main_w = 1.0 - main_x - 0.030

    _rounded_box(ax, main_x, 0.855, main_w, 0.120, face="#FFFFFF", edge="#DBE3EC", lw=0.8, radius=0.010)
    _gui_text(ax, main_x + 0.020, 0.944, "Run footprint workflows", size=11.6, color="#111827")
    _gui_text(
        ax,
        main_x + 0.020,
        0.910,
        "Choose a command, load an example, review YAML, then launch.",
        size=7.1,
        color="#5B6778",
    )

    step_y = 0.735
    step_gap = 0.012
    step_w = (main_w - 3 * step_gap) / 4
    for idx, (title, body) in enumerate(
        [
            ("1. Choose", "sidebar command"),
            ("2. Load", "example YAML"),
            ("3. Review", "editable config"),
            ("4. Inspect", "logs and reports"),
        ]
    ):
        x = main_x + idx * (step_w + step_gap)
        _rounded_box(ax, x, step_y, step_w, 0.092, face="#FFFFFF", edge="#DBE3EC", lw=0.8, radius=0.010)
        _gui_text(ax, x + 0.012, step_y + 0.064, title, size=6.9, color="#111827")
        _gui_text(ax, x + 0.012, step_y + 0.037, body, size=5.7, color="#5B6778")

    left_w = main_w * 0.615
    right_x = main_x + left_w + 0.018
    right_w = main_w - left_w - 0.018

    _rounded_box(ax, main_x, 0.220, left_w, 0.480, face="#FFFFFF", edge="#DBE3EC", lw=0.8, radius=0.010)
    _gui_text(ax, main_x + 0.018, 0.668, "diff-footprints", size=9.0, color="#111827")

    def field_box(x: float, y: float, width: float, label: str, value: str) -> None:
        _gui_text(ax, x, y + 0.052, label, size=5.6, color="#344054")
        _rounded_box(ax, x, y, width, 0.042, face="#F8FAFC", edge="#DBE3EC", lw=0.65, radius=0.007)
        _gui_text(ax, x + 0.009, y + 0.027, value, size=5.4, color="#334155")

    col_gap = 0.016
    col_w = (left_w - 0.052 - col_gap) / 2
    x1 = main_x + 0.020
    x2 = x1 + col_w + col_gap
    fields = [
        ("Footprint bigWigs", "K562_rep1.bw, HepG2_rep1.bw"),
        ("Condition names", "K562 K562 HepG2 HepG2"),
        ("Genome FASTA", "hg38.fa"),
        ("Motif database", "jaspar2026_vertebrates"),
        ("Peaks BED", "merged_peaks.bed"),
        ("Output folder", "results/diff_footprints"),
    ]
    for idx, (label, value) in enumerate(fields):
        row = idx // 2
        col_x = x1 if idx % 2 == 0 else x2
        field_box(col_x, 0.578 - row * 0.112, col_w, label, value)

    _rounded_box(ax, main_x + 0.020, 0.265, left_w - 0.040, 0.090, face="#EAF2FF", edge="#BDD7FF", lw=0.8, radius=0.009)
    _gui_text(ax, main_x + 0.038, 0.328, "Guided workflow", size=7.3, color="#173B73")
    _gui_text(ax, main_x + 0.038, 0.300, "Load example, edit paths, check YAML, run.", size=5.9, color="#173B73")

    _rounded_box(ax, right_x, 0.470, right_w, 0.230, face="#FFFFFF", edge="#DBE3EC", lw=0.8, radius=0.010)
    _gui_text(ax, right_x + 0.016, 0.668, "Run", size=9.0, color="#111827")
    metric_w = (right_w - 0.048) / 3
    for idx, (label, value) in enumerate([("Tool", "Diff"), ("Samples", "6"), ("Report", "HTML")]):
        x = right_x + 0.016 + idx * (metric_w + 0.008)
        _rounded_box(ax, x, 0.603, metric_w, 0.048, face="#F8FAFC", edge="#DBE3EC", lw=0.65, radius=0.007)
        _gui_text(ax, x + 0.008, 0.635, label.upper(), size=4.7, color="#5B6778")
        _gui_text(ax, x + 0.008, 0.616, value, size=5.7, color="#111827")
    _rounded_box(ax, right_x + 0.016, 0.515, right_w - 0.032, 0.065, face="#0F172A", edge="#0F172A", lw=0.8, radius=0.008)
    for idx, line in enumerate(
        [
            "tool: diff-footprints",
            "motif_db: jaspar2026_vertebrates",
            "plot_aggregate: sig",
        ]
    ):
        _gui_text(ax, right_x + 0.030, 0.563 - idx * 0.019, line, size=4.9, color="#DBEAFE")
    _rounded_box(ax, right_x + 0.016, 0.485, right_w - 0.032, 0.030, face="#173B73", edge="#173B73", radius=0.007)
    _gui_text(ax, right_x + right_w / 2, 0.505, "Start run", size=5.8, color="#FFFFFF", ha="center")

    _rounded_box(ax, right_x, 0.220, right_w, 0.220, face="#FFFFFF", edge="#DBE3EC", lw=0.8, radius=0.010)
    _gui_text(ax, right_x + 0.016, 0.408, "Run History", size=8.2, color="#111827")
    for idx, (label, value, color) in enumerate(
        [
            ("Status", "Complete", "#15956B"),
            ("Config", "saved YAML", "#1D4ED8"),
            ("Outputs", "tables + SVG", "#EA580C"),
        ]
    ):
        y0 = 0.360 - idx * 0.050
        ax.add_patch(mpatches.Circle((right_x + 0.030, y0 + 0.012), 0.007, transform=ax.transAxes, facecolor=color, edgecolor=color))
        _gui_text(ax, right_x + 0.045, y0 + 0.025, label, size=5.4, color="#5B6778")
        _gui_text(ax, right_x + 0.045, y0 + 0.009, value, size=5.7, color="#111827")

    _rounded_box(ax, main_x, 0.055, main_w, 0.125, face="#FFFFFF", edge="#DBE3EC", lw=0.8, radius=0.010)
    _gui_text(ax, main_x + 0.018, 0.148, "Outputs stay CLI-compatible", size=8.0, color="#111827")
    for idx, (label, color) in enumerate(
        [
            ("config.yml", "#15956B"),
            ("run.log", "#1D4ED8"),
            ("results.tsv", "#7C3AED"),
            ("report.html", "#EA580C"),
            ("figure.svg", "#334155"),
        ]
    ):
        x = main_x + 0.020 + idx * ((main_w - 0.048) / 5)
        ax.add_patch(mpatches.Circle((x, 0.095), 0.008, transform=ax.transAxes, facecolor=color, edgecolor=color))
        _gui_text(ax, x + 0.014, 0.107, label, size=5.7, color="#111827")


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
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(3.35, 1.65),
        gridspec_kw={"wspace": 0.60},
    )
    plot_metric_bars(axes[0], table, "A", "Runtime")
    plot_metric_bars(axes[1], table, "B", "Peak memory")
    fig.subplots_adjust(left=0.18, right=0.98, top=0.78, bottom=0.28)
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
