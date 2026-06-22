#!/usr/bin/env python3
"""Draw the simplified replicate-aware diff-footprints schematic."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches
from matplotlib.path import Path as MplPath


BLUE = "#1857c9"
BLUE_LIGHT = "#dce9ff"
ORANGE = "#ff6b1a"
ORANGE_LIGHT = "#ffe5d6"
GREEN = "#18823a"
GREEN_LIGHT = "#dcf4e3"
INK = "#141923"
MUTED = "#5b6677"
BORDER = "#d6dde8"
BG = "#ffffff"


def add_round_box(ax, xy, width, height, edge, face="#ffffff", lw=1.8, radius=0.025):
    box = patches.FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        linewidth=lw,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(box)
    return box


def add_arrow(ax, x1, y1, x2, y2, color="#5a6575"):
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(arrowstyle="-|>", lw=2.0, color=color, shrinkA=0, shrinkB=0, mutation_scale=18),
    )


def add_track(ax, x, y, w, h, color, seed=0):
    rng = np.random.default_rng(seed)
    xs = np.linspace(0, 1, 110)
    centers = np.array([0.18, 0.45, 0.72]) + rng.normal(0, 0.025, 3)
    widths = np.array([0.025, 0.035, 0.025])
    vals = 0.08 + sum(np.exp(-0.5 * ((xs - c) / s) ** 2) for c, s in zip(centers, widths))
    vals += 0.035 * rng.normal(size=xs.size)
    vals = np.clip(vals, 0, None)
    vals = vals / vals.max()
    px = x + xs * w
    py = y + vals * h
    verts = [(x, y), *zip(px, py), (x + w, y), (x, y)]
    codes = [MplPath.MOVETO] + [MplPath.LINETO] * len(px) + [MplPath.LINETO, MplPath.CLOSEPOLY]
    ax.add_patch(patches.PathPatch(MplPath(verts, codes), facecolor=color, edgecolor=color, alpha=0.88, lw=0.8))


def add_matrix(ax, x, y, w, h):
    add_round_box(ax, (x, y), w, h, "#c8d1df", "#fbfdff", lw=1.0, radius=0.01)
    rows, cols = 5, 6
    for i in range(1, rows):
        yy = y + h * i / rows
        ax.plot([x, x + w], [yy, yy], color="#d7dee9", lw=0.8)
    for j in range(1, cols):
        xx = x + w * j / cols
        ax.plot([xx, xx], [y, y + h], color="#d7dee9", lw=0.8)
    ax.text(x + w * 0.5, y + h + 0.012, "motif-score matrix", ha="center", va="bottom", fontsize=8.5, weight="bold", color=INK)
    for r in range(rows):
        for c in range(cols):
            color = BLUE if c < 3 else ORANGE
            alpha = 0.25 + 0.12 * ((r + c) % 4)
            ax.add_patch(
                patches.Circle(
                    (x + w * (c + 0.5) / cols, y + h * (rows - r - 0.5) / rows),
                    min(w / cols, h / rows) * 0.17,
                    facecolor=color,
                    edgecolor="none",
                    alpha=alpha,
                )
            )


def add_volcano(ax, x, y, w, h):
    add_round_box(ax, (x, y), w, h, "#c8d1df", "#ffffff", lw=1.0, radius=0.01)
    rng = np.random.default_rng(4)
    xs = rng.normal(0, 0.35, 55)
    ys = rng.gamma(1.3, 0.18, 55) + 0.08
    ax.scatter(x + w * (0.5 + xs * 0.38), y + h * np.clip(ys, 0.06, 0.86), s=7, color="#9aa6b8", alpha=0.75)
    for side, color in [(-1, BLUE), (1, ORANGE)]:
        sx = rng.normal(side * 0.78, 0.08, 10)
        sy = rng.uniform(0.45, 0.9, 10)
        ax.scatter(x + w * (0.5 + sx * 0.38), y + h * sy, s=12, color=color, alpha=0.9)
    ax.plot([x + w * 0.5, x + w * 0.5], [y + h * 0.16, y + h * 0.92], color="#9aa6b8", lw=0.8, ls="--")
    ax.plot([x + w * 0.12, x + w * 0.9], [y + h * 0.32, y + h * 0.32], color="#9aa6b8", lw=0.8, ls="--")
    ax.text(x + w * 0.5, y + 0.018, "volcano evidence", ha="center", va="bottom", fontsize=8, weight="bold", color=INK)


def aggregate_curve(xs, amp=1.0, offset=0.0):
    flank = 0.35 * np.exp(-0.5 * ((np.abs(xs) - 0.45) / 0.18) ** 2)
    center = 0.42 * np.exp(-0.5 * (xs / 0.08) ** 2)
    return offset + amp * (flank - center)


def add_aggregate_plot(ax, x, y, w, h, title="", small=False):
    add_round_box(ax, (x, y), w, h, "#c8d1df", "#ffffff", lw=1.0, radius=0.01)
    px0, py0 = x + w * 0.12, y + h * 0.16
    pw, ph = w * 0.78, h * 0.68
    ax.plot([px0, px0 + pw], [py0 + ph * 0.5, py0 + ph * 0.5], color="#d9e0ea", lw=0.8)
    ax.plot([px0 + pw * 0.5, px0 + pw * 0.5], [py0, py0 + ph], color="#7f8998", lw=0.9)
    xs = np.linspace(-1, 1, 140)
    for i, jitter in enumerate([-0.035, 0.0, 0.03]):
        yv = aggregate_curve(xs, 0.9 + i * 0.07, jitter)
        ax.plot(px0 + (xs + 1) / 2 * pw, py0 + (yv + 0.55) / 1.1 * ph, color=BLUE, alpha=0.28, lw=1.0)
    for i, jitter in enumerate([-0.015, 0.025, 0.055]):
        yv = aggregate_curve(xs, 1.2 + i * 0.05, jitter)
        ax.plot(px0 + (xs + 1) / 2 * pw, py0 + (yv + 0.55) / 1.1 * ph, color=ORANGE, alpha=0.30, lw=1.0)
    ax.plot(px0 + (xs + 1) / 2 * pw, py0 + (aggregate_curve(xs, 1.0, 0.0) + 0.55) / 1.1 * ph, color=BLUE, lw=2.2)
    ax.plot(px0 + (xs + 1) / 2 * pw, py0 + (aggregate_curve(xs, 1.25, 0.035) + 0.55) / 1.1 * ph, color=ORANGE, lw=2.2)
    if title:
        ax.text(x + w * 0.5, y + h - 0.015, title, ha="center", va="top", fontsize=7.8 if small else 9, weight="bold", color=INK)
    if not small:
        ax.text(x + w * 0.5, y + 0.012, "motif center", ha="center", va="bottom", fontsize=7.5, color=MUTED)


def add_browser_mock(ax, x, y, w, h):
    add_round_box(ax, (x, y), w, h, "#aeb8c8", "#ffffff", lw=1.2, radius=0.018)
    ax.add_patch(patches.Rectangle((x, y + h * 0.9), w, h * 0.1, facecolor="#f0f4fa", edgecolor="#d6dde8", lw=0.8))
    for i, color in enumerate(["#fb7185", "#fbbf24", "#34d399"]):
        ax.add_patch(patches.Circle((x + 0.018 + 0.014 * i, y + h * 0.95), 0.0045, color=color))
    ax.text(x + w * 0.52, y + h * 0.948, "interactive diff-footprints report", ha="center", va="center", fontsize=8.2, weight="bold", color=INK)
    add_volcano(ax, x + w * 0.035, y + h * 0.13, w * 0.33, h * 0.68)
    logo_x = x + w * 0.40
    for i, label in enumerate(["SPIB", "CEBPB", "RUNX3"]):
        yy = y + h * (0.69 - i * 0.19)
        add_round_box(ax, (logo_x, yy), w * 0.19, h * 0.135, "#d7dee9", "#fbfdff", lw=0.8, radius=0.008)
        ax.text(logo_x + w * 0.095, yy + h * 0.082, label, ha="center", va="center", fontsize=7.6, weight="bold", color=INK)
        for j, base in enumerate("ACGT"):
            ax.text(logo_x + w * (0.032 + j * 0.034), yy + h * 0.036, base, ha="center", va="center", fontsize=7.5, color=[BLUE, ORANGE, GREEN, "#7c3aed"][j], weight="bold")
    grid_x, grid_y = x + w * 0.63, y + h * 0.16
    cell_w, cell_h = w * 0.155, h * 0.285
    for r in range(2):
        for c in range(2):
            add_aggregate_plot(ax, grid_x + c * w * 0.17, grid_y + (1 - r) * h * 0.32, cell_w, cell_h, ["SPIB", "CEBPB", "RUNX3", "IRF"][r * 2 + c], small=True)


def draw_figure(output_prefix: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "axes.linewidth": 0.8,
        }
    )
    fig = plt.figure(figsize=(16.72, 9.41), dpi=100, facecolor=BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.5, 0.955, "Replicate-aware differential footprint reporting", ha="center", va="center", fontsize=35, weight="black", color=INK)
    ax.text(
        0.5,
        0.914,
        "Q95-scaled signal, replicate-aware TF comparison, customizable aggregate plots, and interactive multi-motif HTML review",
        ha="center",
        va="center",
        fontsize=16,
        style="italic",
        color=MUTED,
    )

    panels = [
        (0.035, 0.205, 0.285, 0.62, BLUE, BLUE_LIGHT, "1", "Replicate-aware comparison"),
        (0.36, 0.205, 0.285, 0.62, GREEN, GREEN_LIGHT, "2", "Normalized aggregate profiles"),
        (0.685, 0.205, 0.28, 0.62, "#3f4a5f", "#f6f8fb", "3", "Interactive HTML output"),
    ]
    for x, y, w, h, edge, face, number, title in panels:
        add_round_box(ax, (x, y), w, h, edge, face, lw=2.0, radius=0.02)
        ax.add_patch(patches.Circle((x + 0.03, y + h - 0.045), 0.019, facecolor=edge, edgecolor=edge))
        ax.text(x + 0.03, y + h - 0.046, number, ha="center", va="center", fontsize=15, color="white", weight="black")
        ax.text(x + 0.06, y + h - 0.047, title, ha="left", va="center", fontsize=14.5, color=edge if number != "3" else INK, weight="black")
    add_arrow(ax, 0.326, 0.515, 0.355, 0.515)
    add_arrow(ax, 0.651, 0.515, 0.68, 0.515)

    # Panel 1
    x, y, w, h = panels[0][:4]
    ax.text(x + w * 0.5, y + h - 0.105, "Repeated condition labels define replicates", ha="center", va="center", fontsize=10.7, color=INK, weight="bold")
    for row, (label, color, light) in enumerate([("K562", BLUE, BLUE_LIGHT), ("HepG2", ORANGE, ORANGE_LIGHT)]):
        yy = y + h * (0.62 - row * 0.23)
        ax.text(x + 0.034, yy + 0.045, label, ha="left", va="center", fontsize=11, color=color, weight="black")
        for i in range(3):
            tx = x + 0.09 + i * 0.058
            add_round_box(ax, (tx, yy), 0.048, 0.083, color, "#ffffff", lw=1.0, radius=0.008)
            add_track(ax, tx + 0.006, yy + 0.018, 0.036, 0.045, color, seed=10 + row * 10 + i)
            ax.text(tx + 0.024, yy + 0.071, f"rep {i + 1}", ha="center", va="center", fontsize=6.8, color=color, weight="bold")
    add_matrix(ax, x + 0.045, y + 0.08, 0.11, 0.15)
    add_arrow(ax, x + 0.165, y + 0.155, x + 0.195, y + 0.155)
    add_volcano(ax, x + 0.205, y + 0.075, 0.095, 0.18)
    ax.text(x + w * 0.5, y + h - 0.137, "Q95-scaled footprint scores", ha="center", va="center", fontsize=10.2, color=MUTED, weight="bold")
    ax.text(x + w * 0.5, y + 0.035, "effect size + p/FDR + replicate support", ha="center", va="center", fontsize=9.8, color=MUTED, weight="bold")

    # Panel 2
    x, y, w, h = panels[1][:4]
    ax.text(x + w * 0.5, y + h - 0.105, "Q95-scaled corrected bigWigs", ha="center", va="center", fontsize=10.7, color=INK, weight="bold")
    bed_x, bed_y = x + 0.035, y + h * 0.58
    add_round_box(ax, (bed_x, bed_y), 0.095, 0.11, GREEN, "#ffffff", lw=1.0, radius=0.008)
    for i in range(5):
        yy = bed_y + 0.085 - i * 0.017
        ax.plot([bed_x + 0.014, bed_x + 0.08], [yy, yy], color="#8aa0b8", lw=1.1)
        ax.add_patch(patches.Rectangle((bed_x + 0.02 + i * 0.006, yy - 0.004), 0.016, 0.008, facecolor=GREEN, edgecolor="none", alpha=0.85))
    ax.text(bed_x + 0.0475, bed_y + 0.017, "*_all.bed sites", ha="center", va="center", fontsize=8.0, color=GREEN, weight="black")
    for i, color in enumerate([BLUE, BLUE, BLUE, ORANGE, ORANGE, ORANGE]):
        add_track(ax, x + 0.16, y + h * (0.67 - i * 0.046), 0.095, 0.028, color, seed=40 + i)
    add_arrow(ax, x + 0.235, y + h * 0.57, x + 0.19, y + h * 0.43, color=GREEN)
    add_aggregate_plot(ax, x + 0.052, y + 0.095, 0.205, 0.24, "motif-centered profiles")
    ax.text(x + w * 0.5, y + 0.055, "significant, top, or all motifs", ha="center", va="center", fontsize=10, color=MUTED, weight="bold")
    ax.text(x + w * 0.5, y + 0.03, "custom colors, samples, conditions, layout", ha="center", va="center", fontsize=9.3, color=MUTED)

    # Panel 3
    x, y, w, h = panels[2][:4]
    ax.text(x + w * 0.5, y + h - 0.105, "Standalone multi-motif HTML review", ha="center", va="center", fontsize=10.7, color=INK, weight="bold")
    add_browser_mock(ax, x + 0.025, y + 0.075, w - 0.05, h * 0.68)
    ax.text(x + w * 0.5, y + 0.055, "volcano + logos + aggregate panels", ha="center", va="center", fontsize=10, color=MUTED, weight="bold")
    ax.text(x + w * 0.5, y + 0.03, "multiple selected motifs + SVG export", ha="center", va="center", fontsize=9.3, color=MUTED)

    # Bottom workflow strip
    strip_y, strip_h = 0.075, 0.085
    add_round_box(ax, (0.09, strip_y), 0.82, strip_h, "#b5bdca", "#fbfcfe", lw=1.2, radius=0.016)
    steps = ["corrected bigWigs", "Q95 scaling", "footprint scoring", "diff-footprints\nnormalization none", "interactive report"]
    xs = np.linspace(0.17, 0.83, len(steps))
    for i, (sx, step) in enumerate(zip(xs, steps)):
        ax.add_patch(patches.Circle((sx, strip_y + strip_h * 0.58), 0.018, facecolor=[BLUE, GREEN, GREEN, ORANGE, "#3f4a5f"][i], edgecolor="white", lw=1.0))
        ax.text(sx, strip_y + strip_h * 0.58, str(i + 1), ha="center", va="center", fontsize=10, weight="black", color="white")
        ax.text(sx, strip_y + strip_h * 0.22, step, ha="center", va="center", fontsize=9.2, color=INK, weight="bold", linespacing=0.92)
        if i < len(steps) - 1:
            add_arrow(ax, sx + 0.032, strip_y + strip_h * 0.58, xs[i + 1] - 0.032, strip_y + strip_h * 0.58, color="#9aa6b8")

    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    for ext in ["png", "svg", "pdf"]:
        fig.savefig(output_prefix.with_suffix(f".{ext}"), dpi=300 if ext == "png" else None, bbox_inches="tight", pad_inches=0.04, facecolor=BG)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("docs/assets/fp-tools-replicate-bindetect"),
        help="Output path prefix; .png, .svg, and .pdf are written.",
    )
    args = parser.parse_args()
    draw_figure(args.output_prefix)
    print(f"Wrote {args.output_prefix.with_suffix('.png')}")
    print(f"Wrote {args.output_prefix.with_suffix('.svg')}")
    print(f"Wrote {args.output_prefix.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
