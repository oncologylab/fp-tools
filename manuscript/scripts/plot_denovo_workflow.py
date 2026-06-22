#!/usr/bin/env python
"""Draw a simplified de novo motif-discovery workflow figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import patches

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
OUT_SVG = REPO_ROOT / "manuscript/figures/fp-tools-de-novo-motif.svg"
BUILD_PDF = REPO_ROOT / "manuscript/build/fp-tools-de-novo-motif.pdf"

BLUE = "#0f4cc8"
BLUE_DARK = "#09379a"
BLUE_LIGHT = "#eef4ff"
GREEN = "#16823a"
GREEN_DARK = "#0b6528"
GREEN_LIGHT = "#eef9f1"
INK = "#111827"
MUTED = "#4b5563"
LINE = "#9aa4b2"
GRAY = "#f8fafc"
BASE_COLORS = {"A": "#2f9e44", "C": "#1971c2", "G": "#f08c00", "T": "#d6336c"}
FONT = "Liberation Sans"


def add_round_box(ax, xy, width, height, edge=BLUE, face="white", lw=1.25, radius=0.75):
    box = patches.FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.35,rounding_size={radius}",
        linewidth=lw,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(box)
    return box


def add_arrow(ax, x1, y1, x2, y2, color="#4b5563", lw=1.25):
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=16, shrinkA=0, shrinkB=0),
    )


def add_motif(ax, x, y, seq, size=12, spacing=0.9):
    for idx, base in enumerate(seq):
        ax.text(
            x + idx * spacing,
            y,
            base,
            ha="center",
            va="center",
            fontsize=size,
            fontweight="bold",
            family=FONT,
            color=BASE_COLORS.get(base, MUTED),
        )


def add_step_badge(ax, x, y, number, title):
    ax.add_patch(patches.Circle((x, y), 0.95, facecolor=BLUE, edgecolor=BLUE_DARK, linewidth=0.9))
    ax.text(x, y - 0.02, str(number), ha="center", va="center", color="white", fontsize=9, fontweight="bold", family=FONT)
    ax.text(x + 1.45, y, title, ha="left", va="center", color=BLUE_DARK, fontsize=9, fontweight="bold", family=FONT)


def add_file_icon(ax, x, y, w=2.6, h=3.3, edge=BLUE_DARK):
    ax.add_patch(patches.Rectangle((x, y), w, h, facecolor="white", edgecolor=edge, linewidth=1.0))
    ax.add_patch(
        patches.Polygon(
            [(x + w * 0.66, y + h), (x + w, y + h), (x + w, y + h * 0.72)],
            facecolor=BLUE_LIGHT,
            edgecolor=edge,
            linewidth=0.8,
        )
    )
    for i in range(3):
        yy = y + h - 0.78 - i * 0.66
        ax.plot([x + 0.38, x + w - 0.38], [yy, yy], color=LINE, linewidth=0.65)


def draw_candidate_box(ax, x, y, w, h):
    add_round_box(ax, (x, y), w, h)
    add_step_badge(ax, x + 2.0, y + h - 1.8, 1, "Candidate intervals")
    ax.text(x + w / 2, y + h - 4.0, "Footprint/BED sites", ha="center", fontsize=9, color=INK, fontweight="bold", family=FONT)

    ax.plot([x + 3.1, x + w - 3.1], [y + 5.9, y + 5.9], color=INK, linewidth=1.0)
    for offset, width in [(3.8, 3.0), (7.5, 4.7), (13.1, 2.9), (16.6, 3.4)]:
        ax.add_patch(patches.Rectangle((x + offset, y + 4.05), width, 0.5, facecolor=BLUE, edgecolor=BLUE_DARK, linewidth=0.65))
    ax.add_patch(patches.Rectangle((x + 9.3, y + 5.45), 4.1, 0.62, facecolor=BLUE, edgecolor=BLUE_DARK, linewidth=0.65))
    ax.plot([x + 11.35, x + 11.35], [y + 3.85, y + 6.65], color=LINE, linewidth=1.1)
    ax.text(x + w / 2, y + 1.75, "chr:start-end", ha="center", fontsize=9, color=INK, fontweight="bold", family=FONT)


def draw_fasta_box(ax, x, y, w, h):
    add_round_box(ax, (x, y), w, h)
    add_step_badge(ax, x + 2.0, y + h - 1.8, 2, "Sequence export")
    ax.text(x + w / 2, y + h - 4.0, "+/-75 bp FASTA", ha="center", fontsize=9, color=INK, fontweight="bold", family=FONT)

    ax.plot([x + 3.2, x + w - 3.2], [y + 5.75, y + 5.75], color=INK, linewidth=1.0)
    ax.add_patch(patches.Rectangle((x + 9.2, y + 5.45), 3.8, 0.58, facecolor=BLUE, edgecolor=BLUE_DARK, linewidth=0.65))
    ax.plot([x + 4.4, x + 18.0], [y + 6.85, y + 6.85], color=INK, linestyle=(0, (2, 1.6)), linewidth=0.8)
    add_file_icon(ax, x + 7.3, y + 1.7, w=2.55, h=3.05)
    ax.text(x + 10.85, y + 3.95, ">cand_0001", ha="left", va="center", fontsize=7.0, family=FONT, color=INK, fontweight="bold")
    add_motif(ax, x + 11.8, y + 2.75, "ACTGAC", size=7.3, spacing=0.5)


def draw_discovery_box(ax, x, y, w, h):
    add_round_box(ax, (x, y), w, h, face=BLUE_LIGHT)
    add_step_badge(ax, x + 2.0, y + h - 1.8, 3, "Motif discovery")
    ax.text(x + w / 2, y + h - 4.0, "MEME / DREME / STREME", ha="center", fontsize=9, color=BLUE_DARK, fontweight="bold", family=FONT)

    add_file_icon(ax, x + 3.0, y + 3.15, w=2.1, h=2.75)
    add_arrow(ax, x + 5.8, y + 4.55, x + 8.4, y + 4.55, color=LINE)
    add_motif(ax, x + 12.7, y + 5.65, "ACTGAC", size=9.6, spacing=0.65)
    add_motif(ax, x + 12.45, y + 3.75, "TGACCA", size=9.6, spacing=0.65)
    ax.text(x + 14.35, y + 1.75, "de novo motifs", ha="center", fontsize=9, color=INK, fontweight="bold", family=FONT)


def draw_output_box(ax, x, y, w, h):
    add_round_box(ax, (x, y), w, h)
    add_step_badge(ax, x + 2.0, y + h - 1.8, 4, "Motif summary")
    ax.text(x + w / 2, y + h - 4.0, "TSV / HTML report", ha="center", fontsize=9, color=BLUE_DARK, fontweight="bold", family=FONT)

    table_x = x + 2.0
    table_y = y + 1.55
    table_w = w - 4.0
    table_h = 5.45
    ax.add_patch(patches.Rectangle((table_x, table_y), table_w, table_h, facecolor=GRAY, edgecolor=LINE, linewidth=0.8))
    for frac in (0.46, 0.74):
        xx = table_x + table_w * frac
        ax.plot([xx, xx], [table_y, table_y + table_h], color=LINE, linewidth=0.6)
    for i in range(1, 4):
        yy = table_y + table_h - i * 1.5
        ax.plot([table_x, table_x + table_w], [yy, yy], color=LINE, linewidth=0.55)
    ax.text(table_x + 2.6, table_y + table_h - 0.85, "motif", ha="center", fontsize=6.0, fontweight="bold", family=FONT)
    ax.text(table_x + table_w * 0.60, table_y + table_h - 0.85, "E/q", ha="center", fontsize=6.0, fontweight="bold", family=FONT)
    ax.text(table_x + table_w * 0.86, table_y + table_h - 0.85, "sites", ha="center", fontsize=6.0, fontweight="bold", family=FONT)
    rows = [("ACTGAC", "1e-12", "1234"), ("TGACCA", "4e-08", "987"), ("CCTGAG", "2e-06", "654")]
    for i, (motif, evalue, sites) in enumerate(rows):
        yy = table_y + table_h - 2.25 - i * 1.5
        add_motif(ax, table_x + 1.35, yy, motif, size=6.7, spacing=0.42)
        ax.text(table_x + table_w * 0.60, yy, evalue, ha="center", va="center", fontsize=6.0, family=FONT)
        ax.text(table_x + table_w * 0.86, yy, sites, ha="center", va="center", fontsize=6.0, family=FONT)


def main() -> int:
    fig, ax = plt.subplots(figsize=(13.25, 3.15))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 22)
    ax.axis("off")

    ax.text(50, 19.85, "De Novo Motif Discovery", ha="center", va="center", fontsize=22, fontweight="bold", color=INK, family=FONT)
    ax.text(
        50,
        18.05,
        "Candidate-centered FASTA and reproducible external motif discovery.",
        ha="center",
        va="center",
        fontsize=9,
        fontweight="bold",
        family=FONT,
        color=INK,
        style="italic",
    )

    panel_y = 3.0
    panel_h = 12.0
    panel_w = 21.4
    xs = [1.7, 26.4, 51.1, 75.8]
    draw_candidate_box(ax, xs[0], panel_y, panel_w, panel_h)
    draw_fasta_box(ax, xs[1], panel_y, panel_w, panel_h)
    draw_discovery_box(ax, xs[2], panel_y, panel_w, panel_h)
    draw_output_box(ax, xs[3], panel_y, panel_w, panel_h)
    for left in xs[:3]:
        add_arrow(ax, left + panel_w + 1.1, panel_y + panel_h / 2, left + 24.0, panel_y + panel_h / 2)

    OUT_SVG.parent.mkdir(parents=True, exist_ok=True)
    BUILD_PDF.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_SVG, bbox_inches="tight")
    fig.savefig(BUILD_PDF, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT_SVG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
