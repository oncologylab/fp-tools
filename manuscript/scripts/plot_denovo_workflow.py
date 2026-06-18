#!/usr/bin/env python
"""Draw the de novo motif-discovery preparation workflow figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import patches

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
OUT_PREFIX = REPO_ROOT / "manuscript/figures/fp-tools-de-novo-motif-hires"

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


def add_round_box(ax, xy, width, height, edge=BLUE, face="white", lw=1.35, radius=0.75):
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


def add_step(ax, x, y, number, title, subtitle):
    circ = patches.Circle((x, y), 1.25, facecolor=BLUE, edgecolor=BLUE_DARK, linewidth=0.9)
    ax.add_patch(circ)
    ax.text(x, y - 0.03, str(number), ha="center", va="center", color="white", fontsize=12, fontweight="bold")
    ax.text(x + 2.1, y + 0.45, title, ha="left", va="center", color=BLUE_DARK, fontsize=10.5, fontweight="bold")
    ax.text(x + 2.1, y - 1.2, subtitle, ha="left", va="center", color=INK, fontsize=7.6)


def add_arrow(ax, x1, y1, x2, y2, color="#4b5563", lw=1.2):
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=14, shrinkA=0, shrinkB=0),
    )


def add_motif(ax, x, y, seq, size=11, spacing=0.8):
    for idx, base in enumerate(seq):
        ax.text(
            x + idx * spacing,
            y,
            base,
            ha="center",
            va="center",
            fontsize=size,
            fontweight="bold",
            family="monospace",
            color=BASE_COLORS.get(base, MUTED),
        )


def add_file_icon(ax, x, y, w=2.4, h=3.0, color=BLUE_DARK):
    ax.add_patch(patches.Rectangle((x, y), w, h, facecolor="white", edgecolor=color, linewidth=1.0))
    ax.add_patch(patches.Polygon([(x + w * 0.65, y + h), (x + w, y + h), (x + w, y + h * 0.72)], facecolor=BLUE_LIGHT, edgecolor=color, linewidth=0.8))
    for i in range(3):
        ax.plot([x + 0.35, x + w - 0.35], [y + h - 0.75 - i * 0.65, y + h - 0.75 - i * 0.65], color=LINE, linewidth=0.6)


def add_database(ax, x, y, w=3.0, h=2.9):
    ax.add_patch(patches.Ellipse((x + w / 2, y + h), w, 0.65, facecolor=GREEN_LIGHT, edgecolor=GREEN_DARK, linewidth=1.0))
    ax.add_patch(patches.Rectangle((x, y), w, h, facecolor=GREEN_LIGHT, edgecolor=GREEN_DARK, linewidth=1.0))
    ax.add_patch(patches.Ellipse((x + w / 2, y), w, 0.65, facecolor=GREEN_LIGHT, edgecolor=GREEN_DARK, linewidth=1.0))
    ax.plot([x, x], [y, y + h], color=GREEN_DARK, linewidth=1.0)
    ax.plot([x + w, x + w], [y, y + h], color=GREEN_DARK, linewidth=1.0)


def draw_candidate_panel(ax):
    add_round_box(ax, (1.5, 9), 18, 35.5, face="white")
    ax.text(10.5, 42.7, "BED-like candidate intervals", ha="center", va="center", color=BLUE_DARK, fontsize=8.2, fontweight="bold")
    ax.plot([4.5, 17.2], [39, 39], color=INK, linewidth=0.9)
    ax.text(9.1, 39.7, "10 kb", fontsize=6.8, ha="center")
    ax.text(15.2, 39.7, "20 kb", fontsize=6.8, ha="center")
    for x, y, w in [(5.0, 36.5, 2.6), (9.1, 36.4, 4.5), (4.2, 34.6, 2.2), (7.8, 34.5, 3.0), (5.0, 29.8, 3.1), (9.0, 29.9, 2.9), (12.8, 29.9, 1.5), (15.4, 29.9, 1.1), (4.5, 28.2, 2.9), (8.9, 28.3, 3.3), (5.9, 26.7, 2.5), (11.0, 26.6, 4.0)]:
        ax.add_patch(patches.Rectangle((x, y), w, 0.48, facecolor=BLUE, edgecolor=BLUE_DARK, linewidth=0.6))
    ax.text(10.5, 32.3, "...", ha="center", fontsize=10)
    ax.plot([6.9, 4.0, 4.0, 17.0, 17.0, 13.0], [26.4, 23.6, 20.0, 20.0, 23.6, 26.4], color=LINE, linestyle="--", linewidth=0.8)
    add_round_box(ax, (3.4, 16.4), 14.2, 6.4, edge=LINE, face=GRAY, lw=0.8, radius=0.45)
    ax.plot([5.5, 15.6], [20.1, 20.1], color=INK, linewidth=1.0)
    ax.add_patch(patches.Rectangle((8.9, 19.7), 3.2, 0.75, facecolor=BLUE, edgecolor=BLUE_DARK, linewidth=0.6))
    ax.text(10.5, 18.9, "Example interval", ha="center", fontsize=7.2, color=BLUE_DARK, fontweight="bold")
    ax.text(10.5, 17.4, "chr1: 12,345-12,360", ha="center", fontsize=7.1)
    ax.text(10.5, 15.1, "Strong footprint candidates may harbor\nshort regulatory motifs.", ha="center", va="top", fontsize=7.1)


def draw_sequence_panel(ax):
    add_round_box(ax, (22.4, 9), 17.3, 35.5, face="white")
    ax.text(31.1, 42.5, "Flank: +/- 75 bp", ha="center", fontsize=8.2, color=BLUE_DARK, fontweight="bold")
    ax.plot([24.3, 37.5], [39.8, 39.8], color=INK, linewidth=0.9)
    ax.plot([31.4, 31.4], [38.7, 36.7], color=LINE, linewidth=1.7)
    add_arrow(ax, 31.4, 36.5, 31.4, 34.5, color=LINE)
    ax.add_patch(patches.Rectangle((29.0, 39.4), 4.0, 0.65, facecolor=BLUE, edgecolor=BLUE_DARK, linewidth=0.6))
    add_round_box(ax, (23.2, 15.2), 15.7, 17.8, edge=BLUE, face=BLUE_LIGHT, lw=1.0, radius=0.55)
    ax.text(31.1, 31.1, "FASTA output", ha="center", fontsize=8.2, color=BLUE_DARK, fontweight="bold")
    rows = [
        (">chr1:12345-12495(+)", "GATCTGACCTAGGCTACGATC"),
        (">chr1:12531-12681(-)", "TGGACATCCTAGCTGACCTAA"),
        (">chr1:7321-7471(+)", "AACCTGACCTAGTTGACCTGC"),
    ]
    y = 28.7
    for header, seq in rows:
        ax.text(24.1, y, header, ha="left", fontsize=6.7, family="monospace", color=INK)
        ax.text(24.1, y - 1.6, seq, ha="left", fontsize=6.7, family="monospace", color=MUTED)
        add_motif(ax, 28.3, y - 1.6, "TGACCTA", size=6.7, spacing=0.38)
        y -= 5.0
    ax.text(31.1, 13.1, "Candidate-centered sequences\nready for external discovery.", ha="center", va="top", fontsize=7.1)


def draw_discovery_panel(ax):
    add_round_box(ax, (42.0, 30.6), 32.8, 13.4, edge=BLUE, face=BLUE_LIGHT, lw=1.0)
    ax.add_patch(patches.Circle((51.0, 42.4), 0.9, facecolor=BLUE, edgecolor=BLUE_DARK, linewidth=0.7))
    ax.text(51.0, 42.35, "A", ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    ax.text(52.8, 42.5, "De novo-only discovery", ha="left", va="center", fontsize=8.2, color=BLUE_DARK, fontweight="bold")
    add_file_icon(ax, 44.0, 35.2, w=2.6, h=3.3)
    ax.text(45.3, 39.5, "Input FASTA", ha="center", fontsize=6.6, color=BLUE_DARK, fontweight="bold")
    add_arrow(ax, 48.0, 36.9, 52.0, 36.9)
    ax.text(56.0, 39.9, "STREME / MEME / DREME", ha="center", fontsize=6.8, color=BLUE_DARK, fontweight="bold")
    add_motif(ax, 55.7, 37.3, "TGACCA", size=11, spacing=0.78)
    ax.text(56.0, 35.2, "discovered motifs", ha="center", fontsize=6.3)
    ax.plot([62.6, 64.2, 64.2, 62.6], [39.4, 39.4, 34.8, 34.8], color=LINE, linewidth=0.9)
    add_arrow(ax, 64.6, 37.1, 67.2, 37.1, color=LINE)
    add_round_box(ax, (67.5, 32.2), 6.3, 8.4, edge=LINE, face="white", lw=0.8)
    ax.text(70.65, 39.2, "Record", ha="center", fontsize=6.7, color=INK, fontweight="bold")
    for i, item in enumerate(["commands", "parameters", "versions", "input FASTA"]):
        ax.text(68.3, 37.8 - i * 1.25, f"- {item}", ha="left", fontsize=6.1)

    add_round_box(ax, (42.0, 13.6), 32.8, 15.7, edge=GREEN, face=GREEN_LIGHT, lw=1.0)
    ax.add_patch(patches.Circle((45.8, 27.4), 0.9, facecolor=GREEN, edgecolor=GREEN_DARK, linewidth=0.7))
    ax.text(45.8, 27.35, "B", ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    ax.text(47.4, 27.5, "Known database + de novo supplement", ha="left", va="center", fontsize=8.0, color=GREEN_DARK, fontweight="bold")
    add_file_icon(ax, 44.0, 20.6, w=2.3, h=3.0)
    ax.text(45.2, 24.4, "Input FASTA", ha="center", fontsize=6.2, color=BLUE_DARK, fontweight="bold")
    add_arrow(ax, 47.2, 22.1, 50.4, 22.1)
    ax.text(54.4, 24.6, "STREME", ha="center", fontsize=6.6, color=BLUE_DARK, fontweight="bold")
    add_motif(ax, 54.2, 22.3, "TGACCA", size=10.5, spacing=0.72)
    add_database(ax, 50.7, 16.6, w=2.7, h=2.4)
    ax.text(52.0, 15.4, "JASPAR2026", ha="center", fontsize=6.2, color=GREEN_DARK, fontweight="bold")
    ax.plot([59.2, 61.0, 61.0, 59.2], [24.0, 24.0, 19.2, 19.2], color=LINE, linewidth=0.9)
    add_arrow(ax, 61.1, 21.7, 64.0, 21.7, color=LINE)
    ax.text(66.2, 24.1, "Tomtom", ha="center", fontsize=6.6, color=GREEN_DARK, fontweight="bold")
    ax.text(66.2, 22.9, "motif comparison", ha="center", fontsize=5.8)
    add_motif(ax, 66.1, 21.0, "ACTGAC", size=10.5, spacing=0.72)
    add_round_box(ax, (63.4, 14.7), 9.7, 3.0, edge=GREEN, face="white", lw=0.8, radius=0.4)
    ax.text(68.2, 17.0, "Optional merged motif set", ha="center", fontsize=5.9, color=GREEN_DARK, fontweight="bold")
    add_motif(ax, 65.6, 15.5, "ACTGAC", size=7.8, spacing=0.55)
    add_motif(ax, 69.9, 15.5, "TGACCA", size=7.8, spacing=0.55)
    add_round_box(ax, (42.0, 8.3), 32.8, 3.8, edge=BLUE, face=GRAY, lw=0.9, radius=0.45)
    ax.text(43.1, 10.7, ">_", ha="left", va="center", fontsize=9, fontweight="bold")
    ax.text(46.0, 10.6, "Saved run script: commands, parameters, versions, inputs", ha="left", va="center", fontsize=6.8, color=BLUE_DARK, fontweight="bold")
    ax.text(46.0, 9.3, "Example run: STREME 5.5.9 + Tomtom 5.5.9 against JASPAR2026", ha="left", va="center", fontsize=6.2, color=INK)


def draw_output_panel(ax):
    add_round_box(ax, (77.0, 12.7), 21.5, 31.4, face="white")
    ax.text(87.8, 42.7, "Motif summary reports", ha="center", fontsize=8.2, color=BLUE_DARK, fontweight="bold")
    ax.text(87.8, 41.1, "TSV/HTML plus optional merged motif set", ha="center", fontsize=6.6)
    for i, (seq, evalue) in enumerate([("ACTGAC", "1.2e-12"), ("TGACCA", "3.4e-08"), ("CCTGAG", "2.1e-06")]):
        add_motif(ax, 82.2, 38.1 - i * 2.3, seq, size=10.5, spacing=0.72)
        ax.text(94.0, 38.1 - i * 2.3, evalue, ha="right", va="center", fontsize=6.4)
    add_round_box(ax, (78.1, 25.3), 19.2, 8.7, edge=LINE, face=GRAY, lw=0.8, radius=0.45)
    ax.text(87.7, 32.6, "Tomtom matches to known motifs", ha="center", fontsize=6.7, color=BLUE_DARK, fontweight="bold")
    headers = ["query", "best match", "q"]
    xs = [80.0, 86.6, 94.7]
    for x, h in zip(xs, headers):
        ax.text(x, 30.8, h, ha="center", fontsize=5.8, color=INK, fontweight="bold")
    rows = [("ACTGAC", "JUNB", "1.2e-7"), ("TGACCA", "AP-1", "4.6e-6"), ("CCTGAG", "NF-kB", "1.3e-5")]
    for i, (query, match, qval) in enumerate(rows):
        y = 29.3 - i * 1.55
        ax.text(xs[0], y, query, ha="center", fontsize=5.8, family="monospace")
        ax.text(xs[1], y, match, ha="center", fontsize=5.8)
        ax.text(xs[2], y, qval, ha="center", fontsize=5.8)
    add_round_box(ax, (78.1, 18.2), 19.2, 5.6, edge=BLUE, face=BLUE_LIGHT, lw=0.8, radius=0.45)
    ax.text(87.7, 22.5, "Run summary", ha="center", fontsize=6.8, color=BLUE_DARK, fontweight="bold")
    for i, item in enumerate(["candidate FASTA", "tool versions", "command plan", "motif_summary.tsv"]):
        ax.text(79.0 + (i // 2) * 9.0, 21.2 - (i % 2) * 1.55, f"check  {item}", ha="left", fontsize=5.9, color=INK)
    add_round_box(ax, (81.0, 7.4), 14.5, 3.9, edge=GREEN, face=GREEN_LIGHT, lw=1.0, radius=0.5)
    ax.text(88.2, 9.35, "Feeds into differential\nfootprint analysis", ha="center", va="center", fontsize=7.0, color=GREEN_DARK, fontweight="bold")
    add_arrow(ax, 87.8, 12.7, 87.8, 11.4, color=GREEN_DARK, lw=1.4)


def draw_legend(ax):
    add_round_box(ax, (3.2, 1.6), 73.5, 4.0, edge=LINE, face="white", lw=0.8, radius=0.45)
    ax.add_patch(patches.Rectangle((5.0, 3.25), 2.8, 0.35, facecolor=BLUE, edgecolor=BLUE_DARK, linewidth=0.6))
    ax.text(8.7, 3.42, "Candidate interval", fontsize=5.9, va="center")
    add_motif(ax, 19.8, 3.42, "ACTGAC", size=7.5, spacing=0.52)
    ax.text(24.4, 3.42, "Motif pattern", fontsize=5.9, va="center")
    add_database(ax, 33.0, 2.6, w=1.9, h=1.5)
    ax.text(35.8, 3.42, "Known database", fontsize=5.9, va="center")
    add_file_icon(ax, 47.5, 2.45, w=1.5, h=1.9)
    ax.text(49.7, 3.42, "FASTA input", fontsize=5.9, va="center")
    ax.text(57.2, 3.42, ">_", fontsize=7.5, va="center", fontweight="bold")
    ax.text(59.2, 3.42, "Saved command plan and provenance", fontsize=5.9, va="center")


def main() -> int:
    fig, ax = plt.subplots(figsize=(13.25, 7.45))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 55)
    ax.axis("off")

    ax.text(50, 52.6, "De novo motif discovery preparation", ha="center", va="center", fontsize=25, fontweight="bold", color=INK)
    ax.text(
        50,
        49.9,
        "Export candidate-centered sequences and prepare reproducible external motif-discovery runs.",
        ha="center",
        va="center",
        fontsize=11.5,
        color=INK,
        style="italic",
    )
    add_step(ax, 5.0, 46.5, 1, "Candidate sites", "Footprint-derived candidate intervals")
    add_step(ax, 25.8, 46.5, 2, "Sequence export", "Extract candidate-centered FASTA")
    add_step(ax, 47.0, 46.5, 3, "External motif discovery", "Established tools + reproducible records")
    add_step(ax, 83.0, 46.5, 4, "Output", "Motif summaries for downstream analysis")

    draw_candidate_panel(ax)
    draw_sequence_panel(ax)
    draw_discovery_panel(ax)
    draw_output_panel(ax)
    draw_legend(ax)
    add_arrow(ax, 19.8, 29.6, 21.9, 29.6)
    add_arrow(ax, 39.8, 29.6, 41.5, 29.6)
    add_arrow(ax, 75.2, 29.6, 76.7, 29.6)

    OUT_PREFIX.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf", ".svg"):
        fig.savefig(OUT_PREFIX.with_suffix(suffix), dpi=450, bbox_inches="tight")
    fig.savefig((REPO_ROOT / "manuscript/figures/fp-tools-de-novo-motif.png"), dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT_PREFIX.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
