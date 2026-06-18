#!/usr/bin/env python
"""Draw the pseudobulk fragment workflow figure."""

from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


BLUE = "#0B45D8"
INK = "#111827"
MUTED = "#4B5563"
GREEN = "#059669"
ORANGE = "#EA580C"
PURPLE = "#6D28D9"
GRAY = "#6B7280"
LIGHT_BLUE = "#EFF6FF"
LIGHT_GRAY = "#F9FAFB"


def add_box(ax, x, y, w, h, edge=BLUE, face="white", lw=0.8, radius=0.7):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.25,rounding_size={radius}",
        linewidth=lw,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(box)
    return box


def add_text(ax, x, y, text, size=6.5, color=INK, weight="normal", ha="left", va="center", width=None, **kwargs):
    if width:
        text = "\n".join(textwrap.wrap(text, width=width))
    return ax.text(x, y, text, fontsize=size, color=color, weight=weight, ha=ha, va=va, **kwargs)


def add_arrow(ax, x1, y1, x2, y2, color=INK, lw=0.8):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=7, lw=lw, color=color))


def add_cell_cluster(ax, x, y, color, n=5, scale=0.55):
    offsets = [(0, 0), (1.0, 0.15), (2.0, 0), (0.55, -0.85), (1.55, -0.8), (2.55, -0.7)]
    for dx, dy in offsets[:n]:
        circ = Circle((x + dx * scale, y + dy * scale), 0.33 * scale, facecolor=color, edgecolor="white", lw=0.45)
        ax.add_patch(circ)
        ax.add_patch(Circle((x + dx * scale, y + dy * scale), 0.20 * scale, facecolor="white", edgecolor=color, lw=0.45, alpha=0.8))


def add_file_icon(ax, x, y, label, color=BLUE):
    ax.add_patch(Rectangle((x, y), 1.5, 2.0, facecolor="white", edgecolor=INK, lw=0.7))
    ax.plot([x + 0.35, x + 1.15], [y + 0.65, y + 0.65], color=color, lw=0.8)
    ax.plot([x + 0.35, x + 1.05], [y + 1.05, y + 1.05], color=color, lw=0.8)
    add_text(ax, x + 0.75, y - 0.35, label, size=4.1, ha="center", va="top", width=9)


def add_track_icon(ax, x, y, label, color=BLUE):
    heights = [0.4, 1.0, 1.7, 0.8, 1.25, 0.55]
    for i, h in enumerate(heights):
        ax.add_patch(Rectangle((x + i * 0.22, y), 0.12, h, facecolor=color, edgecolor=color, lw=0))
    ax.plot([x - 0.1, x + 1.45], [y, y], color=INK, lw=0.45)
    add_text(ax, x + 0.65, y - 0.35, label, size=4.1, ha="center", va="top", width=10)


def add_tiny_umap(ax, x, y):
    colors = [(BLUE, 0.15, 0.25), (ORANGE, 1.1, 0.6), (GREEN, 0.65, 1.25)]
    for color, cx, cy in colors:
        for i in range(9):
            ax.add_patch(Circle((x + cx + (i % 3) * 0.16, y + cy + (i // 3) * 0.13), 0.045, facecolor=color, edgecolor=color, lw=0))


def add_step_header(ax, x, y, number, title, subtitle):
    ax.add_patch(Circle((x, y), 1.35, facecolor="white", edgecolor=BLUE, lw=1.0))
    add_text(ax, x, y, str(number), size=9.5, color=BLUE, weight="bold", ha="center")
    add_text(ax, x + 2.2, y + 0.1, title, size=7.0, color=BLUE, weight="bold", width=13)


def draw_workflow(out_prefix: Path) -> None:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.6, 4.65))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 62)
    ax.axis("off")

    add_text(ax, 50, 59.4, "Pseudobulk fragments", size=17.5, weight="bold", ha="center")
    add_text(
        ax,
        50,
        56.8,
        "Group public 10x Genomics PBMC single-cell ATAC fragments into corrected pseudobulk footprint outputs.",
        size=6.5,
        color=INK,
        ha="center",
    )

    starts = [2, 26.5, 51, 75.5]
    panel_w = 22
    panel_y = 10
    panel_h = 38.5
    headers = [
        ("Single-cell\nATAC input", ""),
        ("Group by\ncell label", ""),
        ("Aggregate\nand QC", ""),
        ("Corrected\noutputs", ""),
    ]
    for i, (x, (title, subtitle)) in enumerate(zip(starts, headers, strict=True), start=1):
        add_step_header(ax, x + 2.0, 52.8, i, title, subtitle)
        add_box(ax, x, panel_y, panel_w, panel_h, edge=BLUE, face="white", lw=0.75, radius=0.85)
        if i < 4:
            add_arrow(ax, x + panel_w + 0.6, 52.8, starts[i] - 1.2, 52.8)

    # Panel 1: input tables.
    x = starts[0]
    add_text(ax, x + panel_w / 2, 46.2, "Fragments.tsv.gz", size=6.4, color=BLUE, weight="bold", ha="center")
    table_x, table_y = x + 2.0, 32.0
    col_w = [4.1, 4.2, 4.2, 5.6]
    headers1 = ["chrom", "start", "end", "barcode"]
    rows1 = [
        ["chr1", "100.6k", "100.7k", "AAAC..."],
        ["chr1", "105.2k", "105.4k", "AAAC..."],
        ["chr2", "150.1k", "150.3k", "AAGT..."],
        ["chr3", "75.9k", "76.1k", "CTTG..."],
    ]
    for j, header in enumerate(headers1):
        ax.add_patch(Rectangle((table_x + sum(col_w[:j]), table_y + 8.0), col_w[j], 2.3, facecolor=LIGHT_GRAY, edgecolor="#9CA3AF", lw=0.35))
        add_text(ax, table_x + sum(col_w[:j]) + col_w[j] / 2, table_y + 9.15, header, size=4.4, weight="bold", ha="center")
    for r, row in enumerate(rows1):
        for j, value in enumerate(row):
            ax.add_patch(Rectangle((table_x + sum(col_w[:j]), table_y + 5.7 - r * 2.1), col_w[j], 2.1, facecolor="white", edgecolor="#D1D5DB", lw=0.3))
            add_text(ax, table_x + sum(col_w[:j]) + col_w[j] / 2, table_y + 6.75 - r * 2.1, value, size=4.1, ha="center")
    add_text(ax, x + panel_w / 2, 26.5, "Metadata", size=6.6, color=BLUE, weight="bold", ha="center")
    meta_rows = [("barcode", "cell_type"), ("AAAC...-1", "T_NK_cell"), ("AAGT...-1", "B_cell"), ("CTTG...-1", "Monocyte")]
    for r, row in enumerate(meta_rows):
        face = LIGHT_GRAY if r == 0 else "white"
        weight = "bold" if r == 0 else "normal"
        ax.add_patch(Rectangle((x + 3.3, 16.4 + (len(meta_rows) - 1 - r) * 2.4), 7.2, 2.4, facecolor=face, edgecolor="#D1D5DB", lw=0.3))
        ax.add_patch(Rectangle((x + 10.5, 16.4 + (len(meta_rows) - 1 - r) * 2.4), 8.0, 2.4, facecolor=face, edgecolor="#D1D5DB", lw=0.3))
        add_text(ax, x + 6.9, 17.6 + (len(meta_rows) - 1 - r) * 2.4, row[0], size=4.8, ha="center", weight=weight)
        add_text(ax, x + 14.5, 17.6 + (len(meta_rows) - 1 - r) * 2.4, row[1], size=4.8, ha="center", weight=weight)

    # Panel 2: grouping.
    x = starts[1]
    groups = [("B_cell", GREEN), ("Monocyte", ORANGE), ("T_NK_cell", PURPLE)]
    display_labels = {"B_cell": "B cell", "Monocyte": "Monocyte", "T_NK_cell": "T/NK"}
    for idx, (group, color) in enumerate(groups):
        gy = 43.0 - idx * 10.0
        add_text(ax, x + 2.0, gy + 1.0, display_labels[group], size=5.7, color=color, weight="bold")
        add_cell_cluster(ax, x + 2.4, gy - 1.0, color, n=5, scale=0.95)
        add_box(ax, x + 9.4, gy - 4.0, 10.4, 7.1, edge=color, face="white", lw=0.6, radius=0.55)
        add_text(ax, x + 10.0, gy + 1.9, "barcodes", size=4.7, color=color, weight="bold")
        add_text(ax, x + 10.0, gy + 0.15, "AAAC...\nAAGT...\nCTTG...", size=4.2, color=INK)
    add_text(ax, x + 2.0, 13.0, "Each barcode contributes fragments to one retained group.", size=4.9, color=MUTED, width=36)

    # Panel 3: aggregation and QC.
    x = starts[2]
    for idx, (group, color) in enumerate(groups):
        gy = 43.2 - idx * 8.2
        add_cell_cluster(ax, x + 2.0, gy, color, n=4, scale=0.65)
        add_arrow(ax, x + 5.1, gy - 0.2, x + 7.1, gy - 0.2, color=INK)
        add_file_icon(ax, x + 7.5, gy - 1.3, "fragments", color=color)
        add_arrow(ax, x + 10.0, gy - 0.2, x + 12.0, gy - 0.2, color=INK)
        add_file_icon(ax, x + 12.4, gy - 1.3, "pseudo-BAM", color=color)
        add_text(ax, x + 16.0, gy + 0.3, display_labels[group], size=4.7, color=color, weight="bold", width=9)
    add_box(ax, x + 1.6, 12.6, 8.8, 7.0, edge="#9CA3AF", face=LIGHT_GRAY, lw=0.45, radius=0.4)
    add_text(ax, x + 2.2, 18.0, "QC filters", size=5.8, color=BLUE, weight="bold")
    add_text(ax, x + 2.2, 16.1, "min cells 300\nmin fragments 50k\nfilter low-depth", size=3.8, color=INK)
    add_box(ax, x + 11.1, 12.6, 9.0, 7.0, edge="#9CA3AF", face=LIGHT_GRAY, lw=0.45, radius=0.4)
    add_text(ax, x + 11.7, 18.0, "Options", size=5.8, color=BLUE, weight="bold")
    add_text(ax, x + 11.7, 16.1, "tabix index\nchrom filters\nbarcode cleanup", size=3.8, color=INK)

    # Panel 4: outputs.
    x = starts[3]
    add_text(ax, x + panel_w / 2, 46.2, "Manifest / QC summary", size=6.8, color=BLUE, weight="bold", ha="center")
    table_x, table_y = x + 1.3, 38.0
    headers4 = ["group", "status", "cells", "outputs"]
    widths4 = [5.2, 4.0, 4.2, 7.5]
    for j, header in enumerate(headers4):
        ax.add_patch(Rectangle((table_x + sum(widths4[:j]), table_y + 4.8), widths4[j], 2.1, facecolor=LIGHT_GRAY, edgecolor="#9CA3AF", lw=0.3))
        add_text(ax, table_x + sum(widths4[:j]) + widths4[j] / 2, table_y + 5.85, header, size=4.4, weight="bold", ha="center")
    rows4 = [("B_cell", "pass", "1,102", "bam/bw/tsv"), ("Monocyte", "pass", "1,045", "bam/bw/tsv"), ("other", "filtered", "112", "-")]
    for r, row in enumerate(rows4):
        for j, value in enumerate(row):
            ax.add_patch(Rectangle((table_x + sum(widths4[:j]), table_y + 2.7 - r * 2.1), widths4[j], 2.1, facecolor="white", edgecolor="#D1D5DB", lw=0.3))
            color = GREEN if value == "pass" else (GRAY if value == "filtered" else INK)
            add_text(ax, table_x + sum(widths4[:j]) + widths4[j] / 2, table_y + 3.75 - r * 2.1, value, size=4.4, ha="center", color=color, weight="bold" if j == 1 else "normal")
    add_box(ax, x + 1.2, 14.0, 20.0, 19.0, edge=BLUE, face=LIGHT_BLUE, lw=0.55, radius=0.55)
    add_text(ax, x + 11.2, 31.0, "Primary fp-tools outputs", size=6.3, color=BLUE, weight="bold", ha="center")
    output_items = [
        ("corrected\ncut-site bw", "track"),
        ("footprint-score\nbigWig", "track"),
        ("diff.\nreport", "file"),
        ("volcano +\nKNN UMAP", "umap"),
    ]
    for i, (label, kind) in enumerate(output_items):
        ox = x + 3.2 + (i % 2) * 9.0
        oy = 23.5 - (i // 2) * 6.8
        add_box(ax, ox - 0.8, oy - 1.4, 7.0, 5.4, edge="#BFDBFE", face="white", lw=0.45, radius=0.35)
        if kind == "track":
            add_track_icon(ax, ox + 0.4, oy + 0.6, label, color=BLUE)
        elif kind == "file":
            add_file_icon(ax, ox + 1.1, oy + 0.3, label, color=BLUE)
        else:
            add_tiny_umap(ax, ox + 0.8, oy + 0.6)
            add_text(ax, ox + 2.0, oy - 0.35, label, size=4.1, ha="center", va="top", width=10)

    # Legend.
    add_box(ax, 18, 2.1, 64, 4.6, edge="#9CA3AF", face="white", lw=0.55, radius=0.45)
    add_cell_cluster(ax, 21.0, 4.6, "#60A5FA", n=1, scale=0.9)
    add_text(ax, 24.0, 4.25, "single cell barcode", size=5.8)
    ax.plot([39, 45], [4.45, 4.45], color=BLUE, lw=1.8)
    add_text(ax, 46.5, 4.25, "pseudobulk fragment", size=5.8)
    add_track_icon(ax, 63.0, 4.0, "", color=BLUE)
    add_text(ax, 66.3, 4.25, "corrected footprint tracks", size=5.8)

    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    repo = Path(__file__).resolve().parents[2]
    manuscript_prefix = repo / "manuscript" / "figures" / "fp-tools-pseudo-bulk-hires"
    docs_prefix = repo / "docs" / "assets" / "fp-tools-pseudo-bulk"
    draw_workflow(manuscript_prefix)
    draw_workflow(docs_prefix)
    # Keep the legacy manuscript PNG name in sync for README/PDF previews.
    (repo / "manuscript" / "figures" / "fp-tools-pseudo-bulk.png").write_bytes(manuscript_prefix.with_suffix(".png").read_bytes())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
