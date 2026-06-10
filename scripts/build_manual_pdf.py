#!/usr/bin/env python3

from __future__ import annotations

import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from textwrap import wrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch, Rectangle


PAGE_WIDTH = 8.27
PAGE_HEIGHT = 11.69
LEFT = 0.07
RIGHT = 0.93
TOP = 0.93
BOTTOM = 0.06
TEXT_WIDTH = RIGHT - LEFT

FONT_SANS = "DejaVu Sans"
FONT_MONO = "DejaVu Sans Mono"

ACCENT = "#1f5aa6"
TEXT = "#1e293b"
MUTED = "#475569"
CODE_BG = "#f1f5f9"
HEADER_BG = "#e8f0fb"
RULE = "#cbd5e1"


@dataclass
class Block:
    kind: str
    text: str


def parse_blocks(markdown: str) -> list[Block]:
    blocks: list[Block] = []
    lines = markdown.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if not stripped:
            i += 1
            continue

        if stripped.startswith("```"):
            code_lines = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("```"):
                code_lines.append(lines[i].rstrip())
                i += 1
            blocks.append(Block("code", "\n".join(code_lines)))
            i += 1
            continue

        if stripped.startswith("# "):
            blocks.append(Block("h1", stripped[2:].strip()))
            i += 1
            continue

        if stripped.startswith("## "):
            blocks.append(Block("h2", stripped[3:].strip()))
            i += 1
            continue

        if stripped.startswith("### "):
            blocks.append(Block("h3", stripped[4:].strip()))
            i += 1
            continue

        if re.match(r"^\d+\.\s", stripped):
            items = []
            while i < len(lines) and re.match(r"^\d+\.\s", lines[i].strip()):
                items.append(lines[i].strip())
                i += 1
            blocks.append(Block("olist", "\n".join(items)))
            continue

        if stripped.startswith("- "):
            items = []
            while i < len(lines) and lines[i].strip().startswith("- "):
                items.append(lines[i].strip()[2:])
                i += 1
            blocks.append(Block("ulist", "\n".join(items)))
            continue

        para = [stripped]
        i += 1
        while i < len(lines):
            nxt = lines[i].strip()
            if not nxt or nxt.startswith("#") or nxt.startswith("```") or nxt.startswith("- ") or re.match(r"^\d+\.\s", nxt):
                break
            para.append(nxt)
            i += 1
        blocks.append(Block("p", " ".join(para)))

    return blocks


def sanitize_inline(text: str) -> str:
    return text.replace("`", "")


def wrapped_lines(text: str, width: int) -> list[str]:
    return wrap(
        sanitize_inline(text),
        width=width,
        replace_whitespace=False,
        drop_whitespace=False,
    ) or [""]


def estimate_height(block: Block) -> float:
    if block.kind == "h1":
        return 0.10
    if block.kind == "h2":
        return 0.055
    if block.kind == "h3":
        return 0.042
    if block.kind == "p":
        return 0.015 * len(wrapped_lines(block.text, 82)) + 0.014
    if block.kind == "ulist":
        lines = sum(len(wrapped_lines(item, 76)) for item in block.text.splitlines())
        return 0.015 * lines + 0.016
    if block.kind == "olist":
        lines = 0
        for item in block.text.splitlines():
            body = re.sub(r"^\d+\.\s*", "", item)
            lines += len(wrapped_lines(body, 74))
        return 0.015 * lines + 0.016
    if block.kind == "code":
        code_lines = block.text.splitlines() or [""]
        wrapped = 0
        for line in code_lines:
            wrapped += max(1, math.ceil(max(1, len(line)) / 68))
        return 0.0145 * wrapped + 0.05
    return 0.03


def split_pages(blocks: list[Block]) -> list[list[Block]]:
    pages: list[list[Block]] = [[]]
    y = TOP
    for block in blocks:
        h = estimate_height(block)
        if pages[-1] and y - h < BOTTOM:
            pages.append([])
            y = TOP
        pages[-1].append(block)
        y -= h
    return pages


def draw_header(ax, page_no: int, total_pages: int, title: str) -> None:
    ax.add_patch(Rectangle((0, 0.94), 1, 0.06, color=ACCENT, transform=ax.transAxes))
    ax.text(LEFT, 0.967, title, color="white", fontsize=18, fontweight="bold", family=FONT_SANS, va="center")
    ax.text(RIGHT, 0.967, f"{page_no}/{total_pages}", color="white", fontsize=10, family=FONT_SANS, va="center", ha="right")
    ax.plot([LEFT, RIGHT], [0.94, 0.94], color=RULE, lw=0.8)


def draw_footer(ax) -> None:
    ax.plot([LEFT, RIGHT], [BOTTOM - 0.012, BOTTOM - 0.012], color=RULE, lw=0.8)
    ax.text(LEFT, BOTTOM - 0.028, "fp-tools manual", fontsize=8.5, color=MUTED, family=FONT_SANS, va="top")


def draw_paragraph(ax, y: float, text: str) -> float:
    for line in wrapped_lines(text, 82):
        ax.text(LEFT, y, line, fontsize=10.5, color=TEXT, family=FONT_SANS, va="top")
        y -= 0.015
    return y - 0.004


def draw_list(ax, y: float, items: list[str], ordered: bool) -> float:
    for idx, item in enumerate(items, start=1):
        prefix = f"{idx}." if ordered else u"\u2022"
        body = re.sub(r"^\d+\.\s*", "", item) if ordered else item
        lines = wrapped_lines(body, 74)
        ax.text(LEFT, y, prefix, fontsize=10.5, color=ACCENT, family=FONT_SANS, va="top", fontweight="bold")
        for line_idx, line in enumerate(lines):
            ax.text(LEFT + 0.03, y - 0.015 * line_idx, line, fontsize=10.5, color=TEXT, family=FONT_SANS, va="top")
        y -= 0.015 * len(lines) + 0.004
    return y - 0.004


def draw_code(ax, y: float, text: str) -> float:
    code_lines = text.splitlines() or [""]
    rendered: list[str] = []
    for line in code_lines:
        rendered.extend(wrap(line, width=68, replace_whitespace=False, drop_whitespace=False) or [""])
    height = 0.0145 * len(rendered) + 0.028
    box_y = y - height + 0.006
    ax.add_patch(
        FancyBboxPatch(
            (LEFT, box_y),
            TEXT_WIDTH,
            height,
            boxstyle="round,pad=0.008,rounding_size=0.01",
            linewidth=0.8,
            edgecolor=RULE,
            facecolor=CODE_BG,
            transform=ax.transAxes,
        )
    )
    text_y = y - 0.008
    for line in rendered:
        ax.text(LEFT + 0.015, text_y, line, fontsize=9.2, color=TEXT, family=FONT_MONO, va="top")
        text_y -= 0.0145
    return box_y - 0.012


def draw_block(ax, y: float, block: Block) -> float:
    if block.kind == "h1":
        ax.add_patch(
            FancyBboxPatch(
                (LEFT, y - 0.055),
                TEXT_WIDTH,
                0.05,
                boxstyle="round,pad=0.01,rounding_size=0.015",
                linewidth=0,
                facecolor=HEADER_BG,
                transform=ax.transAxes,
            )
        )
        ax.text(LEFT + 0.015, y - 0.016, sanitize_inline(block.text), fontsize=22, color=ACCENT, family=FONT_SANS, fontweight="bold", va="top")
        return y - 0.075

    if block.kind == "h2":
        ax.text(LEFT, y, sanitize_inline(block.text), fontsize=15, color=ACCENT, family=FONT_SANS, fontweight="bold", va="top")
        ax.plot([LEFT, RIGHT], [y - 0.022, y - 0.022], color=RULE, lw=0.8)
        return y - 0.038

    if block.kind == "h3":
        ax.text(LEFT, y, sanitize_inline(block.text), fontsize=12.2, color=TEXT, family=FONT_SANS, fontweight="bold", va="top")
        return y - 0.028

    if block.kind == "p":
        return draw_paragraph(ax, y, block.text)

    if block.kind == "ulist":
        return draw_list(ax, y, block.text.splitlines(), ordered=False)

    if block.kind == "olist":
        return draw_list(ax, y, block.text.splitlines(), ordered=True)

    if block.kind == "code":
        return draw_code(ax, y, block.text)

    return y


def build_pdf(src: Path, out: Path) -> None:
    blocks = parse_blocks(src.read_text())
    pages = split_pages(blocks)

    with PdfPages(out) as pdf:
        total = len(pages)
        title = sanitize_inline(next((b.text for b in blocks if b.kind == "h1"), src.stem))
        for page_no, blocks_on_page in enumerate(pages, start=1):
            fig = plt.figure(figsize=(PAGE_WIDTH, PAGE_HEIGHT))
            ax = fig.add_axes([0, 0, 1, 1])
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis("off")

            draw_header(ax, page_no, total, title)
            y = TOP
            for block in blocks_on_page:
                y = draw_block(ax, y, block)

            draw_footer(ax)
            pdf.savefig(fig)
            plt.close(fig)


def main() -> int:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("MANUAL.md")
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("MANUAL.pdf")
    build_pdf(src, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
