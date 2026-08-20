#!/usr/bin/env python3
"""Exercise shared report plot controls in a real browser."""

from __future__ import annotations

import argparse
import re
import xml.etree.ElementTree as ET
from pathlib import Path

from playwright.sync_api import sync_playwright


def downloaded_svg(page, selector: str) -> str:
    with page.expect_download() as pending:
        page.locator(selector).click()
    return Path(pending.value.path()).read_text(encoding="utf-8")


def svg_elements(svg: str, local_name: str) -> list[ET.Element]:
    root = ET.fromstring(svg)
    return [node for node in root.iter() if node.tag.rsplit("}", 1)[-1] == local_name]


def audit_svg_exports(page, report: Path, expected_comparison: str) -> None:
    rank_svg = downloaded_svg(page, "#download-rank")
    if expected_comparison not in rank_svg:
        raise AssertionError(f"{report}: waterfall export lacks comparison name")

    volcano_svg = downloaded_svg(page, "#download-volcano")
    if "paint-order" in volcano_svg:
        raise AssertionError(f"{report}: volcano export contains browser-only paint-order")
    labels = [
        node
        for node in svg_elements(volcano_svg, "text")
        if "volcano-user-label" in node.get("class", "").split()
    ]
    if not labels or not any("JUN" in "".join(node.itertext()).upper() for node in labels):
        raise AssertionError(f"{report}: downloaded volcano lacks requested TF labels")
    for label in labels:
        expected = {
            "fill": "#111827",
            "stroke": "none",
            "font-family": "Helvetica,Arial,sans-serif",
        }
        if "Arial,Helvetica" in label.get("font-family", ""):
            expected["font-family"] = "Arial,Helvetica,sans-serif"
        for attribute, value in expected.items():
            if label.get(attribute) != value:
                raise AssertionError(
                    f"{report}: TF label lacks Illustrator-safe {attribute}={value!r}"
                )
    white_rects = [
        node
        for node in svg_elements(volcano_svg, "rect")
        if node.get("fill", "").lower() in {"#fff", "#ffffff", "white"}
    ]
    if len(white_rects) > 1:
        raise AssertionError(
            f"{report}: downloaded volcano has {len(white_rects)} opaque white backgrounds"
        )

    combined_svg = downloaded_svg(page, "#download-panel")
    if expected_comparison not in combined_svg:
        raise AssertionError(f"{report}: combined export lacks comparison names")
    root = ET.fromstring(combined_svg)
    cells = [
        node for node in root.iter() if node.get("data-export-cell") is not None
    ]
    for cell in cells:
        source_width = float(cell.get("data-source-width", "nan"))
        source_height = float(cell.get("data-source-height", "nan"))
        cell_width = float(cell.get("data-cell-width", "nan"))
        cell_height = float(cell.get("data-cell-height", "nan"))
        scale = float(cell.get("data-scale", "nan"))
        if source_width * scale > cell_width + 1e-6 or source_height * scale > cell_height + 1e-6:
            raise AssertionError(f"{report}: exported subplot exceeds its grid cell")
        direct_white = [
            child
            for child in list(cell)
            if child.tag.rsplit("}", 1)[-1] == "rect"
            and child.get("fill", "").lower() in {"#fff", "#ffffff", "white"}
        ]
        if direct_white:
            raise AssertionError(f"{report}: exported subplot retains an opaque background")


def audit_report(page, report: Path, screenshot: Path | None = None) -> None:
    errors: list[str] = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.on(
        "console",
        lambda message: errors.append(message.text) if message.type == "error" else None,
    )
    page.goto(report.resolve().as_uri(), wait_until="load")
    page.wait_for_selector("#rank-sort-toggle")
    page.wait_for_function(
        "document.querySelectorAll('svg.rank-svg,#rank-chart .rank-bar').length > 0"
    )

    comparison_control = page.locator(
        "[data-comparison-slot],#comparison-selector"
    ).first
    expected_comparison = comparison_control.locator("option:checked").inner_text()
    volcano_text = page.locator("svg.volcano-svg,#chart").first.text_content() or ""
    if expected_comparison not in volcano_text:
        raise AssertionError(
            f"{report}: exported volcano lacks comparison label "
            f"{expected_comparison!r}"
        )

    toggle = page.locator("#rank-sort-toggle")
    if toggle.is_checked():
        raise AssertionError(f"{report}: waterfall must default to differential score")
    initial_fills = page.locator(".rank-bar").evaluate_all(
        "nodes => nodes.slice(0, 20).map(node => node.getAttribute('fill'))"
    )
    directional = page.locator(".rank-bar").evaluate_all(
        """nodes => nodes.slice(0, 40).map(node => ({
          fill: node.getAttribute('fill'),
          title: node.querySelector('title')?.textContent || ''
        }))"""
    )
    def effect(item: dict[str, str]) -> float | None:
        match = re.search(
            r"(?:Differential footprint score|ΔFP)\s+(-?[0-9.]+)",
            item["title"],
        )
        return float(match.group(1)) if match else None

    positive_red = any(
        effect(item) is not None
        and effect(item) > 0
        and int(item["fill"][1:3], 16) > int(item["fill"][5:7], 16)
        for item in directional
        if item["fill"] and item["fill"].startswith("#")
    )
    negative_blue = any(
        effect(item) is not None
        and effect(item) < 0
        and int(item["fill"][5:7], 16) > int(item["fill"][1:3], 16)
        for item in directional
        if item["fill"] and item["fill"].startswith("#")
    )
    if not (positive_red and negative_blue):
        raise AssertionError(
            f"{report}: score-ranked bars do not preserve red/blue direction"
        )

    volcano = page.locator("svg.volcano-svg,#chart").first
    initial_box = volcano.bounding_box()
    if not initial_box or abs(initial_box["width"] - initial_box["height"]) > 2:
        raise AssertionError(f"{report}: volcano viewport is not square: {initial_box}")
    page.locator("#rank-rows").fill("40")
    page.locator("#rank-rows").dispatch_event("input")
    page.wait_for_timeout(50)
    count_control = page.locator("#panel-count,#plot-count")
    if count_control.count() and count_control.first.is_visible():
        options = count_control.first.locator("option").evaluate_all(
            "nodes => nodes.map(node => node.value)"
        )
        count_control.first.select_option(options[-1])
        page.wait_for_timeout(50)
    resized_box = page.locator("svg.volcano-svg,#chart").first.bounding_box()
    if not resized_box or any(
        abs(resized_box[key] - initial_box[key]) > 2 for key in ("width", "height")
    ):
        raise AssertionError(
            f"{report}: volcano resized after row/panel controls: "
            f"{initial_box} -> {resized_box}"
        )
    toggle.check()
    page.wait_for_timeout(50)
    rank_text = page.locator("svg.rank-svg,#rank-chart").all_text_contents()
    if not rank_text or not all("Signed -log10(p-value)" in text or "Signed −log10(p-value)" in text for text in rank_text):
        raise AssertionError(f"{report}: significance-ranked axis is missing")
    significance_fills = page.locator(".rank-bar").evaluate_all(
        "nodes => nodes.slice(0, 20).map(node => node.getAttribute('fill'))"
    )
    if initial_fills == significance_fills:
        raise AssertionError(f"{report}: reciprocal waterfall color encoding did not change")
    titles = page.locator(".rank-bar title").all_text_contents()
    if not titles or not all("log10(p-value)" in title for title in titles[:10]):
        raise AssertionError(f"{report}: waterfall tooltips do not retain both metrics")

    if page.locator("#volcano-highlight").count():
        page.locator("#volcano-highlight").select_option("none")
    else:
        page.locator("#motif-select").select_option("")
    page.wait_for_timeout(50)
    if page.locator("svg.volcano-svg .pt.selected,#chart .pt.selected").count():
        raise AssertionError(f"{report}: volcano still has selected-point highlights")

    page.locator("#volcano-labels").fill("JUN")
    page.wait_for_timeout(50)
    labels = page.locator(".volcano-user-label")
    if labels.count() < 1:
        raise AssertionError(f"{report}: TF-interest labels were not rendered")
    if not any("JUN" in text.upper() for text in labels.all_text_contents()):
        raise AssertionError(f"{report}: rendered TF labels do not match JUN")
    if "volcano-user-label" not in page.locator("svg.volcano-svg,#chart").first.evaluate(
        "node => new XMLSerializer().serializeToString(node)"
    ):
        raise AssertionError(f"{report}: TF labels are absent from serialized SVG")

    audit_svg_exports(page, report, expected_comparison)

    toggle.uncheck()
    page.wait_for_timeout(50)
    if screenshot:
        screenshot.parent.mkdir(parents=True, exist_ok=True)
        page.screenshot(path=str(screenshot), full_page=True)
    if errors:
        raise AssertionError(f"{report}: browser errors: {' | '.join(errors)}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("html", nargs="+", type=Path)
    parser.add_argument("--screenshot-dir", type=Path)
    args = parser.parse_args()
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            for report in args.html:
                page = browser.new_page(viewport={"width": 1800, "height": 1050})
                screenshot = (
                    args.screenshot_dir / f"{report.stem}_plot_controls.png"
                    if args.screenshot_dir
                    else None
                )
                audit_report(page, report, screenshot)
                page.close()
                print(f"PASS {report}")
        finally:
            browser.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
