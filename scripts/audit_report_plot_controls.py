#!/usr/bin/env python3
"""Exercise shared report plot controls in a real browser."""

from __future__ import annotations

import argparse
import functools
import http.server
import re
import tempfile
import threading
import xml.etree.ElementTree as ET
from contextlib import contextmanager
from pathlib import Path
from urllib.parse import quote

from playwright.sync_api import expect, sync_playwright

from fp_tools.tools.review_multi_comparisons import write_review_html


def fixture_payload(include_aggregates: bool) -> dict:
    aggregate_motifs = []
    if include_aggregates:
        aggregate_motifs = [
            {
                "prefix": "JUN_M1",
                "name": "JUN",
                "motif_id": "M1",
                "conditions": [
                    {
                        "name": "KD",
                        "samples": [
                            {"name": "KD_1", "profile": [0.1, 0.2, 0.1]},
                            {"name": "KD_2", "profile": [0.12, 0.23, 0.11]},
                        ],
                    },
                    {
                        "name": "P",
                        "samples": [
                            {"name": "P_1", "profile": [0.2, 0.1, 0.2]},
                            {"name": "P_2", "profile": [0.21, 0.12, 0.22]},
                        ],
                    },
                ],
            }
        ]
    return {
        "conditions": ["KD", "P"],
        "groups": ["KD_up", "P_up", "n.s."],
        "colors": {"KD_up": "#dc2626", "P_up": "#2563eb", "n.s.": "#8a94a6"},
        "change_label": "Differential footprint score",
        "points": [
            {
                "prefix": "JUN_M1",
                "name": "JUN",
                "motif_id": "M1",
                "group": "KD_up",
                "change": 0.8,
                "pvalue": 1e-10,
                "fdr": 1e-8,
                "neglog10p": 10.0,
            },
            {
                "prefix": "FOS_M2",
                "name": "FOS",
                "motif_id": "M2",
                "group": "P_up",
                "change": -0.8,
                "pvalue": 1e-12,
                "fdr": 1e-10,
                "neglog10p": 12.0,
            },
        ],
        "motif_matrices": {},
        "logos": {},
        "aggregate": {"x": [-1, 0, 1], "motifs": aggregate_motifs},
    }


def write_fixtures(output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    reports = []
    for name, include_aggregates in (
        ("with-aggregates.html", True),
        ("aggregate-free.html", False),
    ):
        report = output_dir / name
        write_review_html(
            {
                "schema": "fp-tools.review-multi-comparisons.v1",
                "title": "Plot-control browser audit",
                "comparisons": [
                    {"label": "KD vs P", "payload": fixture_payload(include_aggregates)}
                ],
            },
            report,
        )
        reports.append(report)
    return reports


class QuietRequestHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format: str, *args: object) -> None:
        return


@contextmanager
def serve_report(report: Path):
    handler = functools.partial(QuietRequestHandler, directory=str(report.parent))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        port = server.server_address[1]
        yield f"http://127.0.0.1:{port}/{quote(report.name)}"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


def click_control(page, locator) -> None:
    target = locator
    if locator.get_attribute("role") == "switch":
        target = locator.locator("xpath=ancestor::label[1]")
    target.evaluate("element => element.scrollIntoView({block: 'center', inline: 'nearest'})")
    expect(target).to_be_visible()
    box = locator.bounding_box()
    if not box:
        raise AssertionError("Report control has no clickable browser box")
    x = box["x"] + box["width"] / 2
    y = box["y"] + box["height"] / 2
    hit = locator.evaluate(
        """(element, point) => {
          const target = document.elementFromPoint(point.x, point.y);
          return Boolean(target && (target === element || element.contains(target) || target.closest('label')?.contains(element)));
        }""",
        {"x": x, "y": y},
    )
    if not hit:
        raise AssertionError("Report control center is clipped or covered")
    page.mouse.click(x, y)


def downloaded_svg(page, selector: str) -> str:
    control = page.locator(selector)
    with page.expect_download() as pending:
        click_control(page, control)
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
    expected_cells = 3 if page.locator("#download-aggregate").is_visible() else 2
    if len(cells) != expected_cells:
        raise AssertionError(
            f"{report}: combined export contains {len(cells)} audited cells; "
            f"expected {expected_cells}"
        )
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


def audit_report(page, report: Path, url: str, screenshot: Path | None = None) -> None:
    errors: list[str] = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.on(
        "console",
        lambda message: errors.append(message.text) if message.type == "error" else None,
    )
    page.goto(url, wait_until="load")
    page.wait_for_selector("#rank-sort-toggle")
    page.wait_for_function(
        "document.querySelectorAll('svg.rank-svg,#rank-chart .rank-bar').length > 0"
    )

    comparison_control = page.locator("#comparison-selector")
    if page.locator("#comparison-selector-control").is_visible():
        comparison_control.locator("option").first.wait_for(state="attached", timeout=60_000)
        expected_comparison = comparison_control.locator("option:checked").inner_text()
    else:
        first = page.locator("#condition-1 option:checked").inner_text()
        second = page.locator("#condition-2 option:checked").inner_text()
        expected_comparison = f"{first} vs {second}"
    volcano_text = page.locator("svg.volcano-svg,#chart").first.text_content() or ""
    if expected_comparison not in volcano_text:
        raise AssertionError(
            f"{report}: exported volcano lacks comparison label "
            f"{expected_comparison!r}"
        )

    toggle = page.locator("#rank-sort-toggle")
    click_box = toggle.bounding_box()
    if not click_box or click_box["height"] < 12:
        raise AssertionError(f"{report}: rank switch is clipped: {click_box}")
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
    click_control(page, toggle)
    expect(toggle).to_be_checked()
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

    click_control(page, toggle)
    expect(toggle).not_to_be_checked()
    if screenshot:
        screenshot.parent.mkdir(parents=True, exist_ok=True)
        page.screenshot(path=str(screenshot), full_page=True)
    if errors:
        raise AssertionError(f"{report}: browser errors: {' | '.join(errors)}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("html", nargs="*", type=Path)
    parser.add_argument(
        "--fixture-dir",
        type=Path,
        help="Generate and audit deterministic aggregate and aggregate-free reports here.",
    )
    parser.add_argument("--screenshot-dir", type=Path)
    args = parser.parse_args()
    reports = [path.resolve() for path in args.html]
    temporary_fixture = None
    if args.fixture_dir:
        reports.extend(write_fixtures(args.fixture_dir.resolve()))
    elif not reports:
        temporary_fixture = tempfile.TemporaryDirectory(prefix="fp-tools-report-audit-")
        reports.extend(write_fixtures(Path(temporary_fixture.name)))
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        try:
            for report in reports:
                page = browser.new_page(viewport={"width": 1800, "height": 1050})
                screenshot = (
                    args.screenshot_dir / f"{report.stem}_plot_controls.png"
                    if args.screenshot_dir
                    else None
                )
                with serve_report(report) as url:
                    audit_report(page, report, url, screenshot)
                page.close()
                print(f"PASS {report}")
        finally:
            browser.close()
    if temporary_fixture is not None:
        temporary_fixture.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
