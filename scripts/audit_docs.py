#!/usr/bin/env python3
"""Audit the built documentation with a real browser."""

from __future__ import annotations

import argparse
import contextlib
import http.server
import re
import socketserver
import tempfile
import threading
from pathlib import Path

from playwright.sync_api import TimeoutError as PlaywrightTimeoutError
from playwright.sync_api import sync_playwright


COMMAND_PAGES = (
    "prepare-atac",
    "bulk-footprinting",
    "atac-correct",
    "call-footprints",
    "match-motifs",
    "diff-footprints",
    "normalize-bigwig",
    "plot-aggregate",
    "review-multi-comparisons",
    "run-yaml-workflow",
    "fp-tools-gui",
    "discover-motifs",
    "summarize-motifs",
    "pseudobulk-fragments",
    "find-signature-fp",
    "sc-footprinting",
)
GET_STARTED_PAGES = (
    "",
    "get-started/installation/",
    "get-started/tool-overview/",
    "get-started/workflows/bulk-atac-seq/",
    "get-started/workflows/single-cell/",
    "get-started/workflows/de-novo-motif-discovery/",
    *(f"get-started/commands/{command}/" for command in COMMAND_PAGES),
    "get-started/output-examples/bulk-atac-seq/",
    "get-started/output-examples/single-cell-atac-seq/",
    "get-started/output-examples/region-set-comparison/",
)
CORE_COMMAND_SEQUENCE = (
    "get-started/commands/atac-correct/",
    "get-started/commands/call-footprints/",
    "get-started/commands/match-motifs/",
    "get-started/commands/diff-footprints/",
    "get-started/commands/normalize-bigwig/",
)
PAGES = (
    *GET_STARTED_PAGES,
    "api/",
    "gui/",
    "reports/",
    "demos/reports/diff_footprints_K562_HepG2.html",
    "demos/reports/region_set_HepG2_HNF4A_FOXA2/",
    "demos/gui/fp-tools-gui-static-demo.html",
    "ENCODE-Cancer-Cell-lines-Footprinting/",
)
VIEWPORTS = ((1440, 1000), (1280, 720), (390, 844))
REPORT_IFRAME_PAGES = {
    "reports/",
    "get-started/output-examples/bulk-atac-seq/",
    "get-started/output-examples/region-set-comparison/",
}
STANDALONE_DEMO_PAGES = {
    "demos/reports/diff_footprints_K562_HepG2.html",
    "demos/reports/region_set_HepG2_HNF4A_FOXA2/",
    "demos/gui/fp-tools-gui-static-demo.html",
    "ENCODE-Cancer-Cell-lines-Footprinting/",
}
DOCUMENTATION_RETURN_PAGES = {
    "demos/gui/fp-tools-gui-static-demo.html",
    "demos/reports/region_set_HepG2_HNF4A_FOXA2/",
    "ENCODE-Cancer-Cell-lines-Footprinting/",
}
DOCUMENTATION_RETURN_PAGES = {
    "demos/gui/fp-tools-gui-static-demo.html",
    "ENCODE-Cancer-Cell-lines-Footprinting/",
}
class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format: str, *args) -> None:
        return


@contextlib.contextmanager
def serve(directory: Path):
    def handler(*args, **kwargs):
        return QuietHandler(*args, directory=str(directory), **kwargs)

    with socketserver.TCPServer(("127.0.0.1", 0), handler) as server:
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            yield f"http://127.0.0.1:{server.server_address[1]}/"
        finally:
            server.shutdown()
            thread.join()


def check_expanded_aggregate_layout(surface, label: str, failures: list[str]) -> None:
    """Confirm that opening report controls cannot stack aggregate tiles."""
    details = surface.locator("#options")
    details.evaluate("element => { element.open = true; }")
    surface.wait_for_timeout(100)
    layout = surface.evaluate(
        """() => {
          const grid = document.querySelector("#aggregate-grid");
          const gridBox = grid.getBoundingClientRect();
          const boxes = [...grid.querySelectorAll(".aggregate-tile")].map(tile => {
            const box = tile.getBoundingClientRect();
            return {left: box.left, top: box.top, right: box.right, bottom: box.bottom};
          });
          const overlaps = [];
          for (let first = 0; first < boxes.length; first += 1) {
            for (let second = first + 1; second < boxes.length; second += 1) {
              const a = boxes[first];
              const b = boxes[second];
              const overlapX = Math.min(a.right, b.right) - Math.max(a.left, b.left);
              const overlapY = Math.min(a.bottom, b.bottom) - Math.max(a.top, b.top);
              if (overlapX > 1 && overlapY > 1) overlaps.push([first, second]);
            }
          }
          const outside = boxes
            .map((box, index) => ({box, index}))
            .filter(({box}) =>
              box.left < gridBox.left - 1 || box.top < gridBox.top - 1 ||
              box.right > gridBox.right + 1 || box.bottom > gridBox.bottom + 1
            )
            .map(({index}) => index);
          return {overlaps, outside};
        }"""
    )
    details.evaluate("element => { element.open = false; }")
    if layout["overlaps"]:
        failures.append(
            f"{label}: expanded controls overlap aggregate tiles {layout['overlaps']}"
        )
    if layout["outside"]:
        failures.append(
            f"{label}: expanded controls push aggregate tiles outside their grid "
            f"{layout['outside']}"
        )


def check_grouped_aggregate_legend(surface, label: str, failures: list[str]) -> None:
    """Confirm grouped legend geometry and fixed swatch width."""
    details = surface.locator("#options")
    details.evaluate("element => { element.open = true; }")
    width_control = surface.locator("[data-sample-width]").first
    width_control.fill("0.3")
    width_control.dispatch_event("input")
    surface.wait_for_timeout(100)
    metrics = surface.evaluate(
        """() => {
          const legend = document.querySelector("#aggregate-legend");
          const grid = document.querySelector("#aggregate-grid");
          const legendBox = legend.getBoundingClientRect();
          const gridBox = grid.getBoundingClientRect();
          return {
            groups: legend.querySelectorAll(".legend-group").length,
            widths: [...legend.querySelectorAll(".legend-line")]
              .map(line => getComputedStyle(line).borderTopWidth),
            thinPlotLines: [...document.querySelectorAll(".aggregate-panel path")]
              .some(path => Number(path.getAttribute("stroke-width")) === 0.3),
            separated: legendBox.bottom <= gridBox.top + 1
          };
        }"""
    )
    if metrics["groups"] != 2:
        failures.append(f"{label}: expected two aggregate legend groups")
    if not metrics["widths"] or set(metrics["widths"]) != {"3px"}:
        failures.append(
            f"{label}: aggregate legend widths are not fixed at 3px "
            f"({metrics['widths']})"
        )
    if not metrics["thinPlotLines"]:
        failures.append(f"{label}: aggregate curve width control did not reach 0.3")
    if not metrics["separated"]:
        failures.append(f"{label}: grouped legend overlaps the aggregate grid")
    width_control.fill("2")
    width_control.dispatch_event("input")


def check_wide_command_layout(page, label: str, failures: list[str]) -> None:
    """Confirm command guides use their available desktop content width."""
    metrics = page.evaluate(
        """() => {
          const content = document.querySelector(".md-typeset");
          const selectors = {
            paragraph: ".md-typeset > p",
            list: ".md-typeset > ul, .md-typeset > ol",
            command: ".md-typeset .highlight",
            table: ".md-typeset .md-typeset__table"
          };
          const contentWidth = content?.getBoundingClientRect().width || 0;
          const widths = {};
          for (const [name, selector] of Object.entries(selectors)) {
            const element = document.querySelector(selector);
            widths[name] = element?.getBoundingClientRect().width || null;
          }
          return {contentWidth, widths};
        }"""
    )
    content_width = metrics["contentWidth"]
    if not content_width:
        failures.append(f"{label}: command-guide content width is unavailable")
        return
    for name, width in metrics["widths"].items():
        if width is not None and width < content_width * 0.9:
            failures.append(
                f"{label}: {name} uses only {width:.0f}px of "
                f"{content_width:.0f}px available content width"
            )


def audit_embedded_review(browser, failures: list[str]) -> None:
    """Exercise the portable repeated-comparison report directly from disk."""
    from fp_tools.tools.static_comparison_browser import write_embedded_static_browser

    def payload(effect: float, aggregate: bool = False) -> dict:
        result = {
            "title": "Differential footprint report",
            "report_label": "Biological-replicate empirical-Bayes comparison",
            "conditions": ["shHIF2", "P"],
            "colors": {
                "shHIF2_up": "#dc2626",
                "P_up": "#2563eb",
                "n.s.": "#8a94a6",
            },
            "change_label": "Differential footprint score",
            "points": [
                {
                    "prefix": "TF1",
                    "name": "TF1",
                    "motif_id": "M1",
                    "group": "shHIF2_up",
                    "change": effect,
                    "pvalue": 1e-6,
                    "fdr": 1e-4,
                    "neglog10p": 6.0,
                },
                {
                    "prefix": "TF2",
                    "name": "TF2",
                    "motif_id": "M2",
                    "group": "P_up",
                    "change": -0.3,
                    "pvalue": 1e-5,
                    "fdr": 1e-3,
                    "neglog10p": 5.0,
                },
            ],
            "motif_matrices": {
                "TF1": [[10, 0], [0, 10], [0, 0], [0, 0]],
                "TF2": [[0, 0], [10, 0], [0, 10], [0, 0]],
            },
            "logos": {},
            "aggregate": {"motifs": []},
        }
        if aggregate:
            result["aggregate"] = {
                "x": [-1, 1],
                "motifs": [
                    {
                        "prefix": "TF1",
                        "name": "TF1",
                        "motif_id": "M1",
                        "n_sites": 10,
                        "conditions": [
                            {
                                "name": "shHIF2",
                                "samples": [
                                    {"name": "shHIF2_rep1", "profile": [0.1, 0.2]}
                                ],
                            },
                            {
                                "name": "P",
                                "samples": [
                                    {"name": "P_rep1", "profile": [0.2, 0.1]}
                                ],
                            },
                        ],
                    }
                ],
            }
        return result

    labels = ["shHIF2 replicate 1", "shHIF2 replicate 2", "shHIF2 replicate 3"]
    review = {
        "schema": "fp-tools.review-multi-comparisons.v1",
        "title": "Kidney comparison review",
        "comparisons": [
            {"label": label, "payload": payload(0.4 + index)}
            for index, label in enumerate(labels)
        ],
    }
    aggregate_review = {
        "schema": "fp-tools.review-multi-comparisons.v1",
        "title": "Aggregate comparison review",
        "comparisons": [{"label": "Aggregate", "payload": payload(0.4, True)}],
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        standalone = write_embedded_static_browser(review, root / "standalone.html")
        aggregate = write_embedded_static_browser(
            aggregate_review, root / "aggregate.html"
        )
        for width, height in VIEWPORTS:
            page = browser.new_page(
                viewport={"width": width, "height": height}, accept_downloads=True
            )
            console_errors: list[str] = []
            requests: list[str] = []
            page.on(
                "console",
                lambda message, errors=console_errors: (
                    errors.append(message.text) if message.type == "error" else None
                ),
            )
            page.on("request", lambda request, seen=requests: seen.append(request.url))
            page.goto(standalone.as_uri(), wait_until="load", timeout=60_000)
            page.wait_for_function(
                "document.querySelector('#status').textContent.includes('motifs')",
                timeout=60_000,
            )
            label = f"standalone repeated comparison at {width}x{height}"
            options = page.locator("#comparison-selector option").all_text_contents()
            if options != labels:
                failures.append(f"{label}: comparison order changed {options}")
            if not page.locator("#comparison-selector-control").is_visible():
                failures.append(f"{label}: Comparison selector is hidden")
            if page.locator("#condition-selector-controls").is_visible():
                failures.append(f"{label}: condition selectors are visible")
            if page.locator(".aggregate-card").is_visible():
                failures.append(f"{label}: aggregate card reserves layout space")
            if page.locator(".options-samples").is_visible():
                failures.append(f"{label}: sample controls are visible")
            if page.locator(".selected-motif").count() != 1:
                failures.append(f"{label}: expected one selected-motif card")
            page.locator("#comparison-selector").select_option("2")
            page.wait_for_function(
                "document.querySelector('#report-title').textContent.includes('replicate 3')"
            )
            first_point = page.locator(".pt").first
            prefix = first_point.get_attribute("data-prefix")
            first_point.dispatch_event("click")
            if page.locator(".panel-tf").input_value() != prefix:
                failures.append(f"{label}: volcano selection did not update motif")
            first_bar = page.locator(".rank-bar").first
            prefix = first_bar.get_attribute("data-prefix")
            first_bar.dispatch_event("click")
            if page.locator(".panel-tf").input_value() != prefix:
                failures.append(f"{label}: ranked selection did not update motif")
            metrics = page.evaluate(
                """() => ({
                  scrollWidth: document.documentElement.scrollWidth,
                  clientWidth: document.documentElement.clientWidth
                })"""
            )
            if metrics["scrollWidth"] > metrics["clientWidth"]:
                failures.append(
                    f"{label}: horizontal overflow "
                    f"{metrics['scrollWidth']}>{metrics['clientWidth']}"
                )
            if console_errors:
                failures.append(f"{label}: console errors {console_errors}")
            external = [url for url in requests if url != standalone.as_uri()]
            if external:
                failures.append(f"{label}: unexpected network requests {external}")
            if width == 1440 and height == 1000:
                for button in (
                    "download-logo",
                    "download-rank",
                    "download-volcano",
                    "download-panel",
                    "download-tsv",
                    "download-all",
                ):
                    with page.expect_download(timeout=30_000):
                        page.locator(f"#{button}").dispatch_event("click")
            page.close()

        page = browser.new_page(viewport={"width": 1440, "height": 1000})
        page.goto(aggregate.as_uri(), wait_until="load", timeout=60_000)
        page.locator(".aggregate-panel").first.wait_for(
            state="visible", timeout=60_000
        )
        if not page.locator(".aggregate-card").is_visible():
            failures.append("aggregate-capable standalone: aggregate card is hidden")
        if not page.locator(".options-samples").is_visible():
            failures.append("aggregate-capable standalone: sample controls are hidden")
        page.close()


def audit(site_dir: Path) -> None:
    failures: list[str] = []
    with serve(site_dir) as base_url, sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        for relative in PAGES:
            for width, height in VIEWPORTS:
                for scheme in ("light",):
                    page = browser.new_page(viewport={"width": width, "height": height})
                    page.emulate_media(color_scheme="dark")
                    page.route(
                        re.compile(r"^https?://(?!127\.0\.0\.1)"),
                        lambda route: route.fulfill(status=204, body=""),
                    )
                    console_errors: list[str] = []
                    failed_requests: list[tuple[str, str]] = []
                    successful_requests: set[str] = set()
                    page.on(
                        "console",
                        lambda message, errors=console_errors: (
                            errors.append(message.text) if message.type == "error" else None
                        ),
                    )
                    page.on(
                        "requestfailed",
                        lambda request, errors=failed_requests: errors.append(
                            (request.url, request.failure or "unknown failure")
                        ),
                    )
                    page.on(
                        "response",
                        lambda response, successes=successful_requests: (
                            successes.add(response.url) if response.ok else None
                        ),
                    )
                    response = page.goto(
                        base_url + relative, wait_until="networkidle", timeout=60_000
                    )
                    if relative.endswith("diff_footprints_K562_HepG2.html"):
                        page.locator(".selected-motif").first.wait_for(
                            state="visible", timeout=60_000
                        )
                    if relative == "demos/reports/region_set_HepG2_HNF4A_FOXA2/":
                        page.locator(".aggregate-panel").first.wait_for(
                            state="visible", timeout=60_000
                        )
                    if relative == "ENCODE-Cancer-Cell-lines-Footprinting/":
                        page.locator(".aggregate-panel").first.wait_for(
                            state="visible", timeout=60_000
                        )
                    embedded_frame = None
                    if relative in REPORT_IFRAME_PAGES:
                        selector = "iframe.fp-report-demo"
                        iframe = page.locator(selector)
                        iframe.wait_for(state="attached", timeout=60_000)
                        iframe.evaluate("element => element.scrollIntoView()")
                        handle = iframe.element_handle()
                        embedded_frame = handle.content_frame() if handle else None
                        if embedded_frame is not None:
                            try:
                                embedded_frame.locator(".aggregate-panel").first.wait_for(
                                    state="visible", timeout=60_000
                                )
                            except PlaywrightTimeoutError:
                                failures.append(
                                    f"{relative}: embedded aggregate report timed out"
                                )
                                embedded_frame = None
                    if relative == "gui/":
                        selector = "iframe.fp-gui-demo"
                        iframe = page.locator(selector)
                        iframe.wait_for(state="attached", timeout=60_000)
                        iframe.evaluate("element => element.scrollIntoView()")
                        handle = iframe.element_handle()
                        embedded_frame = handle.content_frame() if handle else None
                        if embedded_frame is not None:
                            try:
                                embedded_frame.locator('[data-page="home"]').wait_for(
                                    state="attached", timeout=60_000
                                )
                            except PlaywrightTimeoutError:
                                failures.append(f"{relative}: embedded GUI timed out")
                                embedded_frame = None
                    label = f"{relative or '/'} at {width}x{height} ({scheme})"
                    if response is None or response.status != 200:
                        failures.append(
                            f"{label}: HTTP {getattr(response, 'status', None)}"
                        )
                    metrics = page.evaluate(
                        """() => ({
                          h1: document.querySelectorAll("h1").length,
                          scrollWidth: document.documentElement.scrollWidth,
                          clientWidth: document.documentElement.clientWidth,
                          scrollHeight: document.documentElement.scrollHeight,
                          clientHeight: document.documentElement.clientHeight,
                          brokenImages: [...document.images].filter(
                            image => !image.complete || image.naturalWidth === 0
                          ).map(image => image.src),
                          duplicateIds: [...document.querySelectorAll("[id]")]
                            .map(node => node.id)
                            .filter((id, index, values) => values.indexOf(id) !== index),
                          colorScheme: getComputedStyle(document.documentElement).colorScheme,
                          mkdocsScheme: document.body.getAttribute("data-md-color-scheme"),
                          paletteControls: document.querySelectorAll('label[for^="__palette"]').length
                        })"""
                    )
                    if "light" not in metrics["colorScheme"].split():
                        failures.append(
                            f"{label}: page does not force a light color scheme "
                            f"under a dark system preference ({metrics['colorScheme']})"
                        )
                    if relative not in STANDALONE_DEMO_PAGES:
                        if metrics["mkdocsScheme"] != "default":
                            failures.append(
                                f"{label}: unexpected MkDocs palette {metrics['mkdocsScheme']}"
                            )
                        if metrics["paletteControls"]:
                            failures.append(f"{label}: theme toggle should be absent")
                    if metrics["h1"] != 1:
                        failures.append(f"{label}: expected one h1, found {metrics['h1']}")
                    if metrics["scrollWidth"] > metrics["clientWidth"]:
                        failures.append(
                            f"{label}: horizontal overflow "
                            f"{metrics['scrollWidth']}>{metrics['clientWidth']}"
                        )
                    if (
                        relative == "ENCODE-Cancer-Cell-lines-Footprinting/"
                        and width >= 1050
                        and metrics["scrollHeight"] > metrics["clientHeight"]
                    ):
                        failures.append(
                            f"{label}: desktop document scroll "
                            f"{metrics['scrollHeight']}>{metrics['clientHeight']}"
                        )
                    if metrics["brokenImages"]:
                        failures.append(
                            f"{label}: broken images {metrics['brokenImages']}"
                        )
                    if metrics["duplicateIds"]:
                        failures.append(
                            f"{label}: duplicate IDs {metrics['duplicateIds']}"
                        )
                    if relative in DOCUMENTATION_RETURN_PAGES:
                        documentation_link = page.locator("a.documentation-link")
                        if documentation_link.count() != 1:
                            failures.append(
                                f"{label}: documentation return link is missing"
                            )
                        else:
                            href = documentation_link.evaluate("element => element.href") or ""
                            target = documentation_link.get_attribute("target") or ""
                            if href != base_url:
                                failures.append(
                                    f"{label}: documentation return link has invalid target {href}"
                                )
                            if target != "_top":
                                failures.append(
                                    f"{label}: documentation return link does not escape an iframe"
                                )
                    if relative in GET_STARTED_PAGES and width >= 1280:
                        primary = page.locator(".md-sidebar--primary")
                        if not primary.is_visible():
                            failures.append(
                                f"{label}: Get Started navigation is not visible"
                            )
                        secondary = page.locator(".md-sidebar--secondary")
                        if secondary.is_visible():
                            failures.append(f"{label}: page TOC should be hidden")
                        previous_count = page.locator(".md-footer__link--prev").count()
                        next_count = page.locator(".md-footer__link--next").count()
                        if relative in CORE_COMMAND_SEQUENCE:
                            index = CORE_COMMAND_SEQUENCE.index(relative)
                            expected_previous = 1
                            expected_next = int(index < len(CORE_COMMAND_SEQUENCE) - 1)
                            if previous_count != expected_previous:
                                failures.append(
                                    f"{label}: expected {expected_previous} core previous link, found {previous_count}"
                                )
                            if next_count != expected_next:
                                failures.append(
                                    f"{label}: expected {expected_next} core next link, found {next_count}"
                                )
                        elif previous_count or next_count:
                            failures.append(
                                f"{label}: Previous/Next banner must be limited to Core analysis"
                            )
                        family = page.locator("body").evaluate(
                            "element => getComputedStyle(element).fontFamily"
                        )
                        if "Helvetica" not in family and "Arial" not in family:
                            failures.append(
                                f"{label}: unexpected body font family {family}"
                            )
                    if (
                        relative.startswith("get-started/commands/")
                        and width >= 1280
                        and scheme == "light"
                    ):
                        check_wide_command_layout(page, label, failures)
                    if width >= 1280 and relative not in STANDALONE_DEMO_PAGES:
                        header = page.locator(".md-header__inner")
                        tabs = page.locator(".fp-header-tabs")
                        if not tabs.is_visible():
                            failures.append(f"{label}: global navigation is not in the header")
                        else:
                            header_box = header.bounding_box()
                            tabs_box = tabs.bounding_box()
                            if (
                                header_box is None
                                or tabs_box is None
                                or tabs_box["y"] < header_box["y"] - 1
                                or tabs_box["y"] + tabs_box["height"]
                                > header_box["y"] + header_box["height"] + 1
                            ):
                                failures.append(
                                    f"{label}: global navigation is outside the header row"
                                )
                        if page.locator(".md-container > .md-tabs:visible").count():
                            failures.append(f"{label}: separate global navigation row is visible")
                        title_link = page.locator("a.md-header__title")
                        if title_link.count() != 1:
                            failures.append(f"{label}: fp-tools header title is not a home link")
                        else:
                            title_href = title_link.get_attribute("href") or ""
                            if not title_href.endswith("/"):
                                failures.append(
                                    f"{label}: fp-tools header title has invalid home link "
                                    f"{title_href}"
                                )
                        source = page.locator(".md-header__source")
                        if source.is_visible():
                            header_box = header.bounding_box()
                            source_box = source.bounding_box()
                            if (
                                header_box is None
                                or source_box is None
                                or abs(
                                    source_box["x"]
                                    + source_box["width"]
                                    - header_box["x"]
                                    - header_box["width"]
                                )
                                > 8
                            ):
                                failures.append(
                                    f"{label}: repository controls are not right-aligned"
                                )
                        expected_tab = (
                            "API Reference"
                            if relative == "api/"
                            else "GUI Demo"
                            if relative == "gui/"
                            else "Output Demo"
                            if relative == "reports/"
                            else "Get Started"
                        )
                        active_tabs = page.locator(
                            ".fp-header-tabs .md-tabs__item--active .md-tabs__link"
                        )
                        if active_tabs.count() != 1:
                            failures.append(
                                f"{label}: expected one active global navigation item"
                            )
                        elif active_tabs.first.inner_text().strip() != expected_tab:
                            failures.append(
                                f"{label}: active global navigation item is not "
                                f"{expected_tab}"
                            )
                    if relative in {"api/", "gui/", "reports/"} and width >= 1280:
                        primary = page.locator(".md-sidebar--primary")
                        if primary.is_visible():
                            failures.append(f"{label}: redundant primary navigation is visible")
                    if relative == "api/" and width >= 1280:
                        secondary = page.locator(".md-sidebar--secondary")
                        content = page.locator(".md-content")
                        if not secondary.is_visible():
                            failures.append(f"{label}: API table of contents is not visible")
                        elif secondary.bounding_box()["x"] >= content.bounding_box()["x"]:
                            failures.append(f"{label}: API table of contents is not on the left")
                    if relative in {"gui/", "reports/"} and width >= 1280:
                        secondary = page.locator(".md-sidebar--secondary")
                        if secondary.is_visible():
                            failures.append(f"{label}: redundant page TOC is visible")
                    if console_errors:
                        failures.append(f"{label}: console errors {console_errors}")
                    meaningful_failures = [
                        f"{url} ({reason})"
                        for url, reason in failed_requests
                        if reason != "net::ERR_ABORTED" or url not in successful_requests
                    ]
                    if meaningful_failures:
                        failures.append(
                            f"{label}: failed requests {meaningful_failures}"
                        )

                    if relative in REPORT_IFRAME_PAGES | {"gui/"}:
                        if embedded_frame is None:
                            failures.append(f"{label}: interactive iframe did not load")
                        else:
                            frame_metrics = embedded_frame.evaluate(
                                """() => ({
                                  scrollWidth: document.documentElement.scrollWidth,
                                  clientWidth: document.documentElement.clientWidth,
                                  brokenImages: [...document.images].filter(
                                    image => !image.complete || image.naturalWidth === 0
                                  ).map(image => image.src)
                                })"""
                            )
                            if frame_metrics["scrollWidth"] > frame_metrics["clientWidth"]:
                                failures.append(
                                    f"{label}: iframe horizontal overflow "
                                    f"{frame_metrics['scrollWidth']}>"
                                    f"{frame_metrics['clientWidth']}"
                                )
                            if frame_metrics["brokenImages"]:
                                failures.append(
                                    f"{label}: iframe broken images "
                                    f"{frame_metrics['brokenImages']}"
                                )

                    if relative == "reports/" and embedded_frame is not None:
                        check_expanded_aggregate_layout(embedded_frame, label, failures)
                        first = embedded_frame.locator("#condition-1")
                        second = embedded_frame.locator("#condition-2")
                        if first.input_value() != "K562" or second.input_value() != "HepG2":
                            failures.append(f"{label}: unexpected embedded comparison")
                        if embedded_frame.locator(".selected-motif").count() != 4:
                            failures.append(
                                f"{label}: expected four embedded motif cards"
                            )
                        if embedded_frame.locator(".aggregate-panel").count() != 4:
                            failures.append(
                                f"{label}: expected four embedded aggregate panels"
                            )
                        if width == 1440 and height == 1000 and scheme == "light":
                            first.select_option("A549")
                            second.select_option("HCT116")
                            embedded_frame.locator("#title-cond1").filter(
                                has_text="A549"
                            ).wait_for(timeout=60_000)
                            embedded_frame.locator("#title-cond2").filter(
                                has_text="HCT116"
                            ).wait_for(timeout=60_000)

                    if relative == "gui/" and embedded_frame is not None:
                        if width == 1440 and height == 1000 and scheme == "light":
                            route = embedded_frame.locator(
                                '[data-page="diff-footprints"]'
                            )
                            route.evaluate("element => element.click()")
                            embedded_frame.wait_for_timeout(50)
                            if embedded_frame.url.rsplit("#", 1)[-1] != "diff-footprints":
                                failures.append(
                                    f"{label}: embedded GUI route did not open"
                                )
                            if route.get_attribute("aria-current") != "page":
                                failures.append(
                                    f"{label}: embedded GUI route lacks aria-current"
                                )

                    if relative in REPORT_IFRAME_PAGES | {"gui/"}:
                        frame_state = page.locator("iframe.fp-live-demo").evaluate(
                            """frame => ({
                              ready: frame.contentDocument?.readyState,
                              title: frame.contentDocument?.title || "",
                              contentLength: frame.contentDocument?.body?.textContent?.length || 0
                            })"""
                        )
                        if (
                            frame_state["ready"] != "complete"
                            or not frame_state["title"]
                            or frame_state["contentLength"] < 100
                        ):
                            failures.append(
                                f"{label}: embedded demo did not load {frame_state}"
                            )

                    if relative.endswith("fp-tools-gui-static-demo.html"):
                        for route in ("diff-footprints", "run-history", "home"):
                            page.locator(f'[data-page="{route}"]').evaluate(
                                "(element) => element.click()"
                            )
                            page.wait_for_timeout(50)
                            if page.url.rsplit("#", 1)[-1] != route:
                                failures.append(
                                    f"{label}: GUI route {route} did not open"
                                )
                            current = page.locator(
                                f'[data-page="{route}"]'
                            ).get_attribute("aria-current")
                            if current != "page":
                                failures.append(
                                    f"{label}: GUI route {route} lacks aria-current"
                                )
                    if relative == "ENCODE-Cancer-Cell-lines-Footprinting/":
                        check_expanded_aggregate_layout(page, label, failures)
                        check_grouped_aggregate_legend(page, label, failures)
                        first = page.locator("#condition-1")
                        second = page.locator("#condition-2")
                        if first.input_value() != "K562" or second.input_value() != "HepG2":
                            failures.append(f"{label}: unexpected default comparison")
                        first.select_option("A549")
                        second.select_option("HCT116")
                        page.wait_for_timeout(3_000)
                        if "A549 vs HCT116" not in page.title():
                            failures.append(f"{label}: dropdown comparison did not load")
                        first.select_option("HCT116")
                        page.wait_for_timeout(3_000)
                        if second.input_value() != "A549":
                            failures.append(f"{label}: duplicate selection did not reverse direction")
                        if page.locator(".selected-motif").count() != 4:
                            failures.append(f"{label}: expected four selected motif cards")
                        if page.locator(".aggregate-panel").count() != 4:
                            failures.append(f"{label}: expected four aggregate panels")
                        if width == 1440 and height == 1000:
                            with page.expect_download(timeout=30_000) as download_info:
                                page.locator("#download-volcano").evaluate(
                                    "element => element.click()"
                                )
                            if not download_info.value.suggested_filename.endswith("_volcano.svg"):
                                failures.append(f"{label}: unexpected SVG export filename")
                            with page.expect_download(timeout=30_000) as aggregate_download:
                                page.locator("#download-aggregate").evaluate(
                                    "element => element.click()"
                                )
                            aggregate_path = aggregate_download.value.path()
                            aggregate_svg = Path(aggregate_path).read_text(encoding="utf-8")
                            legend_match = re.search(
                                r'<g class="aggregate-export-legend">(.*?)</g>',
                                aggregate_svg,
                                flags=re.DOTALL,
                            )
                            if legend_match is None:
                                failures.append(f"{label}: exported aggregate legend is missing")
                            elif 'stroke-width="3"' not in legend_match.group(1):
                                failures.append(
                                    f"{label}: exported aggregate legend is not 3 units thick"
                                )
                    if relative == "demos/reports/region_set_HepG2_HNF4A_FOXA2/":
                        check_expanded_aggregate_layout(page, label, failures)
                        check_grouped_aggregate_legend(page, label, failures)
                        first = page.locator("#condition-1")
                        second = page.locator("#condition-2")
                        if (
                            first.input_value() != "HNF4A + FOXA2"
                            or second.input_value() != "No HNF4A/FOXA2"
                        ):
                            failures.append(f"{label}: unexpected region-set default comparison")
                        if page.locator(".selected-motif").count() != 8:
                            failures.append(f"{label}: expected eight selected motif cards")
                        if page.locator(".aggregate-panel").count() != 8:
                            failures.append(f"{label}: expected eight aggregate panels")
                        first.select_option("HNF4A only")
                        page.wait_for_timeout(100)
                        second.select_option("FOXA2 only")
                        page.wait_for_timeout(1_000)
                        if "HNF4A only vs FOXA2 only" not in page.title():
                            failures.append(f"{label}: region-set dropdown comparison did not load")
                    if relative == "reports/" and embedded_frame is not None:
                        check_grouped_aggregate_legend(
                            embedded_frame, f"{label} embedded report", failures
                        )
                    if (
                        relative == ""
                        and width == 1440
                        and height == 1000
                        and scheme == "light"
                    ):
                        output_tab = page.locator(
                            ".fp-header-tabs .md-tabs__link",
                            has_text="Output Demo",
                        )
                        output_tab.evaluate("element => element.click()")
                        page.wait_for_url(re.compile(r"/reports/$"), timeout=60_000)
                        page.wait_for_function(
                            """() => {
                              const active = document.querySelector(
                                ".fp-header-tabs .md-tabs__item--active .md-tabs__link"
                              );
                              return active?.textContent.trim() === "Output Demo" &&
                                active.getAttribute("aria-current") === "page";
                            }""",
                            timeout=60_000,
                        )
                    page.close()
        audit_embedded_review(browser, failures)
        browser.close()
    if failures:
        raise SystemExit("Documentation browser audit failed:\n- " + "\n- ".join(failures))
    print(
        f"OK: audited {len(PAGES)} pages at {len(VIEWPORTS)} viewport sizes "
        "while emulating a dark system preference"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-dir", type=Path, default=Path("site"))
    args = parser.parse_args()
    audit(args.site_dir.resolve())


if __name__ == "__main__":
    main()
