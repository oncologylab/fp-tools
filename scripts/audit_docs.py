#!/usr/bin/env python3
"""Audit the built documentation with a real browser."""

from __future__ import annotations

import argparse
import contextlib
import http.server
import re
import socketserver
import threading
from pathlib import Path

from playwright.sync_api import sync_playwright


PAGES = (
    "",
    "api/",
    "gui/",
    "reports/",
    "demos/reports/diff_footprints_K562_HepG2.html",
    "demos/gui/fp-tools-gui-static-demo.html",
    "ENCODE-Cancer-Cell-lines-Footprinting/",
)
VIEWPORTS = ((1440, 1000), (1280, 720), (390, 844))


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


def audit(site_dir: Path) -> None:
    failures: list[str] = []
    with serve(site_dir) as base_url, sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        for relative in PAGES:
            schemes = ("light", "dark") if relative in {"", "api/", "gui/", "reports/"} else ("light",)
            for width, height in VIEWPORTS:
                for scheme in schemes:
                    page = browser.new_page(viewport={"width": width, "height": height})
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
                    if scheme == "dark":
                        page.locator('label[for^="__palette"]:visible').first.evaluate(
                            "element => element.click()"
                        )
                        page.wait_for_timeout(100)
                    if relative.endswith("diff_footprints_K562_HepG2.html"):
                        page.locator(".selected-motif").first.wait_for(
                            state="visible", timeout=60_000
                        )
                    if relative == "ENCODE-Cancer-Cell-lines-Footprinting/":
                        page.locator(".aggregate-panel").first.wait_for(
                            state="visible", timeout=60_000
                        )
                    embedded_frame = None
                    if relative == "reports/":
                        selector = "iframe.fp-report-demo"
                        iframe = page.locator(selector)
                        iframe.wait_for(state="attached", timeout=60_000)
                        iframe.evaluate("element => element.scrollIntoView()")
                        handle = iframe.element_handle()
                        embedded_frame = handle.content_frame() if handle else None
                        if embedded_frame is not None:
                            embedded_frame.locator(".aggregate-panel").first.wait_for(
                                state="visible", timeout=60_000
                            )
                    if relative == "gui/":
                        selector = "iframe.fp-gui-demo"
                        iframe = page.locator(selector)
                        iframe.wait_for(state="attached", timeout=60_000)
                        iframe.evaluate("element => element.scrollIntoView()")
                        handle = iframe.element_handle()
                        embedded_frame = handle.content_frame() if handle else None
                        if embedded_frame is not None:
                            embedded_frame.locator('[data-page="home"]').wait_for(
                                state="attached", timeout=60_000
                            )
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
                            .filter((id, index, values) => values.indexOf(id) !== index)
                        })"""
                    )
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

                    if relative in {"reports/", "gui/"}:
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

                    if relative in {"gui/", "reports/"}:
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
                    page.close()
        browser.close()
    if failures:
        raise SystemExit("Documentation browser audit failed:\n- " + "\n- ".join(failures))
    print(
        f"OK: audited {len(PAGES)} pages at {len(VIEWPORTS)} viewport sizes "
        "in supported color schemes"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-dir", type=Path, default=Path("site"))
    args = parser.parse_args()
    audit(args.site_dir.resolve())


if __name__ == "__main__":
    main()
