#!/usr/bin/env python3
"""Audit the built documentation with a real browser."""

from __future__ import annotations

import argparse
import contextlib
import http.server
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
)
VIEWPORTS = ((1440, 1000), (390, 844))


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
                    console_errors: list[str] = []
                    failed_requests: list[str] = []
                    page.on(
                        "console",
                        lambda message, errors=console_errors: (
                            errors.append(message.text) if message.type == "error" else None
                        ),
                    )
                    page.on(
                        "requestfailed",
                        lambda request, errors=failed_requests: errors.append(request.url),
                    )
                    response = page.goto(
                        base_url + relative, wait_until="networkidle", timeout=60_000
                    )
                    if scheme == "dark":
                        page.locator('label[for^="__palette"]:visible').first.click()
                        page.wait_for_timeout(100)
                    if relative.endswith("diff_footprints_K562_HepG2.html"):
                        page.locator(".selected-motif").first.wait_for(
                            state="visible", timeout=60_000
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
                    if failed_requests:
                        failures.append(f"{label}: failed requests {failed_requests}")

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
