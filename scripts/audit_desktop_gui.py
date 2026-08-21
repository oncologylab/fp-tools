#!/usr/bin/env python3
"""Browser-level responsive audit for a frozen fp-tools desktop GUI."""

from __future__ import annotations

import argparse
import os
import socket
import signal
import subprocess
import tempfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

from playwright.sync_api import sync_playwright


SIDEBAR_FIRST_PAGES = (
    ("Overview", "Home"),
    ("Core workflow", "atac-correct"),
    ("Workflow and interface", "bulk-footprinting"),
    ("Signals and reports", "normalize-bigwig"),
    ("De Novo Motif Discovery", "discover-motifs"),
    ("Single-cell ATAC-seq", "pseudobulk-fragments"),
    ("Configuration", "Config"),
)
MINIMUM_SIDEBAR_GAP_PX = 2.0


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as handle:
        handle.bind(("127.0.0.1", 0))
        return int(handle.getsockname()[1])


def _wait_for_health(url: str, process: subprocess.Popen[str], timeout: float) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise RuntimeError("Desktop GUI exited before browser audit")
        try:
            with urllib.request.urlopen(url + "/_stcore/health", timeout=2) as response:
                if response.status == 200:
                    return
        except Exception:
            time.sleep(0.5)
    raise RuntimeError("Desktop GUI health check timed out")


def _stop(process: subprocess.Popen[str]) -> None:
    if os.name == "nt":
        subprocess.run(
            ["taskkill", "/PID", str(process.pid), "/T", "/F"],
            capture_output=True,
            check=False,
        )
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
    try:
        process.wait(timeout=20)
    except subprocess.TimeoutExpired:
        if os.name == "nt":
            process.kill()
        else:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        process.wait(timeout=20)


def _write_valid_bulk_fixture(root: Path) -> tuple[Path, Path, Path]:
    """Create existence-only inputs for the GUI validation state transition."""

    fixture = root / "valid-bulk-inputs"
    fixture.mkdir()
    genome = fixture / "genome.fa"
    peaks = fixture / "peaks.bed"
    for path in (genome, peaks):
        path.write_text("fixture\n", encoding="utf-8")

    sample_rows = ["sample\tcondition\tbam\tpeaks"]
    for sample, condition in (("sample_a", "A"), ("sample_b", "B")):
        bam = fixture / f"{sample}.bam"
        bai = fixture / f"{sample}.bam.bai"
        bam.write_bytes(b"fixture")
        bai.write_bytes(b"fixture")
        sample_rows.append(f"{sample}\t{condition}\t{bam}\t{peaks}")

    samples = fixture / "samples.tsv"
    samples.write_text("\n".join(sample_rows) + "\n", encoding="utf-8")
    comparisons = fixture / "comparisons.tsv"
    comparisons.write_text(
        "comparison\tcond1\tcond2\nA_vs_B\tA\tB\n",
        encoding="utf-8",
    )
    return samples, comparisons, genome


def _audit_page(
    browser,
    base_url: str,
    page_name: str,
    width: int,
    height: int,
    output: Path,
    capture_screenshots: bool,
) -> None:
    context = browser.new_context(viewport={"width": width, "height": height})
    page = context.new_page()
    page.goto(
        f"{base_url}/?page={urllib.parse.quote(page_name)}",
        wait_until="domcontentloaded",
    )
    page.locator(".fp-page-heading h1", has_text=page_name).wait_for(timeout=60_000)
    summary = page.locator(".fp-run-summary").first
    summary.wait_for(timeout=30_000)
    metrics = summary.evaluate(
        """element => ({
          columns: getComputedStyle(element).gridTemplateColumns.split(' ').filter(Boolean).length,
          wordBreaks: [...element.querySelectorAll('strong')].map(node => getComputedStyle(node).wordBreak),
          overflow: document.documentElement.scrollWidth > document.documentElement.clientWidth
        })"""
    )
    expected_columns = 3 if width > 1500 else (2 if width > 900 else 1)
    if metrics["columns"] != expected_columns:
        raise RuntimeError(
            f"{page_name} at {width}px uses {metrics['columns']} run-summary columns; "
            f"expected {expected_columns}"
        )
    if set(metrics["wordBreaks"]) - {"normal"}:
        raise RuntimeError(f"{page_name} at {width}px permits mid-word summary breaks")
    if metrics["overflow"]:
        raise RuntimeError(f"{page_name} at {width}px has horizontal overflow")
    if capture_screenshots:
        page.screenshot(
            path=str(output / f"{page_name}-{width}.png"),
            full_page=False,
            animations="disabled",
            timeout=30_000,
        )
    print(f"Audited {page_name} at {width}x{height}", flush=True)
    context.close()


def _audit_sidebar_navigation(
    browser,
    base_url: str,
    width: int,
    height: int,
    device_scale_factor: float,
    output: Path,
    capture_screenshots: bool,
) -> None:
    """Assert that every group label clears its first button, including active rows."""

    context = browser.new_context(
        viewport={"width": width, "height": height},
        device_scale_factor=device_scale_factor,
    )
    page = context.new_page()
    try:
        page.goto(f"{base_url}/?page=Home", wait_until="domcontentloaded")
        page.locator(".fp-nav-group").first.wait_for(timeout=60_000)
        page.wait_for_function(
            "() => document.querySelectorAll('.fp-nav-group').length === 7",
            timeout=60_000,
        )
        for active_group, page_name in SIDEBAR_FIRST_PAGES:
            print(
                f"Checking active sidebar row {page_name} at {width}px/DPR "
                f"{device_scale_factor}",
                flush=True,
            )
            page.locator(".fp-nav-group").evaluate_all(
                """(headings, activeGroup) => headings.forEach(heading => {
                  const container = heading.closest('[data-testid="stElementContainer"]');
                  let sibling = container ? container.nextElementSibling : null;
                  let button = null;
                  while (sibling && !button) {
                    button = sibling.querySelector('button');
                    sibling = sibling.nextElementSibling;
                  }
                  if (button) button.disabled = heading.textContent.trim() === activeGroup;
                })""",
                active_group,
            )
            page.wait_for_timeout(50)
            measurements = page.locator(".fp-nav-group").evaluate_all(
                """headings => headings.map(heading => {
                  const container = heading.closest('[data-testid="stElementContainer"]');
                  let sibling = container ? container.nextElementSibling : null;
                  let button = null;
                  while (sibling && !button) {
                    button = sibling.querySelector('button');
                    sibling = sibling.nextElementSibling;
                  }
                  const headingBox = heading.getBoundingClientRect();
                  const buttonBox = button ? button.getBoundingClientRect() : null;
                  return {
                    group: heading.textContent.trim(),
                    headingBottom: headingBox.bottom,
                    buttonTop: buttonBox ? buttonBox.top : null,
                    gap: buttonBox ? buttonBox.top - headingBox.bottom : null,
                    buttonText: button ? button.textContent.trim() : null,
                    buttonDisabled: button ? button.disabled : null,
                    clipped: heading.scrollWidth > heading.clientWidth + 1 ||
                      heading.scrollHeight > heading.clientHeight + 1
                  };
                })"""
            )
            failures = []
            for measurement in measurements:
                if measurement["buttonTop"] is None:
                    failures.append(f"{measurement['group']}: first navigation row was not found")
                elif measurement["gap"] < MINIMUM_SIDEBAR_GAP_PX:
                    failures.append(
                        f"{measurement['group']}: heading-to-row gap is "
                        f"{measurement['gap']:.2f}px; expected at least "
                        f"{MINIMUM_SIDEBAR_GAP_PX:.2f}px"
                    )
                if measurement["clipped"]:
                    failures.append(f"{measurement['group']}: heading text is clipped")
            active = next(
                (value for value in measurements if value["group"] == active_group),
                None,
            )
            if active is None:
                failures.append(f"{active_group}: active group heading was not found")
            elif not active["buttonDisabled"]:
                failures.append(
                    f"{active_group}: first row {active['buttonText']!r} is not active"
                )
            elif active["buttonText"] != page_name:
                failures.append(
                    f"{active_group}: expected first row {page_name!r}, found "
                    f"{active['buttonText']!r}"
                )
            if failures:
                screenshot = output / f"sidebar-{page_name}-{width}-failure.png"
                screenshot_note = "screenshot skipped"
                if capture_screenshots:
                    screenshot_note = f"screenshot: {screenshot}"
                    try:
                        page.screenshot(
                            path=str(screenshot),
                            full_page=False,
                            animations="disabled",
                            timeout=10_000,
                        )
                    except Exception as exc:
                        screenshot_note = f"screenshot capture failed: {exc}"
                raise RuntimeError(
                    f"Sidebar audit failed for {page_name} at {width}px/DPR "
                    f"{device_scale_factor}: {'; '.join(failures)}; {screenshot_note}"
                )
        print(
            f"Audited every sidebar group at {width}x{height}, DPR {device_scale_factor}",
            flush=True,
        )
    finally:
        context.close()


def _audit_validation_layout(
    browser,
    base_url: str,
    page_name: str,
    width: int,
    height: int,
    device_scale_factor: float,
    output: Path,
    capture_screenshots: bool,
) -> None:
    """Keep long validation paths and YAML previews inside the Run column."""

    context = browser.new_context(
        viewport={"width": width, "height": height},
        device_scale_factor=device_scale_factor,
    )
    page = context.new_page()
    try:
        initial_page = "sc-footprinting" if page_name == "Config" else page_name
        page.goto(
            f"{base_url}/?page={urllib.parse.quote(initial_page)}",
            wait_until="domcontentloaded",
        )
        page.locator(".fp-page-heading h1", has_text=initial_page).wait_for(timeout=60_000)
        page.locator(".fp-validation-errors").first.wait_for(timeout=30_000)
        if page_name == "Config":
            page.get_by_role("button", name="Config", exact=True).click(force=True)
            page.locator(".fp-page-heading h1", has_text="Config").wait_for(timeout=60_000)
            page.locator(".fp-validation-errors").first.wait_for(timeout=30_000)

        errors = page.locator(".fp-validation-errors").first
        run_column = errors.locator(
            "xpath=ancestor::*[@data-testid='stColumn'][1]"
        )
        main = page.locator("[data-testid='stMain']")
        measurements = errors.evaluate(
            """element => {
              const runColumn = element.closest('[data-testid="stColumn"]');
              const main = element.closest('[data-testid="stMain"]');
              const columnBox = runColumn.getBoundingClientRect();
              return {
                mainOverflow: main.scrollWidth - main.clientWidth,
                listOverflow: element.scrollWidth - element.clientWidth,
                outsideLeft: columnBox.left - element.getBoundingClientRect().left,
                outsideRight: element.getBoundingClientRect().right - columnBox.right,
                itemOverflow: [...element.querySelectorAll('li')].map(
                  item => item.scrollWidth - item.clientWidth
                ),
                itemOutside: [...element.querySelectorAll('li')].map(item => {
                  const box = item.getBoundingClientRect();
                  return Math.max(columnBox.left - box.left, box.right - columnBox.right);
                })
              };
            }"""
        )
        failures = []
        if measurements["mainOverflow"] > 1:
            failures.append(f"main pane overflows by {measurements['mainOverflow']:.1f}px")
        if measurements["listOverflow"] > 1:
            failures.append(f"validation list overflows by {measurements['listOverflow']:.1f}px")
        if max(measurements["outsideLeft"], measurements["outsideRight"]) > 1:
            failures.append("validation list extends outside the Run column")
        if any(value > 1 for value in measurements["itemOverflow"]):
            failures.append("a validation message has horizontal overflow")
        if any(value > 1 for value in measurements["itemOutside"]):
            failures.append("a validation message extends outside the Run column")
        if not page.get_by_role("button", name="Start run", exact=True).is_disabled():
            failures.append("Start run is enabled for an invalid configuration")

        page.get_by_text("Preview runnable YAML", exact=True).evaluate(
            "element => { element.closest('details').open = true; }"
        )
        code = run_column.locator("[data-testid='stCode'], [data-testid='stCodeBlock']").first
        code.wait_for(timeout=30_000)
        code_metrics = code.evaluate(
            """element => {
              const main = element.closest('[data-testid="stMain"]');
              const column = element.closest('[data-testid="stColumn"]');
              const box = element.getBoundingClientRect();
              const columnBox = column.getBoundingClientRect();
              const pre = element.querySelector('pre');
              return {
                mainOverflow: main.scrollWidth - main.clientWidth,
                outside: Math.max(columnBox.left - box.left, box.right - columnBox.right),
                localOverflow: pre ? pre.scrollWidth - pre.clientWidth : 0,
                preOverflowX: pre ? getComputedStyle(pre).overflowX : ''
              };
            }"""
        )
        if code_metrics["mainOverflow"] > 1:
            failures.append("YAML preview creates main-pane overflow")
        if code_metrics["outside"] > 1:
            failures.append("YAML preview extends outside the Run column")
        if code_metrics["localOverflow"] > 1 and code_metrics["preOverflowX"] not in {
            "auto",
            "scroll",
        }:
            failures.append("wide YAML does not scroll within its own code block")

        if failures:
            screenshot = output / f"validation-{page_name}-{width}-failure.png"
            if capture_screenshots:
                page.screenshot(
                    path=str(screenshot),
                    full_page=False,
                    animations="disabled",
                    timeout=30_000,
                )
            raise RuntimeError(
                f"Validation-layout audit failed for {page_name} at {width}px: "
                + "; ".join(failures)
            )
        print(
            f"Audited validation layout for {page_name} at {width}x{height}, "
            f"DPR {device_scale_factor}",
            flush=True,
        )
    finally:
        context.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable")
    parser.add_argument("--output-dir", default="desktop-gui-audit")
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument(
        "--skip-screenshots",
        action="store_true",
        help="Run geometry checks without writing screenshots (useful on headless hosts).",
    )
    args = parser.parse_args()

    executable = Path(args.executable).resolve()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    port = _free_port()
    base_url = f"http://127.0.0.1:{port}"
    with tempfile.TemporaryDirectory(prefix="fp-tools-desktop-browser-") as workdir:
        workdir_path = Path(workdir)
        run_dir = workdir_path / "runs"
        runtime_cache = workdir_path / "runtime-cache"
        process = subprocess.Popen(
            [
                str(executable),
                "--fp-tools-internal-gui-server",
                "--port",
                str(port),
                "--run-dir",
                str(run_dir),
                "--no-browser",
            ],
            cwd=workdir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env={
                **os.environ,
                "FP_TOOLS_RUNTIME_CACHE": str(runtime_cache),
                "LOCALAPPDATA": str(workdir_path / "local-app-data"),
                "STREAMLIT_BROWSER_GATHER_USAGE_STATS": "false",
            },
            start_new_session=os.name != "nt",
        )
        try:
            _wait_for_health(base_url, process, args.timeout)
            with sync_playwright() as playwright:
                browser = playwright.chromium.launch(
                    headless=True,
                    args=["--disable-dev-shm-usage"],
                )
                try:
                    for page_name in ("plot-aggregate", "bulk-footprinting", "review-multi-comparisons"):
                        for width, height in ((1920, 1080), (1280, 720), (390, 844)):
                            _audit_page(
                                browser,
                                base_url,
                                page_name,
                                width,
                                height,
                                output,
                                not args.skip_screenshots,
                            )

                    for width, height, scale in (
                        (1920, 1080, 2.0),
                        (1280, 720, 1.25),
                    ):
                        _audit_sidebar_navigation(
                            browser,
                            base_url,
                            width,
                            height,
                            scale,
                            output,
                            not args.skip_screenshots,
                        )

                    for width, height, scale in (
                        (1920, 1080, 2.0),
                        (1280, 720, 1.25),
                    ):
                        for page_name in (
                            "sc-footprinting",
                            "pseudobulk-fragments",
                            "find-signature-fp",
                            "Config",
                        ):
                            _audit_validation_layout(
                                browser,
                                base_url,
                                page_name,
                                width,
                                height,
                                scale,
                                output,
                                not args.skip_screenshots,
                            )

                    context = browser.new_context(viewport={"width": 1280, "height": 720})
                    page = context.new_page()
                    page.goto(f"{base_url}/?page=bulk-footprinting", wait_until="domcontentloaded")
                    page.locator(".fp-page-heading h1", has_text="bulk-footprinting").wait_for(
                        timeout=60_000
                    )
                    loader = page.locator("details", has_text="Load bulk-footprinting config").first
                    loader.wait_for(timeout=30_000)
                    loader.evaluate("element => { element.open = true; }")
                    page.get_by_text("Example YAML", exact=True).wait_for(timeout=30_000)
                    page.get_by_role("button", name="Load example", exact=True).wait_for(timeout=30_000)
                    page.get_by_text("Config needs fixes before launch.", exact=True).wait_for(timeout=30_000)
                    start_button = page.get_by_role("button", name="Start run", exact=True)
                    if not start_button.is_disabled():
                        raise RuntimeError("Start run remains enabled for an invalid bulk configuration")
                    if run_dir.exists() and any(run_dir.iterdir()):
                        raise RuntimeError("Invalid bulk configuration created a GUI run")
                    if runtime_cache.exists() and any(runtime_cache.iterdir()):
                        raise RuntimeError("Invalid bulk configuration started managed-runtime preparation")
                    if page.get_by_text("Reads Table", exact=True).count():
                        raise RuntimeError("Desktop GUI unexpectedly exposes a FASTQ reads-table field")

                    samples, comparisons, genome = _write_valid_bulk_fixture(workdir_path)
                    page.get_by_label("Sample Table", exact=True).fill(str(samples))
                    page.get_by_label("Comparison Table", exact=True).fill(str(comparisons))
                    page.get_by_label("Genome", exact=True).fill(str(genome))
                    page.get_by_role("button", name="Update page config", exact=True).click(force=True)
                    page.get_by_text("Config is ready to run.", exact=True).wait_for(timeout=30_000)
                    if start_button.is_disabled():
                        raise RuntimeError("Start run remains disabled after the bulk configuration becomes valid")
                    context.close()
                finally:
                    browser.close()
        finally:
            _stop(process)
    if args.skip_screenshots:
        print("Desktop GUI browser audit passed (screenshots skipped on this headless host)")
    else:
        print(f"Desktop GUI browser audit passed; screenshots: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
