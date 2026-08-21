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

import yaml
from playwright.sync_api import expect, sync_playwright


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
GUI_RUN_PAGES = (
    "atac-correct",
    "call-footprints",
    "match-motifs",
    "diff-footprints",
    "normalize-bigwig",
    "plot-aggregate",
    "bulk-footprinting",
    "review-multi-comparisons",
    "discover-motifs",
    "summarize-motifs",
    "pseudobulk-fragments",
    "find-signature-fp",
    "sc-footprinting",
    "Config",
)
SIMPLE_GUI_PAGES = ("Home", "Run History")


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


def _write_loaded_bulk_config(root: Path) -> tuple[Path, dict[str, object]]:
    samples, comparisons, genome = _write_valid_bulk_fixture(root)
    motif = root / "motifs" / "MA0139.1.jaspar"
    motif.parent.mkdir()
    motif.write_text(">MA0139.1 CTCF\nA [1]\nC [0]\nG [0]\nT [0]\n", encoding="utf-8")
    values: dict[str, object] = {
        "sample_table": str(samples),
        "comparison_table": str(comparisons),
        "genome": str(genome),
        "outdir": str(root / "bulk project with spaces"),
        "cores": 2,
        "normalization": "sample-quantile",
        "plot_aggregate": "top",
        "review_format": "bundle",
        "motifs": [str(motif)],
        "dry_run": True,
    }
    config = root / "loaded-bulk.yml"
    config.write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "run_mode": "single",
                "defaults": {},
                "samples": [
                    {
                        "sample_id": "loaded_bulk",
                        "tool": "bulk-footprinting",
                        **values,
                    }
                ],
                "comparisons": [],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return config, values


def _write_uploaded_diff_config(root: Path) -> tuple[Path, dict[str, object]]:
    values: dict[str, object] = {
        "comparison_axis": "regions",
        "motifs": [
            str(root / "motif one.jaspar"),
            str(root / "motif two.jaspar"),
        ],
        "motif_db": "jaspar2026_vertebrates",
        "signals": [str(root / "replicate 1.bw"), str(root / "replicate 2.bw")],
        "sample_names": ["replicate_1", "replicate_2"],
        "cond_names": ["K562", "K562"],
        "regions": [str(root / "CTCF bound.bed"), str(root / "matched control.bed")],
        "region_labels": ["CTCF bound", "matched control"],
        "region_strata_column": 4,
        "genome": str(root / "genome.fa"),
        "outdir": str(root / "region comparison"),
        "cores": 3,
        "skip_excel": True,
    }
    config = root / "uploaded-diff.yml"
    config.write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "run_mode": "single",
                "defaults": {},
                "samples": [
                    {
                        "sample_id": "uploaded_regions",
                        "tool": "diff-footprints",
                        **values,
                    }
                ],
                "comparisons": [],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return config, values


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
    console_errors: list[str] = []
    page.on(
        "console",
        lambda message: console_errors.append(message.text)
        if message.type == "error"
        else None,
    )
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
    if page.locator("[data-testid='stException']").count():
        raise RuntimeError(f"{page_name} at {width}px rendered a Streamlit exception")
    label_locator = page.locator(
        "[data-testid='stMain'] [data-testid='stWidgetLabel']"
    )
    label_locator.first.wait_for(state="visible", timeout=30_000)
    visible_labels = label_locator.evaluate_all(
        """elements => elements.filter(element => {
          const text = (element.innerText || element.textContent || '').trim();
          const style = getComputedStyle(element);
          const box = element.getBoundingClientRect();
          return text && box.width > 0 && box.height > 0 && style.display !== 'none';
        }).map(element => {
          const style = getComputedStyle(element);
          const match = style.color.match(/[\\d.]+/g) || [];
          const rgb = match.slice(0, 3).map(Number);
          return {
            text: (element.innerText || element.textContent || '').trim().slice(0, 80),
            color: style.color,
            opacity: Number(style.opacity),
            visibility: style.visibility,
            nearlyWhite: rgb.length === 3 && rgb.every(value => value > 235)
          };
        })"""
    )
    if not visible_labels:
        raise RuntimeError(f"{page_name} at {width}px has no visible form labels")
    unreadable_labels = [
        label
        for label in visible_labels
        if label["opacity"] < 0.95
        or label["visibility"] != "visible"
        or label["nearlyWhite"]
        or label["color"] == "rgba(0, 0, 0, 0)"
    ]
    if unreadable_labels:
        raise RuntimeError(
            f"{page_name} at {width}px has unreadable form labels: "
            f"{unreadable_labels[:5]}"
        )
    field_metrics = page.locator(
        "[data-testid='stMain'] input, [data-testid='stMain'] textarea, "
        "[data-testid='stMain'] button, [data-testid='stMain'] label"
    ).evaluate_all(
        """elements => elements.filter(element => {
          const style = getComputedStyle(element);
          const box = element.getBoundingClientRect();
          return style.display !== 'none' && style.visibility !== 'hidden' && box.width > 0;
        }).map(element => {
          const box = element.getBoundingClientRect();
          return {
            tag: element.tagName,
            text: (element.innerText || element.textContent || '').trim().slice(0, 80),
            outside: Math.max(0, -box.left, box.right - window.innerWidth),
            clipped: element.scrollWidth - element.clientWidth
          };
        })"""
    )
    bad_fields = [
        metric
        for metric in field_metrics
        if metric["outside"] > 1 or (metric["tag"] != "INPUT" and metric["clipped"] > 2)
    ]
    if bad_fields:
        raise RuntimeError(
            f"{page_name} at {width}px has clipped/out-of-viewport controls: "
            f"{bad_fields[:5]}"
        )
    if console_errors:
        raise RuntimeError(
            f"{page_name} at {width}px logged browser errors: {console_errors[:3]}"
        )
    if capture_screenshots and width in {1280, 390}:
        page.screenshot(
            path=str(output / f"{page_name}-{width}.png"),
            full_page=False,
            animations="disabled",
            timeout=30_000,
        )
    print(f"Audited {page_name} at {width}x{height}", flush=True)
    context.close()


def _audit_simple_page(
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
    try:
        page.goto(
            f"{base_url}/?page={urllib.parse.quote(page_name)}",
            wait_until="domcontentloaded",
        )
        if page_name == "Home":
            page.locator(".fp-hero h1").wait_for(timeout=60_000)
        else:
            page.locator(".fp-page-heading h1", has_text=page_name).wait_for(timeout=60_000)
        metrics = page.evaluate(
            """() => ({
              overflow: document.documentElement.scrollWidth - document.documentElement.clientWidth,
              exceptions: document.querySelectorAll('[data-testid="stException"]').length
            })"""
        )
        if metrics["overflow"] > 1 or metrics["exceptions"]:
            raise RuntimeError(f"{page_name} at {width}px is not rendered cleanly: {metrics}")
        if capture_screenshots and width in {1280, 390}:
            page.screenshot(
                path=str(output / f"{page_name.replace(' ', '-')}-{width}.png"),
                full_page=False,
                animations="disabled",
                timeout=30_000,
            )
    finally:
        context.close()


def _audit_compact_sidebar(
    browser,
    base_url: str,
    output: Path,
    capture_screenshots: bool,
) -> None:
    context = browser.new_context(viewport={"width": 1280, "height": 720})
    page = context.new_page()
    try:
        page.goto(f"{base_url}/?page=Home", wait_until="domcontentloaded")
        page.locator(".fp-sidebar-brand-title", has_text="fp-tools").wait_for(timeout=60_000)
        sidebar = page.locator("[data-testid='stSidebar']")
        text = sidebar.inner_text()
        if "Command-first footprinting workflows" in text:
            raise RuntimeError("Sidebar still shows the persistent branding tagline")
        if "Current run directory" in text:
            raise RuntimeError("Collapsed Workspace exposes the current run directory")
        overview_top = page.locator(".fp-nav-group", has_text="Overview").bounding_box()
        if not overview_top or overview_top["y"] > 175:
            raise RuntimeError(f"Sidebar navigation starts too low: {overview_top}")
        workspace = page.locator("details", has_text="Workspace").first
        workspace.evaluate("element => { element.open = true; }")
        page.get_by_text("Current run directory", exact=True).wait_for(timeout=30_000)
        path = page.locator(".fp-workspace-path")
        path_text = path.inner_text().strip()
        if "very long" not in path_text and "very\\long" not in path_text:
            raise RuntimeError(f"Workspace does not expose the complete test path: {path_text}")
        path_metrics = path.evaluate(
            """element => {
              const box = element.getBoundingClientRect();
              const sidebar = element.closest('[data-testid="stSidebar"]').getBoundingClientRect();
              return {
                outside: Math.max(sidebar.left - box.left, box.right - sidebar.right),
                overflow: element.scrollWidth - element.clientWidth
              };
            }"""
        )
        if path_metrics["outside"] > 1 or path_metrics["overflow"] > 1:
            raise RuntimeError(f"Workspace path overflows the sidebar: {path_metrics}")
        if capture_screenshots:
            page.screenshot(
                path=str(output / "sidebar-workspace-expanded-1280.png"),
                full_page=False,
                animations="disabled",
                timeout=30_000,
            )
    finally:
        context.close()


def _audit_mobile_sidebar_navigation(
    browser,
    base_url: str,
    output: Path,
    capture_screenshots: bool,
) -> None:
    context = browser.new_context(viewport={"width": 390, "height": 844})
    page = context.new_page()
    try:
        page.goto(f"{base_url}/?page=Home", wait_until="domcontentloaded")
        page.locator(".fp-hero h1").wait_for(timeout=60_000)
        open_button = page.get_by_role("button", name="Open navigation")
        open_button.wait_for(state="visible", timeout=30_000)
        box = open_button.bounding_box()
        if not box or min(box["width"], box["height"]) < 40:
            raise RuntimeError(f"Mobile navigation control is too small: {box}")
        open_button.click()
        page.wait_for_function(
            """() => {
              const sidebar = document.querySelector('[data-testid="stSidebar"]');
              return sidebar && sidebar.getBoundingClientRect().left >= -1;
            }""",
            timeout=30_000,
        )
        close_button = page.get_by_role("button", name="Close navigation")
        close_button.wait_for(state="visible", timeout=30_000)
        target = page.get_by_role("button", name="bulk-footprinting", exact=True)
        target.scroll_into_view_if_needed()
        target.click()
        page.locator(".fp-page-heading h1", has_text="bulk-footprinting").wait_for(
            timeout=60_000
        )
        reopened = page.get_by_role("button", name="Open navigation")
        reopened.wait_for(state="visible", timeout=30_000)
        reopened.click()
        page.get_by_role("button", name="Close navigation").click()
        reopened.wait_for(state="visible", timeout=30_000)
        metrics = page.evaluate(
            """() => ({
              overflow: document.documentElement.scrollWidth - document.documentElement.clientWidth,
              exceptions: document.querySelectorAll('[data-testid="stException"]').length
            })"""
        )
        if metrics["overflow"] > 1 or metrics["exceptions"]:
            raise RuntimeError(f"Mobile navigation left a broken page: {metrics}")
        if capture_screenshots:
            reopened.click()
            page.screenshot(
                path=str(output / "mobile-sidebar-open-390.png"),
                full_page=False,
                animations="disabled",
                timeout=30_000,
            )
    finally:
        context.close()


def _assert_control_value(page, label: str, expected: object) -> None:
    control = page.get_by_label(label, exact=True)
    control.wait_for(state="visible", timeout=30_000)
    if isinstance(expected, bool):
        if expected:
            expect(control).to_be_checked(timeout=30_000)
        else:
            expect(control).not_to_be_checked(timeout=30_000)
        actual = control.is_checked()
    else:
        try:
            expect(control).to_have_value(str(expected), timeout=30_000)
            actual = control.input_value()
        except Exception:
            actual = control.inner_text().strip()
    if str(actual) != str(expected):
        raise RuntimeError(
            f"Control {label!r} shows {actual!r}; expected {expected!r}"
        )


def _audit_loaded_config_sync(
    browser,
    base_url: str,
    workdir: Path,
    output: Path,
    capture_screenshots: bool,
) -> None:
    """Loaded files, examples, and uploads must replace visible widget state."""

    context = browser.new_context(viewport={"width": 1440, "height": 960})
    page = context.new_page()
    try:
        bulk_path, values = _write_loaded_bulk_config(workdir)
        page.goto(f"{base_url}/?page=bulk-footprinting", wait_until="domcontentloaded")
        page.locator(".fp-page-heading h1", has_text="bulk-footprinting").wait_for(
            timeout=60_000
        )
        loader = page.locator("details", has_text="Load bulk-footprinting config").first
        loader.evaluate("element => { element.open = true; }")
        page.get_by_label("Config path", exact=True).fill(str(bulk_path))
        page.get_by_role("button", name="Load YAML from path", exact=True).click()
        page.get_by_label("Sample Table", exact=True).wait_for(timeout=30_000)
        for label, key in (
            ("Sample Table", "sample_table"),
            ("Comparison Table", "comparison_table"),
            ("Genome", "genome"),
            ("Outdir", "outdir"),
            ("Cores", "cores"),
            ("Normalization", "normalization"),
            ("Plot Aggregate", "plot_aggregate"),
            ("Review Format", "review_format"),
        ):
            _assert_control_value(page, label, values[key])
        _assert_control_value(page, "Motifs", str(values["motifs"][0]))
        _assert_control_value(page, "Dry Run", True)

        new_outdir = str(workdir / "changed output only")
        page.get_by_label("Outdir", exact=True).fill(new_outdir)
        page.get_by_role("button", name="Update page config", exact=True).click()
        page.get_by_label("Outdir", exact=True).wait_for(timeout=30_000)
        _assert_control_value(page, "Outdir", new_outdir)
        for label, key in (
            ("Sample Table", "sample_table"),
            ("Comparison Table", "comparison_table"),
            ("Genome", "genome"),
            ("Cores", "cores"),
            ("Normalization", "normalization"),
            ("Plot Aggregate", "plot_aggregate"),
            ("Review Format", "review_format"),
        ):
            _assert_control_value(page, label, values[key])

        page.goto(f"{base_url}/?page=normalize-bigwig", wait_until="domcontentloaded")
        page.locator(".fp-page-heading h1", has_text="normalize-bigwig").wait_for(
            timeout=60_000
        )
        normalizer_loader = page.locator(
            "details", has_text="Load normalize-bigwig config"
        ).first
        normalizer_loader.evaluate("element => { element.open = true; }")
        example_select = page.get_by_label("Example YAML", exact=True)
        example_select.click()
        page.get_by_text("normalize_bigwig_single.yml", exact=True).click()
        page.get_by_role("button", name="Load example", exact=True).click()
        page.get_by_label("Background", exact=True).wait_for(timeout=30_000)
        if not page.get_by_label("Background", exact=True).input_value().strip():
            raise RuntimeError("Example YAML did not refresh normalize-bigwig fields")

        diff_path, diff_values = _write_uploaded_diff_config(workdir)
        page.goto(f"{base_url}/?page=diff-footprints", wait_until="domcontentloaded")
        page.locator(".fp-page-heading h1", has_text="diff-footprints").wait_for(
            timeout=60_000
        )
        diff_loader = page.locator("details", has_text="Load diff-footprints config").first
        diff_loader.evaluate("element => { element.open = true; }")
        page.locator("input[type='file']").set_input_files(str(diff_path))
        page.get_by_role("button", name="Apply uploaded YAML", exact=True).click()
        page.get_by_label("Comparison axis", exact=True).wait_for(timeout=30_000)
        _assert_control_value(page, "Comparison axis", diff_values["comparison_axis"])
        _assert_control_value(page, "Motifs", "\n".join(diff_values["motifs"]))
        _assert_control_value(page, "Signals", "\n".join(diff_values["signals"]))
        _assert_control_value(
            page,
            "Condition names",
            "\n".join(diff_values["cond_names"]),
        )
        _assert_control_value(
            page,
            "Region-set BED files",
            "\n".join(diff_values["regions"]),
        )
        _assert_control_value(
            page,
            "Region labels",
            ",".join(diff_values["region_labels"]),
        )
        _assert_control_value(
            page,
            "Matching-stratum BED column (0 = none)",
            diff_values["region_strata_column"],
        )
        _assert_control_value(page, "Cores", diff_values["cores"])
        _assert_control_value(page, "Skip Excel", True)
        if capture_screenshots:
            page.screenshot(
                path=str(output / "loaded-config-sync-diff-footprints.png"),
                full_page=False,
                animations="disabled",
                timeout=30_000,
            )
    finally:
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
        run_dir = (
            workdir_path
            / "very long workspace path"
            / "nested analysis project"
            / "fp-tools runs"
        )
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
                    for page_name in GUI_RUN_PAGES:
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

                    for page_name in SIMPLE_GUI_PAGES:
                        for width, height in ((1280, 720), (390, 844)):
                            _audit_simple_page(
                                browser,
                                base_url,
                                page_name,
                                width,
                                height,
                                output,
                                not args.skip_screenshots,
                            )

                    _audit_compact_sidebar(
                        browser,
                        base_url,
                        output,
                        not args.skip_screenshots,
                    )
                    _audit_mobile_sidebar_navigation(
                        browser,
                        base_url,
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

                    _audit_loaded_config_sync(
                        browser,
                        base_url,
                        workdir_path,
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
