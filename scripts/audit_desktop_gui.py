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
            [str(executable), "--no-browser", "--port", str(port), "--run-dir", str(run_dir)],
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
