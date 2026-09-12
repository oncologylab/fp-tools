from copy import deepcopy

import pytest

from fp_tools.tools.plot_aggregate_batch import merge_payloads, write_html


def legacy_payload():
    return {
        "motif_matrices": {"TF1": [[10, 0], [0, 10], [0, 0], [0, 0]],
                           "unused": [[1], [1], [1], [1]]},
        "logos": {},
        "aggregate": {"x": [-1, 1], "motifs": [
            {"prefix": "TF1", "conditions": [
                {"name": "A", "profile": [0.1, 0.2], "samples": [
                    {"name": "A1", "profile": [0.1, 0.2]}]}]}]},
    }


def test_conversion_and_merge_preserve_selected_matrices(tmp_path):
    source = legacy_payload()
    original = deepcopy(source)
    converted = merge_payloads([source])
    assert converted["motif_matrices"] == {"TF1": source["motif_matrices"]["TF1"]}
    assert source == original
    second = deepcopy(converted)
    second["motif_matrices"] = {}
    merged = merge_payloads([converted, second])
    assert merged["motif_matrices"] == converted["motif_matrices"]
    assert merged["x"] == [-1, 1]
    assert [s["profile"] for s in merged["motifs"][0]["series"]] == [[0.1, 0.2]] * 4
    out = tmp_path / "converted.html"
    write_html(merged, out, "Converted")
    html = out.read_text()
    assert "function motifLogoSvg" in html
    assert "logoDataUri(prefix)" in html


def test_conversion_preserves_prerendered_logo_fallback():
    source = legacy_payload()
    source["logos"] = {"TF1": {"png": "data:image/png;base64,fixture"}}
    source["motif_matrices"] = {}
    converted = merge_payloads([source])
    assert converted["logos"] == source["logos"]


@pytest.mark.parametrize("matrix", [None, [], [[1], [2]], [[1], [2, 3], [0], [0]],
                                    [[float("nan")], [0], [0], [0]],
                                    [[-1], [0], [0], [0]], "invalid"])
def test_malformed_matrix_does_not_break_profiles(matrix):
    source = legacy_payload()
    source["motif_matrices"] = {"TF1": matrix}
    converted = merge_payloads([source])
    assert converted["motif_matrices"] == {}
    assert converted["motifs"][0]["series"][0]["profile"] == [0.1, 0.2]


@pytest.fixture(scope="module")
def browser():
    playwright = pytest.importorskip("playwright.sync_api")
    with playwright.sync_playwright() as runtime:
        try:
            instance = runtime.chromium.launch()
        except playwright.Error as exc:
            if "Executable doesn't exist" in str(exc):
                pytest.skip("Playwright Chromium is not installed")
            raise
        yield instance
        instance.close()


@pytest.mark.parametrize("kind", ["matrix", "png", "svg", "missing", "malformed"])
def test_browser_logos_exports_and_profiles(tmp_path, browser, kind):
    source = legacy_payload()
    if kind != "matrix":
        source["motif_matrices"] = {"TF1": [[1], []]} if kind == "malformed" else {}
    if kind == "png":
        source["logos"] = {"TF1": {"png": "data:image/png;base64,"
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+jRZkAAAAASUVORK5CYII="}}
    if kind == "svg":
        source["logos"] = {"TF1": {"svg": '<svg xmlns="http://www.w3.org/2000/svg"><text>A</text></svg>'}}
    out = tmp_path / "converted.html"
    write_html(source, out, "Converted", default_layout="1x1")
    page = browser.new_page(accept_downloads=True)
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
    try:
        page.goto(out.as_uri())
        page.wait_for_selector(".motif-col .aggregate-panel")
        assert page.evaluate("payload.motifs[0].series[0].profile") == [0.1, 0.2]
        if kind in {"matrix", "svg"}:
            assert page.locator(".motif-logo svg").count() > 0
            assert page.locator(".motif-logo svg").first.evaluate(
                "svg => svg.getBoundingClientRect().height <= svg.parentElement.clientHeight"
            )
        elif kind == "png":
            assert page.locator(".motif-logo img").first.evaluate("img => img.complete && img.naturalWidth > 0")
        else:
            assert "Motif logo unavailable" in page.locator(".motif-logo").first.inner_text()
        with page.expect_download() as download:
            page.click("#download-logo")
        from pathlib import Path
        svg = Path(download.value.path()).read_text()
        if kind in {"matrix", "svg", "png"}:
            assert "<image " in svg
            assert "Motif logo unavailable" not in svg
        else:
            assert "Motif logo unavailable" in svg
        assert not errors
    finally:
        page.close()
