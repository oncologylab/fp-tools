"""The desktop audit must edit forms without an implicit extra submission."""

import pytest

playwright = pytest.importorskip("playwright.sync_api")

from scripts.audit_desktop_gui import _select_example, _submit_text_control


def test_editing_text_waits_for_the_explicit_submit_button():
    with playwright.sync_playwright() as runtime:
        try:
            browser = runtime.chromium.launch()
        except playwright.Error as error:
            if "Executable doesn't exist" in str(error):
                pytest.skip("Playwright Chromium is not installed")
            raise
        try:
            page = browser.new_page()
            page.set_content('''
                <script>window.submissions = 0;</script>
                <form onsubmit="event.preventDefault(); window.submissions++">
                  <label>Output directory<input name="outdir"></label>
                  <button type="submit">Update page config</button>
                </form>
            ''')
            _submit_text_control(page, "Output directory", "changed output")
            assert page.get_by_label("Output directory").input_value() == "changed output"
            assert page.evaluate("window.submissions") == 0
            page.get_by_role("button", name="Update page config").click()
            assert page.evaluate("window.submissions") == 1
        finally:
            browser.close()


def test_example_selection_opens_a_click_triggered_dropdown():
    with playwright.sync_playwright() as runtime:
        try:
            browser = runtime.chromium.launch()
        except playwright.Error as error:
            if "Executable doesn't exist" in str(error):
                pytest.skip("Playwright Chromium is not installed")
            raise
        try:
            page = browser.new_page()
            page.set_default_timeout(1000)
            page.set_content('''
                <input aria-label="Example YAML" role="combobox"
                  onclick="document.querySelector('[role=listbox]').hidden=false">
                <div role="listbox" hidden>
                  <div role="option" onclick="window.selected=true;
                    document.querySelector('input').value=this.textContent;
                    this.parentElement.hidden=true">normalize_bigwig_single.yml</div>
                </div>
            ''')
            _select_example(page, "normalize_bigwig_single.yml")
            assert page.evaluate("window.selected") is True
            assert page.get_by_label("Example YAML").input_value() == "normalize_bigwig_single.yml"
        finally:
            browser.close()
