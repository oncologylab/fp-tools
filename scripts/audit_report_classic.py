#!/usr/bin/env python3
"""Exercise the classic layout, sidebar controls, exports, and payload identity."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from playwright.sync_api import expect, sync_playwright
from audit_report_axes import fixtures, change_number, open_editor
from audit_report_plot_controls import serve_report, downloaded_svg


def ready(page):
    expect(page.locator('#classic-view')).to_have_attribute('data-ready', 'true', timeout=90000)


def audit(browser, report, output, actual=False):
    page = browser.new_page(viewport={'width':1920, 'height':1080})
    errors = []
    page.on('pageerror', lambda error: errors.append(str(error)))
    with serve_report(report) as url:
        page.goto(url, wait_until='load')
        page.wait_for_function('state.payload !== null', timeout=90000)
        original = page.evaluate('JSON.stringify(state.payload)')
        if page.locator('#report-layout').input_value() != 'classic':
            page.locator('#report-layout').select_option('classic')
        ready(page)
        root = page.locator('#classic-view')
        count = page.evaluate('state.metadata.comparisons.length')
        expect(root.locator('.comparison-card')).to_have_count(min(count, 8))
        bounds = page.evaluate('''() => {
            const get=s=>document.querySelector(s).getBoundingClientRect();
            const panel=get('#classic-controls'),logo=get('#classic-motif-logo'),plots=get('#classic-comparison-grid');
            return {below:panel.top>=logo.bottom,left:panel.right<=plots.left};
        }''')
        assert bounds == {'below':True, 'left':True}, bounds
        assert root.locator('.aggregate-panel path[data-sample]').count() > 0
        assert all(float(v) == .7 for v in root.locator('.aggregate-panel path[data-sample]').evaluate_all("nodes=>nodes.map(n=>n.getAttribute('stroke-width'))"))
        initial = output/f'{report.parent.name}-{report.stem}-initial.png'
        page.screenshot(path=str(initial), full_page=True)

        # One shared sidebar controls each panel independently, even duplicate comparisons.
        select = page.locator('#classic-edit-panel')
        select.evaluate("el=>el.closest('details').open=true")
        minimum = page.locator('#classic-ranges input[aria-label="Volcano X minimum"]')
        open_editor(minimum)
        change_number(minimum, -.123)
        ready(page)
        domains=root.locator('.volcano-svg').evaluate_all('nodes=>nodes.map(n=>JSON.parse(n.dataset.xRange))')
        assert domains[0][0] == -.123 and domains[1][0] != -.123
        change_number(minimum, 1e10)
        expect(minimum).to_have_attribute('aria-invalid','true')
        assert json.loads(root.locator('.volcano-svg').first.get_attribute('data-x-range'))[0] == -.123
        root.locator('[data-comparison-slot="1"]').select_option('0');ready(page)
        domains=root.locator('.volcano-svg').evaluate_all('nodes=>nodes.map(n=>JSON.parse(n.dataset.xRange))')
        assert domains[0][0] == -.123 and domains[1][0] != -.123
        open_editor(minimum)
        page.get_by_role('button',name='Reset Volcano X',exact=True).click()
        ready(page)
        for kind in ['volcanoX','volcanoY','rankX']:
            page.locator('#classic-shared-'+kind).check()
            ready(page)
        for selector,attribute in [('.volcano-svg','xRange'),('.volcano-svg','yRange'),('.rank-svg','xRange')]:
            domains=root.locator(selector).evaluate_all('(nodes,attr)=>nodes.map(n=>n.dataset[attr])',attribute)
            assert len(set(domains))==1,domains
        if not actual:
            page.locator('#classic-panel-count').select_option('4');ready(page)
            assert json.loads(root.locator('.volcano-svg').first.get_attribute('data-x-range'))[1] >= 4
        # Editing a slider changes the curve transform, never the observations.
        ymin=page.locator('#classic-ranges input[aria-label="Aggregate Y minimum"]')
        open_editor(ymin)
        change_number(ymin,-2);ready(page)
        assert json.loads(root.locator('.aggregate-panel').first.get_attribute('data-y-range'))[0]==-2
        slider=page.locator('#classic-ranges input[aria-label="Aggregate Y maximum slider"]')
        before=root.locator('.aggregate-panel').first.get_attribute('data-y-range')
        slider.focus();slider.press('ArrowRight');ready(page)
        assert root.locator('.aggregate-panel').first.get_attribute('data-y-range')!=before
        for button in ['rank','volcano','aggregate','panel','logo']:
            exported=downloaded_svg(page,'#classic-download-'+button)
            ET.fromstring(exported)
            if button!='logo':assert 'classic-clip-' in exported

        page.locator('#classic-rank-rows').fill('1')
        page.locator('#classic-rank-rows').dispatch_event('input');ready(page)
        assert root.locator('.comparison-card').first.locator('.rank-bar').count()<=1
        if not actual:
            prefix=page.locator('#classic-motif-select').input_value()
            page.locator('#classic-motif-select').select_option('fail');ready(page)
            expect(root.locator('.aggregate-panel').first).to_contain_text('No aggregate profile')
            page.locator('#classic-motif-select').select_option(prefix);ready(page)
            page.locator('#classic-filter-metric').select_option('pvalue')
            page.locator('#classic-apply').click();ready(page)
            page.locator('#classic-rank-rows').fill('20');page.locator('#classic-rank-rows').dispatch_event('input');ready(page)
            assert root.locator('.rank-bar[data-prefix="fail"]').count()>0
            page.locator('#classic-filter-metric').select_option('fdr');page.locator('#classic-apply').click();ready(page)
            assert root.locator('.rank-bar[data-prefix="fail"]').count()==0
        change_number(page.locator('#classic-delta'),1e9)
        page.locator('#classic-apply').click();ready(page)
        expect(root.locator('.rank-bar')).to_have_count(0)
        expect(root.locator('.pt.selected')).to_have_count(0)
        page.locator('#classic-hide').check();ready(page)
        expect(root.locator('.pt')).to_have_count(0)
        change_number(page.locator('#classic-delta'),.1);page.locator('#classic-apply').click();ready(page)
        page.locator('#classic-hide').uncheck();ready(page)
        if actual:
            for i in range(count):
                root.locator('[data-comparison-slot="0"]').select_option(str(i));ready(page)
                expect(root.locator('.volcano-svg').first.locator('.pt').first).to_be_attached()
                expect(root.locator('.aggregate-panel').first).to_be_attached()
        assert page.evaluate('JSON.stringify(state.payload)') == original
        # Switching views keeps the existing renderer usable.
        page.locator('#classic-layout').select_option('single')
        expect(page.locator('#chart .pt').first).to_be_attached()
        expect(root).to_be_hidden()
        page.locator('#report-layout').select_option('classic');ready(page)
        for width,height in [(1440,1000),(390,844)]:
            page.set_viewport_size({'width':width,'height':height})
            assert page.evaluate('document.documentElement.scrollWidth<=window.innerWidth'),width
            page.screenshot(path=str(output/f'{report.parent.name}-{report.stem}-{width}.png'),full_page=True)
        assert not errors,errors
    page.close()
    print(f'PASS classic layout, all controls, SVGs, view switching and payload: {report}',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('reports',type=Path,nargs='*')
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    generated,_=fixtures(args.output_dir/'fixtures')
    with sync_playwright() as pw:
        browser=pw.chromium.launch()
        for report in [generated[0],generated[2]]:audit(browser,report,args.output_dir)
        for report in args.reports:audit(browser,report,args.output_dir,actual=True)
        browser.close()


if __name__=='__main__':main()
