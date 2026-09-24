#!/usr/bin/env python3
"""Exercise the classic layout, sidebar controls, exports, and payload identity."""
from __future__ import annotations

import argparse
from contextlib import nullcontext
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from playwright.sync_api import expect, sync_playwright
from audit_report_axes import fixtures, change_number, open_editor
from audit_report_plot_controls import serve_report, downloaded_svg


def ready(page):
    expect(page.locator('#classic-view')).to_have_attribute('data-ready', 'true', timeout=90000)


AXES = [
    ('Volcano X', '.volcano-svg', 'xRange', [-.123, .456]),
    ('Volcano Y', '.volcano-svg', 'yRange', [.01, 2.4]),
    ('Waterfall X', '.rank-svg', 'xRange', [-.321, .654]),
    ('Aggregate Y', '.aggregate-panel', 'yRange', [-2, 2]),
]


def domains(root, selector, attribute):
    return root.locator(selector).evaluate_all(
        '(nodes,attr)=>nodes.map(n=>JSON.parse(n.dataset[attr]))', attribute)


def set_limits(page, label, limits):
    minimum = page.locator(f'#classic-ranges input[aria-label="{label} minimum"]')
    maximum = page.locator(f'#classic-ranges input[aria-label="{label} maximum"]')
    open_editor(minimum)
    minimum.fill(str(limits[0])); maximum.fill(str(limits[1]))
    maximum.dispatch_event('change'); ready(page)


def auto_scale(page, label):
    control = page.locator(f'#classic-ranges input[aria-label="{label} minimum"]')
    open_editor(control)
    page.get_by_role('button', name='Auto scale '+label, exact=True).click()
    ready(page)


def audit(browser, report, output, actual=False, local_file=False):
    page = browser.new_page(viewport={'width':1920, 'height':1080})
    errors = []
    page.on('pageerror', lambda error: errors.append(str(error)))
    page.on('console', lambda message: errors.append(message.text) if message.type == 'error' else None)
    page.on('requestfailed', lambda request: errors.append(request.url+': '+str(request.failure)))
    page.on('response', lambda response: errors.append(f'HTTP {response.status}: {response.url}') if response.status >= 400 else None)
    with (nullcontext(report.resolve().as_uri()) if local_file else serve_report(report)) as url:
        page.goto(url, wait_until='load')
        page.wait_for_function('state.payload !== null', timeout=90000)
        original = page.evaluate('JSON.stringify(state.payload)')
        if page.locator('#report-layout').input_value() != 'classic':
            page.locator('#report-layout').select_option('classic')
        ready(page)
        root = page.locator('#classic-view')
        count = page.evaluate('state.metadata.comparisons.length')
        expect(root.locator('.comparison-card')).to_have_count(min(count, 8))
        expect(page.locator('#classic-edit-panel')).to_have_value('all')
        expect(page.locator('#classic-edit-panel option')).to_have_count(count + 1)
        expect(root.locator('[id^="classic-shared-"]')).to_have_count(0)
        if not actual:
            # Last (largest) comparison is never initially displayed.
            assert count > 8
            assert domains(root,'.aggregate-panel','yRange')[0][1] >= .23 * count
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

        # All edits are manual group edits, not merely matching autoscale checkboxes.
        select = page.locator('#classic-edit-panel')
        select.evaluate("el=>el.closest('details').open=true")
        for label, selector, attribute, limits in AXES:
            set_limits(page, label, limits)
            assert all(value == limits for value in domains(root, selector, attribute)), label
        # A hidden comparison can be edited, and All must clear that exception too.
        if count > 8:
            select.select_option(str(count-1))
            set_limits(page,'Volcano X',[-.777,.888])
            assert all(value == AXES[0][3] for value in domains(root,'.volcano-svg','xRange'))
            select.select_option('all')
            set_limits(page,'Volcano X',AXES[0][3])
        # A comparison opened later inherits all four limits.
        root.locator('[data-comparison-slot="0"]').select_option(str(count-1));ready(page)
        for label, selector, attribute, limits in AXES:
            assert domains(root, selector, attribute)[0] == limits, label
        root.locator('[data-comparison-slot="0"]').select_option('0');ready(page)
        # Individual changes follow comparison identity, including duplicate displays.
        select.select_option('0')
        for label, selector, attribute, limits in AXES:
            local = [limits[0], limits[1] + .25]
            set_limits(page, label, local)
            current = domains(root, selector, attribute)
            assert current[0] == local and current[1] == limits, label
        root.locator('[data-comparison-slot="1"]').select_option('0');ready(page)
        for label, selector, attribute, limits in AXES:
            assert domains(root, selector, attribute)[:2] == [[limits[0], limits[1]+.25]] * 2
        # Invalid and empty inputs preserve the last valid range.
        minimum = page.locator('#classic-ranges input[aria-label="Volcano X minimum"]')
        open_editor(minimum)
        change_number(minimum, 1e10)
        expect(minimum).to_have_attribute('aria-invalid','true')
        minimum.fill('');minimum.dispatch_event('change')
        expect(minimum).to_have_attribute('aria-invalid','true')
        assert domains(root,'.volcano-svg','xRange')[0] == [-.123,.706]
        select.select_option('all')
        # All replaces earlier exceptions only for the edited axis.
        set_limits(page, 'Volcano X', [-1,1])
        assert all(value == [-1,1] for value in domains(root,'.volcano-svg','xRange'))
        assert domains(root,'.volcano-svg','yRange')[0] == [.01,2.65]
        for label, selector, attribute, limits in AXES:
            set_limits(page,label,limits)
            assert all(value == limits for value in domains(root,selector,attribute))
            slider=page.locator(f'#classic-ranges input[aria-label="{label} maximum slider"]')
            slider.focus();slider.press('ArrowRight');ready(page)
            values=domains(root,selector,attribute)
            assert all(value == values[0] for value in values) and values[0] != limits
        # Autoscale is scoped: individual first, then All clears the exception.
        page.locator('#classic-ranges').evaluate('''box => {
            const low=box.querySelector('[aria-label="Volcano X minimum"]');
            const high=box.querySelector('[aria-label="Volcano X maximum"]');
            low.value='-.123';
            for(let i=1;i<=10;i++){high.value=i;high.dispatchEvent(new Event('change'));}
        }''')
        ready(page)
        assert all(value == [-.123,10] for value in domains(root,'.volcano-svg','xRange'))
        select.select_option('0')
        for label, selector, attribute, limits in AXES:
            auto_scale(page,label)
            values=domains(root,selector,attribute)
            assert values[0] == values[1] and values[0] != values[2],label
        select.select_option('all')
        for label, selector, attribute, limits in AXES:
            auto_scale(page,label)
            values=domains(root,selector,attribute)
            assert all(value == values[0] for value in values),label
        if not actual:
            page.locator('#classic-panel-count').select_option('4');ready(page)
            assert json.loads(root.locator('.volcano-svg').first.get_attribute('data-x-range'))[1] >= 4
            # The largest profile starts hidden, including in lazy-loaded bundles.
            before=domains(root,'.aggregate-panel','yRange')[0]
            root.locator('[data-comparison-slot="0"]').select_option(str(count-1));ready(page)
            assert domains(root,'.aggregate-panel','yRange')[0] == before
            assert 'Range clips data' not in root.locator('.aggregate-tile').first.inner_text()
            root.locator('[data-comparison-slot="0"]').select_option('0');ready(page)
        for label, selector, attribute, limits in AXES:
            set_limits(page,label,limits)
        assert all(root.locator('.rank-svg').evaluate_all('''nodes=>nodes.map(svg=>{
            const note=[...svg.querySelectorAll('text')].find(n=>n.textContent.includes('Range clips data'));
            if(!note)return true;
            const bounds=note.getBBox();
            return [...svg.querySelectorAll('.tick')].every(t=>bounds.y+bounds.height<t.getBBox().y);
        })''')), 'Waterfall clipping notice overlaps axis ticks'
        for button in ['rank','volcano','aggregate','panel','logo']:
            exported=downloaded_svg(page,'#classic-download-'+button)
            ET.fromstring(exported)
            if button!='logo':assert 'classic-clip-' in exported
            selectors={'rank':'.rank-svg','volcano':'.volcano-svg','aggregate':'.aggregate-panel',
                       'panel':'.rank-svg,.volcano-svg,.aggregate-panel'}
            if button in selectors:
                for markup in root.locator(selectors[button]).evaluate_all('nodes=>nodes.map(n=>n.innerHTML)'):
                    assert markup in exported, f'{button} export differs from displayed plots'
            (output/f'{report.parent.name}-{report.stem}-{button}.svg').write_text(exported)

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
                for label, selector, attribute, limits in AXES:
                    assert domains(root,selector,attribute)[0] == limits,(i,label)
        assert page.evaluate('JSON.stringify(state.payload)') == original
        # Switching views keeps the existing renderer usable.
        page.locator('#classic-layout').select_option('single')
        expect(page.locator('#chart .pt').first).to_be_attached()
        expect(root).to_be_hidden()
        page.locator('#report-layout').select_option('classic');ready(page)
        for label, selector, attribute, limits in AXES:
            assert all(value == limits for value in domains(root,selector,attribute))
        # Leave expanded controls visible for the visual review, with readable autoscaling.
        for label, selector, attribute, limits in AXES:
            auto_scale(page,label)
        page.locator('#classic-ranges details').evaluate_all('nodes=>nodes.forEach(n=>n.open=true)')
        page.locator('#classic-controls').screenshot(path=str(output/f'{report.parent.name}-{report.stem}-controls.png'))
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
    generated,_=fixtures(args.output_dir/'fixtures', comparisons=9)
    with sync_playwright() as pw:
        browser=pw.chromium.launch()
        for report in [generated[0],generated[2]]:audit(browser,report,args.output_dir)
        for report in args.reports:audit(browser,report,args.output_dir,actual=True)
        browser.close()


if __name__=='__main__':main()
