#!/usr/bin/env python3
"""Browser regression for axis controls, display filters and immutable payloads."""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

from playwright.sync_api import expect, sync_playwright
from audit_report_plot_controls import fixture_payload, serve_report, downloaded_svg
from fp_tools.tools.static_comparison_browser import build_static_browser, write_embedded_static_browser
from fp_tools.tools.plot_aggregate_batch import _ensure_batch_payload, write_html


def fixtures(root, comparisons=5):
    root.mkdir(parents=True, exist_ok=True)
    payloads = []
    for i in range(comparisons):
        p = fixture_payload(True)
        p['title'] = f'Axis fixture {i}'
        p['conditions'] = [f'Dose{i}', 'Control']
        for motif in p['aggregate']['motifs']:
            for c, name in zip(motif['conditions'], p['conditions']):
                c['name'] = name
                for s in c['samples']:
                    s['profile'] = [v * (i+1) for v in s['profile']]
        p['points'][0]['change'] *= i+1
        p['points'].extend([
            dict(prefix='fail', name='FAIL', motif_id='F', change=.4, pvalue=.001, fdr=.5, neglog10p=3, group='KD_up'),
            dict(prefix='missing', name='MISSING', motif_id='M', change=.7, pvalue=.002, fdr=None, neglog10p=2.7, group='KD_up'),
        ])
        # The bundle's logo directory must contain its referenced test image.
        p['logos']['JUN_M1'] = {'png': 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAQAAAABCAIAAAB2XpiaAAAAFUlEQVR4nGOQbA/hzfv7dT7DHVNXACE6BVdtHw6KAAAAAElFTkSuQmCC'}
        payloads.append(p)
    review = dict(schema='fp-tools.review-multi-comparisons.v1', title='Axis audit',
                  comparisons=[dict(label=f'Dose {i}',payload=p) for i,p in enumerate(payloads)])
    standalone = root/'review.html'; write_embedded_static_browser(review,standalone)
    single=root/'individual.html'; write_embedded_static_browser(dict(review,comparisons=review['comparisons'][:1]),single,source_payload=payloads[0])
    bundle=build_static_browser(payloads, root/'bundle','Axis audit')
    batch=root/'aggregate.html'; write_html(_ensure_batch_payload(payloads[0]), batch,'Aggregate axes')
    return [standalone,single,bundle],batch


def change_number(locator, value):
    locator.fill(str(value)); locator.dispatch_event('change')


def open_editor(locator):
    locator.evaluate("el => el.closest('details').open = true")


def audit_report(browser, report, screenshots, exhaustive=False):
    page=browser.new_page(viewport={'width':1800,'height':1100})
    errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
    with serve_report(report) as url:
        page.goto(url,wait_until='load');expect(page.locator('#chart .pt').first).to_be_attached(timeout=90000)
        if page.locator('#report-layout').input_value()!='single':
            page.locator('#report-layout').select_option('single')
        original=page.evaluate('JSON.stringify(state.payload)')
        total=page.locator('#chart .pt').count()
        page.locator('#hide-failing').check()
        passing=page.locator('#chart .pt').count()
        assert passing<=total
        change_number(page.locator('#filter-delta'),1e9);page.locator('#apply-filters').click()
        expect(page.locator('#rank-chart')).to_contain_text('No motifs meet')
        expect(page.locator('#chart .pt')).to_have_count(0)
        change_number(page.locator('#filter-delta'),.1);page.locator('#apply-filters').click()
        page.locator('#hide-failing').uncheck()
        ymax=page.locator('#main-axis-controls input[aria-label="Volcano Y maximum slider"]');open_editor(ymax)
        before_y=json.loads(page.locator('#chart').get_attribute('data-y-range'))
        ymax.focus();ymax.press('ArrowRight')
        assert json.loads(page.locator('#chart').get_attribute('data-y-range'))[1]>before_y[1]
        page.get_by_role('button',name='Reset Volcano Y',exact=True).click()
        assert json.loads(page.locator('#chart').get_attribute('data-y-range'))==before_y
        xmin=page.locator('#main-axis-controls input[aria-label="Volcano X minimum"]');open_editor(xmin)
        change_number(xmin,-.25)
        assert json.loads(page.locator('#chart').get_attribute('data-x-range'))[0]==-.25
        previous=page.locator('#chart').get_attribute('data-x-range')
        change_number(xmin,1e20);assert page.locator('#chart').get_attribute('data-x-range')==previous
        expect(xmin).to_have_attribute('aria-invalid','true')
        change_number(xmin,-.25)
        svg=downloaded_svg(page,'#download-volcano');assert 'data-x-range="[-0.25,' in svg and 'clip-path=' in svg
        expect(page.locator('#aggregate-grid .aggregate-panel').first).to_be_attached()
        paths=page.locator('#aggregate-grid .aggregate-panel path').first.get_attribute('d')
        amin=page.locator('#aggregate-grid input[aria-label="Aggregate Y minimum"]').first;open_editor(amin)
        change_number(amin,-2)
        assert json.loads(page.locator('#aggregate-grid .aggregate-panel').first.get_attribute('data-y-range'))[0]==-2
        assert page.locator('#aggregate-grid .aggregate-panel path').first.get_attribute('d')!=paths
        page.locator('#shared-aggregateY').check()
        expect(page.locator('#shared-aggregate-controls')).to_be_visible()
        assert page.evaluate('JSON.stringify(state.payload)')==original
        # Filtering by raw p-value must admit the fixture's FDR-failing point.
        if not exhaustive:
            page.locator('#filter-metric').select_option('pvalue');page.locator('#apply-filters').click()
            assert page.locator('#rank-chart .rank-bar[data-prefix="fail"]').count()==1
            page.locator('#filter-metric').select_option('fdr');page.locator('#apply-filters').click()
            assert page.locator('#rank-chart .rank-bar[data-prefix="fail"]').count()==0
        count=page.evaluate('state.metadata.comparisons.length')
        if count>1:
            page.locator('#report-layout').select_option('side')
            expect(page.locator('.side-volcano')).to_have_count(min(4,count),timeout=90000)
            for kind in ['volcanoX','volcanoY','rankX']:
                page.locator(f'#shared-{kind}').check()
                expect(page.locator(f'#shared-{kind}')).to_be_enabled(timeout=90000)
            page.wait_for_function("document.querySelectorAll('.side-aggregate .aggregate-panel').length > 0")
            # Remove the earlier custom range so every volcano inherits shared X.
            page.evaluate("view.manual.clear(); redrawRanges()")
            domains=page.locator('.side-volcano').evaluate_all("nodes=>nodes.map(n=>n.dataset.xRange)")
            assert len(set(domains))==1,domains
            if count == 5 and not exhaustive:
                # The fifth comparison is hidden, but must determine the shared limit.
                assert json.loads(domains[0]) == [-4.2, 4.2], domains
            local=page.locator('.comparison-card').first.locator('input[aria-label="Volcano X minimum"]');open_editor(local)
            change_number(local,-.123)
            domains=page.locator('.side-volcano').evaluate_all("nodes=>nodes.map(n=>n.dataset.xRange)")
            assert json.loads(domains[0])[0]==-.123 and json.loads(domains[1])[0]!=-.123
            if not exhaustive:
                page.locator('.comparison-card > select').nth(1).select_option('0')
                expect(page.locator('#download-panel')).to_be_enabled()
                domains=page.locator('.side-volcano').evaluate_all("nodes=>nodes.map(n=>n.dataset.xRange)")
                assert json.loads(domains[0])[0]==-.123 and json.loads(domains[1])[0]!=-.123
            exported=downloaded_svg(page,'#download-panel');assert 'clip-path' in exported
            if exhaustive:
                # Visit all actual comparisons, not just the initial four slots.
                page.locator('#report-layout').select_option('single')
                for i in range(count):
                    page.locator('#comparison-selector').select_option(str(i))
                    page.wait_for_function('(i)=>state.comparisonIndex===i && document.querySelector("#chart .pt")',arg=i)
                    expect(page.locator('#aggregate-grid .aggregate-panel').first).to_be_attached(timeout=90000)
                    assert 'Could not' not in page.locator('#status').inner_text()
                page.locator('#report-layout').select_option('side')
                expect(page.locator('.side-volcano')).to_have_count(min(4,count),timeout=90000)
        if screenshots:
            screenshots.mkdir(parents=True,exist_ok=True)
            page.screenshot(path=str(screenshots/f'{report.parent.name}-{report.stem}-desktop.png'),full_page=True)
        page.set_viewport_size({'width':390,'height':844})
        if screenshots:page.screenshot(path=str(screenshots/f'{report.parent.name}-{report.stem}-mobile.png'),full_page=True)
        assert not errors, errors
    page.close(); print(f'PASS axis/filter/layout/export/payload browser checks: {report}',flush=True)


def audit_batch(browser, report):
    page=browser.new_page(viewport={'width':1600,'height':1000});errors=[]
    page.on('pageerror',lambda e:errors.append(str(e)))
    with serve_report(report) as url:
        page.goto(url);expect(page.locator('.aggregate-panel').first).to_be_attached()
        original=page.evaluate('JSON.stringify(payload)')
        editor=page.locator('input[aria-label="Aggregate Y minimum"]').first;open_editor(editor);change_number(editor,-3)
        assert json.loads(page.locator('.aggregate-panel').first.get_attribute('data-y-range'))[0]==-3
        page.locator('#batch-shared-y').check()
        assert 'clip-path' in downloaded_svg(page,'#download-grid')
        assert original==page.evaluate('JSON.stringify(payload)')
        assert not errors,errors
    page.close();print('PASS standalone aggregate sliders, sharing, export and payload',flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,required=True);parser.add_argument('reports',type=Path,nargs='*')
    args=parser.parse_args();reports,batch=fixtures(args.output_dir/'fixtures')
    with sync_playwright() as p:
        browser=p.chromium.launch()
        for report in reports:audit_report(browser,report,args.output_dir/'screenshots')
        audit_batch(browser,batch)
        for report in args.reports:audit_report(browser,report,args.output_dir/'screenshots',exhaustive=True)
        browser.close()


if __name__=='__main__':main()
