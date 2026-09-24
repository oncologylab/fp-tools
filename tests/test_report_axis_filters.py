"""Numerical display contracts: filtering never changes the scientific payload."""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

CONTROLS = Path(__file__).parents[1] / 'src/fp_tools/resources/static_browser/plot_controls.js'


def js(expression):
    if not shutil.which('node'):
        pytest.skip('Node.js required for browser helper contracts')
    result = subprocess.run(['node', '-e',
        'const c=require(process.argv[1]); console.log(JSON.stringify(' + expression + '));',
        str(CONTROLS)], check=True, capture_output=True, text=True)
    return json.loads(result.stdout)


def test_filter_uses_both_inclusive_thresholds_and_rejects_missing_values():
    assert js("[{}, {change:.1,fdr:.05}, {change:-.1,fdr:.05}, {change:.09,fdr:0}, {change:.2,fdr:null}, {change:.2,fdr:''}, {change:.2,fdr:2}, {change:.2,fdr:0}].map(p=>c.passesDisplayFilter(p,{metric:'fdr',alpha:.05,delta:.1}))") == [False, True, True, False, False, False, False, True]
    assert js("c.passesDisplayFilter({change:.2,pvalue:.01,fdr:.2},{metric:'pvalue',alpha:.05,delta:.1})")


def test_ranges_require_finite_ordered_endpoints():
    assert js("[[0,1],[1,1],[2,1],[null,1],['',1],[-2,3]].map(c.validRange)") == [True, False, False, False, False, True]
    assert js("c.autoDomain([-2,1], 'symmetric')") == [-2.1, 2.1]
    assert js("c.autoDomain([0,0], 'positive')") == [0, 1]
    assert js("c.autoDomain([1,40], 'positive')") == [0, 42]


def test_ranking_fills_one_sided_slots_and_honors_one_row():
    expression = "c.rankMotifs([{prefix:'a',effect:-.2},{prefix:'b',effect:-.5},{prefix:'c',effect:-.3}], 'effect', 2)"
    assert [p['prefix'] for p in js(expression)['negative']] == ['b', 'c']
    assert len(js(expression.replace(', 2)', ', 1)'))['negative']) == 1


def test_filter_does_not_mutate_input():
    assert js("(()=>{const p={change:.2,pvalue:.01,fdr:.2,group:'original'}; const before=JSON.stringify(p); c.passesDisplayFilter(p,{metric:'fdr',alpha:.05,delta:.1}); return before===JSON.stringify(p)})()")


def test_rebuild_preserves_complete_payload_and_original_file(tmp_path):
    from fp_tools.tools.review_multi_comparisons import _compressed_json_b64
    from scripts.rebuild_report_html import rebuild, sha256, read_report_payload
    payload = {'title':'Original', 'conditions':['A','B'], 'groups':['n.s.'],
               'points':[{'prefix':'M','change':.2,'pvalue':.04,'fdr':.1,'group':'n.s.'}],
               'aggregate':{'x':[-1,0,1],'motifs':[]}, 'motif_matrices':{},
               'results_tsv':'original\tcolumns\n1\t2\n'}
    source=tmp_path/'old.html'; output=tmp_path/'new.html'
    source.write_text(f'<script>const reportPayloadB64="{_compressed_json_b64(payload)}";</script>')
    before=sha256(source)
    receipt=rebuild(source,output)
    assert receipt['source_sha256']==before==sha256(source)
    assert read_report_payload(output)==payload
    assert receipt['payload_equal'] is True


def test_bundle_and_standalone_ship_same_options(tmp_path):
    from fp_tools.tools.static_comparison_browser import write_embedded_static_browser
    p={'conditions':['A','B'],'points':[{'prefix':'a'}],'aggregate':{'motifs':[]}}
    review={'schema':'fp-tools.review-multi-comparisons.v1','comparisons':[{'payload':p}]}
    target=tmp_path/'report.html'
    write_embedded_static_browser(review,target)
    text=target.read_text()
    assert 'const view = ' in text and 'axisEditor' in text
    assert '<script src=' not in text
    with pytest.raises(ValueError,match='default_view'):
        write_embedded_static_browser(review,target,default_view='unknown')
    from fp_tools.tools.review_multi_comparisons import _compressed_json_b64
    with pytest.raises(ValueError,match='encoded_payload'):
        write_embedded_static_browser(review,target,encoded_payload=_compressed_json_b64({'different':True}))


def test_classic_is_an_optional_embedded_view_with_unchanged_payload(tmp_path):
    from fp_tools.tools.static_comparison_browser import write_embedded_static_browser
    from scripts.rebuild_report_html import read_report_payload
    payload = {'conditions':['A','B'], 'points':[{'prefix':'M','change':.2,'pvalue':.01,'fdr':.02}]}
    review = {'schema':'fp-tools.review-multi-comparisons.v1','comparisons':[{'payload':payload}]}
    path = tmp_path/'classic.html'
    write_embedded_static_browser(review,path,default_view='classic')
    text = path.read_text()
    assert 'defaultView:"classic"' in text
    assert 'createClassicView' in text
    assert '<script src=' not in text
    assert read_report_payload(path) == review
