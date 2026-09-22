#!/usr/bin/env python3
"""Matched historical-region construction isolates control output-layout drift."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
import build_all_site_control_dwm as builder
from review_all_site_control_dwm import compare_expectations
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import canonical_id,checked_path

ROOT=builder.ROOT
SOURCE=ROOT/'benchmarks/results/all_site_control_dwm_20260922_v1'
EXPANDED_REVIEW=ROOT/'benchmarks/results/all_site_control_dwm_review_20260922_v1'
PROTOCOL=dict(evidence='previously_examined_replicate2_region_layout_diagnostic',
    change='historical output regions only; current implementation, input files, model and normalization held fixed',
    expectation_rtol=1e-5,expectation_atol=1e-7,candidates_scored=False,
    thresholds_changed=False,promotion_eligible=False)


def historical_command(command,out):
    result=list(command)
    for option,value in [('--regions-out',result[result.index('--peaks')+1]),
                         ('--outdir',str(out/'tracks')),('--prefix','NakedDNA_rep2_historical')]:
        if result.count(option)!=1:raise ValueError('unique command option required')
        result[result.index(option)+1]=value
    return result


def freeze(out):
    if out.exists():raise ValueError('fresh directory required')
    builder.execute(SOURCE,verify=True)
    old=json.loads((SOURCE/'configuration.freeze.json').read_text())
    expanded=json.loads((EXPANDED_REVIEW/'review_manifest.json').read_text())
    inputs=old['inputs']+[record(SOURCE/'configuration.freeze.json'),record(SOURCE/'tracks_manifest.json'),
        record(EXPANDED_REVIEW/'review_manifest.json')]+list(expanded['outputs'].values())
    for item in inputs:checked_path(item,ROOT)
    value=dict(protocol=PROTOCOL,inputs=inputs,code=old['code']+[record(Path(__file__)),record(Path(compare_expectations.__code__.co_filename))],
        command=historical_command(old['command'],out),original_model=old['original_model'],
        expanded=expanded['outputs'],runtime=builder.audit.screen.prior.builder.extension.runtime())
    value['configuration_id']=canonical_id(value);out.mkdir(parents=True);write_immutable(out/'configuration.freeze.json',value)


def preflight(out):
    ext=builder.audit.screen.prior.builder.extension
    config=ext.read_document(out/'configuration.freeze.json','configuration_id')
    if config['protocol']!=PROTOCOL or config['runtime']!=ext.runtime():raise ValueError('protocol or runtime changed')
    for item in config['inputs']+config['code']:checked_path(item,ROOT)
    return config


def construct(out):
    config=preflight(out);target=out/'tracks_manifest.json'
    if target.exists():
        builder.audit.screen.prior.builder.base.verify_stage(target,config['configuration_id'],ROOT);return
    if (out/'tracks').exists():raise ValueError('partial construction preserved')
    subprocess.run(config['command'],check=True,cwd=ROOT)
    prefix=out/'tracks/NakedDNA_rep2_historical'
    old=builder.AtacBias().from_pickle(str(checked_path(config['original_model'],ROOT)))
    new=builder.AtacBias().from_pickle(str(prefix)+'_AtacBias.pickle')
    if new.stype!=old.stype or new.correction_factor!=1.:raise ValueError('model settings changed')
    for strand in old.bias:
        for attr in ['counts','bg_counts']:np.testing.assert_array_equal(getattr(old.bias[strand],attr),getattr(new.bias[strand],attr))
    result=dict(configuration_id=config['configuration_id'],outputs={k:record(Path(str(prefix)+v)) for k,v in
        [('observed','_uncorrected.bw'),('expected','_expected.bw'),('bias_model','_AtacBias.pickle')]},model_refitted=False,candidates_scored=False)
    result['stage_id']=canonical_id(result);write_immutable(target,result)


def review(out,verify=False):
    config=preflight(out);base=builder.audit.screen.prior.builder.base
    tracks=base.verify_stage(out/'tracks_manifest.json',config['configuration_id'],ROOT)
    target=out/'layout_manifest.json'
    if target.exists():
        base.verify_stage(target,config['configuration_id'],ROOT);print('Verified layout isolation diagnostic',flush=True);return
    if verify or (out/'sites.tsv').exists():raise ValueError('missing or partial review')
    old_config=json.loads((builder.SOURCE/'configuration.freeze.json').read_text())
    safety=builder.audit.screen.prior.safety;rows=[]
    for cell in ['HepG2','K562']:
        entry=old_config['screen']['prior']['direct']['source_configuration']['cells'][cell]['naked_dwm']
        sites,old_observed,old_expected,_,_=safety.load_artifact(entry,builder.DATA,cell,control=True)
        observed,ov=safety.extract_prediction_profiles(sites,checked_path(tracks['outputs']['observed'],ROOT),100,require_dense=False)
        expected,ev=safety.extract_prediction_profiles(sites,checked_path(tracks['outputs']['expected'],ROOT),100,require_dense=False)
        with np.load(checked_path(config['expanded'][cell+'_profiles'],ROOT),allow_pickle=False) as archive:
            np.testing.assert_array_equal(archive['site_hash'],safety.site_hashes(sites));expanded=archive['expected']
        positive,old_match,old_error=compare_expectations(old_expected,expected)
        current_positive,expanded_match,expanded_error=compare_expectations(expected,expanded)
        observed_positive,observed_match,_=compare_expectations(old_observed,observed)
        observed_match |= ~observed_positive & (observed.sum(1)==0)
        table=sites[safety.SITE_KEYS].copy()
        for name,value in dict(coordinate_valid=ov&ev,old_expected_positive=positive,current_expected_positive=current_positive,
            old_expected_shape_match=old_match,old_expected_max_error=old_error,
            expanded_expected_shape_match=expanded_match,expanded_expected_max_error=expanded_error,
            archived_observed_shape_match=observed_match).items():table[name]=value
        rows.append(table)
    data=pd.concat(rows,ignore_index=True)
    if len(data)!=4200:raise ValueError('site denominator changed')
    positive=data.old_expected_positive
    result=dict(configuration_id=config['configuration_id'],sites=4200,old_positive_sites=int(positive.sum()),
        historical_shape_failures=int((~data.loc[positive,'old_expected_shape_match']).sum()),
        historical_support_changes=int((positive!=data.current_expected_positive).sum()),
        historical_observed_shape_failures=int((~data.archived_observed_shape_match).sum()),
        expanded_shape_failures=int((~data.loc[data.current_expected_positive,'expanded_expected_shape_match']).sum()),
        candidates_scored=False,promotion_eligible=False)
    path=out/'sites.tsv';data.to_csv(path,sep='\t',index=False);result['outputs']={'sites':record(path)}
    result['stage_id']=canonical_id(result);write_immutable(target,result);print(json.dumps({k:v for k,v in result.items() if isinstance(v,(int,bool))}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--stage',choices=['freeze','construct','review','verify'],required=True)
    a=p.parse_args();out=a.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'):p.error('use research results')
    if a.stage=='freeze':freeze(out)
    else:
        with builder.audit.screen.prior.builder.extension.exclusive_run(out):
            if a.stage=='construct':construct(out)
            else:review(out,a.stage=='verify')
