#!/usr/bin/env python3
"""Validate new control cuts and old-overlap expectation shape before rescreen."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
import build_all_site_control_dwm as builder
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import canonical_id,checked_path

ROOT=builder.ROOT


def compare_expectations(old,new):
    a,b=np.asarray(old,float),np.asarray(new,float)
    if a.shape!=b.shape or a.ndim!=2 or not np.isfinite(a).all() or not np.isfinite(b).all() or (a<0).any() or (b<0).any():
        raise ValueError('aligned finite nonnegative expectations required')
    at,bt=a.sum(1),b.sum(1)
    ap=np.divide(a,at[:,None],out=np.zeros_like(a),where=at[:,None]>0)
    bp=np.divide(b,bt[:,None],out=np.zeros_like(b),where=bt[:,None]>0)
    old_positive=at>0
    match=(bt>0)&np.isclose(ap,bp,rtol=1e-5,atol=1e-7).all(1)
    return old_positive,match,np.max(abs(ap-bp),axis=1)


def review(source,out):
    if out.exists():raise ValueError('fresh review directory required')
    builder.execute(source,verify=True)
    result=json.loads((source/'tracks_manifest.json').read_text())
    config=json.loads((source/'configuration.freeze.json').read_text())
    original=json.loads((builder.SOURCE/'configuration.freeze.json').read_text())
    rows=[];arrays={};safety=builder.audit.screen.prior.safety
    for item in config['items']:
        cell=item['cell'];entry=original['screen']['prior']['direct']['source_configuration']['cells'][cell]['naked_dwm']
        sites,old_observed,old_expected,_,_=safety.load_artifact(entry,builder.DATA,cell,control=True)
        fresh_sites,fresh=builder.audit.screen.prior.parent.read_profiles(original['screen']['prior']['direct']['corrected_stages'][cell]['naked'],ROOT)
        order=pd.Index(safety.site_hashes(fresh_sites)).get_indexer(safety.site_hashes(sites))
        if (order<0).any():raise ValueError('fresh site identity missing')
        observed,ov=safety.extract_prediction_profiles(sites,checked_path(result['outputs']['observed'],ROOT),100,require_dense=False)
        expected,ev=safety.extract_prediction_profiles(sites,checked_path(result['outputs']['expected'],ROOT),100,require_dense=False)
        old_positive,match,difference=compare_expectations(old_expected,expected)
        fresh_match=np.equal(observed,fresh['observed'][order]).all(1)
        old_comparable,_=builder.audit.screen.prior.parent.coverage.comparable_profiles(observed,old_observed,old_expected)
        table=sites[safety.SITE_KEYS].copy()
        for name,value in dict(coordinate_valid=ov&ev,fresh_observed_exact=fresh_match,
            fresh_valid=fresh['valid'][order],old_expected_positive=old_positive,new_expected_positive=expected.sum(1)>0,
            old_expected_shape_match=match,expected_proportion_max_error=difference,
            old_observed_shape_comparable=old_comparable,observed_total=observed.sum(1),expected_total=expected.sum(1)).items():table[name]=value
        rows.append(table);arrays[cell]=dict(observed=observed,expected=expected,valid=ov&ev,
            fresh_observed_exact=fresh_match,site_hash=safety.site_hashes(sites))
    table=pd.concat(rows,ignore_index=True)
    if len(table)!=4200:raise ValueError('control denominator changed')
    summary=table.groupby(['cell','tf']).agg(sites=('tf','size'),coordinate_valid=('coordinate_valid','sum'),
        fresh_observed_exact=('fresh_observed_exact','sum'),fresh_valid=('fresh_valid','sum'),
        old_positive=('old_expected_positive','sum'),new_positive=('new_expected_positive','sum'),
        old_shape_match=('old_expected_shape_match','sum'),old_observed_comparable=('old_observed_shape_comparable','sum')).reset_index()
    overlap=table.old_expected_positive
    qualified=bool(table.coordinate_valid.all() and table.fresh_valid.all() and table.fresh_observed_exact.all()
                   and table.loc[overlap,'old_expected_shape_match'].all())
    out.mkdir(parents=True);outputs={}
    for name,data in [('sites',table),('summary',summary)]:
        p=out/(name+'.tsv');data.to_csv(p,sep='\t',index=False);outputs[name]=record(p)
    for cell,values in arrays.items():
        p=out/(cell+'.npz');np.savez_compressed(p,**values);outputs[cell+'_profiles']=record(p)
    result=dict(source_configuration=record(source/'configuration.freeze.json'),source_result=record(source/'tracks_manifest.json'),
        code=record(Path(__file__)),outputs=outputs,old_positive_sites=int(overlap.sum()),
        old_expected_shape_failures=int((~table.loc[overlap,'old_expected_shape_match']).sum()),
        fresh_observed_failures=int((~table.fresh_observed_exact).sum()),qualified_for_rescreen=qualified,
        candidates_scored=False,promotion_eligible=False)
    result['stage_id']=canonical_id(result);write_immutable(out/'review_manifest.json',result)
    print(json.dumps({k:v for k,v in result.items() if k in ['old_positive_sites','old_expected_shape_failures','fresh_observed_failures','qualified_for_rescreen']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args();source=a.source.resolve();out=a.outdir.resolve()
    if not source.is_relative_to(ROOT/'benchmarks/results') or not out.is_relative_to(ROOT/'benchmarks/results'):p.error('use research results')
    review(source,out)
