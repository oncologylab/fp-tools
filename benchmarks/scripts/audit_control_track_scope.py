#!/usr/bin/env python3
"""Audit archived replicate-2 track coverage without scoring a candidate."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pyBigWig

sys.path.insert(0,str(Path(__file__).resolve().parent))
import screen_stable_geometry_data_root as screen
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import canonical_id,checked_path

ROOT=screen.ROOT
SOURCE=ROOT/'benchmarks/results/stable_geometry_screen_20260922_v2'
BED=screen.DATA_ROOT/'benchmarks/results/footprint_functional_v1/naked_dna_10m_s2026/ctcf_exact_geometry_safety/ctcf_naked_sites.bed'
PROTOCOL=dict(evidence='previously_examined_replicate2_coverage_only',width=201,
    distinguish=['absent_track_entries','recorded_zero_only','positive_expectation'],
    candidates_scored=False,thresholds_changed=False,final_control_read=False,promotion_eligible=False)


def window_support(values):
    v=np.asarray(values,float)
    if v.shape!=(201,) or np.isinf(v).any() or np.any(v[np.isfinite(v)]<0):
        raise ValueError('201 nonnegative track values or absent NaNs required')
    finite=np.isfinite(v);positive=finite & (v>0)
    return dict(finite_bases=int(finite.sum()),positive_bases=int(positive.sum()),
        category='absent_track_entries' if not finite.any() else 'positive_expectation' if positive.any() else 'recorded_zero_only')


def freeze(out):
    if out.exists():raise ValueError('fresh directory required')
    config,_=screen.preflight(SOURCE)
    documents={};inputs=config['inputs']+[record(SOURCE/'configuration.freeze.json'),record(BED)]
    for cell in ['HepG2','K562']:
        entry=config['prior']['direct']['source_configuration']['cells'][cell]['naked_dwm']
        path=checked_path(entry,screen.DATA_ROOT);doc=json.loads(path.read_text());documents[cell]=doc;inputs.append(record(path))
        for key in ['observed','expected']:
            inputs.append(record(checked_path(dict(path=doc['metadata'][key],sha256=doc['metadata'][key+'_sha256']),screen.DATA_ROOT)))
    value=dict(protocol=PROTOCOL,screen=config,documents=documents,bed=record(BED),inputs=inputs,
        code=config['code']+[record(Path(__file__))],runtime=screen.prior.builder.extension.runtime())
    value['configuration_id']=canonical_id(value);out.mkdir(parents=True);write_immutable(out/'configuration.freeze.json',value)


def execute(out,verify=False):
    config=screen.prior.builder.extension.read_document(out/'configuration.freeze.json','configuration_id')
    if config['protocol']!=PROTOCOL or config['runtime']!=screen.prior.builder.extension.runtime():raise ValueError('scope or runtime changed')
    for r in config['inputs']+config['code']:checked_path(r,ROOT)
    target=out/'coverage_manifest.json'
    if target.exists():
        screen.prior.builder.base.verify_stage(target,config['configuration_id'],ROOT);print('Verified control-track scope audit',flush=True);return
    if verify or (out/'sites.tsv').exists():raise ValueError('missing or partial output')
    bed=pd.read_csv(config['bed']['path'],sep='\t',header=None,usecols=[0,1,2],names=['chrom','start','end'])
    anchors=set(bed.itertuples(index=False,name=None));rows=[]
    for cell in ['HepG2','K562']:
        entry=config['screen']['prior']['direct']['source_configuration']['cells'][cell]['naked_dwm']
        sites,observed,expected,valid,meta=screen.prior.safety.load_artifact(entry,screen.DATA_ROOT,cell,control=True)
        fresh,arrays=screen.prior.parent.read_profiles(config['screen']['prior']['direct']['corrected_stages'][cell]['naked'],ROOT)
        order=pd.Index(screen.prior.safety.site_hashes(fresh)).get_indexer(screen.prior.safety.site_hashes(sites))
        if (order<0).any():raise ValueError('fresh control alignment missing')
        actual=arrays['observed'][order]
        comparable,disagreement=screen.prior.parent.coverage.comparable_profiles(actual,observed,expected)
        expected_path=checked_path(dict(path=meta['expected'],sha256=meta['expected_sha256']),screen.DATA_ROOT)
        observed_path=checked_path(dict(path=meta['observed'],sha256=meta['observed_sha256']),screen.DATA_ROOT)
        with pyBigWig.open(str(expected_path)) as exp,pyBigWig.open(str(observed_path)) as obs:
            for i,s in enumerate(sites.itertuples(index=False)):
                center=(int(s.TFBS_start)+int(s.TFBS_end))//2
                ev=np.asarray(exp.values(s.TFBS_chr,center-100,center+101,numpy=True),float)
                ov=np.asarray(obs.values(s.TFBS_chr,center-100,center+101,numpy=True),float)
                if s.TFBS_strand=='-':ev=ev[::-1];ov=ov[::-1]
                np.testing.assert_array_equal(np.nan_to_num(ev),expected[i]);np.testing.assert_array_equal(np.nan_to_num(ov),observed[i])
                row={key:getattr(s,key) for key in screen.prior.safety.SITE_KEYS}
                row.update(window_support(ev),archived_coordinate_valid=bool(valid[i]),observed_finite_bases=int(np.isfinite(ov).sum()),
                    archived_observed_total=float(observed[i].sum()),fresh_observed_total=float(actual[i].sum()),
                    expected_total=float(expected[i].sum()),comparable=bool(comparable[i]),shape_disagreement=bool(disagreement[i]),
                    in_ctcf_reference_bed=(s.TFBS_chr,int(s.TFBS_start),int(s.TFBS_end)) in anchors)
                rows.append(row)
    table=pd.DataFrame(rows);summaries=[]
    for (cell,tf),g in table.groupby(['cell','tf']):
        summaries.append(dict(cell=cell,tf=tf,sites=len(g),coordinate_valid=int(g.archived_coordinate_valid.sum()),
            absent_entries=int(g.category.eq('absent_track_entries').sum()),recorded_zero_only=int(g.category.eq('recorded_zero_only').sum()),
            positive_expectation=int(g.category.eq('positive_expectation').sum()),fully_recorded=int(g.finite_bases.eq(201).sum()),
            fresh_informative=int(g.fresh_observed_total.gt(0).sum()),comparable=int(g.comparable.sum()),
            shape_disagreement=int(g.shape_disagreement.sum()),in_ctcf_reference_bed=int(g.in_ctcf_reference_bed.sum())))
    if len(table)!=4200 or len(summaries)!=21:raise ValueError('control denominator changed')
    outputs={}
    for name,data in [('sites',table),('summary',pd.DataFrame(summaries))]:
        p=out/(name+'.tsv');data.to_csv(p,sep='\t',index=False);outputs[name]=record(p)
    result=dict(configuration_id=config['configuration_id'],outputs=outputs,coverage_only=True,
        candidates_scored=False,promotion_eligible=False)
    result['stage_id']=canonical_id(result);write_immutable(target,result);print('Audited all 4200 sites and 21 tasks',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',type=Path,required=True)
    m=p.add_mutually_exclusive_group();m.add_argument('--freeze-only',action='store_true');m.add_argument('--verify-only',action='store_true')
    a=p.parse_args();out=a.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'):p.error('use research results')
    if a.freeze_only:freeze(out)
    else:
        with screen.prior.builder.extension.exclusive_run(out):execute(out,a.verify_only)
