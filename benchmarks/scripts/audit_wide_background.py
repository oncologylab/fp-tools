#!/usr/bin/env python3
"""Verify archived flank summaries and compare radii on common sequence support."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
import build_wide_training_profiles as source
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import checked_path, canonical_id

ROOT=source.ROOT


def prefix_diagnostic(counts,log_bias,flank,outer):
    counts,log_bias=np.asarray(counts,float),np.asarray(log_bias,float)
    if counts.ndim!=2 or counts.shape!=log_bias.shape or counts.shape[1]!=2*flank+1 or not 50<outer<=flank:
        raise ValueError('invalid diagnostic inputs')
    if not np.isfinite(counts).all() or (counts<0).any() or np.isinf(log_bias).any():
        raise ValueError('invalid counts or exposure')
    def prefix(values): return np.pad(np.cumsum(values,axis=1),((0,0),(1,0)))
    def interval(values,start,stop): return values[:,stop]-values[:,start]
    left,right=(flank-outer,flank-49),(flank+50,flank+outer+1)
    finite=np.isfinite(log_bias);n=prefix(finite.astype(int))
    available=(interval(n,*left)+interval(n,*right))==2*(outer-49)
    c=prefix(counts);lc,rc=interval(c,*left),interval(c,*right)
    discrepancy=np.full(len(counts),np.nan)
    # Sum the two exposures independently of the production helper. Work in
    # log space so arbitrary offsets and extremely small exposures stay safe.
    from scipy.special import logsumexp
    values=log_bias[available]
    lexp=logsumexp(values[:,left[0]:left[1]],axis=1)
    rexp=logsumexp(values[:,right[0]:right[1]],axis=1)
    discrepancy[available]=np.abs((np.log(lc[available]+.5)-np.log(rc[available]+.5)+rexp-lexp)/np.log(2))
    return dict(available=available,left_count=lc,right_count=rc,disagreement=discrepancy)


def summary(diagnostic,selected):
    left=diagnostic['left_count'][selected];right=diagnostic['right_count'][selected]
    values=diagnostic['disagreement'][selected]
    return dict(supported=int(selected.sum()),zero_background=int(((left+right)==0).sum()),
        zero_left=int((left==0).sum()),zero_right=int((right==0).sum()),
        median_background_count=float(np.median(left+right)) if len(left) else np.nan,
        median_absolute_log2_rate_disagreement=float(np.median(values)) if len(left) else np.nan,
        p90_absolute_log2_rate_disagreement=float(np.quantile(values,.9)) if len(left) else np.nan)


def freeze(directory,out):
    if out.exists(): raise ValueError('fresh directory required')
    configuration=directory/'configuration.freeze.json'
    config=source.old.parent.extension.read_document(configuration,'configuration_id')
    source.old.parent.base.verify_stage(directory/'profiles_manifest.json',config['configuration_id'],ROOT)
    inputs=[record(configuration),record(directory/'profiles_manifest.json')]
    stages={}
    for cell in ('HepG2','K562'):
        path=directory/(cell+'.json')
        stage=source.old.parent.base.verify_stage(path,config['configuration_id'],ROOT)
        inputs += [record(path)]+list(stage['outputs'].values());stages[cell]=stage
    value=dict(inputs=inputs,stages=stages,code=config['code']+[record(Path(__file__))],
        runtime=source.old.parent.extension.runtime(),flank=config['protocol']['extraction_flank'],
        radii=config['protocol']['background_outer_grid'],labels_read=False,promotion_eligible=False,
        purpose='independent diagnostic verification and common support comparison; no window selection')
    value['configuration_id']=canonical_id(value)
    out.mkdir(parents=True);write_immutable(out/'configuration.freeze.json',value)


def execute(out,verify=False):
    config=source.old.parent.extension.read_document(out/'configuration.freeze.json','configuration_id')
    if config['runtime']!=source.old.parent.extension.runtime(): raise ValueError('runtime changed')
    for entry in config['inputs']+config['code']: checked_path(entry,ROOT)
    manifest=out/'audit_manifest.json'
    if manifest.exists():
        source.old.parent.base.verify_stage(manifest,config['configuration_id'],ROOT)
        print('Verified completed wider-background audit',flush=True);return
    if verify or (out/'common_support.tsv').exists(): raise ValueError('missing or partial output')
    rows=[];verified=0
    for cell,stage in config['stages'].items():
        sites=pd.read_csv(stage['outputs']['sites']['path'],sep='\t')
        reported=pd.read_csv(stage['outputs']['diagnostics']['path'],sep='\t',float_precision='round_trip')
        with np.load(stage['outputs']['profiles']['path'],allow_pickle=False) as archive:
            np.testing.assert_array_equal(source.old.strands.site_hashes(sites),archive['site_hash'])
            valid=archive['bam_valid']
            for strand in ('plus','minus'):
                counts,bias=archive[strand+'_observed'],archive[strand+'_log_bias']
                diagnostics={radius:prefix_diagnostic(counts,bias,config['flank'],radius) for radius in config['radii']}
                common=valid & np.logical_and.reduce([d['available'] for d in diagnostics.values()])
                for radius,diagnostic in diagnostics.items():
                    for tf in sorted(sites.tf.unique()):
                        task=sites.tf.eq(tf).to_numpy();selected=task & valid & diagnostic['available']
                        previous=reported.loc[reported.tf.eq(tf)&reported.strand.eq(strand)&reported.outer.eq(radius)]
                        if len(previous)!=1: raise ValueError('missing or duplicate task row')
                        previous=previous.iloc[0]
                        if previous.sites!=task.sum(): raise ValueError('task denominator changed')
                        calculated=summary(diagnostic,selected)
                        for key,value in calculated.items():
                            np.testing.assert_allclose(value,previous[key],rtol=1e-9,atol=1e-10,equal_nan=True)
                        verified+=1
                        rows.append(dict(cell=cell,tf=tf,strand=strand,outer=radius,sites=int(task.sum()),
                            **summary(diagnostic,task&common)))
    if verified!=126: raise ValueError('full 21-task strand/radius denominator changed')
    path=out/'common_support.tsv';pd.DataFrame(rows).to_csv(path,sep='\t',index=False)
    result=dict(configuration_id=config['configuration_id'],outputs=dict(common_support=record(path)),
        verified_summary_rows=verified,labels_read=False,promotion_eligible=False,window_selected=False)
    result['stage_id']=canonical_id(result);write_immutable(manifest,result)
    print('Independently verified 126 rows; retained every task on common sequence support',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path)
    parser.add_argument('--outdir',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--freeze-only',action='store_true');mode.add_argument('--verify-only',action='store_true')
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'): parser.error('use research results')
    if args.freeze_only:
        if args.source is None: parser.error('--source required to freeze')
        freeze(args.source.resolve(),out)
    else:
        with source.old.parent.extension.exclusive_run(out): execute(out,args.verify_only)
