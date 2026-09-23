"""Portable fit-free reanalysis of frozen effects, calls and metrics, without BAMs.

Write only to a new destination. This is the run entry point; verify.py never
calls it or silently regenerates missing evidence.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from shrinkage import posterior
from verify import read_table,metrics


def calculate(args):
    key,method,c,n,p,rho=args
    if method.endswith('count_tail'):
        from conditional_protection_score import depletion_score
        score=round(depletion_score(*n,sum(c),*p),10)
    else:score=posterior(*c,*n,*p,rho)['mean']
    return key,score


def execute(root,out,workers):
    if out.exists():raise ValueError('new output directory required')
    out.mkdir(parents=True)
    original=read_table(root/'tables/per_site_scores.tsv.gz')
    unique=original.drop_duplicates(['site_id','library','depth','method'])
    lookup={};jobs=[]
    for path in sorted((root/'tables').glob('counts.*.tsv.gz')):
        frame=read_table(path)
        for r in frame.itertuples():lookup[(r.site_id,r.library,r.depth)]=r
    result={}
    for row in unique.itertuples():
        key=(row.site_id,row.library,row.depth,row.method)
        if not np.isfinite(row.score):result[key]=np.nan;continue
        if row.method=='archived_DWM_residual':continue
        counts=lookup[key[:3]];arm=row.method.split('_')[0]
        jobs.append((key,row.method,[int(counts.plus_center),int(counts.minus_center)],
                     [int(counts.plus_total),int(counts.minus_total)],
                     [getattr(counts,arm+'_plus_p'),getattr(counts,arm+'_minus_p')],float(row.rho)))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for key,value in pool.map(calculate,jobs,chunksize=64):result[key]=value
    replay=read_table(root/'tables/replay_sites.tsv.gz');positions=np.arange(-100,101)
    weights=((abs(positions)>=12)&(abs(positions)<=40))/58-(abs(positions)<=7)/15
    with np.load(root/'arrays/archived_replay.npz',allow_pickle=False) as z:
        values=z['DWM_residual']@weights;values[z['raw'].sum(1)==0]=0.
    for row,value in zip(replay.itertuples(),values):result[(row.site_id,row.cell+'_rep1','archived_10M','archived_DWM_residual')]=value
    original['recomputed_score']=[result[(r.site_id,r.library,r.depth,r.method)] for r in original.itertuples()]
    original['recomputed_call']=original.call_eligible.eq(True)&~original.prior_only.eq(True)&(original.recomputed_score>original.threshold)
    original.to_csv(out/'recomputed_scores.tsv.gz',sep='\t',index=False)
    difference=np.abs(original.recomputed_score-original.score)
    print('Maximum frozen-score discrepancy:',difference.max())
    disagreements=int(original.recomputed_call.ne(original.call.fillna(False)).sum())
    print('Call disagreements:',disagreements)
    if difference.max()>1e-8 or disagreements or not np.array_equal(np.isfinite(original.recomputed_score),np.isfinite(original.score)):
        raise AssertionError('portable replay differs from frozen scores/calls; outputs retained for inspection')
    rows=[]
    reference=read_table(root/'tables/fold_metrics.tsv')
    columns=['cohort','library','depth','cell','tf','method','fold']
    expected_rows={tuple(getattr(r,c) for c in columns):r.original_sites for r in reference.itertuples()}
    if len(expected_rows)!=len(reference):raise AssertionError('duplicate reference fold keys')
    for keys,frame in original.groupby(columns):
        finite=frame[np.isfinite(frame.recomputed_score)]
        complete=keys in expected_rows and len(finite)==int(expected_rows[keys])
        auc,ap=metrics(finite.chip_label,finite.recomputed_score) if complete else (np.nan,np.nan)
        rows.append(dict(zip(columns,keys),auroc=auc,auprc=ap,finite_sites=len(finite),complete=complete))
    pd.DataFrame(rows).to_csv(out/'recomputed_fold_metrics.tsv',sep='\t',index=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--outdir',type=Path,required=True);p.add_argument('--workers',type=int,default=4)
    a=p.parse_args();execute(a.root.resolve(),a.outdir.resolve(),a.workers)
