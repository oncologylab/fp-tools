"""Pre-score cohort construction, matching and covariate-only diagnostics."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from sklearn.linear_model import LogisticRegression
from access import digest
from prepare import ROOT,dump
sys.path.insert(0,str(ROOT/'benchmarks/scripts'))
from caliper_matching import caliper_assignment


def smd(x,y):
    a,b=x[y==1],x[y==0]
    if min(len(a),len(b))<2:return np.full(x.shape[1],np.nan)
    variance=((len(a)-1)*a.var(0,ddof=1)+(len(b)-1)*b.var(0,ddof=1))/(len(a)+len(b)-2)
    difference=a.mean(0)-b.mean(0)
    return np.divide(difference,np.sqrt(variance),out=np.where(difference==0,0.,np.inf),where=variance>0)


def match(x,y,scale):
    p,n=np.flatnonzero(y==1),np.flatnonzero(y==0)
    if not len(p) or not len(n):return p[:0],n[:0],np.array([])
    d=cdist(x[p]/scale,x[n]/scale,'sqeuclidean');pi,ni=caliper_assignment(d,.5)
    distance=np.sqrt(d[pi,ni]);p,n=p[pi],n[ni]
    order=np.lexsort((np.arange(len(p)),distance))
    for k in range(len(p),1,-1):
        ix=order[:k]
        if np.max(abs(smd(x[np.r_[p[ix],n[ix]]],np.repeat([1,0],k))))<=.1:
            return p[ix],n[ix],distance[ix]
    return p[:0],n[:0],distance[:0]


def execute(out):
    target=out/'provenance/cohorts.freeze.json'
    if target.exists():raise ValueError('cohort attempt already exists')
    pool=pd.read_csv(out/'tables/pool.tsv.gz',sep='\t')
    tasks=pd.read_csv(out/'tables/task_ledger.tsv',sep='\t')
    rows=[];balance=[];baseline=[];folds=[]
    for task in tasks.itertuples():
        selected=pool[pool.cell.eq(task.cell)&pool.tf.eq(task.tf)&pool.supported].copy()
        x=np.column_stack([selected.motif_score,np.log1p(selected.central_201),np.log1p(selected.outer_101_500)])
        y=selected.chip_label.to_numpy(int)
        for chrom in [f'chr{i}' for i in range(1,16)]:
            train=np.flatnonzero(selected.TFBS_chr.ne(chrom));test=np.flatnonzero(selected.TFBS_chr.eq(chrom))
            scale=x[train].std(0) if len(train) else np.ones(3);scale=np.where(scale>0,scale,1.)
            pi,ni,dist=match(x[test],y[test],scale)
            tp,tn,td=match(x[train],y[train],scale)
            fitting=train[np.r_[tp,tn]];held=test[np.r_[pi,ni]]
            if task.cell=='HepG2' and task.tf=='HMGA1':held=held[:0];pi=pi[:0];ni=ni[:0];dist=dist[:0]
            ids=[f'{task.cell}:{task.tf}:{chrom}:{i}' for i in range(len(pi))]
            if len(held):
                frame=selected.iloc[held].copy();frame['pair_id']=ids+ids;frame['match_distance']=np.r_[dist,dist]
                frame['cohort']='nested_training';rows.append(frame)
                if len(fitting) and np.unique(y[fitting]).size==2:
                    mean=x[fitting].mean(0)
                    model=LogisticRegression(C=1.,solver='lbfgs',max_iter=1000,random_state=2026).fit((x[fitting]-mean)/scale,y[fitting])
                    prediction=model.decision_function((x[held]-mean)/scale)
                    baseline.extend(dict(site_id=sid,score=float(score),fold=chrom,fitting_sites=len(fitting)) for sid,score in zip(frame.site_id,prediction))
                else:
                    baseline.extend(dict(site_id=sid,score=np.nan,fold=chrom,fitting_sites=len(fitting)) for sid in frame.site_id)
            for phase,ix in [('before',test),('after',held),('fitting_matches',fitting)]:
                values=smd(x[ix],y[ix])
                balance.append(dict(cell=task.cell,tf=task.tf,fold=chrom,phase=phase,sites=len(ix),
                    positives=int(y[ix].sum()),motif_smd=values[0],central_smd=values[1],outer_smd=values[2]))
            folds.append(dict(cell=task.cell,tf=task.tf,fold=chrom,held_pairs=len(pi),fitting_pairs=len(tp),
                              status='AVAILABLE' if len(pi) else 'UNAVAILABLE',reason='' if len(pi) else 'no balanced held pairs',
                              scale_motif=scale[0],scale_central=scale[1],scale_outer=scale[2]))
    nested=pd.concat(rows,ignore_index=True)
    if nested.site_id.duplicated().any():raise ValueError('held site appears in multiple folds')
    # Whole-chromosome exclusion means every overlapping/shared locus is
    # excluded from all fitting tasks when that chromosome is evaluated.
    nested.to_csv(out/'tables/nested_sites.tsv.gz',sep='\t',index=False)
    pd.DataFrame(balance).to_csv(out/'tables/cohort_balance.tsv',sep='\t',index=False)
    pd.DataFrame(baseline).to_csv(out/'tables/covariate_baseline.tsv.gz',sep='\t',index=False)
    fd=pd.DataFrame(folds);fd.to_csv(out/'tables/cohort_folds.tsv',sep='\t',index=False)
    counts=nested.groupby(['cell','tf']).size().div(2)
    tasks['nested_pairs']=[int(counts.get((r.cell,r.tf),0)) for r in tasks.itertuples()]
    tasks['nested_gate_eligible_before_scoring']=tasks.nested_pairs>=200
    tasks['cohort_status']=np.where(tasks.nested_pairs>0,'AVAILABLE','UNAVAILABLE')
    tasks['missingness_reason']=np.where(tasks.nested_pairs>0,'','no balanced held chromosome pairs; retained in denominator')
    tasks.to_csv(out/'tables/task_ledger.tsv',sep='\t',index=False)
    dump(target,dict(code_sha256=digest(__file__),pool_sha256=digest(out/'tables/pool.tsv.gz'),
                    nested_sha256=digest(out/'tables/nested_sites.tsv.gz'),task_ledger_sha256=digest(out/'tables/task_ledger.tsv'),
                    nested_sites=len(nested),original_pool_sites=len(pool),task_denominator=21,
                    scoring_started=False,folds='held chromosome excluded across all tasks/contexts',
                    conditioning='held labels used only by frozen matched evaluation sampling, never fitting scaling or nuisance',
                    fragment_groups='overlapping windows recorded; shared-fragment membership awaits complete-pair extraction'))
    print('Frozen held-chromosome cohort:',len(nested),'sites; all21 task ledger retained',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);a=p.parse_args();execute(a.outdir.resolve())
