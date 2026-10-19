"""Read-only independent count/tail/rank checks plus declared aggregation replay.

Rank and AP calculations below do not call the scientific scorer or sklearn.
The aggregation replay deliberately reuses the frozen weighting contract; it is
not an independent validation of cohort construction or full-refit uncertainty.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import binom, rankdata
from chromatin_metrics import screen_pairs, arm_means


def verify(root, output=None):
    sub=root/'chromatin'
    read=lambda name:pd.read_csv(sub/'tables'/name,sep='\t',float_precision='round_trip',low_memory=False)
    scores=read('per_site_scores.tsv.gz');sites=read('site_ledger.tsv.gz')
    assert len(sites)==6020 and sites.site_id.is_unique
    assert sites.TFBS_chr.isin([f'chr{i}' for i in range(1,16)]).all()
    prediction_rows=0
    with np.load(sub/'arrays/sequence_contexts.npz',allow_pickle=False) as contexts:
        for arm,l2 in [('nuclear_primary',.01),('nuclear_ablation',.0001)]:
            with np.load(sub/f'arrays/{arm}.logits.npz',allow_pickle=False) as z:
                np.testing.assert_array_equal(z['site_id'],contexts['site_id'])
                for chromosome in [f'chr{i}' for i in range(1,16)]:
                    ix=np.flatnonzero(contexts['chromosome']==chromosome)
                    if not len(ix):continue
                    encoded=np.array([[{'A':0,'C':1,'G':2,'T':3}[b] for b in s] for s in contexts['sequence'][ix]])
                    windows=np.lib.stride_tricks.sliding_window_view(encoded,21,axis=1)[:,31:232,:]
                    with np.load(root/f'models/{chromosome}.l2_{l2}.ensemble.npz',allow_pickle=False) as model:
                        predictions=[]
                        for w in [windows,3-windows[:,:,::-1]]:
                            v=sum(model['main'][j,w[:,:,j]] for j in range(21))
                            for d in [1,5,9]:
                                v+=sum(model[f'pair_{d}'][j,w[:,:,j],w[:,:,j+d]] for j in range(21-d))
                            predictions.append(v)
                    p,m=predictions;reverse=contexts['motif_strand'][ix]=='-'
                    a,b=p.copy(),m.copy();a[reverse]=m[reverse,::-1];b[reverse]=p[reverse,::-1]
                    np.testing.assert_allclose(np.stack([a,b],axis=1),z['logits'][ix],atol=1e-12,rtol=1e-12)
                    prediction_rows+=len(ix)
    offset=np.arange(-100,101);center=abs(offset)<=7
    domain=center|((abs(offset)>=12)&(abs(offset)<=40))
    checked=0
    for (library,depth),group in scores.groupby(['library','depth']):
        with np.load(sub/f'arrays/{library}.{depth}.npz',allow_pickle=False) as z:
            ids=z['site_id'];observed=z['observed']
        for arm,tab in group.groupby('arm'):
            tab=tab.set_index('site_id').loc[ids]
            for si,strand in enumerate(['plus','minus']):
                np.testing.assert_array_equal(observed[:,si,center].sum(1),tab[strand+'_center'])
                np.testing.assert_array_equal(observed[:,si,domain].sum(1),tab[strand+'_total'])
            if not arm.startswith('nuclear_'):continue
            with np.load(sub/f'arrays/{arm}.logits.npz',allow_pickle=False) as z:
                ix=pd.Index(z['site_id']).get_indexer(ids);assert (ix>=0).all();log=z['logits'][ix]
            q=np.exp(logsumexp(log[:,:,center],axis=2)-logsumexp(log[:,:,domain],axis=2))
            np.testing.assert_allclose(q,tab[['q_plus','q_minus']],rtol=1e-12,atol=1e-14)
            for row in tab.itertuples():
                n,m,c=int(row.plus_total),int(row.minus_total),int(row.plus_center+row.minus_center)
                if not np.isfinite([row.q_plus,row.q_minus]).all():
                    assert row.status=='UNAVAILABLE';continue
                # Independent truncated convolution in log space; neither
                # historical depletion_score nor its implementation is imported.
                x=np.arange(min(n,c)+1)
                lp=logsumexp(binom.logpmf(x,n,row.q_plus)+binom.logcdf(c-x,m,row.q_minus))
                value=round(min(300.,max(0.,-lp/np.log(10))),10)
                ceiling=round(min(300.,max(0.,-(n*np.log1p(-row.q_plus)+m*np.log1p(-row.q_minus))/np.log(10))),10)
                np.testing.assert_allclose([value,ceiling],[row.score,row.ceiling],rtol=0,atol=2e-10)
                status=('NO_INFORMATION' if n+m==0 else 'ABSTAIN_UNREACHABLE' if ceiling<=1.3010299957 else 'CALL' if value>1.3010299957 else 'NO_CALL')
                assert status==row.status and (status=='CALL')==row.call
                checked+=1
    folds=read('fold_metrics.tsv');cov=read('covariates.tsv.gz')
    index={k:g for k,g in scores.groupby(['library','depth','tf','arm','TFBS_chr'])}
    numerical=0
    for r in folds.itertuples():
        if r.status!='AVAILABLE':continue
        key=(r.library,r.depth,r.tf,'LOG21_count_tail' if r.arm=='covariates_only' else r.arm,r.fold)
        g=index[key]
        if r.arm=='covariates_only':
            v=cov[cov.library.eq(r.library)&cov.depth.eq(r.depth)].set_index('site_id').loc[g.site_id,'score'].to_numpy()
        else:v=g.score.to_numpy()
        y=g.chip_label.to_numpy(int);p=y.sum();n=len(y)-p
        auc=(rankdata(v)[y==1].sum()-p*(p+1)/2)/(p*n)
        order=np.argsort(-v,kind='stable');ys=y[order];vs=v[order]
        ends=np.r_[np.flatnonzero(vs[:-1]!=vs[1:]),len(vs)-1]
        tp=np.cumsum(ys)[ends];ap=np.sum(np.diff(np.r_[0,tp])/p*tp/(ends+1))
        np.testing.assert_allclose([auc,ap],[r.auroc,r.auprc],rtol=1e-12,atol=1e-12)
        numerical+=1
    tasks=read('task_ledger.tsv');summaries=[];deltas=[]
    for row in read('screen_summary.tsv').itertuples():
        d,s=screen_pairs(folds,tasks,row.depth,row.method,row.reference)
        for field in ['mean_auroc_gain','mean_relative_auprc_gain','bootstrap_low','bootstrap_high','valid_draws']:
            np.testing.assert_allclose(s[field],getattr(row,field),atol=1e-12,rtol=1e-12)
        summaries.append(s);deltas.extend(d)
    if output:
        output.mkdir(parents=True,exist_ok=False)
        pd.DataFrame(summaries).to_csv(output/'screen_summary.tsv',sep='\t',index=False)
        pd.DataFrame(deltas).to_csv(output/'paired_deltas.tsv',sep='\t',index=False)
        arm_means(folds,tasks).to_csv(output/'task_metrics.tsv',sep='\t',index=False)
    return dict(independent_sequence_prediction_sites=prediction_rows,independent_new_tail_and_call_rows=checked,independent_rank_AP_rows=numerical,
        replayed_paired_screens=len(summaries),cohort_sites=len(sites),status='VERIFIED',
        limitation='aggregation uses frozen contract; no independent raw extraction or full-refit uncertainty claim')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--output',type=Path)
    a=p.parse_args();print(json.dumps(verify(a.root,a.output),indent=2))
