"""Frozen fitting-fold covariate baseline and complete stratified diagnostics."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from scipy.stats import spearmanr
from cohorts import match
from verify import metrics,read_table
from access import digest
from prepare import dump


def execute(out):
    dump(out/'provenance/confounds.freeze.json',dict(code_sha256=digest(__file__),
        baseline='logistic C1 lbfgs max_iter1000; motif score, log1p whole201 cuts, log1p library complete depth; no central shape',
        fitting='whole held chromosome excluded; frozen anchor-covariate matched training sample; fit-only centering and scaling',
        strata='frozen regional count and expected-count bins plus fitting-chromosome motif/accessibility quartiles',
        role='diagnostic, never an additional candidate or rescue classifier'))
    pool=read_table(out/'tables/pool.tsv.gz');held=read_table(out/'tables/nested_sites.tsv.gz')
    all_scores=read_table(out/'tables/per_site_scores.tsv.gz');all_scores=all_scores[all_scores.cohort.eq('nested')]
    predictions=[];rows=[];balance=[]
    for path in sorted((out/'tables').glob('counts.*.tsv.gz')):
        counts=read_table(path);library=counts.library.iloc[0];depth=counts.depth.iloc[0];cell=library.split('_')[0]
        sites=pool[pool.cell.eq(cell)].reset_index(drop=True)
        array_path=out/(f'arrays/{cell}.replay_pool.npz' if depth=='archived_10M' else f'arrays/{library}.{depth}.npz')
        covariate_path=out/'tables/pool_library_covariates.tsv.gz'
        if covariate_path.exists():
            covariates=read_table(covariate_path)
            covariates=covariates[covariates.library.eq(library)&covariates.depth.eq(depth)].set_index('site_id').loc[sites.site_id]
            whole=covariates.actual_whole201.to_numpy()
            if covariates.depth_covariate_constant.nunique()!=1:raise ValueError('ambiguous library depth covariate')
            fragment_depth=float(covariates.depth_covariate_constant.iloc[0])
        else:
            with np.load(array_path,allow_pickle=False) as z:
                np.testing.assert_array_equal(z['site_id'],sites.site_id.to_numpy(str))
                whole=z['plus_observed'].sum(1)+z['minus_observed'].sum(1)
                fragment_depth=float(z['achieved_fragments']) if 'achieved_fragments' in z else 10000000.
        sites['actual_whole201']=whole
        for tf,task in sites[sites.supported].groupby('tf'):
            anchor=np.column_stack([task.motif_score,np.log1p(task.central_201),np.log1p(task.outer_101_500)])
            features=np.column_stack([task.motif_score,np.log1p(task.actual_whole201),np.full(len(task),np.log1p(fragment_depth))])
            labels=task.chip_label.to_numpy(int)
            for chrom in [f'chr{i}' for i in range(1,16)]:
                train=np.flatnonzero(task.TFBS_chr.ne(chrom));test=np.flatnonzero(task.TFBS_chr.eq(chrom)&task.site_id.isin(held.site_id))
                scale=anchor[train].std(0);scale=np.where(scale>0,scale,1.)
                a,b,d=match(anchor[train],labels[train],scale);fit=train[np.r_[a,b]]
                if not len(test) or not len(fit):continue
                mean=features[fit].mean(0);sd=features[fit].std(0);sd=np.where(sd>0,sd,1.)
                model=LogisticRegression(C=1.,solver='lbfgs',max_iter=1000,random_state=2026).fit((features[fit]-mean)/sd,labels[fit])
                pred=model.decision_function((features[test]-mean)/sd)
                for i,value in zip(test,pred):predictions.append(dict(site_id=task.iloc[i].site_id,library=library,depth=depth,cell=cell,tf=tf,fold=chrom,score=value,chip_label=int(labels[i]),fitting_sites=len(fit)))
                observed=all_scores[all_scores.library.eq(library)&all_scores.depth.eq(depth)&all_scores.tf.eq(tf)&all_scores.fold.eq(chrom)].copy()
                current=task.iloc[test][['site_id','motif_score','actual_whole201']].merge(counts[['site_id','count_stratum','expected_stratum']],on='site_id',validate='one_to_one')
                for name,column in [('motif',0),('accessibility',1)]:
                    cuts=np.quantile(features[fit,column],[.25,.5,.75])
                    current[name+'_quartile']=np.searchsorted(cuts,features[test,column],side='right')+1
                for field,levels in [('all',['all']),('count_stratum',['0','1','2-4','5-9','10-19','20-49','50-99','100+']),
                                     ('expected_stratum',['<0.1','0.1-<0.5','0.5-<1','1-<3','3+']),('motif_quartile',[1,2,3,4]),('accessibility_quartile',[1,2,3,4])]:
                    for level in levels:
                        ids=current.site_id if field=='all' else current.loc[current[field].eq(level),'site_id']
                        for method,group in observed.groupby('method'):
                            g=group[group.site_id.isin(ids)&np.isfinite(group.score)]
                            auc,ap=metrics(g.chip_label,g.score)
                            rows.append(dict(library=library,depth=depth,cell=cell,tf=tf,fold=chrom,method=method,
                                stratum_type=field,stratum=level,original_sites=len(ids),finite_sites=len(g),positives=int(g.chip_label.sum()),
                                auroc=auc,auprc=ap,status='AVAILABLE' if np.isfinite(auc) else 'UNDERPOWERED'))
        print('Confounding diagnostics complete',library,depth,flush=True)
    pred=pd.DataFrame(predictions);pred.to_csv(out/'tables/covariate_baseline_all_libraries.tsv.gz',sep='\t',index=False)
    summary=[]
    for keys,group in pred.groupby(['library','depth','cell','tf','fold']):
        auc,ap=metrics(group.chip_label,group.score);summary.append(dict(zip(['library','depth','cell','tf','fold'],keys),auroc=auc,auprc=ap,sites=len(group)))
    pd.DataFrame(summary).to_csv(out/'tables/covariate_baseline_metrics.tsv',sep='\t',index=False)
    pd.DataFrame(rows).to_csv(out/'tables/stratified_metrics.tsv.gz',sep='\t',index=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);a=p.parse_args();execute(a.outdir.resolve())
