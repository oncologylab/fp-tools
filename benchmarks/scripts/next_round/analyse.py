"""Read-only scientific summaries and plots from completed per-site outputs."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.special import logsumexp
from prepare import ROOT,dump
from access import digest
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'benchmarks/scripts'))
from fp_tools.tools.functional_footprints import standardized_functional_separation
from fp_tools.tools.parametric_bias import calibrated_residuals
from fixed_shape_normalization import normalize_fixed


def read_table(path):
    return pd.read_csv(path,sep='\t',float_precision='round_trip',low_memory=False)

X=np.arange(-100,101)
WEIGHTS=((abs(X)>=12)&(abs(X)<=40))/58-(abs(X)<=7)/15


def wilson(k,n):
    if not n:return np.nan,np.nan
    z=norm.ppf(.975);p=k/n;d=1+z*z/n
    mid=(p+z*z/(2*n))/d;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return mid-half,mid+half


def metrics(out):
    task=read_table(out/'tables/task_ledger.tsv')
    cohorts={k:read_table(out/f'tables/{v}_sites.tsv.gz') for k,v in [('replay','replay'),('nested','nested')]}
    inputs=sorted((out/'tables').glob('scores.*.tsv.gz'))
    rows=[];all_scores=[];run_ledger=[]
    for path in inputs:
        scores=read_table(path);library=scores.library.iloc[0];depth=scores.depth.iloc[0];cell=library.split('_')[0]
        for cohort,sites in cohorts.items():
            current=scores[scores.site_id.isin(sites.site_id)].copy()
            if cohort=='replay' and depth=='archived_10M':
                subset=sites[sites.cell.eq(cell)].copy();indices=np.flatnonzero(sites.cell.eq(cell))
                with np.load(out/'arrays/archived_replay.npz',allow_pickle=False) as z:
                    values=z['DWM_residual'][indices]@WEIGHTS;values[z['raw'][indices].sum(1)==0]=0.
                archived=subset[['site_id','cell','tf','TFBS_chr','chip_label']].copy();archived['score']=values
                archived['method']='archived_DWM_residual';archived['status']='AVAILABLE';archived['library']=library;archived['depth']=depth
                archived['fold']=archived.TFBS_chr
                archived['call']=False;archived['call_eligible']=False;archived['exclusion_reason']='ranking only; archived count-tail policy is a separate arm'
                current=pd.concat([current,archived],ignore_index=True)
            current['cohort']=cohort;all_scores.append(current)
            for t in task[task.cell.eq(cell)].itertuples():
                expected=sites[sites.cell.eq(cell)&sites.tf.eq(t.tf)]
                for method in ['archived_DWM_residual','uniform_count_tail','LOG21_count_tail','LOG21_binomial','LOG21_beta_binomial','DWM_binomial','DWM_beta_binomial']:
                    group=current[current.tf.eq(t.tf)&current.method.eq(method)]
                    for chrom in [f'chr{i}' for i in range(1,16)]:
                        original=expected[expected.TFBS_chr.eq(chrom)]
                        subset=group[group.TFBS_chr.eq(chrom)]
                        finite=subset[np.isfinite(subset.score)]
                        available=len(finite)==len(original) and len(original)>0 and finite.chip_label.nunique()==2
                        row=dict(cohort=cohort,library=library,depth=depth,cell=cell,tf=t.tf,method=method,fold=chrom,
                            original_sites=len(original),finite_sites=len(finite),positives=int(finite.chip_label.sum()),
                            prevalence=float(finite.chip_label.mean()) if len(finite) else np.nan,
                            auroc=roc_auc_score(finite.chip_label,finite.score) if available else np.nan,
                            auprc=average_precision_score(finite.chip_label,finite.score) if available else np.nan,
                            status='AVAILABLE' if available else 'INCOMPLETE',
                            reason='' if available else ('no balanced held pairs' if not len(original) else 'missing arm/support or unidentified model'))
                        rows.append(row)
                    run_ledger.append(dict(cohort=cohort,library=library,depth=depth,cell=cell,tf=t.tf,method=method,
                        original_sites=len(expected),finite_sites=int(np.isfinite(group.score).sum()),
                        status='AVAILABLE' if len(expected)>0 and len(group)==len(expected) and np.isfinite(group.score).all() else 'INCOMPLETE'))
    if not rows:raise ValueError('no completed score outputs')
    folds=pd.DataFrame(rows);folds.to_csv(out/'tables/fold_metrics.tsv',sep='\t',index=False)
    pd.concat(all_scores,ignore_index=True).to_csv(out/'tables/per_site_scores.tsv.gz',sep='\t',index=False)
    pd.DataFrame(run_ledger).to_csv(out/'tables/run_ledger.tsv',sep='\t',index=False)
    summary=folds.groupby(['cohort','library','depth','cell','tf','method'],sort=True)[['auroc','auprc']].mean().reset_index()
    summary.to_csv(out/'tables/task_metrics.tsv',sep='\t',index=False)
    deltas=[];bootstrap=[]
    weights=np.random.default_rng(2026).multinomial(15,np.repeat(1/15,15),size=2000)
    for (cohort,library,depth),group in folds.groupby(['cohort','library','depth']):
        for reference in ['archived_DWM_residual','LOG21_count_tail','LOG21_binomial','DWM_binomial']:
            base=group[group.method.eq(reference)]
            for method in group.method.unique():
                if method==reference:continue
                merged=group[group.method.eq(method)].merge(base,on=['cell','tf','fold'],suffixes=('','_ref'))
                merged['auroc_delta']=merged.auroc-merged.auroc_ref
                merged['auprc_delta']=merged.auprc-merged.auprc_ref
                for (cell,tf),task_rows in merged.groupby(['cell','tf']):
                    common=task_rows[np.isfinite(task_rows.auroc_delta)&np.isfinite(task_rows.auprc_delta)]
                    base_ap=common.auprc_ref.mean()
                    deltas.append(dict(cohort=cohort,library=library,depth=depth,cell=cell,tf=tf,method=method,reference=reference,
                        nonempty_blocks=len(common),auroc_gain=common.auroc_delta.mean(),auprc_gain=common.auprc_delta.mean(),
                        relative_auprc_gain=common.auprc_delta.mean()/base_ap if base_ap>0 else np.nan,
                        scope='paired complete-fold common support; see full coverage ledger'))
                # Equal chromosome then equal task, never pooled scores/AP.
                estimates=[];blocks=[]
                for (cell,tf),tr in merged.groupby(['cell','tf']):
                    a=tr.set_index('fold').auroc_delta.reindex([f'chr{i}' for i in range(1,16)]).to_numpy()
                    mask=np.isfinite(a)
                    denom=weights@mask.astype(float)
                    values=np.divide(weights@np.nan_to_num(a),denom,out=np.full(2000,np.nan),where=denom>0)
                    if mask.any():estimates.append(values);blocks.append(int(mask.sum()))
                if estimates:
                    sample=np.mean(np.stack(estimates),axis=0);valid=np.isfinite(sample)
                    lower,upper=np.quantile(sample[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan]
                else:valid=np.zeros(2000,bool);lower=upper=np.nan
                bootstrap.append(dict(cohort=cohort,library=library,depth=depth,method=method,reference=reference,
                    defined_tasks=len(estimates),original_task_denominator=21,valid_draws=int(valid.sum()),
                    low=lower,high=upper,minimum_nonempty_blocks=min(blocks) if blocks else 0,
                    principal_qualified=False,reason='per-library fixed-prediction diagnostic; missing tasks/complete full-depth screen not waived'))
    pd.DataFrame(deltas).to_csv(out/'tables/paired_deltas.tsv',sep='\t',index=False)
    pd.DataFrame(bootstrap).to_csv(out/'tables/paired_bootstrap.tsv',sep='\t',index=False)


def calibration(out):
    rows=[]
    for path in sorted((out/'tables').glob('control_scores.*.tsv.gz')):
        frame=read_table(path)
        for (depth,cell,tf,method),group in frame.groupby(['depth','cell','tf','method']):
            for stratum_column in ['all','count_stratum','expected_stratum']:
                values=['all'] if stratum_column=='all' else (['0','1','2-4','5-9','10-19','20-49','50-99','100+'] if stratum_column=='count_stratum' else ['<0.1','0.1-<0.5','0.5-<1','1-<3','3+'])
                for value in values:
                    subset=group if value=='all' else group[group[stratum_column].eq(value)]
                    n=len(subset);informative=int((subset.total>0).sum());eligible=int(subset.call_eligible.sum());k=int(subset.call.sum())
                    lo,hi=wilson(k,n);ilo,ihi=wilson(k,eligible)
                    rows.append(dict(depth=depth,cell=cell,tf=tf,method=method,stratum_type=stratum_column,stratum=value,
                        supported_control_sites=n,informative_sites=informative,call_eligible_sites=eligible,calls=k,
                        rate=k/n if n else np.nan,wilson_low=lo,wilson_high=hi,
                        informative_rate=k/eligible if eligible else np.nan,informative_wilson_high=ihi,
                        interpretation='examined development data; Wilson is descriptive under overlapping sites'))
    pd.DataFrame(rows).to_csv(out/'tables/calibration.tsv',sep='\t',index=False)


def shape(out):
    pool=read_table(out/'tables/pool.tsv.gz')
    nested=read_table(out/'tables/nested_sites.tsv.gz')
    tasks=read_table(out/'tables/task_ledger.tsv')
    metrics_rows=[];identity=[]
    handle=gzip.open(out/'tables/shape_curves.tsv.gz','wt')
    writer=csv.DictWriter(handle,delimiter='\t',fieldnames=['library','depth','cell','tf','representation','label','position','value','aggregation','strand','units','stratum_type','stratum','sites','weight_sum'])
    writer.writeheader()
    class CurveSink:
        def append(self,row):writer.writerow(row)
    curve_rows=CurveSink()
    sources=sorted((out/'tables').glob('counts.*.tsv.gz'))
    with PdfPages(out/'figures/all_task_shapes.pdf') as pdf:
        for source in sources:
            frame=read_table(source);library=frame.library.iloc[0];depth=frame.depth.iloc[0];cell=library.split('_')[0]
            allsites=pool[pool.cell.eq(cell)].reset_index(drop=True);selected=allsites.site_id.isin(nested.site_id).to_numpy()
            sites=allsites[selected].reset_index(drop=True)
            path=out/(f'arrays/{cell}.replay_pool.npz' if depth=='archived_10M' else f'arrays/{library}.{depth}.npz')
            with np.load(path,allow_pickle=False) as z:observed=np.stack([z[s+'_observed'][selected] for s in ['plus','minus']],axis=1)
            with np.load(out/f'arrays/{library}.propensities.npz',allow_pickle=False) as z:
                logits=np.stack([z['LOG21_'+s][selected] for s in ['plus','minus']],axis=1)
            expected=np.exp(logits-logsumexp(logits,axis=2,keepdims=True))*observed.sum(2,keepdims=True)
            residual=calibrated_residuals(observed.sum(1),expected.sum(1),'deviance',dispersion=0. if cell=='HepG2' else 1.6144318781609228)
            panels={'raw':observed.sum(1),'LOG21_expected':expected.sum(1),'LOG21_residual':residual}
            celltasks=tasks[tasks.cell.eq(cell)]
            fig,axes=plt.subplots(len(celltasks),3,figsize=(13,2*len(celltasks)),squeeze=False)
            for row,t in enumerate(celltasks.itertuples()):
                mask=sites.tf.eq(t.tf).to_numpy();labels=sites.loc[mask,'chip_label'].to_numpy(int)
                for col,(arm,values) in enumerate(panels.items()):
                    ax=axes[row,col];ax.axvspan(-7,7,color='0.9');ax.set_title(f'{cell} {t.tf}: {arm}',fontsize=9)
                    if mask.sum()<4 or len(np.unique(labels))<2:
                        metrics_rows.append(dict(library=library,depth=depth,cell=cell,tf=t.tf,representation=arm,cohort='nested',
                            sites=int(mask.sum()),positive_sites=int(labels.sum()),status='UNAVAILABLE',reason='no balanced held pairs'))
                        ax.text(.5,.5,'UNAVAILABLE',transform=ax.transAxes,ha='center');continue
                    data=values[mask];normal=normalize_fixed(data,X,1.)
                    delta=normal[labels==1].mean(0)-normal[labels==0].mean(0)
                    center=delta[abs(X)<=5].mean();left=delta[(X>=-40)&(X<=-15)].mean();right=delta[(X>=15)&(X<=40)].mean()
                    metrics_rows.append(dict(library=library,depth=depth,cell=cell,tf=t.tf,representation=arm,cohort='nested',
                        sites=len(data),positive_sites=int(labels.sum()),center_contrast=center,left_shoulder=left,right_shoulder=right,
                        left_minus_center=left-center,right_minus_center=right-center,
                        canonical_protection=bool(center<0 and left>center and right>center),
                        standardized_functional_separation=standardized_functional_separation(normal,labels,X),
                        units='fixed outside-flank standardized profile for shape statistics',status='AVAILABLE'))
                    for label,color in [(0,'#777777'),(1,'#1675ae')]:
                        mean=data[labels==label].mean(0);ax.plot(X,mean,color=color,lw=.8,label=f'label {label}')
                        for position,value in zip(X,mean):curve_rows.append(dict(library=library,depth=depth,cell=cell,tf=t.tf,representation=arm,
                            label=label,position=position,value=value,aggregation='site_equal',strand='combined',
                            units='cuts per site per base' if arm!='LOG21_residual' else 'signed deviance residual'))
                    # All strand and count-weighted source curves are retained,
                    # including empty strata. Plot selection never uses scores.
                    taskframe=frame.set_index('site_id').reindex(sites.loc[mask,'site_id'])
                    cuts=observed[mask].sum(axis=(1,2))
                    groups=[('all','all',np.ones(len(labels),bool))]
                    for field,levels in [('count_stratum',['0','1','2-4','5-9','10-19','20-49','50-99','100+']),
                                         ('expected_stratum',['<0.1','0.1-<0.5','0.5-<1','1-<3','3+'])]:
                        groups.extend((field,level,taskframe[field].eq(level).to_numpy()) for level in levels)
                    for strand in ['plus','minus','combined']:
                        if strand=='combined':strand_data=data
                        else:
                            si=0 if strand=='plus' else 1
                            strand_data=(observed[mask,si] if arm=='raw' else expected[mask,si] if arm=='LOG21_expected' else
                                calibrated_residuals(observed[mask,si],expected[mask,si],'deviance',dispersion=0. if cell=='HepG2' else 1.6144318781609228))
                        for field,level,in_stratum in groups:
                            for label in [0,1]:
                                chosen=in_stratum & (labels==label)
                                for aggregation in ['site_equal','count_weighted']:
                                    weights=np.ones(chosen.sum()) if aggregation=='site_equal' else cuts[chosen]
                                    values_mean=np.average(strand_data[chosen],axis=0,weights=weights) if weights.sum()>0 else np.full(201,np.nan)
                                    for position,value in zip(X,values_mean):curve_rows.append(dict(library=library,depth=depth,cell=cell,tf=t.tf,representation=arm,
                                        label=label,position=position,value=value,aggregation=aggregation,strand=strand,
                                        stratum_type=field,stratum=level,sites=int(chosen.sum()),weight_sum=float(weights.sum()),
                                        units='cuts per site per base' if arm!='LOG21_residual' else 'signed deviance residual'))
                    ax.set_xlabel('base pairs from motif centre',fontsize=8);ax.tick_params(labelsize=7)
                identity.append(dict(library=library,depth=depth,cell=cell,tf=t.tf,
                    detector_comparison='LOG21_binomial versus LOG21_beta_binomial',same_sites=True,
                    observed_max_difference=0.,expected_max_difference=0.,residual_max_difference=0.,
                    detector_only_shape_gain=0.,qualification='algebraic identity on identical observations, not a new biological result'))
            fig.suptitle(f'{library} / {depth}: held cohort, no score selection',fontsize=12);fig.tight_layout(rect=[0,0,1,.985]);pdf.savefig(fig);plt.close(fig)
    pd.DataFrame(metrics_rows).to_csv(out/'tables/shape_metrics.tsv',sep='\t',index=False)
    handle.close()
    pd.DataFrame(identity).to_csv(out/'tables/detector_shape_identity.tsv',sep='\t',index=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);p.add_argument('--stage',choices=['metrics','calibration','shape'],required=True)
    a=p.parse_args();globals()[a.stage](a.outdir.resolve())
