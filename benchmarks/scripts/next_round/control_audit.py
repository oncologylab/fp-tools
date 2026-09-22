"""Separate null mean error from excess variance on all development controls."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from verify import read_table
from prepare import dump
from access import digest


def execute(out):
    dump(out/'provenance/control_mean_variance.freeze.json',dict(code_sha256=digest(__file__),
        readout='observed-minus-expected central/left/right fractions, independent-cut variance ratio, block mean-error interval',
        bootstrap_draws=2000,seed=2026,block='chromosome',interpretation='examined control development; no mean refit or threshold change',
        strata=['all','count','expected_count','GC_fixed_quartiles','native_DWM_positive_support']))
    positions=np.arange(-100,101);masks=[abs(positions)<=7,(positions>=-40)&(positions<=-12),(positions>=12)&(positions<=40)]
    selected=np.logical_or.reduce(masks)
    weights=np.random.default_rng(2026).multinomial(15,np.repeat(1/15,15),size=2000)
    probabilities={};rows=[];per_site=[]
    with np.load(out/'arrays/NakedDNA_rep1.propensities.npz',allow_pickle=False) as z:
        gc=z['gc']
        for arm in ['LOG21','DWM']:
            probabilities[arm]=np.stack([np.column_stack([np.exp(logsumexp(z[arm+'_'+strand][:,m],axis=1)-logsumexp(z[arm+'_'+strand][:,selected],axis=1)) for m in masks]) for strand in ['plus','minus']],axis=1)
    probabilities['uniform']=np.broadcast_to(np.asarray([15,29,29])/73,(4200,2,3))
    support=read_table(out/'tables/dwm_support.v3.tsv')
    for depth in ['10M','full']:
        control=read_table(out/f'tables/control_counts.{depth}.tsv.gz')
        counts=np.stack([control[[strand+'_'+region for region in ['center','left','right']]].to_numpy(int) for strand in ['plus','minus']],axis=1)
        total=counts.sum(2)
        control['GC_stratum']=pd.cut(gc,[-.01,.25,.5,.75,1.01],labels=['0-.25','.25-.5','.5-.75','.75-1']).astype(str)
        control['DWM_support']=np.where(support.full_expected_sum.to_numpy()>0,'positive_native_expectation','valid_zero_native_expectation')
        for arm,prob in probabilities.items():
            for strand,si in [('plus',0),('minus',1),('combined',None)]:
                for region,ri in [('center',0),('left',1),('right',2)]:
                    observed=counts[:,:,ri].sum(1) if si is None else counts[:,si,ri]
                    expected=(total*prob[:,:,ri]).sum(1) if si is None else total[:,si]*prob[:,si,ri]
                    denominator=total.sum(1) if si is None else total[:,si]
                    variance=(total*prob[:,:,ri]*(1-prob[:,:,ri])).sum(1) if si is None else total[:,si]*prob[:,si,ri]*(1-prob[:,si,ri])
                    for (cell,tf),group in control.groupby(['cell','tf']):
                        indices=group.index[group.eligible]
                        strata=[('all','all',indices)]
                        for field in ['count_stratum','expected_stratum','GC_stratum','DWM_support']:
                            levels={'count_stratum':['0','1','2-4','5-9','10-19','20-49','50-99','100+'],
                                    'expected_stratum':['<0.1','0.1-<0.5','0.5-<1','1-<3','3+'],
                                    'GC_stratum':['0-.25','.25-.5','.5-.75','.75-1'],
                                    'DWM_support':['positive_native_expectation','valid_zero_native_expectation']}[field]
                            strata.extend((field,level,indices[group.loc[indices,field].eq(level)]) for level in levels)
                        for field,level,ix in strata:
                            block=pd.DataFrame(dict(chromosome=control.loc[ix,'TFBS_chr'],error=observed[ix]-expected[ix],total=denominator[ix])).groupby('chromosome')[['error','total']].sum().reindex([f'chr{i}' for i in range(1,16)],fill_value=0)
                            den=weights@block.total.to_numpy();sample=np.divide(weights@block.error.to_numpy(),den,out=np.full(2000,np.nan),where=den>0)
                            valid=np.isfinite(sample);bounds=np.quantile(sample[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan]
                            n=denominator[ix].sum();v=variance[ix].sum()
                            rows.append(dict(depth=depth,cell=cell,tf=tf,arm=arm,strand=strand,region=region,stratum_type=field,stratum=level,
                                original_task_sites=len(group),eligible_sites=len(ix),informative_sites=int((denominator[ix]>0).sum()),cuts=int(n),
                                observed=float(observed[ix].sum()),expected=float(expected[ix].sum()),
                                mean_error_fraction=float((observed[ix]-expected[ix]).sum()/n) if n else np.nan,
                                mean_error_low=bounds[0],mean_error_high=bounds[1],valid_block_draws=int(valid.sum()),
                                independent_variance_ratio=float(((observed[ix]-expected[ix])**2).sum()/v) if v else np.nan,
                                biological_replication=1,mean_refitted=False))
        print('Null mean/variance audit complete',depth,flush=True)
    result=pd.DataFrame(rows);result.to_csv(out/'tables/control_mean_variance.tsv.gz',sep='\t',index=False)
    tasks=read_table(out/'tables/task_ledger.tsv');fig,axes=plt.subplots(1,2,figsize=(10,9))
    for ax,depth in zip(axes,['10M','full']):
        subset=result[result.depth.eq(depth)&result.strand.eq('combined')&result.region.eq('center')&result.stratum_type.eq('all')]
        matrix=subset.pivot(index=['cell','tf'],columns='arm',values='mean_error_fraction').reindex(pd.MultiIndex.from_frame(tasks[['cell','tf']]))
        im=ax.imshow(matrix.to_numpy(),aspect='auto',cmap='coolwarm',vmin=-.1,vmax=.1)
        ax.set_xticks(np.arange(len(matrix.columns)),matrix.columns,rotation=30);ax.set_yticks(np.arange(len(tasks)),tasks.cell+' '+tasks.tf,fontsize=7)
        ax.set_title(depth+' complete fragments');fig.colorbar(im,ax=ax,label='observed minus expected central fraction')
    fig.tight_layout();fig.savefig(out/'figures/all_task_control_calibration.pdf');plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);a=p.parse_args();execute(a.outdir.resolve())
