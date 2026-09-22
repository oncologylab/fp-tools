#!/usr/bin/env python3
"""Frozen chromosome-block sensitivity for training scores and signed shape."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
import evaluate_radius_full_pool as source
from fixed_shape_normalization import normalize_fixed
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import canonical_id,checked_path

ROOT=source.ROOT
SOURCE=ROOT/'benchmarks/results/radius_full_pool_comparison_20260922_v1'
PROTOCOL=dict(schema='fp-tools-radius-training-chromosome-bootstrap-v1',replicates=2000,
    chromosomes=[f'chr{i}' for i in range(1,16)],seeds={'HepG2':2026,'K562':2027},
    resampling='15 chromosome draws with replacement; identical multiplicities across all methods and tasks within each cell',
    matching='original fixed cohort; no rematching; chromosome reweighting may change pair and covariate balance',
    metrics=['AUROC gain over DWM','relative AUPRC gain over DWM','signed center contrast','left shoulder','right shoulder','canonical direction fraction'],
    intervals='pointwise percentile 2.5 and 97.5; not calibrated simultaneous functional bands',
    invalid='retain missing-class resamples as NaN, report counts, no redrawing',
    candidate_selection=False,promotion_eligible=False,principal_holdout_evidence=False,
    multiplicity_adjusted=False,training_only=True)


def weighted_ranks(scores,labels,blocks,draws):
    s,y,b=np.asarray(scores,float),np.asarray(labels),np.asarray(blocks)
    draws=np.asarray(draws,float)
    if (s.ndim!=1 or y.shape!=s.shape or b.shape!=s.shape or draws.ndim!=2
            or not np.isfinite(s).all() or not np.isin(y,[0,1]).all() or not np.issubdtype(b.dtype,np.integer)
            or (b<0).any() or (b>=draws.shape[1]).any() or not np.isfinite(draws).all() or (draws<0).any()):
        raise ValueError('aligned finite scores and valid block multiplicities required')
    _,group=np.unique(s,return_inverse=True);width=int(group.max())+1 if len(group) else 0
    pos=np.zeros((draws.shape[1],width));neg=np.zeros_like(pos)
    np.add.at(pos,(b[y==1],group[y==1]),1);np.add.at(neg,(b[y==0],group[y==0]),1)
    p,n=draws@pos,draws@neg;pt,nt=p.sum(1),n.sum(1);valid=(pt>0)&(nt>0)
    roc=np.full(len(draws),np.nan);ap=np.full(len(draws),np.nan)
    roc[valid]=(p[valid]*(np.cumsum(n[valid],axis=1)-.5*n[valid])).sum(1)/(pt[valid]*nt[valid])
    cp=np.cumsum(p[valid,::-1],axis=1);total=cp+np.cumsum(n[valid,::-1],axis=1)
    precision=np.divide(cp,total,out=np.zeros_like(cp),where=total>0)
    ap[valid]=(precision*p[valid,::-1]).sum(1)/pt[valid]
    return roc,ap


def weighted_shape(values,labels,blocks,draws):
    labels=np.asarray(labels);blocks=np.asarray(blocks)
    results=[];totals=[]
    for label in (1,0):
        count=np.bincount(blocks[labels==label],minlength=draws.shape[1])
        sums=np.zeros((draws.shape[1],values.shape[1]));np.add.at(sums,blocks[labels==label],values[labels==label])
        total=draws@count;totals.append(total)
        results.append(np.divide(draws@sums,total[:,None],out=np.full((len(draws),values.shape[1]),np.nan),where=total[:,None]>0))
    return results[0]-results[1]


def freeze(out):
    if out.exists(): raise ValueError('fresh directory required')
    source.execute(SOURCE,verify=True)
    old=json.loads((SOURCE/'configuration.freeze.json').read_text());result=json.loads((SOURCE/'comparison_manifest.json').read_text())
    config=dict(protocol=PROTOCOL,inputs=old['inputs']+[record(SOURCE/'configuration.freeze.json'),
        record(SOURCE/'comparison_manifest.json')]+list(result['outputs'].values()),outputs_source=result['outputs'],
        code=old['code']+[record(Path(__file__))],runtime=source.cohorts.prior.wide.old.parent.extension.runtime(),promotion_eligible=False)
    for entry in config['inputs']+config['code']: checked_path(entry,ROOT)
    config['configuration_id']=canonical_id(config);out.mkdir(parents=True);write_immutable(out/'configuration.freeze.json',config)


def execute(out,verify=False):
    config=source.cohorts.prior.wide.old.parent.extension.read_document(out/'configuration.freeze.json','configuration_id')
    if config['protocol']!=PROTOCOL or config['runtime']!=source.cohorts.prior.wide.old.parent.extension.runtime():
        raise ValueError('protocol or runtime changed')
    for entry in config['inputs']+config['code']: checked_path(entry,ROOT)
    target=out/'bootstrap_manifest.json'
    if target.exists():
        source.cohorts.prior.wide.old.parent.base.verify_stage(target,config['configuration_id'],ROOT)
        print('Verified training chromosome-block sensitivity',flush=True);return
    if verify or list(out.glob('*.tsv')) or list(out.glob('*.npz')): raise ValueError('missing or partial output')
    paths=config['outputs_source'];read=lambda name:pd.read_csv(paths[name]['path'],sep='\t',float_precision='round_trip')
    sites,metrics,scores=read('sites'),read('metrics'),read('scores')
    if not sites.TFBS_chr.isin(PROTOCOL['chromosomes']).all(): raise ValueError('nontraining row')
    with np.load(paths['panels']['path'],allow_pickle=False) as archive: panels={arm:archive[arm] for arm in source.ARMS}
    draws={cell:np.random.default_rng(seed).multinomial(15,np.full(15,1/15),size=2000) for cell,seed in PROTOCOL['seeds'].items()}
    blocks_map={chrom:i for i,chrom in enumerate(PROTOCOL['chromosomes'])};rows=[];saved={};x=np.arange(-100,101)
    for task in metrics[['cell','tf','matched_pairs','status']].drop_duplicates().itertuples(index=False):
        ix=np.flatnonzero((sites.cell.eq(task.cell)&sites.tf.eq(task.tf)).to_numpy());y=sites.iloc[ix].chip_label.to_numpy(int)
        blocks=sites.iloc[ix].TFBS_chr.map(blocks_map).to_numpy(int);weights=draws[task.cell]
        ranks={};shape={}
        for arm in source.ARMS:
            if not len(ix):
                ranks[arm]=(np.full(2000,np.nan),np.full(2000,np.nan));shape[arm]=np.full((2000,3),np.nan);continue
            frame=scores.loc[scores.cell.eq(task.cell)&scores.tf.eq(task.tf)&scores.arm.eq(arm)]
            pd.testing.assert_frame_equal(frame[source.cohorts.prior.previous.KEYS+['chip_label']].reset_index(drop=True),
                sites.iloc[ix][source.cohorts.prior.previous.KEYS+['chip_label']].reset_index(drop=True))
            ranks[arm]=weighted_ranks(frame.score,y,blocks,weights)
            norm=normalize_fixed(panels[arm][ix],x,1.)
            summaries=np.column_stack([norm[:,abs(x)<=5].mean(1),norm[:,(x>=-40)&(x<=-15)].mean(1),norm[:,(x>=15)&(x<=40)].mean(1)])
            shape[arm]=weighted_shape(summaries,y,blocks,weights)
        base_roc,base_ap=ranks['DWM_expected']
        for arm in source.ARMS:
            roc,ap=ranks[arm];delta=shape[arm];valid=np.isfinite(roc)&np.isfinite(ap)&np.isfinite(delta).all(1)
            values=np.column_stack([roc-base_roc,np.divide(ap-base_ap,base_ap,out=np.full(2000,np.nan),where=base_ap>0),delta])
            row=dict(cell=task.cell,tf=task.tf,arm=arm,matched_pairs=task.matched_pairs,status=task.status,
                valid_resamples=int(valid.sum()),invalid_resamples=int((~valid).sum()),observed_chromosomes=len(np.unique(blocks)))
            for j,name in enumerate(['auroc_gain','relative_auprc_gain','center','left','right']):
                finite=values[:,j][np.isfinite(values[:,j])]
                low,high=np.quantile(finite,[.025,.975]) if len(finite) else [np.nan,np.nan]
                row[name+'_lower'],row[name+'_upper']=low,high
            row['canonical_fraction']=float(((delta[valid,0]<0)&(delta[valid,1]>delta[valid,0])&(delta[valid,2]>delta[valid,0])).mean()) if valid.any() else np.nan
            rows.append(row);saved[f'{task.cell}.{task.tf}.{arm}']=values
        print('Bootstrapped:',task.cell,task.tf,flush=True)
    if len(rows)!=147: raise ValueError('task denominator changed')
    outputs={}
    path=out/'intervals.tsv';pd.DataFrame(rows).to_csv(path,sep='\t',index=False);outputs['intervals']=record(path)
    path=out/'draws.npz';np.savez_compressed(path,**draws);outputs['draws']=record(path)
    path=out/'resamples.npz';np.savez_compressed(path,**saved);outputs['resamples']=record(path)
    result=dict(configuration_id=config['configuration_id'],outputs=outputs,training_only=True,
        principal_holdout_evidence=False,uncertainty_calibrated=False,promotion_eligible=False)
    result['stage_id']=canonical_id(result);write_immutable(target,result)
    print('Completed all 147 training uncertainty rows',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--outdir',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--freeze-only',action='store_true');mode.add_argument('--verify-only',action='store_true')
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'): parser.error('use research results')
    if args.freeze_only: freeze(out)
    else:
        with source.cohorts.prior.wide.old.parent.extension.exclusive_run(out): execute(out,args.verify_only)
