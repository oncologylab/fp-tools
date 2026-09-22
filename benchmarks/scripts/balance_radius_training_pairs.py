#!/usr/bin/env python3
"""Freeze a covariate-only, complete-pair secondary training cohort."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_radius_accessibility as parent
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import checked_path,canonical_id

ROOT=parent.ROOT
SOURCE=ROOT/'benchmarks/results/radius_accessibility_audit_20260922_v1'
PROTOCOL=dict(schema='fp-tools-radius-secondary-balance-v1',training_only=True,
    covariates=['motif_score','log1p_central_201','log1p_outer_101_500'],
    ordering='Euclidean paired covariate differences divided by within-class pooled standard deviations; constant scales use 1; stable pair-ID ties',
    selection='largest prefix with absolute SMD <= 0.1 in all three covariates; minimum two pairs; no rematching or new caliper',
    preserve_original_if_balanced=True,minimum_pairs_for_power=200,task_denominator=21,
    score_or_shape_selection=False,primary_result_replaced=False,promotion_eligible=False,
    limitation='Mean balance does not ensure distributional balance or identify causal confounding.')


def select_pairs(positive,negative,pair_ids):
    p,n=np.asarray(positive,float),np.asarray(negative,float);ids=np.asarray(pair_ids)
    if (p.ndim!=2 or p.shape!=n.shape or p.shape[1]!=3 or ids.shape!=(len(p),)
            or len(np.unique(ids))!=len(ids) or not np.isfinite(p).all() or not np.isfinite(n).all()):
        raise ValueError('finite three-covariate complete pairs with unique IDs required')
    if len(p)<2: return np.array([],int),np.full(3,np.nan),np.zeros(len(p))
    scale=np.sqrt((p.var(axis=0)+n.var(axis=0))/2);scale=np.where(scale>0,scale,1.)
    distance=np.linalg.norm((p-n)/scale,axis=1);order=np.lexsort((ids,distance))
    for count in range(len(order),1,-1):
        ix=order[:count]
        smd=parent.prior.balance.smd(np.concatenate([p[ix],n[ix]]),np.repeat([1,0],count))
        if np.isfinite(smd).all() and np.max(np.abs(smd))<=.1:
            return ix,smd,distance
    return np.array([],int),np.full(3,np.nan),distance


def freeze(out):
    if out.exists(): raise ValueError('fresh output directory required')
    parent.execute(SOURCE,verify=True)
    original=json.loads((SOURCE/'configuration.freeze.json').read_text())
    result=json.loads((SOURCE/'audit_manifest.json').read_text())
    old_path=parent.prior.previous.OLD;old=json.loads(old_path.read_text())
    config=dict(protocol=PROTOCOL,inputs=original['inputs']+[record(SOURCE/'configuration.freeze.json'),
        record(SOURCE/'audit_manifest.json'),result['outputs']['features'],record(old_path)],
        features=result['outputs']['features'],old_configuration=record(old_path),
        code=original['code']+[record(Path(__file__))],runtime=parent.prior.wide.old.parent.extension.runtime(),
        scores_read=False,promotion_eligible=False)
    for entry in config['inputs']+config['code']: checked_path(entry,ROOT)
    config['configuration_id']=canonical_id(config);out.mkdir(parents=True)
    write_immutable(out/'configuration.freeze.json',config)
    print('Frozen secondary covariate-only pair selection',flush=True)


def execute(out,verify=False):
    config=parent.prior.wide.old.parent.extension.read_document(out/'configuration.freeze.json','configuration_id')
    if config['protocol']!=PROTOCOL or config['runtime']!=parent.prior.wide.old.parent.extension.runtime():
        raise ValueError('protocol or runtime changed')
    for entry in config['inputs']+config['code']: checked_path(entry,ROOT)
    manifest=out/'cohorts_manifest.json'
    if manifest.exists():
        parent.prior.wide.old.parent.base.verify_stage(manifest,config['configuration_id'],ROOT)
        print('Verified secondary balanced cohorts',flush=True);return
    if verify or list(out.glob('*.tsv')): raise ValueError('missing or partial output')
    features=pd.read_csv(config['features']['path'],sep='\t',float_precision='round_trip')
    old=json.loads(checked_path(config['old_configuration'],ROOT).read_text())
    cohort=json.loads(checked_path(old['cohort_manifest'],ROOT).read_text())
    selected=[];diagnostics=[];all_pairs=[]
    keys=parent.prior.previous.KEYS
    for entry in old['protocol']['profile_inputs']:
        cell=entry['cell'];stage=json.loads(checked_path(entry,ROOT).read_text())
        # Read coordinates and motif score only: no source tuning labels.
        coordinates=pd.read_csv(checked_path(stage['outputs']['sites'],ROOT),sep='\t',usecols=keys+['motif_score'])
        pairs=pd.read_csv(checked_path(cohort['outputs'][cell+'_pairs'],ROOT),sep='\t')
        pairs=pairs.loc[pairs.split.eq('train')].reset_index(drop=True)
        sites=coordinates.iloc[pairs.site_index].reset_index(drop=True)
        if not sites.TFBS_chr.isin([f'chr{i}' for i in range(1,16)]).all(): raise ValueError('nontraining row')
        order=parent.prior.previous.align_indices(sites,features)
        aligned=features.iloc[order].reset_index(drop=True)
        np.testing.assert_array_equal(pairs.chip_label,aligned.chip_label)
        sites['chip_label']=pairs.chip_label;sites['pair_id']=pairs.pair_id
        sites['central_201']=aligned.central_201;sites['outer_101_500']=aligned.outer_101_500
        for tf,group in sites.groupby('tf',sort=True):
            positive=group.loc[group.chip_label.eq(1)].sort_values('pair_id')
            negative=group.loc[group.chip_label.eq(0)].sort_values('pair_id')
            np.testing.assert_array_equal(positive.pair_id,negative.pair_id)
            def covariates(table): return np.column_stack([table.motif_score,np.log1p(table.central_201),np.log1p(table.outer_101_500)])
            p,n=covariates(positive),covariates(negative)
            ix,smd,distance=select_pairs(p,n,positive.pair_id.to_numpy())
            before=parent.prior.balance.smd(np.concatenate([p,n]),np.repeat([1,0],len(p)))
            diagnostics.append(dict(cell=cell,tf=tf,original_pairs=len(p),retained_pairs=len(ix),
                status='eligible' if len(ix)>=200 else ('underpowered' if len(ix) else 'no_balanced_prefix'),
                before_motif_smd=before[0],before_central_smd=before[1],before_outer_smd=before[2],
                after_motif_smd=smd[0],after_central_smd=smd[1],after_outer_smd=smd[2]))
            kept=set(positive.iloc[ix].pair_id)
            audit=group.copy();audit['selected']=audit.pair_id.isin(kept);all_pairs.append(audit)
            selected.append(pd.concat([positive.iloc[ix],negative.iloc[ix]],ignore_index=True))
    if len(diagnostics)!=21: raise ValueError('task denominator changed')
    outputs={}
    for name,table in [('diagnostics',pd.DataFrame(diagnostics)),('selected_sites',pd.concat(selected)),('all_sites',pd.concat(all_pairs))]:
        path=out/(name+'.tsv');table.to_csv(path,sep='\t',index=False);outputs[name]=record(path)
    result=dict(configuration_id=config['configuration_id'],outputs=outputs,scores_read=False,
        primary_result_replaced=False,promotion_eligible=False)
    result['stage_id']=canonical_id(result);write_immutable(manifest,result)
    print('Prepared all 21 secondary cohorts without reading footprint scores',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--outdir',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--freeze-only',action='store_true');mode.add_argument('--verify-only',action='store_true')
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'): parser.error('use research results')
    if args.freeze_only: freeze(out)
    else:
        with parent.prior.wide.old.parent.extension.exclusive_run(out): execute(out,args.verify_only)
