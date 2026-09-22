#!/usr/bin/env python3
"""Three-covariate matching from the complete supported training pool."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
sys.path.insert(0,str(Path(__file__).resolve().parent))
import compare_radius_training as prior
import caliper_matching
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import canonical_id,checked_path

ROOT=prior.ROOT
PROTOCOL=dict(schema='fp-tools-radius-full-pool-cohorts-v1',training_chromosomes=[f'chr{i}' for i in range(1,16)],
    covariates=['motif_score','log1p_central_201','log1p_outer_101_500'],
    scaling='standard deviation over all supported positive and negative sites per task; constant scales use 1',
    assignment='maximum cardinality inside Euclidean caliper 0.5, then minimum total squared distance',
    positive_subsampling=False,balance='largest prefix sorted by matched distance then pair ID with all absolute SMD <= 0.1',
    minimum_pairs=200,task_denominator=21,score_selection=False,shape_selection=False,
    evidence='additional_training_cohort_diagnostic',primary_replaced=False,promotion_eligible=False,
    support='old corrected and DWM expectation support plus alignment-valid and complete finite wider sequence context',
    label_scope='existing source tables include examined tuning rows; only chr1-15 selected; no tuning evaluation')


def match(features,labels):
    x,y=np.asarray(features,float),np.asarray(labels)
    if x.ndim!=2 or x.shape[1]!=3 or y.shape!=(len(x),) or not np.isfinite(x).all() or not np.isin(y,[0,1]).all():
        raise ValueError('finite three-covariate binary-labelled rows required')
    p,n=np.flatnonzero(y==1),np.flatnonzero(y==0)
    if not len(p) or not len(n): return p[:0],n[:0],np.array([]),np.full(3,np.nan),0
    scale=x.std(axis=0);scale=np.where(scale>0,scale,1.)
    distance=cdist(x[p]/scale,x[n]/scale,metric='sqeuclidean')
    pi,ni=caliper_matching.caliper_assignment(distance,.5)
    length=np.sqrt(distance[pi,ni]);p,n=p[pi],n[ni]
    order=np.lexsort((np.arange(len(p)),length))
    for count in range(len(p),1,-1):
        ix=order[:count]
        smd=prior.balance.smd(x[np.r_[p[ix],n[ix]]],np.repeat([1,0],count))
        if np.isfinite(smd).all() and np.max(abs(smd))<=.1:
            return p[ix],n[ix],length[ix],smd,len(p)
    return p[:0],n[:0],length[:0],np.full(3,np.nan),len(p)


def freeze(out):
    if out.exists(): raise ValueError('fresh directory required')
    prior.wide.execute(prior.SOURCE,verify=True)
    old_path=prior.previous.OLD;old=json.loads(old_path.read_text())
    wide_config=json.loads((prior.SOURCE/'configuration.freeze.json').read_text())
    cohort=json.loads(checked_path(old['cohort_manifest'],ROOT).read_text())
    inputs=old['inputs']+wide_config['inputs']+[record(old_path),old['cohort_manifest'],old['cohort_configuration'],
        record(prior.SOURCE/'configuration.freeze.json'),record(prior.SOURCE/'profiles_manifest.json')]
    stages={}
    for cell in ('HepG2','K562'):
        path=prior.SOURCE/(cell+'.json');stage=json.loads(path.read_text());stages[cell]=stage
        inputs += [record(path)]+list(stage['outputs'].values())+[cohort['outputs'][cell+'_dwm']]
    value=dict(protocol=PROTOCOL,inputs=inputs,stages=stages,old_configuration=record(old_path),
        code=old['code']+wide_config['code']+[record(Path(__file__)),record(Path(caliper_matching.__file__))],
        runtime=prior.wide.old.parent.extension.runtime(),scores_read=False,promotion_eligible=False)
    for entry in value['inputs']+value['code']: checked_path(entry,ROOT)
    value['configuration_id']=canonical_id(value);out.mkdir(parents=True)
    write_immutable(out/'configuration.freeze.json',value)
    print('Frozen full-pool three-covariate matching before execution',flush=True)


def execute(out,verify=False):
    config=prior.wide.old.parent.extension.read_document(out/'configuration.freeze.json','configuration_id')
    if config['runtime']!=prior.wide.old.parent.extension.runtime() or config['protocol']!=PROTOCOL:
        raise ValueError('runtime or protocol changed')
    for entry in config['inputs']+config['code']: checked_path(entry,ROOT)
    manifest=out/'cohorts_manifest.json'
    if manifest.exists():
        prior.wide.old.parent.base.verify_stage(manifest,config['configuration_id'],ROOT)
        print('Verified full-pool cohorts',flush=True);return
    if verify or list(out.glob('*.tsv')): raise ValueError('missing or partial output')
    old=json.loads(checked_path(config['old_configuration'],ROOT).read_text())
    cohort=json.loads(checked_path(old['cohort_manifest'],ROOT).read_text())
    all_sites=[];pairs=[];diagnostics=[]
    for entry in old['protocol']['profile_inputs']:
        cell=entry['cell'];stage=config['stages'][cell]
        original,arrays=prior.previous.old.base.source.read_profiles(entry,ROOT)
        sites=pd.read_csv(stage['outputs']['sites']['path'],sep='\t')
        if not sites.TFBS_chr.isin(PROTOCOL['training_chromosomes']).all(): raise ValueError('nontraining source')
        order=prior.previous.align_indices(sites,original)
        sites['old_index']=order;sites['wide_index']=np.arange(len(sites))
        sites['chip_label']=original.iloc[order].chip_label.to_numpy(int)
        sites['motif_score']=original.iloc[order].motif_score.to_numpy(float)
        observed=arrays['observed'][order];expected=arrays['expected'][order]
        with np.load(checked_path(cohort['outputs'][cell+'_dwm'],ROOT),allow_pickle=False) as data:
            np.testing.assert_array_equal(data['site_hash'],arrays['site_hash'])
            dwm=data['expected'][order];dwm_valid=data['support'][order]
        with np.load(stage['outputs']['profiles']['path'],allow_pickle=False) as data:
            np.testing.assert_array_equal(data['site_hash'],prior.wide.old.strands.site_hashes(sites))
            counts=data['plus_observed']+data['minus_observed']
            shared=counts[:,500:701].copy();reverse=sites.TFBS_strand.eq('-').to_numpy();shared[reverse]=shared[reverse,::-1]
            np.testing.assert_array_equal(shared,observed)
            support=data['bam_valid']&np.isfinite(data['plus_log_bias']).all(axis=1)&np.isfinite(data['minus_log_bias']).all(axis=1)
        sites['central_201']=observed.sum(axis=1)
        offset=np.abs(np.arange(-600,601));sites['outer_101_500']=counts[:,(offset>=101)&(offset<=500)].sum(axis=1)
        features=np.column_stack([sites.motif_score,np.log1p(sites.central_201),np.log1p(sites.outer_101_500)])
        support &= arrays['valid'][order]&dwm_valid&np.isfinite(features).all(axis=1)
        for values in (observed,expected,dwm): support &= np.isfinite(values).all(axis=1)&(values>=0).all(axis=1)
        support &= (dwm.sum(axis=1)>0)&((sites.central_201.to_numpy()==0)|(expected.sum(axis=1)>0))
        sites['supported']=support;all_sites.append(sites)
        for tf,task in sites.groupby('tf',sort=True):
            ix=task.index[support[task.index]].to_numpy();labels=sites.iloc[ix].chip_label.to_numpy()
            p,n,length,smd,assigned=match(features[ix],labels)
            status='eligible' if len(p)>=200 else ('underpowered' if len(p) else 'no_balanced_pairs')
            diagnostics.append(dict(cell=cell,tf=tf,original_sites=len(task),supported_sites=len(ix),
                supported_positive=int((labels==1).sum()),supported_negative=int((labels==0).sum()),
                caliper_pairs=assigned,matched_pairs=len(p),status=status,
                motif_smd=smd[0],central_smd=smd[1],outer_smd=smd[2]))
            for indexes in (p,n):
                selected=sites.iloc[ix[indexes]].copy();selected['pair_id']=np.arange(len(p));selected['match_distance']=length
                pairs.append(selected)
        print('Matched full training pool:',cell,flush=True)
    if len(diagnostics)!=21: raise ValueError('task denominator changed')
    outputs={}
    for name,table in [('all_sites',pd.concat(all_sites)),('pairs',pd.concat(pairs)),('diagnostics',pd.DataFrame(diagnostics))]:
        path=out/(name+'.tsv');table.to_csv(path,sep='\t',index=False);outputs[name]=record(path)
    result=dict(configuration_id=config['configuration_id'],outputs=outputs,scores_read=False,
        primary_replaced=False,promotion_eligible=False)
    result['stage_id']=canonical_id(result);write_immutable(manifest,result)
    print('Completed all 21 full-pool training cohorts',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--outdir',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--freeze-only',action='store_true');mode.add_argument('--verify-only',action='store_true')
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'): parser.error('use research results')
    if args.freeze_only: freeze(out)
    else:
        with prior.wide.old.parent.extension.exclusive_run(out): execute(out,args.verify_only)
