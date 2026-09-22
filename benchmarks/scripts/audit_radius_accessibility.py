#!/usr/bin/env python3
"""Frozen training accessibility audit and uniform-sequence spatial comparator."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score,average_precision_score
sys.path.insert(0,str(Path(__file__).resolve().parent))
import compare_radius_training as prior
from fixed_shape_normalization import normalize_fixed
from fp_tools.tools.functional_footprints import standardized_functional_separation
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import checked_path,canonical_id

ROOT=prior.ROOT
SOURCE=ROOT/'benchmarks/results/radius_training_comparison_20260922_v1'
PROTOCOL=ROOT/'benchmarks/manifests/radius_accessibility_audit_v1.json'
FEATURES={'central_201':(0,100),'outer_50_100':(50,100),'outer_50_200':(50,200),
          'outer_50_500':(50,500),'outer_101_500':(101,500)}


def count_features(plus,minus):
    plus,minus=np.asarray(plus,float),np.asarray(minus,float)
    if plus.shape!=minus.shape or plus.ndim!=2 or plus.shape[1]!=1201:
        raise ValueError('aligned 1201-position counts required')
    counts=plus+minus
    if not np.isfinite(counts).all() or (plus<0).any() or (minus<0).any(): raise ValueError('invalid counts')
    x=np.abs(np.arange(-600,601))
    return {name:counts[:,(x>=low)&(x<=high)].sum(axis=1) for name,(low,high) in FEATURES.items()}


def association(values,labels):
    values,labels=np.asarray(values,float),np.asarray(labels,int)
    if values.shape!=labels.shape or not np.isfinite(values).all() or set(labels)!={0,1}:
        raise ValueError('finite feature and both binary classes required')
    return dict(smd=float(prior.balance.smd(values[:,None],labels)[0]),
        auroc_increasing=roc_auc_score(labels,values),auroc_decreasing=roc_auc_score(labels,-values))


def freeze(out):
    if out.exists(): raise ValueError('fresh output directory required')
    prior.execute(SOURCE,verify=True)
    config=json.loads((SOURCE/'configuration.freeze.json').read_text())
    result=json.loads((SOURCE/'comparison_manifest.json').read_text())
    value=dict(protocol=json.loads(PROTOCOL.read_text()),stages=config['stages'],source_outputs=result['outputs'],
        inputs=config['inputs']+[record(PROTOCOL),record(SOURCE/'configuration.freeze.json'),
            record(SOURCE/'comparison_manifest.json')]+list(result['outputs'].values()),
        code=config['code']+[record(Path(__file__))],runtime=prior.wide.old.parent.extension.runtime(),
        training_only=True,promotion_eligible=False)
    for entry in value['inputs']+value['code']: checked_path(entry,ROOT)
    value['configuration_id']=canonical_id(value);out.mkdir(parents=True)
    write_immutable(out/'configuration.freeze.json',value)
    print('Frozen accessibility audit before label associations or uniform scoring',flush=True)


def execute(out,verify=False):
    config=prior.wide.old.parent.extension.read_document(out/'configuration.freeze.json','configuration_id')
    if config['runtime']!=prior.wide.old.parent.extension.runtime() or config['protocol']!=json.loads(PROTOCOL.read_text()):
        raise ValueError('runtime or protocol changed')
    for entry in config['inputs']+config['code']: checked_path(entry,ROOT)
    manifest=out/'audit_manifest.json'
    if manifest.exists():
        prior.wide.old.parent.base.verify_stage(manifest,config['configuration_id'],ROOT)
        print('Verified accessibility audit',flush=True);return
    if verify or list(out.glob('*.tsv')): raise ValueError('missing or partial output')
    scores=pd.read_csv(config['source_outputs']['scores']['path'],sep='\t',float_precision='round_trip')
    sites=scores.loc[scores.arm.eq('raw'),prior.previous.KEYS+['chip_label']].reset_index(drop=True)
    if not sites.TFBS_chr.isin([f'chr{i}' for i in range(1,16)]).all(): raise ValueError('nontraining sites')
    values={key:np.full(len(sites),np.nan) for key in FEATURES};uniform=np.full((len(sites),201),np.nan)
    for cell,stage in config['stages'].items():
        ix=np.flatnonzero(sites.cell.eq(cell));source=pd.read_csv(stage['outputs']['sites']['path'],sep='\t')
        order=prior.previous.align_indices(sites.iloc[ix],source)
        with np.load(stage['outputs']['profiles']['path'],allow_pickle=False) as archive:
            np.testing.assert_array_equal(archive['site_hash'],prior.wide.old.strands.site_hashes(source))
            if not archive['bam_valid'][order].all(): raise ValueError('invalid retained alignment window')
            plus,minus=archive['plus_observed'][order],archive['minus_observed'][order]
        for name,feature in count_features(plus,minus).items(): values[name][ix]=feature
        zero=np.zeros(1201)
        for row,dest in enumerate(ix):
            result=prior.radius_model.score_tracks(plus[row],minus[row],zero,zero,500,
                0. if cell=='HepG2' else 1.6144318781609228)
            residual=result['corrected'][500:701]
            uniform[dest]=residual[::-1] if sites.iloc[dest].TFBS_strand=='-' else residual
        print('Computed uniform-background control:',cell,flush=True)
    if not np.isfinite(uniform).all(): raise ValueError('uniform comparator changed support')
    weights=prior.previous.old.candidates.geometry_weights((7,12,40));uniform_score=uniform@weights
    uniform_score[values['central_201']==0]=0
    frame=sites.copy();frame['arm']='uniform_radius_500';frame['score']=uniform_score
    scores=pd.concat([scores,frame],ignore_index=True)
    rows=[];correlations=[];metrics=[];curves=[];x=np.arange(-100,101)
    for (cell,tf),group in sites.groupby(['cell','tf'],sort=False):
        ix=group.index.to_numpy();y=group.chip_label.to_numpy(int)
        for name,feature in values.items():
            transformed=np.log1p(feature[ix]);rows.append(dict(cell=cell,tf=tf,feature=name,sites=len(ix),
                positive_median_count=float(np.median(feature[ix][y==1])),negative_median_count=float(np.median(feature[ix][y==0])),
                **association(transformed,y)))
            for arm,table in scores.loc[scores.cell.eq(cell)&scores.tf.eq(tf)].groupby('arm',sort=False):
                order=prior.previous.align_indices(group,table);s=table.iloc[order].score.to_numpy()
                rho=float(spearmanr(transformed,s).statistic) if np.ptp(transformed)>0 and np.ptp(s)>0 else np.nan
                correlations.append(dict(cell=cell,tf=tf,feature=name,arm=arm,spearman=rho))
        norm=normalize_fixed(uniform[ix],x,1.);delta=norm[y==1].mean(axis=0)-norm[y==0].mean(axis=0)
        metrics.append(dict(cell=cell,tf=tf,arm='uniform_radius_500',matched_pairs=len(ix)//2,
            auroc=roc_auc_score(y,uniform_score[ix]),auprc=average_precision_score(y,uniform_score[ix]),
            center_contrast=delta[abs(x)<=5].mean(),left_shoulder=delta[(x>=-40)&(x<=-15)].mean(),
            right_shoulder=delta[(x>=15)&(x<=40)].mean(),canonical_protection=prior.previous.old.candidates.protection_direction(delta),
            functional_separation=standardized_functional_separation(norm,y,x)))
        curves.append(pd.DataFrame(dict(cell=cell,tf=tf,arm='uniform_radius_500',position=x,
            positive_mean=norm[y==1].mean(axis=0),negative_mean=norm[y==0].mean(axis=0))))
    if len(metrics)!=21 or len(rows)!=105 or len(correlations)!=735: raise ValueError('task denominator changed')
    outputs={}
    features=sites.copy()
    for name,feature in values.items(): features[name]=feature
    for name,table in [('features',features),('associations',pd.DataFrame(rows)),('correlations',pd.DataFrame(correlations)),
        ('uniform_metrics',pd.DataFrame(metrics)),('uniform_scores',frame),('uniform_curves',pd.concat(curves))]:
        path=out/(name+'.tsv');table.to_csv(path,sep='\t',index=False);outputs[name]=record(path)
    result=dict(configuration_id=config['configuration_id'],outputs=outputs,training_only=True,
        promotion_eligible=False,candidate_selected=False,sites_rematched=False)
    result['stage_id']=canonical_id(result);write_immutable(manifest,result)
    print('Completed all 21 tasks without rematching or candidate selection',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--outdir',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--freeze-only',action='store_true');mode.add_argument('--verify-only',action='store_true')
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'): parser.error('use research results')
    if args.freeze_only: freeze(out)
    else:
        with prior.wide.old.parent.extension.exclusive_run(out): execute(out,args.verify_only)
