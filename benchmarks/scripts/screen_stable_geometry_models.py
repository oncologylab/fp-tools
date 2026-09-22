#!/usr/bin/env python3
"""Exact stable-geometry screen on previously examined naked-DNA replicate 2."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import freeze_stable_geometry_models as models
import screen_fixed_shape_operators as prior
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = models.ROOT
MODELS = ROOT / 'benchmarks/results/stable_geometry_models_20260922_v1'
PRIOR = ROOT / 'benchmarks/results/fixed_shape_operator_screen_20260922_v1'
SUPPORTED = ['raw', 'DWM_expected', 'parametric_direct']
PROTOCOL = dict(evidence='previously_examined_replicate2_development_screen',
    task_count=21, arms=models.source.source.ARMS, available_control_arms=SUPPORTED,
    unavailable_arms='retain ledger rows; matching wider control profiles not prepared',
    calibration_chromosomes=['chr16','chr17','chr18'], calibration_labels='negative_only',
    calibration_support='prior common valid positive observed total and positive DWM/parametric expectations',
    minimum_calibration_negatives=200, target_negative_call_rate=.05,
    threshold_rule='existing conformal_upper_threshold; strict score greater than threshold',
    baseline='stable rounded unnormalized DWM residual conventional geometry (7,12,40)',
    zero_observed='never called; score zero', model_refitting=False,
    final_control_read=False, new_holdout_labels_read=False, promotion_eligible=False)


def available(model):
    return model['status'] != 'unavailable' and model['arm'] in SUPPORTED


def baseline_scores(panels, totals):
    return models.source.stable.stable_scores(panels['DWM_expected'],
        models.source.geometry_weights((7,12,40)), np.asarray(totals)==0)


def freeze(out):
    if out.exists(): raise ValueError('fresh directory required')
    models.execute(MODELS, verify=True)
    old = prior.preflight(ROOT, PRIOR)
    config = json.loads((MODELS/'configuration.freeze.json').read_text())
    result = json.loads((MODELS/'models.freeze.json').read_text())
    value = dict(protocol=PROTOCOL, prior=old, models=result['outputs']['detectors'],
        inputs=config['inputs'] + old['inputs'] + [record(MODELS/'configuration.freeze.json'),
            record(MODELS/'models.freeze.json'),record(PRIOR/'configuration.freeze.json')] + list(result['outputs'].values()),
        code=config['code'] + old['code'] + [record(Path(__file__))], runtime=prior.builder.extension.runtime())
    value['configuration_id']=canonical_id(value); out.mkdir(parents=True)
    write_immutable(out/'configuration.freeze.json',value)


def preflight(out):
    config=prior.builder.extension.read_document(out/'configuration.freeze.json','configuration_id')
    if canonical_id(config['protocol'])!=canonical_id(PROTOCOL) or config['runtime']!=prior.builder.extension.runtime():
        raise ValueError('protocol or runtime changed')
    for item in config['inputs']+config['code']: checked_path(item,ROOT)
    artifact=prior.builder.extension.read_document(checked_path(config['models'],ROOT),'artifact_id')
    fitted=artifact['models']
    if len(fitted)!=147 or len({(m['cell'],m['tf'],m['arm']) for m in fitted})!=147:
        raise ValueError('model denominator changed')
    return config,fitted


def calibrate(out):
    config,fitted=preflight(out);target=out/'thresholds.freeze.json'
    if target.exists():
        prior.builder.base.verify_stage(target,config['configuration_id'],ROOT);print('Verified thresholds; no recalibration',flush=True);return
    if (out/'thresholds.tsv').exists():raise ValueError('partial calibration preserved')
    rows=[]
    for cell in ['HepG2','K562']:
        sites,totals,panels,_,_,common=prior.load_cell(config['prior'],ROOT,cell,False)
        baseline=baseline_scores(panels,totals)
        for model in [m for m in fitted if m['cell']==cell]:
            ix=np.flatnonzero(sites.tf.eq(model['tf']));task=sites.iloc[ix]
            negatives=task.TFBS_chr.isin(PROTOCOL['calibration_chromosomes']).to_numpy() & task.chip_label.eq(0).to_numpy()
            mask=negatives & common[ix] & np.isfinite(baseline[ix])
            score=models.score_model(model,panels[model['arm']][ix],totals[ix]==0) if available(model) else None
            usable=mask & np.isfinite(score) if score is not None else np.zeros(len(ix),bool)
            n=int(usable.sum());eligible=n>=200
            rows.append(dict(cell=cell,tf=model['tf'],arm=model['arm'],model_status=model['status'],
                model_available=model['status']!='unavailable',control_profiles_available=model['arm'] in SUPPORTED,
                calibration_negatives=n,calibration_exclusions=int(negatives.sum())-n,threshold_eligible=eligible,
                candidate_threshold=prior.safety.conformal_upper_threshold(score[usable],.05) if eligible else np.nan,
                baseline_threshold=prior.safety.conformal_upper_threshold(baseline[ix][usable],.05) if eligible else np.nan))
    table=pd.DataFrame(rows)
    if len(table)!=147:raise ValueError('threshold denominator changed')
    path=out/'thresholds.tsv';table.to_csv(path,sep='\t',index=False)
    result=dict(configuration_id=config['configuration_id'],outputs={'thresholds':record(path)},control_profiles_scored=False)
    result['stage_id']=canonical_id(result);write_immutable(target,result)
    print('Frozen all 147 thresholds before control scoring',flush=True)


def evaluate(out,verify=False):
    config,fitted=preflight(out)
    thresholds=prior.builder.base.verify_stage(out/'thresholds.freeze.json',config['configuration_id'],ROOT)
    if thresholds['control_profiles_scored'] is not False:raise ValueError('threshold provenance changed')
    target=out/'screen_manifest.json'
    if target.exists():
        result=prior.builder.base.verify_stage(target,config['configuration_id'],ROOT)
        if result['thresholds']!=record(out/'thresholds.freeze.json'):raise ValueError('threshold identity changed')
        print('Verified exact stable-geometry control screen',flush=True);return
    if verify or list(out.glob('screen_*.tsv')):raise ValueError('missing or partial screen')
    table=pd.read_csv(thresholds['outputs']['thresholds']['path'],sep='\t',float_precision='round_trip').set_index(['cell','tf','arm'])
    rows=[];site_rows=[]
    for cell in ['HepG2','K562']:
        sites,totals,panels,_,valid,comparable=prior.load_cell(config['prior'],ROOT,cell,True)
        baseline=baseline_scores(panels,totals)
        for model in [m for m in fitted if m['cell']==cell]:
            tf,arm=model['tf'],model['arm'];ix=np.flatnonzero(sites.tf.eq(tf));threshold=table.loc[(cell,tf,arm)]
            score=models.score_model(model,panels[arm][ix],totals[ix]==0) if available(model) else np.full(len(ix),np.nan)
            support=valid[arm][ix] if available(model) else np.zeros(len(ix),bool)
            row,calls,base_calls=prior.rate_summary(score,baseline[ix],totals[ix],support,comparable[ix],
                threshold.candidate_threshold,threshold.baseline_threshold,bool(threshold.threshold_eligible))
            status='screened' if threshold.threshold_eligible else 'insufficient_calibration'
            if not threshold.model_available:status='missing_model'
            elif not threshold.control_profiles_available:status='missing_control_profiles'
            row.update(cell=cell,tf=tf,arm=arm,status=status,**threshold.to_dict());rows.append(row)
            frame=sites.iloc[ix][prior.safety.SITE_KEYS].copy();frame['arm']=arm
            for name,value in dict(candidate_score=score,baseline_score=baseline[ix],total_cuts=totals[ix],
                candidate_valid=support & np.isfinite(score),baseline_comparable=comparable[ix],candidate_call=calls,
                baseline_call=base_calls,candidate_threshold=threshold.candidate_threshold,
                baseline_threshold=threshold.baseline_threshold,threshold_eligible=threshold.threshold_eligible).items():frame[name]=value
            site_rows.append(frame)
    if len(rows)!=147:raise ValueError('control denominator changed')
    outputs={}
    for name,data in [('summary',pd.DataFrame(rows)),('sites',pd.concat(site_rows,ignore_index=True))]:
        path=out/f'screen_{name}.tsv';data.to_csv(path,sep='\t',index=False);outputs[name]=record(path)
    result=dict(configuration_id=config['configuration_id'],outputs=outputs,thresholds=record(out/'thresholds.freeze.json'),
        evidence=PROTOCOL['evidence'],promotion_eligible=False,model_refitted=False,final_control_read=False,new_holdout_labels_read=False)
    result['stage_id']=canonical_id(result);write_immutable(target,result);print('Completed all 147 development safety ledger rows',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--outdir',type=Path,required=True)
    parser.add_argument('--stage',choices=['freeze','calibrate','evaluate','verify'],required=True)
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'):parser.error('use research results')
    if args.stage=='freeze':freeze(out)
    else:
        with prior.builder.extension.exclusive_run(out):
            if args.stage=='calibrate':calibrate(out)
            else:evaluate(out,args.stage=='verify')
