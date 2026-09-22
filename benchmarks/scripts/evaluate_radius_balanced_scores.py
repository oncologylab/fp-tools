#!/usr/bin/env python3
"""Recompute ranks on frozen secondary cohorts; no new fit or shape claim."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score
sys.path.insert(0,str(Path(__file__).resolve().parent))
import balance_radius_training_pairs as cohorts
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import canonical_id,checked_path
ROOT=cohorts.ROOT
SOURCE=ROOT/'benchmarks/results/radius_balanced_pairs_20260922_v1'
ARMS=['raw','DWM_expected','parametric_direct','radius_100','radius_200','radius_500','uniform_radius_500']


def freeze(out):
    if out.exists(): raise ValueError('fresh output directory required')
    cohorts.execute(SOURCE,verify=True)
    config=json.loads((SOURCE/'configuration.freeze.json').read_text())
    result=json.loads((SOURCE/'cohorts_manifest.json').read_text())
    primary=cohorts.parent.SOURCE/'comparison_manifest.json'
    control=cohorts.SOURCE/'audit_manifest.json'
    paths=[json.loads(primary.read_text())['outputs']['scores'],json.loads(control.read_text())['outputs']['uniform_scores']]
    value=dict(inputs=config['inputs']+[record(SOURCE/'configuration.freeze.json'),record(SOURCE/'cohorts_manifest.json'),
        record(primary),record(control)]+list(result['outputs'].values())+paths,
        code=config['code']+[record(Path(__file__))],cohort_outputs=result['outputs'],scores=paths,arms=ARMS,
        runtime=cohorts.parent.prior.wide.old.parent.extension.runtime(),training_only=True,
        evidence='secondary_covariate_balanced_rank_diagnostic',promotion_eligible=False,shape_reassessed=False)
    for entry in value['inputs']+value['code']: checked_path(entry,ROOT)
    value['configuration_id']=canonical_id(value);out.mkdir(parents=True)
    write_immutable(out/'configuration.freeze.json',value)


def execute(out,verify=False):
    extension=cohorts.parent.prior.wide.old.parent.extension
    config=extension.read_document(out/'configuration.freeze.json','configuration_id')
    if config['runtime']!=extension.runtime() or config['arms']!=ARMS: raise ValueError('scope changed')
    for entry in config['inputs']+config['code']: checked_path(entry,ROOT)
    manifest=out/'comparison_manifest.json'
    if manifest.exists():
        cohorts.parent.prior.wide.old.parent.base.verify_stage(manifest,config['configuration_id'],ROOT)
        print('Verified secondary rank comparison',flush=True);return
    if verify or (out/'metrics.tsv').exists(): raise ValueError('missing or partial output')
    read=lambda item:pd.read_csv(item['path'],sep='\t',float_precision='round_trip')
    selected=read(config['cohort_outputs']['selected_sites']);diag=read(config['cohort_outputs']['diagnostics'])
    scores=pd.concat([read(item) for item in config['scores']],ignore_index=True)
    keys=cohorts.parent.prior.previous.KEYS+['chip_label'];rows=[]
    if not selected.TFBS_chr.isin([f'chr{i}' for i in range(1,16)]).all(): raise ValueError('nontraining row')
    for task in diag.itertuples(index=False):
        sites=selected.loc[selected.cell.eq(task.cell)&selected.tf.eq(task.tf)]
        for arm in ARMS:
            frame=scores.loc[scores.cell.eq(task.cell)&scores.tf.eq(task.tf)&scores.arm.eq(arm)]
            subset=sites[keys].merge(frame,on=keys,validate='one_to_one')
            if len(subset)!=len(sites) or len(sites)!=2*task.retained_pairs: raise ValueError('support changed')
            row=dict(cell=task.cell,tf=task.tf,arm=arm,original_pairs=task.original_pairs,
                retained_pairs=task.retained_pairs,status=task.status,auroc=np.nan,auprc=np.nan)
            if len(subset):
                if subset.chip_label.sum()!=task.retained_pairs: raise ValueError('class imbalance')
                row.update(auroc=roc_auc_score(subset.chip_label,subset.score),auprc=average_precision_score(subset.chip_label,subset.score))
            rows.append(row)
    if len(rows)!=147: raise ValueError('incomplete task denominator')
    path=out/'metrics.tsv';pd.DataFrame(rows).to_csv(path,sep='\t',index=False)
    result=dict(configuration_id=config['configuration_id'],outputs=dict(metrics=record(path)),
        training_only=True,promotion_eligible=False,shape_reassessed=False,primary_result_replaced=False)
    result['stage_id']=canonical_id(result);write_immutable(manifest,result)
    print('Recomputed all 147 secondary rank rows; no shape or promotion claim',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--outdir',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--freeze-only',action='store_true');mode.add_argument('--verify-only',action='store_true')
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'): parser.error('use research results')
    if args.freeze_only: freeze(out)
    else:
        with cohorts.parent.prior.wide.old.parent.extension.exclusive_run(out): execute(out,args.verify_only)
