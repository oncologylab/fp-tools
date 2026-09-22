#!/usr/bin/env python3
"""Frozen identical-site training comparison of three background radii."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score

sys.path.insert(0,str(Path(__file__).resolve().parent))
import compare_gamma_training as previous
import radius_background as radius_model
import build_wide_training_profiles as wide
import prepare_corrected_shape_cohorts as balance
from fixed_shape_normalization import normalize_fixed
from fp_tools.tools.functional_footprints import standardized_functional_separation
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import checked_path,canonical_id

ROOT=previous.ROOT
PROTOCOL=ROOT/'benchmarks/manifests/radius_training_comparison_v1.json'
SOURCE=ROOT/'benchmarks/results/wide_training_profiles_20260922_v1'
GAMMA=ROOT/'benchmarks/results/gamma_training_comparison_20260922_v2'


def freeze(out):
    if out.exists(): raise ValueError('fresh directory required')
    wide.execute(SOURCE,verify=True)
    source=json.loads((SOURCE/'configuration.freeze.json').read_text())
    old=json.loads(previous.OLD.read_text())
    stages={}
    inputs=source['inputs']+old['inputs']+[record(PROTOCOL),record(previous.OLD),
        old['cohort_manifest'],old['cohort_configuration'],record(SOURCE/'configuration.freeze.json'),
        record(SOURCE/'profiles_manifest.json')]
    for cell in ('HepG2','K562'):
        path=SOURCE/(cell+'.json');stage=json.loads(path.read_text());stages[cell]=stage
        inputs += [record(path)]+list(stage['outputs'].values())
    gamma=json.loads((GAMMA/'comparison_manifest.json').read_text())
    inputs += [record(GAMMA/'comparison_manifest.json'),gamma['outputs']['scores']]
    files=[Path(__file__),Path(previous.__file__),Path(radius_model.__file__),Path(balance.__file__),
           ROOT/'benchmarks/scripts/fixed_shape_normalization.py',ROOT/'benchmarks/scripts/global_parametric_scoring.py',
           ROOT/'src/fp_tools/tools/parametric_bias.py']
    config=dict(protocol=json.loads(PROTOCOL.read_text()),inputs=inputs,stages=stages,
        code=source['code']+old['code']+[record(p) for p in files],runtime=wide.old.parent.extension.runtime(),
        prior_scores=gamma['outputs']['scores'],training_only=True,promotion_eligible=False)
    if config['protocol']['fixed_dispersion']!=old['dispersion']: raise ValueError('dispersion changed')
    for entry in config['inputs']+config['code']: checked_path(entry,ROOT)
    config['configuration_id']=canonical_id(config);out.mkdir(parents=True)
    write_immutable(out/'configuration.freeze.json',config)
    print('Frozen radius comparison before scoring',flush=True)


def execute(out,verify=False):
    config=wide.old.parent.extension.read_document(out/'configuration.freeze.json','configuration_id')
    if config['runtime']!=wide.old.parent.extension.runtime() or config['protocol']!=json.loads(PROTOCOL.read_text()):
        raise ValueError('runtime or protocol changed')
    for entry in config['inputs']+config['code']: checked_path(entry,ROOT)
    target=out/'comparison_manifest.json'
    if target.exists():
        wide.old.parent.base.verify_stage(target,config['configuration_id'],ROOT)
        print('Verified radius training comparison',flush=True);return
    if verify or list(out.glob('*.tsv')): raise ValueError('missing or partial output preserved')
    old=json.loads(previous.OLD.read_text());protocol=config['protocol']
    sites,panels,ledger=previous.old.load_split(old,ROOT,'train')
    if not sites.TFBS_chr.isin(protocol['chromosomes']).all(): raise ValueError('nontraining row')
    cohort=json.loads(checked_path(old['cohort_manifest'],ROOT).read_text());pairs=[]
    extra={f'radius_{r}':np.full((len(sites),201),np.nan) for r in protocol['background_radii']}
    for cell,stage in config['stages'].items():
        ix=np.flatnonzero(sites.cell.eq(cell));source=pd.read_csv(stage['outputs']['sites']['path'],sep='\t')
        order=previous.align_indices(sites.iloc[ix],source)
        pair=pd.read_csv(checked_path(cohort['outputs'][cell+'_pairs'],ROOT),sep='\t')
        pair=pair.loc[pair.split.eq('train')];pairs.append(pair)
        np.testing.assert_array_equal(pair.chip_label,sites.iloc[ix].chip_label)
        with np.load(stage['outputs']['profiles']['path'],allow_pickle=False) as data:
            np.testing.assert_array_equal(data['site_hash'],wide.old.strands.site_hashes(source))
            tracks=previous.load_tracks(data,order);bam_valid=data['bam_valid'][order]
        for row,dest in enumerate(ix):
            if not bam_valid[row]: continue
            reverse=sites.iloc[dest].TFBS_strand=='-'
            observed=(tracks['plus_observed'][row]+tracks['minus_observed'][row])[500:701]
            np.testing.assert_array_equal(observed[::-1] if reverse else observed,panels['raw'][dest])
            for radius in protocol['background_radii']:
                result=radius_model.score_tracks(tracks['plus_observed'][row],tracks['minus_observed'][row],
                    tracks['plus_log_bias'][row],tracks['minus_log_bias'][row],radius,old['dispersion'][cell])
                values=result['corrected'][500:701]
                extra[f'radius_{radius}'][dest]=values[::-1] if reverse else values
        print('Scored fixed radii:',cell,flush=True)
    panels.update(extra)
    keep=previous.complete_pairs(pd.concat(pairs,ignore_index=True),
        np.logical_and.reduce([np.isfinite(v).all(axis=1) for v in panels.values()]))
    archive=pd.read_csv(config['prior_scores']['path'],sep='\t',float_precision='round_trip')
    archive=archive.loc[archive.arm.eq('gamma_fixed_dispersion')].reset_index(drop=True)
    old_scores=archive.iloc[previous.align_indices(sites,archive)].score.to_numpy()
    weights=previous.old.candidates.geometry_weights(protocol['detector_geometry']);x=np.arange(-100,101)
    rows=[];curves=[];scores=[];parity=[]
    for task in ledger.itertuples(index=False):
        selected=(sites.cell.eq(task.cell)&sites.tf.eq(task.tf)).to_numpy() & keep
        ix=np.flatnonzero(selected);y=sites.iloc[ix].chip_label.to_numpy(int)
        if len(ix):
            covariates=np.column_stack([sites.iloc[ix].motif_score.to_numpy(float),np.log1p(panels['raw'][ix].sum(axis=1))])
            smd=balance.smd(covariates,y)
        else: smd=[np.nan,np.nan]
        for arm,values in panels.items():
            row=dict(cell=task.cell,tf=task.tf,arm=arm,original_pairs=task.matched_pairs,
                retained_pairs=len(ix)//2,motif_smd=smd[0],count_smd=smd[1],
                status='underpowered' if len(ix)<400 else ('unbalanced' if np.max(np.abs(smd))>.1 else 'training_diagnostic'))
            if len(np.unique(y))==2:
                score=values[ix]@weights;score[panels['raw'][ix].sum(axis=1)==0]=0
                norm=normalize_fixed(values[ix],x,1.);delta=norm[y==1].mean(axis=0)-norm[y==0].mean(axis=0)
                left=delta[(x>=-40)&(x<=-15)].mean();right=delta[(x>=15)&(x<=40)].mean()
                row.update(auroc=roc_auc_score(y,score),auprc=average_precision_score(y,score),
                    center_contrast=delta[abs(x)<=5].mean(),left_shoulder=left,right_shoulder=right,
                    asymmetry=(right-left)/(abs(right)+abs(left)) if abs(right)+abs(left)>0 else 0,
                    canonical_protection=previous.old.candidates.protection_direction(delta),
                    functional_separation=standardized_functional_separation(norm,y,x))
                if arm=='radius_100':
                    np.testing.assert_allclose(score,old_scores[ix],rtol=1e-9,atol=1e-9)
                    parity.append(dict(cell=task.cell,tf=task.tf,sites=len(ix),max_absolute_score_difference=float(np.max(abs(score-old_scores[ix]))),
                        auroc_difference=row['auroc']-roc_auc_score(y,old_scores[ix]),
                        auprc_difference=row['auprc']-average_precision_score(y,old_scores[ix])))
                curves.append(pd.DataFrame(dict(cell=task.cell,tf=task.tf,arm=arm,position=x,
                    positive_mean=norm[y==1].mean(axis=0),negative_mean=norm[y==0].mean(axis=0))))
                frame=sites.iloc[ix][previous.KEYS+['chip_label']].copy();frame['arm'],frame['score']=arm,score;scores.append(frame)
            rows.append(row)
    if len(rows)!=126: raise ValueError('task denominator changed')
    outputs={}
    for name,table in [('metrics',pd.DataFrame(rows)),('curves',pd.concat(curves)),('scores',pd.concat(scores)),('radius100_parity',pd.DataFrame(parity))]:
        path=out/(name+'.tsv');table.to_csv(path,sep='\t',index=False);outputs[name]=record(path)
    result=dict(configuration_id=config['configuration_id'],outputs=outputs,training_only=True,
        promotion_eligible=False,uncertainty_calibrated=False,candidate_selected=False)
    result['stage_id']=canonical_id(result);write_immutable(target,result)
    print('Completed all 126 radius training rows',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--outdir',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--freeze-only',action='store_true');mode.add_argument('--verify-only',action='store_true')
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'): parser.error('use research results')
    if args.freeze_only: freeze(out)
    else:
        with wide.old.parent.extension.exclusive_run(out): execute(out,args.verify_only)
