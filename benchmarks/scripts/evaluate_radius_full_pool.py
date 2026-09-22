#!/usr/bin/env python3
"""Recompute seven fixed arms and shape on three-covariate training cohorts."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score
sys.path.insert(0,str(Path(__file__).resolve().parent))
import prepare_radius_full_pool as cohorts
from fixed_shape_normalization import normalize_fixed
from fp_tools.tools.functional_footprints import standardized_functional_separation
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import checked_path,canonical_id

ROOT=cohorts.ROOT
SOURCE=ROOT/'benchmarks/results/radius_full_pool_20260922_v1'
ARMS=['raw','DWM_expected','parametric_direct','radius_100','radius_200','radius_500','uniform_radius_500']
PROTOCOL=dict(evidence='full_pool_three_covariate_training_comparison',arms=ARMS,
    geometry=[7,12,40],normalization_fallback=1.,dispersion={'HepG2':0.,'K562':1.6144318781609228},
    support='exact frozen complete pairs; abort on support loss',task_denominator=21,
    refitting=False,candidate_selection=False,uncertainty_calibrated=False,promotion_eligible=False)


def evaluate_panels(sites,panels,diagnostics):
    if set(panels)!=set(ARMS) or any(v.shape!=(len(sites),201) or not np.isfinite(v).all() for v in panels.values()):
        raise ValueError('seven complete finite aligned panels required')
    rows=[];curves=[];scores=[];x=np.arange(-100,101)
    weights=cohorts.prior.previous.old.candidates.geometry_weights(PROTOCOL['geometry'])
    for task in diagnostics.itertuples(index=False):
        ix=np.flatnonzero((sites.cell.eq(task.cell)&sites.tf.eq(task.tf)).to_numpy())
        y=sites.iloc[ix].chip_label.to_numpy(int)
        if len(ix)!=2*task.matched_pairs or y.sum()!=task.matched_pairs: raise ValueError('pair membership changed')
        for arm in ARMS:
            row=dict(cell=task.cell,tf=task.tf,arm=arm,matched_pairs=task.matched_pairs,status=task.status)
            if len(ix):
                values=panels[arm][ix];score=values@weights;score[panels['raw'][ix].sum(axis=1)==0]=0
                norm=normalize_fixed(values,x,1.);delta=norm[y==1].mean(0)-norm[y==0].mean(0)
                left=delta[(x>=-40)&(x<=-15)].mean();right=delta[(x>=15)&(x<=40)].mean()
                row.update(auroc=roc_auc_score(y,score),auprc=average_precision_score(y,score),
                    center_contrast=delta[abs(x)<=5].mean(),left_shoulder=left,right_shoulder=right,
                    asymmetry=(right-left)/(abs(left)+abs(right)) if abs(left)+abs(right)>0 else 0,
                    canonical_protection=cohorts.prior.previous.old.candidates.protection_direction(delta),
                    functional_separation=standardized_functional_separation(norm,y,x))
                curves.append(pd.DataFrame(dict(cell=task.cell,tf=task.tf,arm=arm,position=x,
                    positive_mean=norm[y==1].mean(0),negative_mean=norm[y==0].mean(0))))
                frame=sites.iloc[ix][cohorts.prior.previous.KEYS+['chip_label']].copy();frame['arm'],frame['score']=arm,score;scores.append(frame)
            rows.append(row)
    if len(rows)!=147: raise ValueError('full task denominator changed')
    return pd.DataFrame(rows),pd.concat(curves,ignore_index=True),pd.concat(scores,ignore_index=True)


def freeze(out):
    if out.exists(): raise ValueError('fresh directory required')
    cohorts.execute(SOURCE,verify=True)
    config=json.loads((SOURCE/'configuration.freeze.json').read_text());result=json.loads((SOURCE/'cohorts_manifest.json').read_text())
    paths={entry['path'] for entry in config['code']}
    # Pin imported research helpers as well as the direct scientific code.
    for module in list(sys.modules.values()):
        path=getattr(module,'__file__',None)
        if path:
            path=Path(path).resolve()
            if path.is_relative_to(ROOT) and path.suffix=='.py': paths.add(str(path))
    paths.add(str(Path(__file__).resolve()))
    value=dict(protocol=PROTOCOL,inputs=config['inputs']+[record(SOURCE/'configuration.freeze.json'),
        record(SOURCE/'cohorts_manifest.json')]+list(result['outputs'].values()),
        code=[record(Path(p)) for p in sorted(paths)],cohorts=result['outputs'],source_configuration=record(SOURCE/'configuration.freeze.json'),
        runtime=cohorts.prior.wide.old.parent.extension.runtime(),training_only=True,promotion_eligible=False)
    for entry in value['inputs']+value['code']: checked_path(entry,ROOT)
    value['configuration_id']=canonical_id(value);out.mkdir(parents=True);write_immutable(out/'configuration.freeze.json',value)
    print('Frozen exact-cohort score and shape comparison',flush=True)


def execute(out,verify=False):
    config=cohorts.prior.wide.old.parent.extension.read_document(out/'configuration.freeze.json','configuration_id')
    if config['runtime']!=cohorts.prior.wide.old.parent.extension.runtime() or config['protocol']!=PROTOCOL:
        raise ValueError('runtime or protocol changed')
    for entry in config['inputs']+config['code']: checked_path(entry,ROOT)
    target=out/'comparison_manifest.json'
    if target.exists():
        cohorts.prior.wide.old.parent.base.verify_stage(target,config['configuration_id'],ROOT)
        print('Verified exact-cohort full-pool comparison',flush=True);return
    if verify or list(out.glob('*.tsv')) or list(out.glob('*.npz')): raise ValueError('missing or partial output')
    source=json.loads(checked_path(config['source_configuration'],ROOT).read_text())
    old=json.loads(checked_path(source['old_configuration'],ROOT).read_text())
    old_cohort=json.loads(checked_path(old['cohort_manifest'],ROOT).read_text())
    if old['dispersion']!=PROTOCOL['dispersion']: raise ValueError('dispersion changed')
    sites=pd.read_csv(config['cohorts']['pairs']['path'],sep='\t',float_precision='round_trip')
    diagnostics=pd.read_csv(config['cohorts']['diagnostics']['path'],sep='\t',float_precision='round_trip')
    if not sites.TFBS_chr.isin([f'chr{i}' for i in range(1,16)]).all(): raise ValueError('nontraining row')
    panels={arm:np.full((len(sites),201),np.nan) for arm in ARMS}
    for entry in old['protocol']['profile_inputs']:
        cell=entry['cell'];ix=np.flatnonzero(sites.cell.eq(cell));selected=sites.iloc[ix]
        original,arrays=cohorts.prior.previous.old.base.source.read_profiles(entry,ROOT)
        order=selected.old_index.to_numpy(int)
        pd.testing.assert_frame_equal(original.iloc[order][cohorts.prior.previous.KEYS+['chip_label']].reset_index(drop=True),
            selected[cohorts.prior.previous.KEYS+['chip_label']].reset_index(drop=True),check_dtype=False)
        observed=arrays['observed'][order];panels['raw'][ix]=observed
        with np.load(old_cohort['outputs'][cell+'_dwm']['path'],allow_pickle=False) as archive:
            np.testing.assert_array_equal(archive['site_hash'],arrays['site_hash']);dwm=archive['expected'][order]
        for arm,expected in [('DWM_expected',dwm),('parametric_direct',arrays['expected'][order])]:
            panels[arm][ix]=cohorts.prior.previous.old.base.source.safety.direct_scores(observed,expected,old['dispersion'][cell])[2]
        stage=source['stages'][cell];order=selected.wide_index.to_numpy(int)
        with np.load(stage['outputs']['profiles']['path'],allow_pickle=False) as archive:
            np.testing.assert_array_equal(archive['site_hash'][order],cohorts.prior.wide.old.strands.site_hashes(selected))
            if not archive['bam_valid'][order].all(): raise ValueError('alignment support changed')
            tracks=cohorts.prior.previous.load_tracks(archive,order)
        zero=np.zeros(1201)
        for row,dest in enumerate(ix):
            reverse=sites.iloc[dest].TFBS_strand=='-'
            counts=(tracks['plus_observed'][row]+tracks['minus_observed'][row])[500:701]
            np.testing.assert_array_equal(counts[::-1] if reverse else counts,observed[row])
            for arm,radius in [('radius_100',100),('radius_200',200),('radius_500',500),('uniform_radius_500',500)]:
                uniform=arm=='uniform_radius_500'
                result=cohorts.prior.radius_model.score_tracks(tracks['plus_observed'][row],tracks['minus_observed'][row],
                    zero if uniform else tracks['plus_log_bias'][row],zero if uniform else tracks['minus_log_bias'][row],radius,old['dispersion'][cell])
                values=result['corrected'][500:701];panels[arm][dest]=values[::-1] if reverse else values
        print('Recomputed all seven arms:',cell,flush=True)
    metrics,curves,scores=evaluate_panels(sites,panels,diagnostics);outputs={}
    for name,table in [('metrics',metrics),('curves',curves),('scores',scores),('sites',sites)]:
        path=out/(name+'.tsv');table.to_csv(path,sep='\t',index=False);outputs[name]=record(path)
    path=out/'panels.npz';np.savez_compressed(path,**panels);outputs['panels']=record(path)
    result=dict(configuration_id=config['configuration_id'],outputs=outputs,training_only=True,
        candidate_selected=False,uncertainty_calibrated=False,promotion_eligible=False)
    result['stage_id']=canonical_id(result);write_immutable(target,result)
    print('Completed all 147 score and shape rows',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--outdir',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--freeze-only',action='store_true');mode.add_argument('--verify-only',action='store_true')
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'): parser.error('use research results')
    if args.freeze_only: freeze(out)
    else:
        with cohorts.prior.wide.old.parent.extension.exclusive_run(out): execute(out,args.verify_only)
