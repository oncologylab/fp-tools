"""Prespecified fixed-cohort diagnostic after the successful control screen.

Native production DWM remains explicitly unqualified. This stage cannot promote
a method or substitute its sequence-propensity comparison for native TOBIAS.
"""
import argparse
import gzip
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy.special import logsumexp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from access import deny,digest,dump
from chromatin_metrics import ARMS,fold_metrics,screen_pairs,arm_means
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'benchmarks/scripts'))
from fp_tools.tools.parametric_bias import ConditionalSequenceBiasModel,strand_log_bias,calibrated_residuals
from fp_tools.tools.functional_footprints import standardized_functional_separation
from conditional_protection_score import depletion_score
from fixed_shape_normalization import normalize_fixed
import pysam

LIBRARIES=[f'{c}_rep{r}' for c in ['HepG2','K562'] for r in [1,2,3]]
DEPTHS=['full','10M'];X=np.arange(-100,101)


def freeze(out):
    decision=json.loads((out/'control.decision.json').read_text())
    assert all(decision['gates'][k] for k in ['likelihood','mean','factorial_controls'])
    sub=out/'chromatin';sub.mkdir(exist_ok=False)
    for name in ['tables','arrays','figures','logs','provenance']:(sub/name).mkdir()
    prior=ROOT.parent/'fp-tools-next-round/benchmarks/results/next_round_20260922_v1'
    old=ROOT.parent/'fp-tools-calibration/benchmarks/results/calibration_round_20261019_v1'
    registry=json.loads((old/'inputs.json').read_text())['permitted'];records={}
    keys=['tables/nested_sites.tsv.gz','tables/task_ledger.tsv','tables/covariate_baseline_all_libraries.tsv.gz']
    keys += [f'tables/counts.{lib}.{depth}.tsv.gz' for lib in LIBRARIES for depth in DEPTHS]
    for key in keys:
        r=registry[key];p=Path(r['path']);deny(p,r['identity'])
        assert str(p.resolve())==r['resolved_path'] and digest(p)==r['sha256'];records[key]=r
    extra={'baseline_scores':old/'tables/chromatin_site_scores.tsv.gz',
           'native_qualification':prior/'provenance/dwm_reference_qualification.v3.json',
           'native_layout':prior/'provenance/dwm_layout_isolation_result.json'}
    for lib in LIBRARIES:
        extra[f'arrays/{lib}.propensities.npz']=prior/f'arrays/{lib}.propensities.npz'
        for depth in DEPTHS:extra[f'arrays/{lib}.{depth}.npz']=prior/f'arrays/{lib}.{depth}.npz'
    for key,p in extra.items():
        identity='permitted chr1-15 development derived data: '+key;deny(p,identity)
        records[key]=dict(path=str(p),resolved_path=str(p.resolve()),identity=identity,sha256=digest(p))
    genome=json.loads((out/'inputs.local.json').read_text())['genome'];records['genome']=genome
    dump(sub/'inputs.local.json',records)
    native=json.loads(extra['native_qualification'].read_text());layout=json.loads(extra['native_layout'].read_text())
    dump(sub/'provenance/native_reference_audit.json',dict(latest_qualification=native,layout_isolation=layout,
        decision='INCOMPLETE_NATIVE_STRAND_AND_COMPLETE_PAIR_CONTRACT',
        reuse='do not rerun checksum-valid support stages; 4200 processed, 180 valid zero expectations, 2 layout-dependent shapes',
        comparison_scope='sequence-propensity factorial only; archived native residual available only in historical 8106-site first10M replay'))
    science=[Path(__file__),Path(__file__).with_name('chromatin_metrics.py'),ROOT/'benchmarks/scripts/fixed_shape_normalization.py',
             ROOT/'src/fp_tools/tools/functional_footprints.py',ROOT/'src/fp_tools/tools/parametric_bias.py',ROOT/'benchmarks/scripts/conditional_protection_score.py']
    dump(sub/'EXPERIMENT_PLAN.freeze.json',dict(created_unix=time.time(),
        authorization='conditional Phase4 of plan committed before training; control component passed; native comparator unqualified so diagnostic only',
        policy=dict(cohort='unchanged 6020 nested training sites; no rematching; all21 tasks',depths=DEPTHS,libraries=LIBRARIES,
            primary_l2=.01,ablation_l2=.0001,geometry=[7,12,40],seeds=list(range(2026,2031)),threshold=1.3010299957,
            metrics='unchanged calibration equal-chromosome/library/task functions; 2000 synchronized chromosome draws seed2026',
            minimum_positive_and_negative_sites=200,shape='unchanged outside-flank normalization fallback1; fixed deviance dispersion HepG2=0 K562=1.6144318781609228; original RMS standardized separation ±50',
            shape_geometry=[[5,15,40],[7,12,40]],native_DWM_qualified=False,protected_contents_allowed=False,
            expansion='only original full-depth accuracy/shape/positive-control/stability gates plus qualified native reference; otherwise stop'),
        code={str(p.relative_to(ROOT)):digest(p) for p in science},
        inputs={k:dict(identity=v['identity'],sha256=v['sha256']) for k,v in records.items()},
        ensembles={p.name:digest(p) for p in sorted((out/'models').glob('*.ensemble.npz'))},
        control_decision_sha256=digest(out/'control.decision.json')))


def run(out):
    sub=out/'chromatin';frozen=json.loads((sub/'EXPERIMENT_PLAN.freeze.json').read_text())
    for p,h in frozen['code'].items():assert digest(ROOT/p)==h,p
    for p,h in frozen['ensembles'].items():assert digest(out/'models'/p)==h,p
    registry=json.loads((sub/'inputs.local.json').read_text())
    def access(key):
        r=registry[key];p=Path(r['path']);deny(p,r['identity'])
        assert r['sha256']==frozen['inputs'][key]['sha256'] and r['identity']==frozen['inputs'][key]['identity']
        assert str(p.resolve())==r['resolved_path'] and digest(p)==r['sha256']
        with (sub/'logs/access.jsonl').open('a') as log:log.write(json.dumps(dict(input_id=key,identity=r['identity'],sha256=r['sha256']))+'\n')
        return p
    read=lambda key:pd.read_csv(access(key),sep='\t',float_precision='round_trip',low_memory=False)
    sites=read('tables/nested_sites.tsv.gz')
    assert len(sites)==6020 and sites.site_id.is_unique and sites.TFBS_chr.isin([f'chr{i}' for i in range(1,16)]).all()
    tasks=read('tables/task_ledger.tsv');cov=read('tables/covariate_baseline_all_libraries.tsv.gz');old=read('baseline_scores')
    sites.to_csv(sub/'tables/site_ledger.tsv.gz',sep='\t',index=False)
    cov.to_csv(sub/'tables/covariates.tsv.gz',sep='\t',index=False)
    tasks.to_csv(sub/'tables/task_ledger.tsv',sep='\t',index=False)
    logits={};sequences=[]
    with pysam.FastaFile(str(access('genome'))) as fasta:
        for row in sites.itertuples():sequences.append(fasta.fetch(row.TFBS_chr,int(row.center)-141,int(row.center)+142).upper())
    np.savez_compressed(sub/'arrays/sequence_contexts.npz',sequence=np.array(sequences),site_id=sites.site_id.to_numpy(str),
        chromosome=sites.TFBS_chr.to_numpy(str),motif_strand=sites.TFBS_strand.to_numpy(str))
    for arm,l2 in [('nuclear_primary',.01),('nuclear_ablation',.0001)]:
        values=np.full((len(sites),2,201),np.nan)
        for chrom in [f'chr{i}' for i in range(1,16)]:
            model=ConditionalSequenceBiasModel.load(out/f'models/{chrom}.l2_{l2}.ensemble.npz')
            for i in np.flatnonzero(sites.TFBS_chr.eq(chrom)):
                p,m,valid=strand_log_bias(model,sequences[i],np.arange(41,242))
                if sites.iloc[i].TFBS_strand=='-':p,m=m[::-1],p[::-1]
                values[i]=[p,m]
        logits[arm]=values
        np.savez_compressed(sub/f'arrays/{arm}.logits.npz',site_id=sites.site_id.to_numpy(str),logits=values)
    frames=[];shape_rows=[];regional_rows=[];curve_frames=[];parity=[]
    center=abs(X)<=7;domain=center|((abs(X)>=12)&(abs(X)<=40))
    for library in LIBRARIES:
        cell=library.split('_')[0];global_index=np.flatnonzero(sites.cell.eq(cell));local=sites.iloc[global_index].reset_index(drop=True)
        with np.load(access(f'arrays/{library}.propensities.npz'),allow_pickle=False) as z:
            ix=pd.Index(z['site_id']).get_indexer(local.site_id);assert (ix>=0).all()
            bias={a:np.stack([z[a+'_'+s][ix] for s in ['plus','minus']],axis=1) for a in ['LOG21','DWM']}
        np.savez_compressed(sub/f'arrays/{library}.reference_logits.npz',site_id=local.site_id.to_numpy(str),**bias)
        for depth in DEPTHS:
            counts=read(f'tables/counts.{library}.{depth}.tsv.gz').set_index('site_id').loc[local.site_id].reset_index()
            assert counts.cohort_nested.all()
            for field in ['TFBS_chr','tf','cell','chip_label']:np.testing.assert_array_equal(counts[field],local[field])
            with np.load(access(f'arrays/{library}.{depth}.npz'),allow_pickle=False) as z:
                ix=pd.Index(z['site_id']).get_indexer(local.site_id);assert (ix>=0).all()
                observed=np.stack([z[s+'_observed'][ix] for s in ['plus','minus']],axis=1)
                achieved=int(z['achieved_fragments'])
            np.savez_compressed(sub/f'arrays/{library}.{depth}.npz',site_id=local.site_id.to_numpy(str),observed=observed,achieved_fragments=achieved)
            for si,s in enumerate(['plus','minus']):
                np.testing.assert_array_equal(observed[:,si,center].sum(1),counts[s+'_center'])
                np.testing.assert_array_equal(observed[:,si,domain].sum(1),counts[s+'_total'])
            baseline=old[old.library.eq(library)&old.depth.eq(depth)&old.arm.isin(ARMS)]
            assert set(baseline.site_id)==set(local.site_id)
            for arm in ['uniform_count_tail','LOG21_count_tail','DWM_count_tail','LOG21NC1_count_tail']:
                b=baseline[baseline.arm.eq(arm)].set_index('site_id').loc[local.site_id]
                for field in ['plus_center','minus_center','plus_total','minus_total','chip_label']:
                    np.testing.assert_array_equal(counts[field],b[field])
            frames.append(baseline.copy())
            for arm in ['nuclear_primary','nuclear_ablation']:
                log=logits[arm][global_index]
                probability=np.exp(logsumexp(log[:,:,center],axis=2)-logsumexp(log[:,:,domain],axis=2))
                valid=np.isfinite(log[:,:,domain]).all(axis=(1,2))
                probability[~valid]=np.nan
                result=counts.copy();result['arm']=arm;result['q_plus']=probability[:,0];result['q_minus']=probability[:,1]
                result['pair_id']=local.pair_id.to_numpy();result['fitting_provenance']='nakedDNArep1; held TFBS_chr excluded; fixed five-seed ensemble'
                result['expected_center']=counts.plus_total*probability[:,0]+counts.minus_total*probability[:,1]
                values=[]
                for row,(p,m) in zip(counts.itertuples(),probability):
                    if not np.isfinite([p,m]).all():values.append((np.nan,np.nan));continue
                    values.append(tuple(round(depletion_score(int(row.plus_total),int(row.minus_total),int(c),p,m),10)
                        for c in [row.plus_center+row.minus_center,0]))
                result['score']=[v[0] for v in values];result['ceiling']=[v[1] for v in values]
                result['status']=np.select([~valid,result.total.eq(0),result.ceiling.le(1.3010299957),result.score.gt(1.3010299957)],
                    ['UNAVAILABLE','NO_INFORMATION','ABSTAIN_UNREACHABLE','CALL'],default='NO_CALL')
                result['call']=result.status.eq('CALL');frames.append(result)
            parity.append(dict(library=library,depth=depth,sites=len(local),regional_count_mismatches=0,baseline_count_mismatches=0,achieved_fragments=achieved))
            # Original shape operator: full201 strand-normalized expectations,
            # historical fixed deviance variance, then unchanged flank normalization.
            for arm,log in [('nuclear_primary',logits['nuclear_primary'][global_index]),('nuclear_ablation',logits['nuclear_ablation'][global_index]),
                            ('DWM_propensity',bias['DWM']),('LOG21_propensity',bias['LOG21'])]:
                expected=np.exp(log-logsumexp(log,axis=2,keepdims=True))*observed.sum(2,keepdims=True)
                phi=0. if cell=='HepG2' else 1.6144318781609228
                residual=calibrated_residuals(observed.sum(1),expected.sum(1),'deviance',dispersion=phi)
                for task in tasks[tasks.cell.eq(cell)].itertuples():
                    mask=local.tf.eq(task.tf).to_numpy();label=local.loc[mask,'chip_label'].to_numpy(int)
                    if len(label)<4 or len(np.unique(label))!=2:
                        shape_rows.append(dict(library=library,depth=depth,cell=cell,tf=task.tf,arm=arm,status='UNAVAILABLE',sites=len(label)));continue
                    data=residual[mask]
                    if not np.isfinite(data).all() or not np.isfinite(log[mask]).all():
                        shape_rows.append(dict(library=library,depth=depth,cell=cell,tf=task.tf,arm=arm,status='UNSUPPORTED_SEQUENCE',sites=len(label)));continue
                    normal=normalize_fixed(data,X,1.);delta=normal[label==1].mean(0)-normal[label==0].mean(0)
                    c=delta[abs(X)<=5].mean();left=delta[(X>=-40)&(X<=-15)].mean();right=delta[(X>=15)&(X<=40)].mean()
                    shape_rows.append(dict(library=library,depth=depth,cell=cell,tf=task.tf,arm=arm,status='AVAILABLE',sites=len(label),positive_sites=int(label.sum()),
                        center_contrast=c,left_minus_center=left-c,right_minus_center=right-c,canonical_protection=bool(c<0 and left>c and right>c),
                        standardized_functional_separation=standardized_functional_separation(normal,label,X),
                        scope='same-operator propensity/deviance diagnostic; not native production DWM'))
                    obs=observed[mask];exp=expected[mask];coverage=obs.sum(axis=(1,2))
                    for representation,panel in [('raw',obs),('expected',exp),('residual',None)]:
                        for strand in ['plus','minus','combined']:
                            if representation=='residual':
                                panel_values=data if strand=='combined' else calibrated_residuals(obs[:,0 if strand=='plus' else 1],exp[:,0 if strand=='plus' else 1],'deviance',dispersion=phi)
                            else:panel_values=panel.sum(1) if strand=='combined' else panel[:,0 if strand=='plus' else 1]
                            for label_value in [0,1]:
                                chosen=label==label_value
                                for aggregation in ['site_equal','count_weighted']:
                                    weights=np.ones(chosen.sum()) if aggregation=='site_equal' else coverage[chosen]
                                    mean=np.average(panel_values[chosen],axis=0,weights=weights) if weights.sum()>0 else np.full(201,np.nan)
                                    curve_frames.append(pd.DataFrame(dict(library=library,depth=depth,cell=cell,tf=task.tf,arm=arm,
                                        representation=representation,strand=strand,label=label_value,aggregation=aggregation,
                                        position=X,value=mean,sites=int(chosen.sum()),units='signed deviance residual' if representation=='residual' else 'cuts per site per base')))
                                    for geom,cw,inner,outer in [('canonical',5,15,40),('scoring',7,12,40)]:
                                        mc=mean[abs(X)<=cw].mean();ml=mean[(X>=-outer)&(X<=-inner)].mean();mr=mean[(X>=inner)&(X<=outer)].mean()
                                        regional_rows.append(dict(library=library,depth=depth,cell=cell,tf=task.tf,arm=arm,representation=representation,strand=strand,
                                            label=label_value,aggregation=aggregation,geometry=geom,center=mc,left=ml,right=mr,left_minus_center=ml-mc,right_minus_center=mr-mc,asymmetry=mr-ml))
            print('Fixed-cohort scoring and shapes complete',library,depth,flush=True)
    scored=pd.concat(frames,ignore_index=True);scored.to_csv(sub/'tables/per_site_scores.tsv.gz',sep='\t',index=False)
    pd.DataFrame(parity).to_csv(sub/'tables/observation_parity.tsv',sep='\t',index=False)
    pd.DataFrame(shape_rows).to_csv(sub/'tables/shape_metrics.tsv',sep='\t',index=False)
    pd.DataFrame(regional_rows).to_csv(sub/'tables/regional_shapes.tsv',sep='\t',index=False)
    pd.concat(curve_frames,ignore_index=True).to_csv(sub/'tables/shape_plot_source.tsv.gz',sep='\t',index=False)
    metrics=fold_metrics(scored,cov,tasks);metrics.to_csv(sub/'tables/fold_metrics.tsv',sep='\t',index=False)
    arm_means(metrics,tasks).to_csv(sub/'tables/task_metrics.tsv',sep='\t',index=False)
    deltas=[];summary=[]
    for depth in DEPTHS:
        for method in ['nuclear_primary','nuclear_ablation']:
            for reference in ['DWM_count_tail','LOG21_count_tail','LOG21NC1_count_tail']:
                d,s=screen_pairs(metrics,tasks,depth,method,reference);deltas.extend(d);summary.append(s)
    pd.DataFrame(deltas).to_csv(sub/'tables/paired_deltas.tsv',sep='\t',index=False)
    pd.DataFrame(summary).to_csv(sub/'tables/screen_summary.tsv',sep='\t',index=False)
    dump(sub/'complete.json',dict(scored_rows=len(scored),fold_metric_rows=len(metrics),sites=6020,tasks=21,native_comparator_qualified=False,
        scope='fixed pre-existing training cohort; repeatedly examined development; no locked labels',scientific_promotion=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','run']);p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args();globals()[a.stage](a.outdir.resolve())
