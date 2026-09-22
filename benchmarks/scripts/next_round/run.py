"""Execute frozen counts, cross-fitted nuisance fits and signed-shrinkage arms."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime,timezone
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.special import logsumexp

from access import digest
from prepare import ROOT,dump
from shrinkage import posterior,fit_dispersion,threshold
sys.path.insert(0,str(ROOT/'benchmarks/scripts'))
from conditional_protection_score import depletion_score

X=np.arange(-100,101)
REGIONS={'center':abs(X)<=7,'left':(X>=-40)&(X<=-12),'right':(X>=12)&(X<=40)}
UNION=np.logical_or.reduce(list(REGIONS.values()))
ARMS=['uniform_count_tail','LOG21_count_tail','LOG21_binomial','LOG21_beta_binomial','DWM_binomial','DWM_beta_binomial']


def components(sites,arrays,propensities):
    frame=sites.copy().reset_index(drop=True)
    for strand in ['plus','minus']:
        observed=arrays[strand+'_observed']
        if observed.shape!=(len(frame),201) or not np.isfinite(observed).all() or (observed<0).any() or (observed!=np.floor(observed)).any():raise ValueError('invalid observed strand array')
        for region,mask in REGIONS.items():frame[strand+'_'+region]=observed[:,mask].sum(1).astype(int)
        frame[strand+'_total']=observed[:,UNION].sum(1).astype(int)
        for arm in ['LOG21','DWM']:
            logits=propensities[arm+'_'+strand]
            if logits.shape!=observed.shape:raise ValueError('propensity length mismatch')
            frame[arm+'_'+strand+'_valid']=np.isfinite(logits[:,UNION]).all(1)
            frame[arm+'_'+strand+'_p']=np.exp(logsumexp(logits[:,REGIONS['center']],axis=1)-logsumexp(logits[:,UNION],axis=1))
        frame['uniform_'+strand+'_p']=15/73
    frame['total']=frame.plus_total+frame.minus_total
    frame['LOG21_expected_center']=frame.plus_total*frame.LOG21_plus_p+frame.minus_total*frame.LOG21_minus_p
    frame['count_stratum']=pd.cut(frame.total,[-1,0,1,4,9,19,49,99,np.inf],labels=['0','1','2-4','5-9','10-19','20-49','50-99','100+']).astype(str)
    frame['expected_stratum']=pd.cut(frame.LOG21_expected_center,[-np.inf,.1,.5,1,3,np.inf],right=False,labels=['<0.1','0.1-<0.5','0.5-<1','1-<3','3+']).astype(str)
    frame['no_cuts']=frame.total.eq(0)
    return frame


def posterior_job(values):
    try:
        return posterior(*values)
    except (ValueError,ArithmeticError) as error:
        return {'mean':np.nan,'median':np.nan,'sd':np.nan,'low':np.nan,'high':np.nan,
                'probability_positive':np.nan,'no_information':False,'prior_only':False,
                'variance_reduction':np.nan,'integration_error':np.nan,'failure':str(error)}


def effect_rows(frame,arm,rho,pool):
    args=[(int(r.plus_center),int(r.minus_center),int(r.plus_total),int(r.minus_total),
           float(getattr(r,arm+'_plus_p')),float(getattr(r,arm+'_minus_p')),float(rho)) for r in frame.itertuples()]
    return pd.DataFrame(list(pool.map(posterior_job,args,chunksize=32)))


def control_counts(out,depth):
    sites=pd.read_csv(out/'tables/control_sites.tsv.gz',sep='\t')
    if depth=='archived_full':
        arrays={s+'_observed':np.concatenate([np.load(out/f'arrays/{c}.control_full.npz',allow_pickle=False)[s+'_observed'] for c in ['HepG2','K562']]) for s in ['plus','minus']}
    else:
        with np.load(out/f'arrays/NakedDNA_rep1.{depth}.npz',allow_pickle=False) as z:arrays={k:z[k] for k in ['plus_observed','minus_observed']}
    with np.load(out/'arrays/NakedDNA_rep1.propensities.npz',allow_pickle=False) as z:props={k:z[k] for k in z.files}
    np.testing.assert_array_equal(sites.site_id.to_numpy(str),props['site_id'])
    result=components(sites,arrays,props)
    old=pd.read_csv(out/'tables/archived_control_regions.tsv.gz',sep='\t')
    key=['cell','tf','TFBS_chr','TFBS_start','TFBS_end','TFBS_strand']
    eligible=old.groupby(key,sort=False).eligible.all().reset_index()
    result=result.merge(eligible,on=key,how='left',validate='one_to_one')
    if result.eligible.isna().any():raise ValueError('control exclusion provenance missing')
    for strand in ['plus','minus']:
        previous=old[old.arm.eq('parametric')&old.region.eq('center')&old.strand.eq(strand)]
        joined=result.merge(previous[key+['probability','observed','total']],on=key,validate='one_to_one',suffixes=('','_old'))
        np.testing.assert_allclose(joined['LOG21_'+strand+'_p'],joined.probability,atol=1e-10,rtol=1e-10)
        if depth=='archived_full':
            np.testing.assert_array_equal(joined[strand+'_center'],joined.observed)
            np.testing.assert_array_equal(joined[strand+'_total'],joined.total_old)
    result['depth']=depth;result['library']='NakedDNA_rep1';result['seed']=2026
    result['physical_id']=result.TFBS_chr+':'+result.center.astype(str)
    result['deduplicated_fit_representative']=~result.physical_id.duplicated()
    # Duplicates share one control library, even when entered under two cells.
    result.to_csv(out/f'tables/control_counts.{depth}.tsv.gz',sep='\t',index=False)
    return result


def calibrate(out,depth,workers):
    model_path=out/f'models/calibration.{depth}.json'
    if model_path.exists():raise ValueError('calibration already exists')
    controls=control_counts(out,depth)
    records=[];scored=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for chrom in [f'chr{i}' for i in range(1,16)]:
            for arm in ['LOG21','DWM']:
                valid=controls.eligible & controls[arm+'_plus_valid'] & controls[arm+'_minus_valid']
                train=controls[valid & controls.deduplicated_fit_representative & controls.TFBS_chr.ne(chrom)]
                held=controls[valid & controls.TFBS_chr.eq(chrom)]
                dispersion=fit_dispersion(train[['plus_center','minus_center']],train[['plus_total','minus_total']],train[[arm+'_plus_p',arm+'_minus_p']])
                for method,rho in [(arm+'_binomial',0.),(arm+'_beta_binomial',dispersion['rho'])]:
                    identified=method.endswith('_binomial') and not method.endswith('_beta_binomial') or dispersion['identified']
                    rec=dict(fold=chrom,method=method,control_depth=depth,rho=rho,dispersion=dispersion,
                        fitting_unique_sites=len(train),fitting_informative_sites=int((train.total>0).sum()),
                        status='AVAILABLE' if identified else 'UNIDENTIFIED',threshold=None)
                    if identified:
                        informative=train[train.total>0]
                        fitted=effect_rows(informative,arm,rho,pool)
                        if 'failure' in fitted and fitted.failure.notna().any():raise ArithmeticError('calibration posterior failed')
                        cutoff=threshold(fitted['mean'])
                        rec['threshold']=float(cutoff) if np.isfinite(cutoff) else None
                        result=effect_rows(held,arm,rho,pool)
                        frame=held[['site_id','cell','tf','TFBS_chr','physical_id','total','eligible','deduplicated_fit_representative','count_stratum','expected_stratum']].reset_index(drop=True)
                        frame=pd.concat([frame,result],axis=1)
                        frame['method']=method;frame['fold']=chrom;frame['threshold']=cutoff;frame['rho']=rho
                        frame['call_eligible']=np.isfinite(cutoff)&~frame.prior_only&np.isfinite(frame['mean'])
                        frame['call']=frame.call_eligible & (frame['mean']>cutoff)
                        frame['depth']=depth;scored.append(frame)
                    records.append(rec)
            print('Calibrated held control fold',depth,chrom,flush=True)
    dump(model_path,dict(records=records,working_likelihood=True,final_control=False,
        code_sha256=digest(__file__),input_counts_sha256=digest(out/f'tables/control_counts.{depth}.tsv.gz')))
    pd.concat(scored,ignore_index=True).to_csv(out/f'tables/control_scores.{depth}.tsv.gz',sep='\t',index=False)


def score(out,library,depth,calibration_depth,workers):
    destination=out/f'tables/scores.{library}.{depth}.tsv.gz'
    if destination.exists():raise ValueError('scores already exist')
    cell=library.split('_')[0]
    poolsites=pd.read_csv(out/'tables/pool.tsv.gz',sep='\t');sites=poolsites[poolsites.cell.eq(cell)].reset_index(drop=True)
    nested=pd.read_csv(out/'tables/nested_sites.tsv.gz',sep='\t')
    replay=pd.read_csv(out/'tables/replay_sites.tsv.gz',sep='\t')
    select=sites.site_id.isin(set(nested.site_id)|set(replay.site_id)).to_numpy()
    path=out/(f'arrays/{cell}.replay_pool.npz' if depth=='archived_10M' else f'arrays/{library}.{depth}.npz')
    with np.load(path,allow_pickle=False) as z:
        np.testing.assert_array_equal(z['site_id'],sites.site_id.to_numpy(str))
        arrays={k:z[k][select] for k in ['plus_observed','minus_observed']}
    with np.load(out/f'arrays/{library}.propensities.npz',allow_pickle=False) as z:
        np.testing.assert_array_equal(z['site_id'],sites.site_id.to_numpy(str))
        props={k:z[k][select] for k in ['LOG21_plus','LOG21_minus','DWM_plus','DWM_minus']}
    frame=components(sites.loc[select],arrays,props)
    frame['library']=library;frame['depth']=depth;frame['seed']=2026
    frame['cohort_nested']=frame.site_id.isin(nested.site_id)
    frame['cohort_replay']=frame.site_id.isin(replay.site_id)
    frame.to_csv(out/f'tables/counts.{library}.{depth}.tsv.gz',sep='\t',index=False)
    calibration=json.loads((out/f'models/calibration.{calibration_depth}.json').read_text())
    lookup={(r['method'],r['fold']):r for r in calibration['records']}
    scores=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for chrom,held in frame.groupby('TFBS_chr',sort=False):
            held=held.reset_index(drop=True)
            for method in ARMS:
                arm=method.split('_')[0]
                result=held[['site_id','cell','tf','TFBS_chr','chip_label','total','count_stratum','expected_stratum']].copy()
                result['method']=method;result['library']=library;result['depth']=depth;result['seed']=2026
                result['fold']=chrom;result['geometry']='7,12,40';result['calibration_depth']=calibration_depth
                result['status']='AVAILABLE';result['exclusion_reason']=''
                if method.endswith('count_tail'):
                    values=[];ceilings=[]
                    for r in held.itertuples():
                        pp,pm=getattr(r,arm+'_plus_p'),getattr(r,arm+'_minus_p')
                        if not np.isfinite([pp,pm]).all():values.append(np.nan);ceilings.append(np.nan);continue
                        values.append(round(depletion_score(r.plus_total,r.minus_total,r.plus_center+r.minus_center,pp,pm),10))
                        ceilings.append(round(depletion_score(r.plus_total,r.minus_total,0,pp,pm),10))
                    result['score']=values;result['attainable_old_score_ceiling']=ceilings;result['threshold']=1.3010299957
                    result['prior_only']=held.total.eq(0);result['rho']=0.
                else:
                    rec=lookup[(method,chrom)];result['rho']=rec['rho']
                    result['threshold']=rec['threshold'];result['calibration_provenance']=f'{calibration_depth}:{method}:exclude_{chrom}'
                    if rec['status']=='UNIDENTIFIED':
                        result['score']=np.nan;result['prior_only']=held.total.eq(0);result['status']='UNIDENTIFIED'
                        result['exclusion_reason']='dispersion profile not qualified'
                    else:
                        values=effect_rows(held,arm,rec['rho'],pool)
                        result=pd.concat([result,values],axis=1);result['score']=result['mean']
                        if 'failure' in values:
                            failed=values.failure.notna();result.loc[failed,'status']='NUMERICAL_FAILURE'
                            result.loc[failed,'exclusion_reason']=values.loc[failed,'failure']
                result['threshold']=pd.to_numeric(result.threshold)
                result['call_eligible']=np.isfinite(result.score)&np.isfinite(result.threshold)&~result.prior_only
                result['call']=result.call_eligible & (result.score>result.threshold)
                scores.append(result)
            print('Scored held fold',library,depth,chrom,flush=True)
    pd.concat(scores,ignore_index=True).to_csv(destination,sep='\t',index=False)
    dump(out/f'provenance/score.{library}.{depth}.json',dict(code_sha256=digest(__file__),counts_sha256=digest(out/f'tables/counts.{library}.{depth}.tsv.gz'),
        score_sha256=digest(destination),calibration_depth=calibration_depth,rows=len(frame),completed_utc=datetime.now(timezone.utc).isoformat(),
        evidence='development only',full_panel_screen_complete=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--stage',choices=['calibrate','score'],required=True);p.add_argument('--depth',required=True)
    p.add_argument('--library');p.add_argument('--calibration-depth',default='archived_full');p.add_argument('--workers',type=int,default=12)
    p.add_argument('--attempt',type=int,default=1)
    a=p.parse_args();out=a.outdir.resolve()
    # Per-stage immutable implementation identity is written before outcomes.
    name=f'{a.stage}.{a.library or "control"}.{a.depth}'+(f'.attempt{a.attempt}' if a.attempt!=1 else '')
    freeze=out/f'provenance/{name}.implementation.freeze.json'
    dump(freeze,dict(code={f.name:digest(f) for f in Path(__file__).parent.glob('*.py')},
                    command=sys.argv,created_utc=datetime.now(timezone.utc).isoformat(),scientific_policy_changed=False))
    if a.stage=='calibrate':calibrate(out,a.depth,a.workers)
    else:score(out,a.library,a.depth,a.calibration_depth,a.workers)
