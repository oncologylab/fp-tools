"""Declared known-count fixtures: conditional uncertainty and synthetic power."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.special import expit,logit
from scipy.stats import betabinom
from shrinkage import posterior,threshold
from access import digest
from prepare import dump
from protocol import POLICY


def simulate_one(args):
    scenario_id,n,p,theta,rho_true,shift=args
    rng=np.random.default_rng(np.random.SeedSequence([2026,scenario_id]))
    n=np.asarray(n,int);p=np.asarray(p,float)
    q=expit(logit(p)-theta+np.asarray(shift))
    repetitions=2000
    if rho_true:
        probabilities=rng.beta(q*(1/rho_true-1),(1-q)*(1/rho_true-1),size=(repetitions,2))
    else:probabilities=np.broadcast_to(q,(repetitions,2))
    c=rng.binomial(n,probabilities)
    rows=[];summaries=[]
    for method,rho in [('binomial',0.),('beta_binomial',.1)]:
        # Known-model calibration fixtures use independent generated windows.
        # These are not the biological calling thresholds or safety evidence.
        calibration_p=rng.beta(p*(1/rho-1),(1-p)*(1/rho-1),size=(repetitions,2)) if rho else np.broadcast_to(p,(repetitions,2))
        calibration=rng.binomial(n,calibration_p)
        cutoff=threshold([posterior(*cc,*n,*p,rho)['mean'] for cc in calibration]) if n.sum() else np.nan
        result=[posterior(*cc,*n,*p,rho) for cc in c]
        means=np.array([r['mean'] for r in result]);lower=np.array([r['low'] for r in result]);upper=np.array([r['high'] for r in result])
        calls=(means>cutoff)&bool(n.sum())
        covered=(lower<=theta)&(upper>=theta)
        frame=pd.DataFrame(dict(scenario=scenario_id,draw=np.arange(repetitions),method=method,
            plus_center=c[:,0],minus_center=c[:,1],plus_total=n[0],minus_total=n[1],p_plus=p[0],p_minus=p[1],
            true_theta=theta,generating_rho=rho_true,assumed_rho=rho,mean=means,low=lower,high=upper,
            threshold=cutoff,call=calls,covered=covered,shift_plus=shift[0],shift_minus=shift[1]))
        rows.append(frame)
        summaries.append(dict(scenario=scenario_id,method=method,n_plus=n[0],n_minus=n[1],p_plus=p[0],p_minus=p[1],
            theta=theta,generating_rho=rho_true,assumed_rho=rho,shift_plus=shift[0],shift_minus=shift[1],
            repetitions=repetitions,independent_source_blocks='model-generated independent windows',
            model_correct=bool(rho==rho_true and not any(shift)),calls=int(calls.sum()),call_rate=float(calls.mean()),
            monte_carlo_se=float(np.sqrt(calls.mean()*(1-calls.mean())/repetitions)),coverage=float(covered.mean()),
            mean_error=float((means-theta).mean()),rmse=float(np.sqrt(((means-theta)**2).mean())),
            mean_expected_central=float(n@p),mean_observed_central=float(c.sum(1).mean()),
            calling_policy='separate independent synthetic-null order-statistic fixture, not biological calibrated call',
            status='AVAILABLE' if n.sum() else 'PRIOR_ONLY_NO_CALL'))
    return pd.concat(rows,ignore_index=True),pd.DataFrame(summaries)


def execute(out,workers):
    freeze=out/'provenance/model_simulation.freeze.json'
    scenarios=[]
    for n in POLICY['model_simulation']['n_pairs']:
        for p in POLICY['model_simulation']['p_pairs']:
            for theta in POLICY['model_simulation']['theta']:
                for rho in POLICY['model_simulation']['rho']:
                    scenarios.append((len(scenarios),n,p,theta,rho,[0.,0.]))
            for shift in [[.5,.5],[.5,-.5]]:
                scenarios.append((len(scenarios),n,p,0.,0.,shift))
    dump(freeze,dict(code_sha256=digest(__file__),scenarios=scenarios,repetitions=2000,
        calibration_repetitions=2000,seed=2026,posterior_code_sha256=digest(Path(__file__).with_name('shrinkage.py')),
        empirical_control_evidence=False,workers=workers))
    summaries=[]
    directory=out/'tables/model_simulations';directory.mkdir()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i,(rows,summary) in enumerate(pool.map(simulate_one,scenarios,chunksize=1)):
            rows.to_csv(directory/f'{i:03d}.tsv.gz',sep='\t',index=False)
            summaries.append(summary)
            print('Completed known-count scenario',i+1,'/',len(scenarios),flush=True)
    pd.concat(summaries,ignore_index=True).to_csv(out/'tables/model_power.tsv',sep='\t',index=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);p.add_argument('--workers',type=int,default=4)
    a=p.parse_args();execute(a.outdir.resolve(),a.workers)
