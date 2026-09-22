"""Predeclared prior-width diagnostic, never candidate or threshold selection."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import time
from shrinkage import posterior
from verify import read_table
from prepare import dump
from access import digest
import numpy as np
import pandas as pd


def calculate(args):
    identity,c,n,p,rho,width=args
    result=posterior(*c,*n,*p,rho,prior_sd=width)
    return dict(identity,prior_sd=width,mean=result['mean'],sd=result['sd'],low=result['low'],high=result['high'],
                prior_only=result['prior_only'],variance_reduction=result['variance_reduction'])


def execute(out):
    dump(out/'provenance/prior_sensitivity.freeze.json',dict(code_sha256=digest(__file__),
        prior_sd=[.5,2.],cohort='all nested sites, all six full-depth libraries',
        methods=['LOG21_binomial','LOG21_beta_binomial'],selection_allowed=False,
        calls='not computed; original N(0,1) thresholds do not transfer to changed priors'))
    receipt=out/'logs/initial_panel_queue.attempt2/exit.json'
    while not receipt.exists():time.sleep(20)
    if json.loads(receipt.read_text())['returncode']!=0:raise RuntimeError('initial panel failed')
    ids=set(read_table(out/'tables/nested_sites.tsv.gz').site_id);jobs=[]
    for path in sorted((out/'tables').glob('scores.*.full.tsv.gz')):
        scores=read_table(path);scores=scores[scores.site_id.isin(ids)&scores.method.isin(['LOG21_binomial','LOG21_beta_binomial'])]
        counts=read_table(out/'tables'/path.name.replace('scores.','counts.')).set_index('site_id')
        for r in scores.itertuples():
            if not np.isfinite(r.score):continue
            v=counts.loc[r.site_id]
            identity=dict(site_id=r.site_id,library=r.library,method=r.method,original_mean=r.score,depth='full')
            for width in [.5,2.]:jobs.append((identity,[int(v.plus_center),int(v.minus_center)],
                [int(v.plus_total),int(v.minus_total)],[v.LOG21_plus_p,v.LOG21_minus_p],r.rho,width))
    with ProcessPoolExecutor(max_workers=6) as pool:rows=list(pool.map(calculate,jobs,chunksize=64))
    pd.DataFrame(rows).to_csv(out/'tables/prior_sensitivity.tsv.gz',sep='\t',index=False)
    print('Prior diagnostic rows',len(rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);execute(p.parse_args().outdir.resolve())
