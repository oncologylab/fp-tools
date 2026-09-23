"""Explicit calibration replay from packaged controls, separate from verification."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import numpy as np
import pandas as pd
from shrinkage import posterior,fit_dispersion,threshold
from verify import read_table


def calculate(job):
    record,frame,refit=job;method=record['method'];arm=method.split('_')[0]
    dispersion=fit_dispersion(frame[['plus_center','minus_center']],frame[['plus_total','minus_total']],frame[[arm+'_plus_p',arm+'_minus_p']]) if refit else record['dispersion']
    rho=dispersion['rho'] if 'beta' in method else 0.
    if 'beta' in method and not dispersion['identified']:
        return dict(method=method,fold=record['fold'],status='UNIDENTIFIED'),[]
    fitting=frame[frame.total>0];rows=[]
    for r in fitting.itertuples():
        result=posterior(int(r.plus_center),int(r.minus_center),int(r.plus_total),int(r.minus_total),
            float(getattr(r,arm+'_plus_p')),float(getattr(r,arm+'_minus_p')),rho)
        rows.append(dict(site_id=r.site_id,fold=record['fold'],method=method,score=result['mean']))
    cutoff=threshold([r['score'] for r in rows])
    expected=record['threshold'] if record['threshold'] is not None else np.nan
    return dict(method=method,fold=record['fold'],status='AVAILABLE',rho=rho,original_rho=record['rho'],
                threshold=cutoff,original_threshold=expected,threshold_difference=cutoff-expected,
                fitting_informative_sites=len(rows),refitted_dispersion=refit),rows


def execute(root,out,depth,workers,refit):
    if out.exists():raise ValueError('new output directory required')
    out.mkdir(parents=True)
    controls=read_table(root/f'tables/control_counts.{depth}.tsv.gz')
    records=json.loads((root/f'models/calibration.{depth}.json').read_text())['records'];jobs=[]
    for r in records:
        arm=r['method'].split('_')[0]
        fit=controls[controls.eligible&controls.deduplicated_fit_representative&controls.TFBS_chr.ne(r['fold'])&controls[arm+'_plus_valid']&controls[arm+'_minus_valid']]
        jobs.append((r,fit,refit))
    summaries=[];scores=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for summary,rows in pool.map(calculate,jobs):summaries.append(summary);scores.extend(rows)
    summary=pd.DataFrame(summaries)
    summary.to_csv(out/'calibration_replay.tsv',sep='\t',index=False)
    pd.DataFrame(scores).to_csv(out/'fitting_scores.tsv.gz',sep='\t',index=False)
    print('Calibration records',len(summary),'maximum threshold discrepancy',summary.threshold_difference.abs().max(),flush=True)
    if summary.threshold_difference.abs().max()>1e-8:raise AssertionError('calibration replay differs beyond declared tolerance')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--depth',choices=['10M','full','archived_full'],default='full');p.add_argument('--workers',type=int,default=4)
    p.add_argument('--refit-dispersion',action='store_true',help='re-estimate rho on fitting chromosomes; otherwise replay frozen rho')
    a=p.parse_args();execute(a.root.resolve(),a.outdir.resolve(),a.depth,a.workers,a.refit_dispersion)
