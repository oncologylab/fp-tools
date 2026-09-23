"""Retain complete fitting-pool nuisance covariates and depth-specific balance."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from verify import read_table
from access import digest
from prepare import dump


def execute(out,attempt=1):
    suffix='' if attempt==1 else f'.attempt{attempt}'
    dump(out/f'provenance/pool_covariates{suffix}.freeze.json',dict(code_sha256=digest(__file__),
        purpose='complete actual library count covariates for the prespecified confound baseline and balance audit',
        depth_metadata='archived10M is a nominal historical label, not a verified complete-fragment count; the original constant depth covariate is retained separately and centers out within each library',
        cohort_or_scores_changed=False,coverage='all original pool sites, including those outside the headline union'))
    pool=read_table(out/'tables/pool.tsv.gz');nested=set(read_table(out/'tables/nested_sites.tsv.gz').site_id)
    all_rows=[];balance=[]
    for path in sorted((out/'tables').glob('counts.*.tsv.gz')):
        meta=read_table(path).iloc[0];library,depth=meta.library,meta.depth;cell=meta.cell
        sites=pool[pool.cell.eq(cell)].reset_index(drop=True).copy()
        source=out/(f'arrays/{cell}.replay_pool.npz' if depth=='archived_10M' else f'arrays/{library}.{depth}.npz')
        with np.load(source,allow_pickle=False) as z:
            np.testing.assert_array_equal(z['site_id'],sites.site_id.to_numpy(str))
            plus=z['plus_observed'].sum(1);minus=z['minus_observed'].sum(1)
            achieved=int(z['achieved_fragments']) if 'achieved_fragments' in z else np.nan
            depth_constant=achieved if np.isfinite(achieved) else 10000000
        sites['actual_whole201']=plus+minus;sites['plus_whole201']=plus;sites['minus_whole201']=minus
        sites['library']=library;sites['depth']=depth;sites['complete_library_fragments']=achieved;sites['depth_covariate_constant']=depth_constant;sites['cohort_nested']=sites.site_id.isin(nested)
        sites['depth_identity']='verified complete pairs' if np.isfinite(achieved) else 'nominal historical10M label; actual complete pairs unverified'
        all_rows.append(sites[['site_id','library','depth','actual_whole201','plus_whole201','minus_whole201','complete_library_fragments','depth_covariate_constant','depth_identity','cohort_nested']])
        for tf,group in sites.groupby('tf'):
            for phase,data in [('before',group[group.supported]),('after',group[group.cohort_nested])]:
                for field in ['motif_score','actual_whole201','plus_whole201','minus_whole201']:
                    value=data[field].to_numpy(float)
                    if field!='motif_score':value=np.log1p(value)
                    labels=data.chip_label.to_numpy(bool);a=value[labels];b=value[~labels]
                    difference=a.mean()-b.mean() if len(a) and len(b) else np.nan
                    variance=((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2) if len(a)>1 and len(b)>1 else np.nan
                    smd=difference/np.sqrt(variance) if variance>0 else 0. if difference==0 else np.nan
                    balance.append(dict(library=library,depth=depth,cell=cell,tf=tf,phase=phase,covariate=field,
                        transform='identity' if field=='motif_score' else 'log1p',sites=len(data),positives=len(a),negatives=len(b),
                        positive_mean=a.mean() if len(a) else np.nan,negative_mean=b.mean() if len(b) else np.nan,
                        standardized_mean_difference=smd,status='AVAILABLE' if np.isfinite(smd) else 'UNAVAILABLE'))
    pd.concat(all_rows,ignore_index=True).to_csv(out/'tables/pool_library_covariates.tsv.gz',sep='\t',index=False)
    pd.DataFrame(balance).to_csv(out/'tables/actual_library_covariate_balance.tsv',sep='\t',index=False)
    print('Complete fitting-pool library covariates and balance saved',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);p.add_argument('--attempt',type=int,default=1)
    a=p.parse_args();execute(a.outdir.resolve(),a.attempt)
