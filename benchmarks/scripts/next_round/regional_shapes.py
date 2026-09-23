"""Summarize both prespecified geometries from existing principal curves.

Descriptive reporting only: no fitting, selection, scores or gate changes.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def execute(out):
    keys=['library','depth','cell','tf','representation','label','aggregation','strand','units']
    rows=[]
    for chunk in pd.read_csv(out/'tables/shape_curves.tsv.gz',sep='\t',chunksize=201000,
                             float_precision='round_trip',low_memory=False):
        # Explicit all-site strata include each strand and both weightings.
        # The legacy empty-stratum rows duplicate combined/site-equal panels.
        principal=chunk[chunk.stratum_type.eq('all')]
        for values,group in principal.groupby(keys,dropna=False,sort=False):
            if len(group)!=201 or set(group.position)!=set(range(-100,101)):
                raise AssertionError('principal curve panel split or incomplete')
            x=group.position.to_numpy();y=group.value.to_numpy()
            for name,center,inner,outer in [('canonical',5,15,40),('scoring',7,12,40)]:
                c=y[np.abs(x)<=center].mean()
                left=y[(x>=-outer)&(x<=-inner)].mean()
                right=y[(x>=inner)&(x<=outer)].mean()
                rows.append(dict(zip(keys,values),geometry=name,center_mean=c,
                    left_mean=left,right_mean=right,left_minus_center=left-c,
                    right_minus_center=right-c,right_minus_left=right-left,
                    status='AVAILABLE' if np.isfinite([c,left,right]).all() else 'UNAVAILABLE'))
    pd.DataFrame(rows).to_csv(out/'tables/principal_regional_shapes.tsv',sep='\t',index=False)
    print('Descriptive principal region summaries:',len(rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True)
    execute(p.parse_args().outdir.resolve())
