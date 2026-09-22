"""Replay unchanged historical count-tail policies on all development controls."""
import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from prepare import ROOT,dump
from access import digest
from verify import read_table
sys.path.insert(0,str(ROOT/'benchmarks/scripts'))
from conditional_protection_score import depletion_score
from analyse import wilson


def execute(out):
    dump(out/'provenance/control_old_calls.freeze.json',dict(code_sha256=digest(__file__),
        scorer_sha256=digest(ROOT/'benchmarks/scripts/conditional_protection_score.py'),
        cap=300,decimals=10,threshold=1.3010299957,strict_greater_than=True,
        scope='existing uniform/LOG21 count arms only; not a native historical DWM safety reconstruction'))
    all_rows=[];summaries=[]
    for depth in ['10M','full','archived_full']:
        data=read_table(out/f'tables/control_counts.{depth}.tsv.gz')
        for arm in ['uniform','LOG21']:
            frame=data.copy();scores=[]
            for r in frame.itertuples():
                pp,pm=getattr(r,arm+'_plus_p'),getattr(r,arm+'_minus_p')
                scores.append(round(depletion_score(r.plus_total,r.minus_total,r.plus_center+r.minus_center,pp,pm),10) if np.isfinite([pp,pm]).all() else np.nan)
            frame['score']=scores;frame['method']=arm+'_count_tail';frame['threshold']=1.3010299957
            frame['call_eligible']=frame.eligible&frame.total.gt(0)&np.isfinite(frame.score)
            frame['call']=frame.call_eligible&(frame.score>frame.threshold);all_rows.append(frame)
            for (cell,tf),group in frame.groupby(['cell','tf']):
                n=int(group.eligible.sum());info=int(group.call_eligible.sum());k=int(group.call.sum())
                summaries.append(dict(depth=depth,cell=cell,tf=tf,method=arm+'_count_tail',original_sites=len(group),eligible_sites=n,
                    informative_sites=info,calls=k,eligible_rate=k/n if n else np.nan,informative_rate=k/info if info else np.nan,
                    informative_wilson_high=wilson(k,info)[1],policy='historical cap/rounding/strict threshold unchanged'))
    pd.concat(all_rows,ignore_index=True).to_csv(out/'tables/control_old_scores.tsv.gz',sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(out/'tables/control_old_rates.tsv',sep='\t',index=False)
    print('Old control count calls replayed on all4200sites at each depth; excluded rows retained',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);execute(p.parse_args().outdir.resolve())
