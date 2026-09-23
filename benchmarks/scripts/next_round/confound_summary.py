"""Summarize prespecified covariate diagnostics without fitting another model."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from verify import read_table


def execute(out):
    base=read_table(out/'tables/covariate_baseline_metrics.tsv')
    tasks=base.groupby(['library','depth','cell','tf'])[['auroc','auprc']].mean().reset_index()
    tasks=tasks.groupby(['depth','cell','tf'])[['auroc','auprc']].mean().reset_index()
    tasks['method']='covariates_only';tasks.to_csv(out/'tables/covariate_baseline_task_metrics.tsv',sep='\t',index=False)
    data=read_table(out/'tables/stratified_metrics.tsv.gz');rows=[]
    keys=['library','depth','cell','tf','fold','stratum_type','stratum']
    for method in ['LOG21_binomial','LOG21_beta_binomial']:
        for reference in ['LOG21_count_tail','DWM_binomial']:
            pair=data[data.method.eq(method)].merge(data[data.method.eq(reference)],on=keys,suffixes=('','_reference'))
            pair['auroc_delta']=pair.auroc-pair.auroc_reference;pair['auprc_delta']=pair.auprc-pair.auprc_reference
            for (depth,kind,stratum),group in pair.groupby(['depth','stratum_type','stratum'],dropna=False):
                library=group.groupby(['cell','tf','library'])[['auroc_delta','auprc_delta']].mean()
                task=library.groupby(['cell','tf']).mean()
                rows.append(dict(depth=depth,stratum_type=kind,stratum=stratum,method=method,reference=reference,
                    available_tasks=int(task.auroc_delta.notna().sum()),original_tasks=21,
                    mean_auroc_delta=task.auroc_delta.mean(),mean_auprc_delta=task.auprc_delta.mean(),
                    scope='descriptive common-defined strata; no causal adjustment, no favorable subset used for advancement'))
    pd.DataFrame(rows).to_csv(out/'tables/stratified_summary.tsv',sep='\t',index=False)
    print('Covariate baseline and complete prespecified strata summarized',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);execute(p.parse_args().outdir.resolve())
