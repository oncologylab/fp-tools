#!/usr/bin/env python3
"""Independently recompute all saved radius metrics and render every task."""
import argparse
import json
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score
sys.path.insert(0,str(Path(__file__).resolve().parent))
from compare_gamma_training import ROOT,KEYS
from review_gamma_training import read_table
from audit_chrombpnet_comparison import checked_path,canonical_id
from audit_direct_bias_safety import record,write_immutable


def review(source,out):
    manifest=source/'comparison_manifest.json';result=json.loads(manifest.read_text())
    identity=result.pop('stage_id')
    if canonical_id(result)!=identity or result['promotion_eligible'] or not result['training_only']:
        raise ValueError('invalid source result')
    tables={k:read_table(checked_path(v,ROOT)) for k,v in result['outputs'].items()}
    metrics,curves,scores=tables['metrics'],tables['curves'],tables['scores']
    if len(metrics)!=126 or len(metrics[['cell','tf']].drop_duplicates())!=21:
        raise ValueError('task denominator changed')
    if not scores.TFBS_chr.isin([f'chr{i}' for i in range(1,16)]).all(): raise ValueError('nontraining rows')
    checked=0
    for (cell,tf),task in scores.groupby(['cell','tf'],sort=False):
        arms=list(task.groupby('arm',sort=False));reference=arms[0][1][KEYS+['chip_label']].reset_index(drop=True)
        if len(arms)!=6: raise ValueError('arm missing')
        for arm,frame in arms:
            pd.testing.assert_frame_equal(reference,frame[KEYS+['chip_label']].reset_index(drop=True))
            found=metrics.loc[metrics.cell.eq(cell)&metrics.tf.eq(tf)&metrics.arm.eq(arm)]
            if len(found)!=1: raise ValueError('missing or duplicate metric')
            row=found.iloc[0]
            np.testing.assert_allclose([row.auroc,row.auprc],[roc_auc_score(frame.chip_label,frame.score),
                average_precision_score(frame.chip_label,frame.score)],rtol=0,atol=1e-12)
            if len(frame)!=2*row.retained_pairs or frame.chip_label.sum()!=row.retained_pairs:
                raise ValueError('pair count changed')
            checked+=1
    if checked!=126 or out.exists(): raise ValueError('incomplete review or existing output')
    out.mkdir(parents=True);outputs={}
    tasks=list(metrics[['cell','tf']].drop_duplicates().itertuples(index=False,name=None))
    colors=dict(raw='#999999',DWM_expected='#111111',parametric_direct='#0072B2',
                radius_100='#E69F00',radius_200='#CC79A7',radius_500='#009E73')
    pdf=out/'all_task_contrasts.pdf'
    with PdfPages(pdf) as pages:
        for page,start in enumerate(range(0,21,7),1):
            fig,axes=plt.subplots(4,2,figsize=(11.7,12),constrained_layout=True)
            for ax,(cell,tf) in zip(axes.flat,tasks[start:start+7]):
                subset=curves.loc[curves.cell.eq(cell)&curves.tf.eq(tf)]
                for arm,frame in subset.groupby('arm',sort=False):
                    ax.plot(frame.position,frame.positive_mean-frame.negative_mean,label=arm,color=colors[arm],lw=1)
                count=int(metrics.loc[metrics.cell.eq(cell)&metrics.tf.eq(tf),'retained_pairs'].iloc[0])
                ax.set_title(f'{cell} {tf}: {count} matched pairs');ax.axhline(0,color='#bbbbbb',lw=.5)
                ax.axvspan(-7,7,color='#dddddd',alpha=.35)
                ax.set_xlabel('Position relative to motif center (base pairs)');ax.set_ylabel('Normalized positive − negative')
            axes.flat[-1].axis('off');axes.flat[-1].legend(*axes.flat[0].get_legend_handles_labels(),loc='center')
            fig.suptitle('Training only; fixed geometry; no uncertainty, safety or promotion claim')
            pages.savefig(fig);path=out/f'contrasts_{page}.png';fig.savefig(path,dpi=110);plt.close(fig)
            outputs[path.name]=record(path)
    outputs['pdf']=record(pdf)
    review=dict(source=record(manifest),code=record(Path(__file__)),outputs=outputs,
        verified_metric_rows=checked,identical_site_sets=True,training_only=True,promotion_eligible=False,
        visual_review_complete=False)
    review['stage_id']=canonical_id(review);write_immutable(out/'review_manifest.json',review)
    print('Recomputed all 126 metric rows; rendered all 21 task contrasts',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args()
    if not args.outdir.resolve().is_relative_to(ROOT/'benchmarks/results'): parser.error('use research results')
    review(args.source.resolve(),args.outdir.resolve())
