#!/usr/bin/env python3
"""Independently check saved training scores and render all task contrasts."""
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
from sklearn.metrics import roc_auc_score, average_precision_score

sys.path.insert(0,str(Path(__file__).resolve().parent))
from compare_gamma_training import ROOT, KEYS
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path


def read_table(path):
    """Keep binary-float score ranks unchanged through text serialization."""
    return pd.read_csv(path,sep='\t',float_precision='round_trip')


def review(source,out):
    manifest=source/'comparison_manifest.json'
    result=json.loads(manifest.read_text())
    identity=result.pop('stage_id')
    if canonical_id(result)!=identity or result['promotion_eligible'] is not False:
        raise ValueError('invalid source result')
    # Preserve exact saved binary floats: the default parser can merge or
    # split adjacent scores and thereby change rank metrics at tied counts.
    tables={key:read_table(checked_path(item,ROOT))
            for key,item in result['outputs'].items()}
    metrics,curves,scores=tables['metrics'],tables['curves'],tables['scores']
    if len(metrics)!=105 or metrics[['cell','tf']].drop_duplicates().shape[0]!=21:
        raise ValueError('incomplete task denominator')
    if not scores.TFBS_chr.isin([f'chr{i}' for i in range(1,16)]).all():
        raise ValueError('evaluation outside training chromosomes')
    checked=0
    for (cell,tf),task in scores.groupby(['cell','tf'],sort=False):
        arms=list(task.groupby('arm',sort=False))
        reference=arms[0][1][KEYS+['chip_label']].reset_index(drop=True)
        for arm,frame in arms:
            pd.testing.assert_frame_equal(reference,frame[KEYS+['chip_label']].reset_index(drop=True))
            row=metrics.loc[metrics.cell.eq(cell)&metrics.tf.eq(tf)&metrics.arm.eq(arm)].iloc[0]
            np.testing.assert_allclose([row.auroc,row.auprc],
                [roc_auc_score(frame.chip_label,frame.score),average_precision_score(frame.chip_label,frame.score)],rtol=0,atol=1e-12)
            if len(frame)!=2*row.retained_pairs or frame.chip_label.sum()!=row.retained_pairs:
                raise ValueError('class balance or pair count changed')
            checked+=1
    if out.exists():
        raise ValueError('fresh review directory required')
    out.mkdir(parents=True)
    outputs={}
    tasks=list(metrics[['cell','tf']].drop_duplicates().itertuples(index=False,name=None))
    colors=dict(raw='#777777',DWM_expected='#111111',parametric_direct='#0072B2',
                gamma_fixed_dispersion='#D55E00',gamma_predictive='#009E73')
    pdf=out/'all_task_contrasts.pdf'
    with PdfPages(pdf) as pages:
        for page,start in enumerate(range(0,len(tasks),7),1):
            fig,axes=plt.subplots(4,2,figsize=(11.7,12),constrained_layout=True)
            for ax,(cell,tf) in zip(axes.flat,tasks[start:start+7]):
                subset=curves.loc[curves.cell.eq(cell)&curves.tf.eq(tf)]
                for arm,frame in subset.groupby('arm',sort=False):
                    ax.plot(frame.position,frame.positive_mean-frame.negative_mean,label=arm,color=colors[arm],lw=1)
                row=metrics.loc[metrics.cell.eq(cell)&metrics.tf.eq(tf)].iloc[0]
                ax.set_title(f'{cell} {tf}: {row.retained_pairs} matched pairs')
                ax.axhline(0,color='#aaaaaa',lw=.5);ax.axvspan(-7,7,color='#dddddd',alpha=.35)
                ax.set_xlabel('Position relative to motif center (base pairs)')
                ax.set_ylabel('Normalized positive − negative')
            axes.flat[-1].axis('off')
            axes.flat[-1].legend(*axes.flat[0].get_legend_handles_labels(),loc='center',fontsize=10)
            fig.suptitle('Training chromosomes 1–15; fixed geometry; no uncertainty or promotion claim')
            pages.savefig(fig)
            png=out/f'contrasts_{page}.png';fig.savefig(png,dpi=110);plt.close(fig)
            outputs[png.name]=record(png)
    outputs['pdf']=record(pdf)
    report=dict(source=record(manifest),code=record(Path(__file__)),outputs=outputs,
        score_rows_recomputed=checked,identical_site_sets_verified=True,training_chromosomes_verified=True,
        visual_review_complete=False,promotion_eligible=False)
    report['stage_id']=canonical_id(report)
    write_immutable(out/'review_manifest.json',report)
    print(f'Independently verified {checked} metric rows; rendered all 21 tasks')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args()
    if not args.outdir.resolve().is_relative_to(ROOT/'benchmarks/results'):
        parser.error('use research results')
    review(args.source.resolve(),args.outdir.resolve())
