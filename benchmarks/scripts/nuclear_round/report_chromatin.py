"""Report the frozen failed screen; no selection, fitting or new scoring."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def report(out):
    sub=out/'chromatin';read=lambda name:pd.read_csv(sub/'tables'/name,sep='\t',low_memory=False)
    summary=read('screen_summary.tsv');primary=summary[summary.method.eq('nuclear_primary')&summary.reference.eq('DWM_count_tail')].set_index('depth')
    r=primary.loc['full'];assert r.mean_auroc_gain<.03 and r.mean_relative_auprc_gain<.1 and r.bootstrap_low<=0
    decision=dict(decision='REJECTED_FOR_FOOTPRINT_IMPROVEMENT',control_component='DEVELOPMENT_CRITERION_MET',
        chromatin_screen='FAIL',native_DWM='INCOMPLETE',locked_controls='NOT_TESTED_LOCKED',
        expansion='INTENTIONALLY_NOT_RUN_FAILED_ACCURACY_SCREEN',scientific_promotion=False)
    (out/'final.decision.json').write_text(json.dumps(decision,indent=2)+'\n')
    tasks=read('task_ledger.tsv');tasks['new_chromatin_status']=np.where(tasks.tf.eq('HMGA1')&tasks.cell.eq('HepG2'),'UNAVAILABLE','COMPLETED_DEVELOPMENT_DIAGNOSTIC')
    tasks.to_csv(out/'tables/task_ledger.tsv',sep='\t',index=False)
    for name in ['fold_metrics.tsv','paired_deltas.tsv','shape_metrics.tsv']:
        read(name).to_csv(out/'tables'/name,sep='\t',index=False)
    gates=pd.read_csv(out/'tables/promotion_gates.tsv',sep='\t')
    for key in ['AUROC_gain','relative_AP_gain']:
        gates.loc[gates.gate.eq(key),['status','scope']]=['FAIL','full-depth fixed-cohort DWM-propensity diagnostic; native comparison unqualified']
    gates.loc[gates.scope.eq('not qualified by this control-only round'),'scope']='not qualified by this bounded development round'
    gates.to_csv(out/'tables/promotion_gates.tsv',sep='\t',index=False)
    panel=pd.read_csv(out/'tables/depth_replication.tsv',sep='\t');done=panel.depth.isin(['full','10M'])&panel.seed.eq(2026)
    panel['status']=np.where(done,'COMPLETED_DEVELOPMENT_DIAGNOSTIC','INTENTIONALLY_NOT_RUN_FAILED_SCREEN')
    panel['evidence']=np.where(done,'chromatin/tables; all six libraries','no expansion after failed screen')
    panel.to_csv(out/'tables/depth_replication.tsv',sep='\t',index=False)
    resources=pd.read_csv(out/'tables/resource_results.tsv',sep='\t');resources['status']='INTENTIONALLY_NOT_RUN_FAILED_SCREEN'
    resources.to_csv(out/'tables/resource_results.tsv',sep='\t',index=False)
    runs=pd.read_csv(out/'tables/run_ledger.tsv',sep='\t');runs.loc[runs.stage.eq('chromatin'),'status']='COMPLETED_DEVELOPMENT_DIAGNOSTIC'
    runs.loc[runs.stage.isin(['robustness_panel','full_resource_matrix']),'status']='INTENTIONALLY_NOT_RUN_FAILED_SCREEN'
    runs.loc[runs.stage.eq('locked_controls'),'status']='NOT_TESTED_LOCKED';runs.to_csv(out/'tables/run_ledger.tsv',sep='\t',index=False)
    findings=f'''# Nuclear-bias follow-up: rejected for footprint improvement

The nuclear model improves held-chromosome naked-DNA calibration but fails the
prespecified chromatin accuracy screen. Main, releases, manuscript and historical
scientific results remain unchanged. No locked labels or final controls were opened.

## What ran and changed

150 nuclear LOG21 fits (15 held chromosomes, five seeds, two frozen regularizations),
30 ensembles, all 4,200 development-control windows, and 6,020 fixed training-cohort
sites in all six biological libraries at full depth and 10M. All 21 tasks are retained;
HepG2 HMGA1 is unavailable. The old count-tail detector, geometry, rounding and
strict threshold did not change. The primary model was not replaced by its ablation.
Historical calibration plots were corrected separately without changing 19 numerical tables.

## Findings

1. Nuclear training passes the development control component: 2/21 task mean
   discrepancies; paired held-control log-likelihood gain 101.39 [64.48,145.52]
   versus mitochondrial LOG21 and 73.36 [46.25,103.35] versus intercept-only.
   Zero calls among 3,470 informative task-window rows is not proof of safety or power.
2. Full-depth chromatin AUROC gain is {r.mean_auroc_gain:.6f}, interval
   [{r.bootstrap_low:.6f},{r.bootstrap_high:.6f}], and relative average-precision gain
   {100*r.mean_relative_auprc_gain:.3f}%. Required gains are 0.03 and 10%; both fail.
   At 10M the gains are {primary.loc['10M','mean_auroc_gain']:.6f} and
   {100*primary.loc['10M','mean_relative_auprc_gain']:.3f}%. Comparators are same-detector
   DWM sequence propensities, not native TOBIAS or qualified production correction.
3. The primary model also loses ranking accuracy against frozen mitochondrial and
   intercept-only LOG21. Better naked-DNA likelihood did not translate into better
   footprint discrimination. This contradicts the proposed rescue by this nuclear
   sequence model; it does not establish that all detectors or models must fail.
4. Raw, expected and residual all-task shape tables and plots are provided separately.
   Raw cuts are unchanged. No stronger canonical protection or functional uncertainty
   is claimed from a ranking metric or model-derived expected profile.
5. Native DWM strand/complete-pair qualification remains incomplete. Full depth/seed
   expansion and end-to-end resources were intentionally stopped after the failed
   screen. Locked safety, functional uncertainty and general promotion are unqualified.

## Interpretation and limits

Data are repeatedly examined development observations. The cohort is the previously
constructed nested cohort, not newly independent evidence. Paired intervals condition
on fixed predictions, with 15 chromosome blocks and repeated selection history.
ChIP-negative call rates are not known false-positive rates. Low-count non-detection
remains a mixture of attainable-significance limits and unresolved biological/mean
calibration limitations; the fragment-based synthetic challenge is not binding biology.
Power tables retain missing perturbations and dependence on one biological source.
Mean calibration is improved on nuclear naked DNA, not established in chromatin.
No beta-binomial parameter was added in this round; previously rejected shrinkage
models remain historical. Nuclear-source, preparation, depth and window-selection
differences prevent attributing the change uniquely to DNA source.

Stop this frozen candidate. Do not reopen locked data, relax gates, or add another
detector to find a positive result. A future round needs a distinct falsifiable question.
'''
    (out/'EXECUTIVE_FINDINGS.md').write_text(findings)
    with (out/'NEXT_DEEP_RESEARCH_PROMPT.md').open('a') as f:
        f.write('\nThe conditional chromatin stage completed and failed accuracy: inspect chromatin/ and final.decision.json. Challenge whether better naked-DNA likelihood without chromatin improvement supports stopping this model, whether shape claims remain appropriately limited, and whether cohort, count/confounding and native-reference gaps invalidate stronger conclusions.\n')
    # Every task is displayed, including the unavailable task; no best-task selection.
    shapes=read('shape_metrics.tsv');fig,axes=plt.subplots(7,3,figsize=(14,19))
    for ax,t in zip(axes.flat,tasks.itertuples()):
        g=shapes[shapes.cell.eq(t.cell)&shapes.tf.eq(t.tf)&shapes.depth.eq('full')]
        for arm,color in [('nuclear_primary','tab:blue'),('DWM_propensity','tab:gray')]:
            h=g[g.arm.eq(arm)&g.status.eq('AVAILABLE')]
            ax.scatter([1,2,3][:len(h)],h.center_contrast,color=color,label=arm,s=18)
        ax.axhline(0,color='black',lw=.5);ax.set_title(f'{t.cell} {t.tf}',fontsize=9);ax.set_xticks([1,2,3]);ax.set_xlabel('Biological replicate',fontsize=8)
        if not g.status.eq('AVAILABLE').any():ax.text(.5,.5,'UNAVAILABLE',transform=ax.transAxes,ha='center')
    axes[0,0].legend(fontsize=6);fig.suptitle('Full-depth signed center contrast: negative alone does not establish canonical protection')
    fig.tight_layout(rect=[0,0,1,.98]);fig.savefig(sub/'figures/all_task_center_contrasts.pdf');plt.close(fig)
    # Raw measured cuts, expected cuts and residuals remain separate units.
    source=read('shape_plot_source.tsv.gz')
    source=source[source.arm.eq('nuclear_primary')&source.strand.eq('combined')&source.aggregation.eq('site_equal')&source.depth.eq('full')]
    with PdfPages(sub/'figures/all_task_profiles.pdf') as pdf:
        for t in tasks.itertuples():
            fig,axes=plt.subplots(1,3,figsize=(12,3.2))
            for ax,rep in zip(axes,['raw','expected','residual']):
                for label,color in [(0,'tab:gray'),(1,'tab:blue')]:
                    g=source[source.cell.eq(t.cell)&source.tf.eq(t.tf)&source.representation.eq(rep)&source.label.eq(label)]
                    mean=g.groupby('position').value.mean()
                    ax.plot(mean.index,mean,color=color,label=f'ChIP label {label}')
                ax.axvspan(-5,5,color='gray',alpha=.15);ax.set_title(rep);ax.set_xlabel('Bases from motif center')
                ax.set_ylabel('Signed deviance residual' if rep=='residual' else 'Cuts/site/base')
            axes[0].legend(fontsize=8);fig.suptitle(f'{t.cell} {t.tf}: full depth, equal biological-library means')
            fig.tight_layout();pdf.savefig(fig);plt.close(fig)
    # Descriptive ROC/PR curves retain original fold/library/task averaging.
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'calibration_round'))
    try:from curve_reporting import curves
    except ImportError:
        sys.path.insert(0,str(out/'reporting_code'));from curve_reporting import curves
    scores=read('per_site_scores.tsv.gz');cov=read('covariates.tsv.gz')
    arms=['DWM_count_tail','LOG21_count_tail','nuclear_primary','covariates_only']
    grid,rows,plot,metrics=curves(scores,cov,tasks,arms,'LOG21_count_tail')
    plot.to_csv(sub/'tables/ROC_PR_plot_source.tsv.gz',sep='\t',index=False)
    metrics.to_csv(sub/'tables/ROC_PR_scalar_source.tsv',sep='\t',index=False)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for arm in arms:
        values=[v for k,v in rows.items() if k[2]==arm and v is not None]
        axes[0].plot(grid,np.mean([v['tpr'] for v in values],axis=0),label=arm)
        axes[1].plot(grid,np.mean([v['prec'] for v in values],axis=0),label=arm)
    axes[0].plot([0,1],[0,1],'k:',lw=.6);axes[0].set(xlabel='False positive rate',ylabel='True positive rate',title='ROC: full-depth development')
    axes[1].set(xlabel='Recall',ylabel='Precision',title='PR: descriptive mean curve, not scalar AP')
    axes[0].legend(fontsize=7);fig.tight_layout();fig.savefig(sub/'figures/ROC_PR.png',dpi=180);fig.savefig(sub/'figures/ROC_PR.pdf');plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);report(p.parse_args().outdir)
