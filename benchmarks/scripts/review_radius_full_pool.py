#!/usr/bin/env python3
"""Check scores/curves against saved panels and render all training tasks."""
import argparse,json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score
sys.path.insert(0,str(Path(__file__).resolve().parent))
from evaluate_radius_full_pool import ROOT,ARMS
from compare_gamma_training import KEYS
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import canonical_id,checked_path


def normalized(values):
    x=np.arange(-100,101);outer=abs(x)>=60
    design=np.column_stack([np.ones(outer.sum()),x[outer]])
    coefficients=np.linalg.lstsq(design,values[:,outer].T,rcond=None)[0]
    residual=values-coefficients[0,:,None]-coefficients[1,:,None]*x
    scale=np.sqrt((residual[:,outer]**2).mean(axis=1));scale=np.where(scale>1e-8,scale,1.)
    return np.clip(residual/scale[:,None],-12,12)


def review(source,out):
    manifest=source/'comparison_manifest.json';result=json.loads(manifest.read_text());identity=result.pop('stage_id')
    if canonical_id(result)!=identity or result['promotion_eligible'] or not result['training_only']: raise ValueError('invalid result')
    paths={key:checked_path(value,ROOT) for key,value in result['outputs'].items()}
    read=lambda name:pd.read_csv(paths[name],sep='\t',float_precision='round_trip')
    metrics,curves,scores,sites=[read(name) for name in ('metrics','curves','scores','sites')]
    if len(metrics)!=147 or len(metrics[['cell','tf']].drop_duplicates())!=21: raise ValueError('task denominator changed')
    if not sites.TFBS_chr.isin([f'chr{i}' for i in range(1,16)]).all(): raise ValueError('nontraining sites')
    with np.load(paths['panels'],allow_pickle=False) as archive: panels={arm:archive[arm] for arm in ARMS}
    x=np.arange(-100,101);inner=abs(x)<=7;outer=(abs(x)>=12)&(abs(x)<=40)
    weights=outer.astype(float)/outer.sum()-inner.astype(float)/inner.sum();checked=0;unavailable=0
    for row in metrics.itertuples(index=False):
        ix=np.flatnonzero((sites.cell.eq(row.cell)&sites.tf.eq(row.tf)).to_numpy());y=sites.iloc[ix].chip_label.to_numpy(int)
        saved=scores.loc[scores.cell.eq(row.cell)&scores.tf.eq(row.tf)&scores.arm.eq(row.arm)]
        if len(ix)!=2*row.matched_pairs: raise ValueError('pair count changed')
        if not len(ix):
            if len(saved) or not np.isnan(row.auroc): raise ValueError('unavailable task has scores')
            unavailable+=1;continue
        pd.testing.assert_frame_equal(saved[KEYS+['chip_label']].reset_index(drop=True),sites.iloc[ix][KEYS+['chip_label']].reset_index(drop=True))
        calculated=panels[row.arm][ix]@weights;calculated[panels['raw'][ix].sum(1)==0]=0
        np.testing.assert_allclose(calculated,saved.score,rtol=0,atol=1e-12)
        n=y.sum();roc=(rankdata(saved.score)[y==1].sum()-n*(n+1)/2)/(n*(len(y)-n))
        np.testing.assert_allclose([roc,average_precision_score(y,saved.score)],[row.auroc,row.auprc],rtol=0,atol=1e-12)
        norm=normalized(panels[row.arm][ix]);positive=norm[y==1].mean(0);negative=norm[y==0].mean(0)
        curve=curves.loc[curves.cell.eq(row.cell)&curves.tf.eq(row.tf)&curves.arm.eq(row.arm)].sort_values('position')
        np.testing.assert_array_equal(curve.position,x)
        np.testing.assert_allclose(positive,curve.positive_mean,rtol=1e-9,atol=1e-9)
        np.testing.assert_allclose(negative,curve.negative_mean,rtol=1e-9,atol=1e-9)
        delta=positive-negative;center=delta[abs(x)<=5].mean();left=delta[(x>=-40)&(x<=-15)].mean();right=delta[(x>=15)&(x<=40)].mean()
        np.testing.assert_allclose([center,left,right],[row.center_contrast,row.left_shoulder,row.right_shoulder],atol=1e-9)
        if bool(center<0 and left>center and right>center)!=bool(row.canonical_protection): raise ValueError('canonical direction differs')
        checked+=1
    if out.exists(): raise ValueError('fresh review directory required')
    out.mkdir(parents=True);outputs={};pdf=out/'all_task_contrasts.pdf'
    tasks=list(metrics[['cell','tf']].drop_duplicates().itertuples(index=False,name=None))
    colors=dict(raw='#999999',DWM_expected='#111111',parametric_direct='#0072B2',radius_100='#E69F00',radius_200='#CC79A7',radius_500='#009E73',uniform_radius_500='#D55E00')
    with PdfPages(pdf) as pages:
        for page,start in enumerate(range(0,21,7),1):
            fig,axes=plt.subplots(4,2,figsize=(11.7,12),constrained_layout=True)
            for ax,(cell,tf) in zip(axes.flat,tasks[start:start+7]):
                subset=curves.loc[curves.cell.eq(cell)&curves.tf.eq(tf)]
                for arm,frame in subset.groupby('arm',sort=False):
                    ax.plot(frame.position,frame.positive_mean-frame.negative_mean,label=arm,color=colors[arm],lw=1,
                        ls='--' if arm=='uniform_radius_500' else '-')
                m=metrics.loc[metrics.cell.eq(cell)&metrics.tf.eq(tf)].iloc[0]
                ax.set_title(f'{cell} {tf}: {int(m.matched_pairs)} pairs ({m.status})',fontsize=10)
                if subset.empty: ax.text(.5,.5,'No balanced pairs',ha='center',transform=ax.transAxes)
                ax.axhline(0,color='#bbbbbb',lw=.5);ax.axvspan(-7,7,color='#dddddd',alpha=.35)
                ax.set_xlabel('Position relative to motif center (base pairs)');ax.set_ylabel('Normalized positive − negative')
            axes.flat[-1].axis('off');axes.flat[-1].legend(*axes.flat[0].get_legend_handles_labels(),loc='center')
            fig.suptitle('Three-covariate matched training cohorts; no safety or uncertainty claim')
            pages.savefig(fig);path=out/f'contrasts_{page}.png';fig.savefig(path,dpi=110);plt.close(fig);outputs[path.name]=record(path)
    outputs['pdf']=record(pdf)
    value=dict(source=record(manifest),code=record(Path(__file__)),outputs=outputs,
        verified_score_and_curve_rows=checked,verified_unavailable_rows=unavailable,
        visual_review_complete=False,promotion_eligible=False)
    value['stage_id']=canonical_id(value);write_immutable(out/'review_manifest.json',value)
    print(f'Independently verified {checked} score/curve rows and {unavailable} unavailable rows; rendered all 21 tasks',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args()
    if not args.outdir.resolve().is_relative_to(ROOT/'benchmarks/results'): parser.error('use research results')
    review(args.source.resolve(),args.outdir.resolve())
