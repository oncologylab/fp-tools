"""Render the existing principal curve table with explicit units and legend."""
import argparse
from pathlib import Path
import shutil
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def execute(out):
    keep=[]
    source=out/'tables/shape_curves.tsv.gz'
    if not source.exists():source=out/'tables/shape_curves.tsv.xz'
    for chunk in pd.read_csv(source,sep='\t',chunksize=100000,float_precision='round_trip'):
        keep.append(chunk[chunk.stratum_type.isna()&chunk.strand.eq('combined')&chunk.aggregation.eq('site_equal')])
    curves=pd.concat(keep,ignore_index=True);tasks=pd.read_csv(out/'tables/task_ledger.tsv',sep='\t')
    path=out/'figures/all_task_shapes.pdf';old=out/'historical_context/all_task_shapes_initial_unlabelled.pdf'
    if path.exists() and not old.exists():shutil.copyfile(path,old)
    with PdfPages(path) as pdf:
        for (library,depth),frame in curves.groupby(['library','depth'],sort=True):
            cell=library.split('_')[0];celltasks=tasks[tasks.cell.eq(cell)]
            fig,axes=plt.subplots(len(celltasks),3,figsize=(14,2*len(celltasks)),squeeze=False)
            for row,t in enumerate(celltasks.itertuples()):
                for col,representation in enumerate(['raw','LOG21_expected','LOG21_residual']):
                    ax=axes[row,col];ax.axvspan(-7,7,color='.9')
                    selected=frame[frame.tf.eq(t.tf)&frame.representation.eq(representation)]
                    for label,color,text in [(0,'#777777','label-negative'),(1,'#1675ae','label-positive')]:
                        g=selected[selected.label.eq(label)].sort_values('position')
                        if len(g):ax.plot(g.position,g.value,color=color,lw=.8,label=text)
                    if not len(selected):ax.text(.5,.5,'UNAVAILABLE',transform=ax.transAxes,ha='center')
                    ax.set_title(f'{cell} {t.tf}: {representation}',fontsize=9)
                    ax.set_xlabel('base pairs from motif centre',fontsize=8)
                    ax.set_ylabel('signed deviance residual' if representation=='LOG21_residual' else 'cuts per site per base',fontsize=7)
                    ax.tick_params(labelsize=7);ax.set_xlim(-100,100)
            fig.suptitle(f'{library} / {depth}: complete held cohort; no score selection',fontsize=12,y=.998)
            fig.legend(*axes[0,0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.986),ncol=2,fontsize=8)
            fig.tight_layout(rect=[0,0,1,.975]);pdf.savefig(fig);plt.close(fig)
    print('Rendered14 complete endpoint pages from existing curve table; explicit units/legend',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);execute(p.parse_args().outdir.resolve())
