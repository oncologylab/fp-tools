"""Independent representative figure from packaged observations/coefficients.

This intentionally uses the first task in alphabetical order, never the best
result. It writes only to a new explicit path outside the evidence package.
"""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import logsumexp


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('package',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();root=a.package.resolve();out=a.output.resolve()
    if out.is_relative_to(root) or out.exists():raise ValueError('use a new output path outside the evidence package')
    controls=pd.read_csv(root/'inputs/control_counts.full.tsv.gz',sep='\t')
    cell,tf=sorted(set(zip(controls.cell,controls.tf)))[0]
    indices=np.flatnonzero(controls.cell.eq(cell)&controls.tf.eq(tf));x=np.arange(-100,101)
    domain=(abs(x)<=7)|((abs(x)>=12)&(abs(x)<=40))
    with np.load(root/'arrays/control.l2_0.01.npz',allow_pickle=False) as z:
        observed=z['plus_observed'][indices]+z['minus_observed'][indices]
        expected=np.zeros_like(observed,dtype=float)
        for strand in ['plus','minus']:
            log=z[strand+'_log_bias'][indices]
            total=z[strand+'_observed'][indices][:,domain].sum(1)
            expected+=np.exp(log-logsumexp(log[:,domain],axis=1)[:,None])*total[:,None]
    fig,ax=plt.subplots(figsize=(7,4));ax.plot(x,observed.mean(0),label='Observed cuts');ax.plot(x,expected.mean(0),label='Conditional expected cuts')
    ax.axvspan(-7,7,color='grey',alpha=.15);ax.legend()
    ax.set(xlabel='Position relative to motif center (bases)',ylabel='Cuts per original site',title=f'{cell} {tf}: examined naked-DNA replicate 1')
    fig.tight_layout();out.parent.mkdir(parents=True,exist_ok=True);fig.savefig(out,dpi=180)
    print(f'Regenerated first alphabetical task from all {len(indices)} original windows: {out.name}')
