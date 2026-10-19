"""Independent representative raw-profile reconstruction from packaged arrays."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def replot(root,output):
    sub=root/'chromatin'
    sites=pd.read_csv(sub/'tables/site_ledger.tsv.gz',sep='\t')
    first=sites[['cell','tf']].drop_duplicates().sort_values(['cell','tf']).iloc[0]
    library=first.cell+'_rep1'
    with np.load(sub/f'arrays/{library}.full.npz',allow_pickle=False) as z:
        table=sites.set_index('site_id').loc[z['site_id']];observed=z['observed'].sum(1)
    source=pd.read_csv(sub/'tables/shape_plot_source.tsv.gz',sep='\t')
    source=source[source.library.eq(library)&source.depth.eq('full')&source.tf.eq(first.tf)&source.arm.eq('nuclear_primary')&source.representation.eq('raw')&source.strand.eq('combined')&source.aggregation.eq('site_equal')]
    fig,ax=plt.subplots(figsize=(6,3.5))
    for label in [0,1]:
        mean=observed[table.tf.eq(first.tf)&table.chip_label.eq(label)].mean(0)
        expected=source[source.label.eq(label)].sort_values('position')
        np.testing.assert_allclose(mean,expected.value,rtol=1e-12,atol=1e-12)
        ax.plot(np.arange(-100,101),mean,label=f'ChIP label {label}')
    ax.set(title=f'{library} {first.tf}: raw measured cuts',xlabel='Bases from motif center',ylabel='Cuts/site/base');ax.legend()
    fig.tight_layout();fig.savefig(output,dpi=180);plt.close(fig)
    print('VERIFIED: 402 raw profile positions recomputed from complete primary arrays')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();replot(a.root,a.output)
