"""Recover and verify exact window membership for paired control depth arrays.

The global cutoff was not serialized by the first extraction implementation.
The observed-window fragment membership is nevertheless identifiable from
stable hash order and the complete 201-position 10M array. This records that
limited fact, not a fabricated exact global cutoff.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from access import digest
from prepare import dump


def execute(out):
    with np.load(out/'arrays/control_complete_fragments.npz',allow_pickle=False) as z:
        site=z['site_index'];hashes=z['fragment_hash'];ends=z['endpoints']
    with np.load(out/'arrays/NakedDNA_rep1.10M.npz',allow_pickle=False) as z:
        target=np.stack([z['plus_observed'],z['minus_observed']],axis=1);ids=z['site_id']
    contributions=(abs(ends[:,0])<=100).astype(int)+(abs(ends[:,2])<=100).astype(int)
    unique,inverse=np.unique(hashes,return_inverse=True)
    weights=np.bincount(inverse,weights=contributions).astype(int)
    cumulative=np.cumsum(weights);total=int(target.sum());index=int(np.searchsorted(cumulative,total))
    if index>=len(unique) or cumulative[index]!=total:raise AssertionError('target cuts split a complete fragment')
    cutoff=int(unique[index]);membership=hashes<=cutoff
    result=np.zeros(target.shape,np.uint32)
    for end in [0,2]:
        use=membership&(abs(ends[:,end])<=100)
        np.add.at(result,(site[use],ends[use,end+1],ends[use,end]+100),1)
    np.testing.assert_array_equal(result,target)
    table=pd.DataFrame(dict(site_id=ids,unique_fragments_full=np.bincount(site,minlength=len(ids)),
        unique_fragments_10M=np.bincount(site[membership],minlength=len(ids))))
    table.to_csv(out/'tables/fragment_links.NakedDNA_rep1.tsv.gz',sep='\t',index=False)
    np.savez_compressed(out/'arrays/control_fragment_10M_membership.npz',selected_row=membership)
    dump(out/'provenance/control_fragment_membership.json',dict(code_sha256=digest(__file__),
        array_source_sha256=digest(out/'arrays/control_complete_fragments.npz'),
        exact_201_position_parity=True,selected_window_fragments=int(len(np.unique(hashes[membership]))),
        last_selected_window_hash_uint64=cutoff,
        first_excluded_window_hash_uint64=int(unique[index+1]) if index+1<len(unique) else None,
        exact_global_sampling_cutoff_known=False,
        scope='exact membership of fragments overlapping the4200 windows; global cutoff only bounded between adjacent observed hashes'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);a=p.parse_args();execute(a.outdir.resolve())
