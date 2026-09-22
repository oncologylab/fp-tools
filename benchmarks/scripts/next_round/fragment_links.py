"""Augment qualified arrays with unique-fragment and cross-window dependence.

This additional inventory never replaces or rescores completed observations.
It reuses the generated name-collated cache and checks the original pair count.
"""
import argparse
from array import array
from bisect import bisect_left,bisect_right
from collections import defaultdict,Counter
from hashlib import blake2b
from itertools import groupby
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
from fragments import complete_pair,oracle_cut
from access import digest,deny
from prepare import dump


def execute(out,library):
    inventory=json.loads((out/f'provenance/{library}.fragment_inventory.json').read_text())
    cache=out/f'scratch/{library}/name_collated.bam';deny(cache,library+' permitted generated name cache')
    if not cache.resolve().is_relative_to((out/'scratch').resolve()):raise PermissionError('cache outside generated workspace')
    dump(out/f'provenance/{library}.fragment_links.freeze.json',dict(code_sha256=digest(__file__),
        cache_sha256=digest(cache),prior_inventory_sha256=digest(out/f'provenance/{library}.fragment_inventory.json'),
        purpose='additional unique-fragment/window graph inventory, not rerun of completed scientific scoring',
        numerical_counts_replaced=False))
    sites=pd.read_csv(out/'tables/pool.tsv.gz',sep='\t');sites=sites[sites.cell.eq(library.split('_')[0])].reset_index(drop=True)
    intervals={}
    for chrom,group in sites.groupby('TFBS_chr'):
        g=group.sort_values('center');intervals[chrom]=(g.center.to_numpy(int),g.index.to_numpy(int))
    indexes=array('I');hashes=array('Q');lengths=array('I');all_hashes=array('Q');cross=Counter()
    parent=list(range(len(sites)))
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    with pysam.AlignmentFile(cache,'rb') as handle:
        for name,group in groupby(handle.fetch(until_eof=True),key=lambda r:r.query_name):
            by_rg=defaultdict(list)
            for read in group:by_rg[read.get_tag('RG') if read.has_tag('RG') else ''].append(read)
            for rg,reads in by_rg.items():
                pair,status=complete_pair(reads)
                if pair is None:continue
                identity=int.from_bytes(blake2b(('2026\0'+library+'\0'+rg+'\0'+name).encode(),digest_size=8).digest(),'little')
                all_hashes.append(identity)
                matched=set();chromosomes=set()
                for read in pair:
                    chrom=read.reference_name;chromosomes.add(chrom)
                    if chrom not in intervals:continue
                    cut=oracle_cut(read);centers,rows=intervals[chrom]
                    matched.update(int(i) for i in rows[bisect_left(centers,cut-100):bisect_right(centers,cut+100)])
                cross['cross_chromosome_complete_fragments']+=len(chromosomes)>1
                if len({sites.iloc[i].TFBS_chr for i in matched})>1:cross['cross_chromosome_shared_training_windows']+=1
                length=abs(int(pair[0].template_length)) if len(chromosomes)==1 else 0
                for index in matched:indexes.append(index);hashes.append(identity);lengths.append(length)
                if len(matched)>1:
                    first=next(iter(matched))
                    for index in matched:parent[root(index)]=root(first)
    if len(all_hashes)!=inventory['complete_fragments']:raise AssertionError('pair count changed in independent augmentation')
    h=np.asarray(all_hashes,np.uint64);cutoff=int(np.partition(h,9999999)[9999999]) if len(h)>=10000000 else None
    if cutoff is not None and (h<=cutoff).sum()!=10000000:raise AssertionError('depth hash tie')
    index=np.asarray(indexes,np.uint32);source_hash=np.asarray(hashes,np.uint64);length=np.asarray(lengths,np.uint32)
    table=sites[['site_id','TFBS_chr','cell','tf']].copy()
    table['shared_fragment_component']=[library+':'+str(root(i)) for i in range(len(sites))]
    for depth,selected in [('full',np.ones(len(index),bool)),('10M',source_hash<=cutoff if cutoff is not None else np.zeros(len(index),bool))]:
        counts=np.bincount(index[selected],minlength=len(sites))
        sums=np.bincount(index[selected],weights=length[selected],minlength=len(sites))
        table['unique_fragments_'+depth]=counts
        table['mean_fragment_length_'+depth]=np.divide(sums,counts,out=np.full(len(sites),np.nan),where=counts>0)
        with np.load(out/f'arrays/{library}.{depth}.npz',allow_pickle=False) as z:
            cuts=z['plus_observed'].sum(1)+z['minus_observed'].sum(1)
        if ((cuts<counts)|(cuts>2*counts)).any():raise AssertionError('unique-fragment/cut accounting inconsistent')
        table['cut_assignments_'+depth]=cuts
    table.to_csv(out/f'tables/fragment_links.{library}.tsv.gz',sep='\t',index=False)
    dump(out/f'provenance/{library}.fragment_links.json',dict(cross,complete_fragments=len(h),
        sampling_cutoff_uint64=cutoff,seed=2026,fragment_window_assignments=len(index),
        independent_components=int(table.shared_fragment_component.nunique()),
        cross_chromosome_fit_leakage_qualified=not cross['cross_chromosome_shared_training_windows']))
    print('Unique-fragment augmentation completed',library,len(h),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);p.add_argument('--library',required=True)
    a=p.parse_args();execute(a.outdir.resolve(),a.library)
