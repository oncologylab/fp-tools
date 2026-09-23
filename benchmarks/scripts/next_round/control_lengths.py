"""Recover retained control template lengths; audit means without changing fits."""
import argparse
from bisect import bisect_left,bisect_right
from collections import defaultdict
from hashlib import blake2b
from itertools import groupby
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
from fragments import complete_pair,oracle_cut
from verify import read_table
from prepare import dump
from access import digest,deny


def execute(out):
    cache=out/'scratch/NakedDNA_rep1/name_collated.bam';deny(cache,'permitted generated naked-DNA replicate1 cache')
    if not cache.resolve().is_relative_to((out/'scratch').resolve()):raise PermissionError('external cache')
    dump(out/'provenance/control_lengths.freeze.json',dict(code_sha256=digest(__file__),cache_sha256=digest(cache),
        scope='requested descriptive fragment-length dependence; no new candidate, fit, threshold or cohort change',
        length='absolute recorded template length only when both mates agree and map to same chromosome',
        strata='per-window mean template length: <100,100-<200,200-<400,400+,unavailable',
        policy='fixed bins before this length summary; no outcome-dependent selection'))
    sites=read_table(out/'tables/control_sites.tsv.gz');intervals={}
    for chrom,g in sites.groupby('TFBS_chr'):
        g=g.sort_values('center');intervals[chrom]=(g.center.to_numpy(int),g.index.to_numpy(int))
    cutoff=json.loads((out/'provenance/control_fragment_membership.json').read_text())['last_selected_window_hash_uint64']
    counts=np.zeros((len(sites),2),np.int64);valid=counts.copy();sums=counts.astype(float);pairs=0
    with pysam.AlignmentFile(cache,'rb') as bam:
        for name,group in groupby(bam.fetch(until_eof=True),key=lambda r:r.query_name):
            by_rg=defaultdict(list)
            for r in group:by_rg[r.get_tag('RG') if r.has_tag('RG') else ''].append(r)
            for rg,records in by_rg.items():
                pair,status=complete_pair(records)
                if pair is None:continue
                pairs+=1;matched=set()
                for read in pair:
                    if read.reference_name not in intervals:continue
                    center,index=intervals[read.reference_name];cut=oracle_cut(read)
                    matched.update(int(i) for i in index[bisect_left(center,cut-100):bisect_right(center,cut+100)])
                if not matched:continue
                h=int.from_bytes(blake2b(('2026\0NakedDNA_rep1\0'+rg+'\0'+name).encode(),digest_size=8).digest(),'little')
                a,b=pair;length=abs(a.template_length)
                qualified=a.reference_id==b.reference_id and length==abs(b.template_length) and length>0
                for depth in [0,1] if h<=cutoff else [0]:
                    ix=list(matched);counts[ix,depth]+=1
                    if qualified:valid[ix,depth]+=1;sums[ix,depth]+=length
    original=json.loads((out/'provenance/NakedDNA_rep1.fragment_inventory.json').read_text())
    if pairs!=original['complete_fragments']:raise AssertionError('control complete-pair inventory changed')
    previous=read_table(out/'tables/fragment_links.NakedDNA_rep1.tsv.gz')
    table=sites[['site_id']].copy()
    for j,depth in enumerate(['full','10M']):
        np.testing.assert_array_equal(counts[:,j],previous['unique_fragments_'+depth])
        table['mean_template_length_'+depth]=np.divide(sums[:,j],valid[:,j],out=np.full(len(sites),np.nan),where=valid[:,j]>0)
        table['length_qualified_fragments_'+depth]=valid[:,j];table['length_unavailable_fragments_'+depth]=counts[:,j]-valid[:,j]
    table.to_csv(out/'tables/control_fragment_lengths.tsv.gz',sep='\t',index=False)
    rows=[]
    for depth in ['10M','full']:
        control=read_table(out/f'tables/control_counts.{depth}.tsv.gz').merge(table,on='site_id',validate='one_to_one')
        control['length_stratum']=pd.cut(control['mean_template_length_'+depth],[-np.inf,100,200,400,np.inf],right=False,labels=['<100','100-<200','200-<400','400+']).astype(str)
        control['length_stratum']=control.length_stratum.replace('nan','unavailable')
        for (cell,tf),task in control.groupby(['cell','tf']):
            for arm in ['LOG21','DWM','uniform']:
                for strand in ['plus','minus','combined']:
                    for stratum in ['<100','100-<200','200-<400','400+','unavailable']:
                        g=task[task.eligible&task.length_stratum.eq(stratum)]
                        strands=['plus','minus'] if strand=='combined' else [strand]
                        observed=sum(g[s+'_center'].sum() for s in strands)
                        total=sum(g[s+'_total'].sum() for s in strands)
                        expected=sum((g[s+'_total']*g[arm+'_'+s+'_p']).sum() for s in strands)
                        rows.append(dict(depth=depth,cell=cell,tf=tf,arm=arm,strand=strand,length_stratum=stratum,
                            original_sites=len(task),eligible_sites=len(g),regional_cuts=total,observed_center=observed,expected_center=expected,
                            mean_error_fraction=(observed-expected)/total if total else np.nan,scope='examined control, descriptive mean error; not a new calibrated model'))
    pd.DataFrame(rows).to_csv(out/'tables/control_length_calibration.tsv',sep='\t',index=False)
    print('Control length audit completed',pairs,'complete pairs; exact window fragment-count parity',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);execute(p.parse_args().outdir.resolve())
