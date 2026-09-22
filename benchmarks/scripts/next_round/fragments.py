"""Exact name/read-group pair accounting and pair-preserving regional extraction.

Unrelated alignments are never read as labels. Only training-chromosome windows
are extracted; all chromosome records are counted to qualify whole-library depth.
Temporary collated BAMs are generated outputs and are retained until packaging.
"""
import argparse
from array import array
from bisect import bisect_left,bisect_right
from collections import Counter,defaultdict
from hashlib import blake2b
from itertools import groupby
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import pysam
from access import Access,digest
from prepare import dump,ROOT


def usable(read):
    return not (read.is_unmapped or read.is_secondary or read.is_supplementary or
                read.is_qcfail or read.is_duplicate or read.mapping_quality<30)


def oracle_cut(read):
    """Independent CIGAR implementation of the established soft-clip convention."""
    cigar=read.cigartuples or []
    if read.is_reverse:
        trailing=0
        for operation,length in reversed(cigar):
            if operation==5: continue
            if operation==4: trailing=length
            break
        consumed=sum(length for operation,length in cigar if operation in (0,2,3,7,8))
        return int(read.reference_start)+consumed+trailing-5
    leading=0
    for operation,length in cigar:
        if operation==5: continue
        if operation==4: leading=length
        break
    return int(read.reference_start)-leading+4


def complete_pair(reads):
    good=[r for r in reads if usable(r)]
    if len(good)!=2:
        return None, 'singleton' if len(good)==1 else ('ambiguous' if len(good)>2 else 'no_retained_records')
    a,b=good
    if not (a.is_paired and b.is_paired and a.is_read1 != b.is_read1 and a.is_read2 != b.is_read2
            and a.is_read1 != a.is_read2 and b.is_read1 != b.is_read2):
        return None,'ambiguous'
    if (a.mate_is_unmapped or b.mate_is_unmapped or a.next_reference_id!=b.reference_id
            or b.next_reference_id!=a.reference_id or a.next_reference_start!=b.reference_start
            or b.next_reference_start!=a.reference_start):
        return None,'mate_metadata_mismatch'
    return (a,b),'complete'


def extract(out,library,threads):
    access=Access(out/'provenance/input_access.local.json',out/'logs/access.jsonl')
    bam=access.path(library+'_bam')
    stage=out/'provenance'/f'{library}.fragments.freeze.json'
    if stage.exists(): raise ValueError('attempt exists; preserve outputs and use explicit verification')
    dump(stage,dict(library=library,code_sha256=digest(__file__),input_sha256=digest(bam),
        grouping='library/readgroup/queryname; exactly one retained read1 and read2, reciprocal mate coordinates',
        retain_filter='mapped primary nonsupplementary nonQCfail nonduplicate MAPQ>=30',
        depth='complete retained fragments across all chromosomes; no TLEN-sign representative',
        partial_pairs='excluded from new complete-pair arrays, separately counted',
        subsampling='exact lowest 10000000 blake2b-64(seed2026,library,readgroup,queryname); preserve mates',
        regions='permitted chr1-15; exact201-position motif windows',threads=threads,
        duplicate_query_groups='ambiguous groups excluded; no arbitrary representative',
        created=time.time()))
    start=time.monotonic()
    scratch=out/'scratch'/library;scratch.mkdir(parents=True)
    collated=scratch/'name_collated.bam'
    print('Collating complete fragments:',library,flush=True)
    pysam.collate('-@',str(threads),'-o',str(collated),str(bam),str(scratch/'collate'))
    sites=pd.read_csv(out/'tables'/('control_sites.tsv.gz' if library=='NakedDNA_rep1' else 'pool.tsv.gz'),sep='\t')
    if library!='NakedDNA_rep1': sites=sites[sites.cell.eq(library.split('_')[0])].reset_index(drop=True)
    if not sites.TFBS_chr.isin([f'chr{i}' for i in range(1,16)]).all(): raise ValueError('nontraining extraction')
    intervals={}
    for chrom,group in sites.groupby('TFBS_chr'):
        group=group.sort_values('center')
        intervals[chrom]=(group.center.to_numpy(int),group.index.to_numpy(int))
    # Sparse cut assignments retain fragment IDs. No independent end resampling.
    indices=array('I');strands=array('B');offsets=array('H');hashes=array('Q')
    empirical=[];summary=Counter();fragment_hashes=array('Q')
    sys.path.insert(0,str(ROOT/'src'))
    from fp_tools.tools.parametric_bias import cut_position_from_alignment
    with pysam.AlignmentFile(collated,'rb') as handle:
        lengths=dict(zip(handle.references,handle.lengths))
        for query,group in groupby(handle.fetch(until_eof=True),key=lambda r:r.query_name):
            by_rg=defaultdict(list)
            for read in group:
                summary['records']+=1
                summary['mapped_records']+=int(not read.is_unmapped)
                summary['excluded_records']+=int(not usable(read))
                by_rg[read.get_tag('RG') if read.has_tag('RG') else ''].append(read)
            for rg,reads in by_rg.items():
                pair,status=complete_pair(reads);summary[status+'_groups']+=1
                if pair is None:
                    summary['unpaired_retained_records']+=sum(usable(r) for r in reads)
                    continue
                identity=int.from_bytes(blake2b(('2026\0'+library+'\0'+rg+'\0'+query).encode(),digest_size=8).digest(),'little')
                fragment_hashes.append(identity)
                cuts=[]
                for read in pair:
                    cut=oracle_cut(read)
                    if cut!=cut_position_from_alignment(read,(4,-5)):
                        raise AssertionError('independent CIGAR coordinate oracle differs')
                    chrom=read.reference_name;strand=int(read.is_reverse)
                    cuts.append((chrom,cut,strand))
                    if chrom not in intervals: continue
                    centers,rows=intervals[chrom]
                    lo,hi=bisect_left(centers,cut-100),bisect_right(centers,cut+100)
                    for j in range(lo,hi):
                        index=int(rows[j]);center=int(centers[j])
                        if center-100<0 or center+101>lengths[chrom]:continue
                        offset=cut-center+100
                        reverse=sites.iloc[index].TFBS_strand=='-'
                        indices.append(index);strands.append(1-strand if reverse else strand)
                        offsets.append(200-offset if reverse else offset);hashes.append(identity)
                if library=='NakedDNA_rep1':
                    # Preserve both endpoints even when only one overlaps a
                    # window; deduplicate later by physical window, not TF.
                    matched=set()
                    for chrom,cut,strand in cuts:
                        if chrom not in intervals:continue
                        centers,rows=intervals[chrom]
                        for j in range(bisect_left(centers,cut-100),bisect_right(centers,cut+100)):
                            matched.add(int(rows[j]))
                    for index in matched:
                        row=sites.iloc[index];center=int(row.center);reverse=row.TFBS_strand=='-'
                        encoded=[]
                        for chrom,cut,strand in cuts:
                            offset=cut-center if chrom==row.TFBS_chr else 2**30
                            encoded.extend([-offset if reverse else offset,1-strand if reverse else strand])
                        empirical.append([index,identity,*encoded])
                if summary['complete_groups']%5000000==0:
                    print(library,'complete fragments',summary['complete_groups'],flush=True)
    total=len(fragment_hashes)
    h=np.asarray(fragment_hashes,dtype=np.uint64)
    assigned_h=np.asarray(hashes,dtype=np.uint64)
    i=np.asarray(indices,dtype=np.uint32);s=np.asarray(strands,dtype=np.uint8);x=np.asarray(offsets,dtype=np.uint16)
    thresholds={'full':None}
    if total>=10000000:
        cutoff=int(np.partition(h,9999999)[9999999])
        if int((h<=cutoff).sum())!=10000000:raise ValueError('hash collision at exact depth boundary')
        thresholds['10M']=cutoff
    summary['complete_fragments']=total
    summary['unique_fragment_definition']='exact one retained read1/read2 per library/RG/queryname, ambiguous groups excluded'
    summary['motif_window_cut_assignments']=len(indices)
    summary['coordinate_oracle_comparisons']=2*total
    summary['elapsed_seconds']=time.monotonic()-start
    summary['new_depths_available']=list(thresholds)
    for depth,cutoff in thresholds.items():
        select=np.ones(len(i),bool) if cutoff is None else assigned_h<=cutoff
        counts=np.zeros((len(sites),2,201),np.uint32)
        np.add.at(counts,(i[select],s[select],x[select]),1)
        np.savez_compressed(out/f'arrays/{library}.{depth}.npz',plus_observed=counts[:,0],minus_observed=counts[:,1],site_id=sites.site_id.to_numpy(str),
                            achieved_fragments=total if cutoff is None else 10000000)
    if empirical:
        # Separate integer fields so fragment hashes cannot lose precision.
        np.savez_compressed(out/'arrays/control_complete_fragments.npz',
                            site_index=np.asarray([r[0] for r in empirical],np.uint32),
                            fragment_hash=np.asarray([r[1] for r in empirical],np.uint64),
                            endpoints=np.asarray([r[2:] for r in empirical],np.int64))
    dump(out/f'provenance/{library}.fragment_inventory.json',dict(summary))
    print('Completed pair-accounted extraction',library,total,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--library',required=True,choices=['NakedDNA_rep1']+[c+'_rep'+str(r) for c in ['HepG2','K562'] for r in [1,2,3]])
    p.add_argument('--threads',type=int,default=4)
    a=p.parse_args();extract(a.outdir.resolve(),a.library,a.threads)
