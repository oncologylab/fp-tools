"""Follow-up support audit: sparse bigWig entries are not output support.

Preserves the overly strict initial finite-position audit. Zeros are restored
ONLY inside an explicit native processed-region contract, never outside it.
"""
import argparse
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pysam
import pyBigWig
from access import Access,deny,digest
from prepare import ROOT,DATA,dump


def regions(bed,lengths):
    table=pd.read_csv(bed,sep='\t',header=None,usecols=[0,1,2],names=['chrom','start','end'])
    result={}
    for chrom,group in table.groupby('chrom'):
        if chrom not in lengths:continue
        intervals=[]
        for row in group.sort_values('start').itertuples():
            start,end=max(0,int(row.start)-162),min(lengths[chrom],int(row.end)+162)
            if intervals and start<=intervals[-1][1]:intervals[-1][1]=max(intervals[-1][1],end)
            else:intervals.append([start,end])
        result[chrom]=[(a+62,b-62) for a,b in intervals]
    return result


def execute(out):
    access=Access(out/'provenance/dwm_reference_access.local.json',out/'logs/access.jsonl')
    oldbed=DATA/'benchmarks/results/footprint_functional_v1/naked_dna_10m_s2026/ctcf_exact_geometry_safety/ctcf_naked_sites.bed'
    deny(oldbed,'naked DNA replicate1 historical label-free geometry BED')
    with pysam.AlignmentFile(access.path('NakedDNA_rep1_bam'),'rb') as bam:lengths=dict(zip(bam.references,bam.lengths))
    oldregions=regions(oldbed,lengths);newregions=regions(out/'models/all_control_sites.bed',lengths)
    dump(out/'provenance/dwm_sparse_verification.v3.freeze.json',dict(code_sha256=digest(__file__),oldbed_sha256=digest(oldbed),
        change='initial finite-position and subsequent mixed-unit audits retained; restore sparse zeros only on processed support and compare raw cuts after recorded normalization; no model/threshold/cohort changes',
        sparse_writer='utilities.bigwig_writer emits only np.flatnonzero(signal)',
        output_contract='expand motif intervals162, merge, chromosomeclip, shrink62; matches fixed100 extension +12 flank+50 halfwindow',
        old_contract_status='inferred historical region contract, not sufficient on its own to qualify archived output provenance'))
    sites=pd.read_csv(out/'tables/control_sites.tsv.gz',sep='\t')
    with np.load(out/'arrays/control_production_DWM.npz',allow_pickle=False) as z:panels={k:z[k].copy() for k in z.files if k!='site_id'}
    sys.path.insert(0,str(ROOT/'src'))
    from fp_tools.tools.atacorrect_functions import AtacBias
    old_model=AtacBias().from_pickle(str(access.path('control_old_AtacBias_pickle')))
    correction_factor=float(old_model.correction_factor)
    # The archived writer stored normalized uncorrected cuts; the new run has
    # norm-off. Compare raw-cut units, retaining both earlier audits unchanged.
    panels['old_observed']/=correction_factor
    rows=[]
    for i,row in enumerate(sites.itertuples()):
        start,end=int(row.center)-100,int(row.center)+101
        old_support=any(a<=start and end<=b for a,b in oldregions.get(row.TFBS_chr,[]))
        new_support=any(a<=start and end<=b for a,b in newregions.get(row.TFBS_chr,[]))
        r=dict(site_id=row.site_id,cell=row.cell,tf=row.tf,old_region_supported=old_support,new_region_supported=new_support)
        for key in panels:
            permitted=old_support if key.startswith('old_') else new_support
            r[key+'_recorded_bases']=int(np.isfinite(panels[key][i]).sum())
            if permitted:panels[key][i]=np.nan_to_num(panels[key][i])
            else:panels[key][i]=np.nan
            r[key+'_sum']=float(panels[key][i].sum())
        overlap=old_support and new_support and r['old_expected_sum']>0 and r['original_expected_sum']>0
        r['old_positive_overlap']=overlap
        r['expected_proportion_parity']=False;r['observed_parity']=False;r['max_expected_proportion_difference']=np.nan
        if overlap:
            a=panels['old_expected'][i]/r['old_expected_sum'];b=panels['original_expected'][i]/r['original_expected_sum']
            r['expected_proportion_parity']=bool(np.allclose(a,b,rtol=1e-5,atol=1e-7))
            r['max_expected_proportion_difference']=float(abs(a-b).max())
            r['observed_parity']=bool(np.allclose(panels['old_observed'][i],panels['original_observed'][i],rtol=1e-5,atol=1e-7))
        rows.append(r)
    table=pd.DataFrame(rows);table.to_csv(out/'tables/dwm_support.v3.tsv',sep='\t',index=False)
    np.savez_compressed(out/'arrays/control_production_DWM.v3.npz',**panels,site_id=sites.site_id.to_numpy(str))
    overlap=table[table.old_positive_overlap]
    dump(out/'provenance/dwm_reference_qualification.v3.json',dict(sites=len(table),processed_new=int(table.new_region_supported.sum()),
        original_correction_factor=correction_factor,observed_units='raw cuts after dividing archived uncorrected values by their recorded correction factor',
        old_overlap=len(overlap),expected_shape_mismatches=int((~overlap.expected_proportion_parity).sum()),
        observed_mismatches=int((~overlap.observed_parity).sum()),positive_full_expectations=int((table.full_expected_sum>0).sum()),
        valid_zero_full_expectations=int((table.new_region_supported & table.full_expected_sum.eq(0)).sum()),
        qualified_for_safety=False,reason='combined-track audit only; historic context assumptions and any shape discrepancies remain explicit',
        native_strand_expectation_qualification='INCOMPLETE'))
    print('DWM output support audited',len(table),'sites;',len(overlap),'historical positive overlap;',int((~overlap.expected_proportion_parity).sum()),'shape mismatches',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);a=p.parse_args();execute(a.outdir.resolve())
