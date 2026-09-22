"""Repair all-site production DWM output support without declaring parity."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
import pyBigWig
from access import Access,digest,deny
from prepare import ROOT,DATA,dump


def prepare(out):
    access=Access(out/'provenance/input_access.local.json',out/'logs/access.jsonl')
    sites=pd.read_csv(out/'tables/control_sites.tsv.gz',sep='\t')
    old_bam=DATA/'benchmarks/results/footprint_functional_v1/naked_dna_10m_s2026/NakedDNA_rep1.first10m.bam'
    deny(old_bam,'naked DNA replicate1 original first10M preparation')
    registry=json.loads((out/'provenance/input_access.local.json').read_text())
    for key,path in [('control_original_bam',old_bam),('control_original_bai',Path(str(old_bam)+'.bai'))]:
        registry['inputs'][key]=dict(path=str(path),resolved_path=str(path.resolve()),sha256=digest(path),bytes=path.stat().st_size,
                                    identity='naked DNA replicate1 first10M original preparation; historical overlap parity only')
    dump(out/'provenance/dwm_reference_access.local.json',registry)
    bed=out/'models/all_control_sites.bed'
    sites[['TFBS_chr','TFBS_start','TFBS_end']].drop_duplicates().sort_values(['TFBS_chr','TFBS_start','TFBS_end']).to_csv(bed,sep='\t',header=False,index=False)
    genome=access.path('genome');model=access.path('control_old_AtacBias_pickle')
    jobs=[]
    for depth,bam in [('original',old_bam),('full',access.path('NakedDNA_rep1_bam'))]:
        command=[sys.executable,'-B','-m','fp_tools.cli','--bams',str(bam),'--genome',str(genome),
            '--peaks',str(bed),'--regions-out',str(bed),'--bias-pkl',str(model),'--norm-off',
            '--write-tracks','uncorrected','expected','--skip-qc','--scale-corrected','none',
            '--extend','100','--window','100','--k_flank','12','--read_shift','4','-5',
            '--score_mat','DWM','--cores','4','--split','100','--outdir',str(out/f'scratch/dwm_{depth}'),
            '--prefix','NakedDNA_rep1_all_sites']
        jobs.append(dict(depth=depth,command=command))
    dump(out/'provenance/dwm_reference.freeze.json',dict(code_sha256=digest(__file__),jobs=jobs,
        context_contract='production region merge/extend100, native k_flank12, rolling correction window100 and step10; output layout can change fitted expectations',
        qualification='all4200 rows, explicit support and expected sum; old overlap observed and normalized expected rtol1e-5 atol1e-7',
        layout_parity_required=True,all_site_reference_is_not_safety=True))


def run(out):
    config=json.loads((out/'provenance/dwm_reference.freeze.json').read_text())
    for job in config['jobs']:
        print('Production DWM support construction',job['depth'],flush=True)
        subprocess.run(job['command'],cwd=ROOT,env={**__import__('os').environ,'PYTHONPATH':str(ROOT/'src')},check=True)


def verify(out):
    sites=pd.read_csv(out/'tables/control_sites.tsv.gz',sep='\t')
    access=Access(out/'provenance/dwm_reference_access.local.json',out/'logs/access.jsonl')
    paths={'old_expected':access.path('control_old_expected_bw'),'old_observed':access.path('control_old_uncorrected_bw')}
    for depth in ['original','full']:
        for kind,suffix in [('expected','expected'),('observed','uncorrected')]:
            paths[depth+'_'+kind]=out/f'scratch/dwm_{depth}/NakedDNA_rep1_all_sites_{suffix}.bw'
    panels={}
    for key,path in paths.items():
        values=np.full((len(sites),201),np.nan)
        with pyBigWig.open(str(path)) as track:
            for i,row in enumerate(sites.itertuples()):
                if row.TFBS_chr not in track.chroms() or row.center<100:continue
                value=np.asarray(track.values(row.TFBS_chr,row.center-100,row.center+101))
                values[i]=value[::-1] if row.TFBS_strand=='-' else value
        panels[key]=values
    table=sites[['site_id','cell','tf','TFBS_chr']].copy()
    for key,values in panels.items():
        table[key+'_finite']=np.isfinite(values).all(1)
        table[key+'_sum']=np.nansum(values,axis=1)
    overlap=table.old_expected_finite & (table.old_expected_sum>0)
    old=panels['old_expected'][overlap];new=panels['original_expected'][overlap]
    old=old/old.sum(1,keepdims=True);new=new/new.sum(1,keepdims=True)
    expected_match=np.isclose(old,new,rtol=1e-5,atol=1e-7).all(1)
    observed_match=np.isclose(panels['old_observed'][overlap],panels['original_observed'][overlap],rtol=1e-5,atol=1e-7).all(1)
    table['historical_overlap']=overlap;table['expected_proportion_parity']=False;table['observed_parity']=False
    table.loc[overlap,'expected_proportion_parity']=expected_match;table.loc[overlap,'observed_parity']=observed_match
    table.to_csv(out/'tables/dwm_support.tsv',sep='\t',index=False)
    np.savez_compressed(out/'arrays/control_production_DWM.npz',**panels,site_id=sites.site_id.to_numpy(str))
    dump(out/'provenance/dwm_reference_qualification.json',dict(sites=len(sites),historical_overlap=int(overlap.sum()),
        expected_parity_failures=int((~expected_match).sum()),observed_parity_failures=int((~observed_match).sum()),
        all4200_finite=bool(table.full_expected_finite.all()),positive_expected_sites=int((table.full_expected_sum>0).sum()),
        combined_overlap_qualified=bool(overlap.any() and expected_match.all() and observed_match.all() and table.full_expected_finite.all()),
        qualified_for_safety=False,
        strand_parity='INCOMPLETE: combined production tracks do not qualify separate strand expectations',
        unresolved='output-grid phase and context effects require discrepancy explanation; valid zero expectations still need calling support checks',
        production_counts_are_not_complete_pair_arrays=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','verify']);p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args();globals()[a.stage](a.outdir.resolve())
