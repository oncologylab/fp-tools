"""Identify permitted originals, freeze policy, and prepare replayable arrays.

Only exact paths selected from already-authorized manifests are read. No
recursive data discovery, protected controls, or historic preflight chains.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
from access import Access, DENY, deny, digest
from protocol import POLICY, PROMOTION

ROOT=Path(__file__).resolve().parents[3]
OLD=ROOT.parent/'fp-tools-research'
DATA=ROOT.parent/'fp-tools'
R=OLD/'benchmarks/results'


def dump(path,value):
    with Path(path).open('x') as handle:
        json.dump(value,handle,indent=2,sort_keys=True,allow_nan=False)
        handle.write('\n')


def freeze(out):
    if out.exists(): raise ValueError('fresh attempt required')
    paths={}
    def add(key,path,identity,sha=None):
        path=Path(path);deny(path,identity)
        paths[key]=dict(path=str(path),resolved_path=str(path.resolve()),identity=identity,expected_sha256=sha)
    for key,rel in [('pool','radius_full_pool_20260922_v1/all_sites.tsv'),
                    ('replay_sites','radius_full_pool_comparison_20260922_v1/sites.tsv'),
                    ('replay_panels','radius_full_pool_comparison_20260922_v1/panels.npz'),
                    ('tasks','radius_full_pool_20260922_v1/diagnostics.tsv'),
                    ('control_regions','full_rep1_region_calibration_20260922_v1/calibration/sites.tsv')]:
        add(key,R/rel,'training chr1-15; '+key)
    for cell in ['HepG2','K562']:
        manifest=json.loads((R/f'wide_training_profiles_20260922_v1/{cell}.json').read_text())
        add(cell+'_wide',manifest['outputs']['profiles']['path'],cell+' rep1 training chr1-15 observed cuts and frozen LOG21 propensities',manifest['outputs']['profiles']['sha256'])
        manifest=json.loads((R/'full_rep1_region_calibration_20260922_v1/profiles_manifest.json').read_text())
        add(cell+'_control',manifest['outputs'][cell]['path'],'naked DNA replicate1 '+cell+' locus panel',manifest['outputs'][cell]['sha256'])
    add('resource_protocol',OLD/'benchmarks/manifests/end_to_end_resource_protocol_v1.json','immutable resource protocol')
    model=DATA/'benchmarks/results/frozen_parametric_factorization_v1_run2/enzyme_control/mitochondrial/ensembles/pooled.mitochondrial.cross_cell_pooled.loglinear21.shift_4_-5.l2_0p0001.depth_full.seed_ensemble.npz'
    add('LOG21',model,'mitochondrial five-seed bias model',POLICY['bias']['sha256'])
    add('genome',DATA/'data/public/raw/genome/hg38.fa','hg38 sequence', '5be01555d98347fdb3714dc84c6f77c9d8bc774adcf32c6f7a8fa06f5baf5e51')
    add('genome_index',DATA/'data/public/raw/genome/hg38.fa.fai','hg38 sequence index')
    full=DATA/'benchmarks/results/frozen_parametric_factorization_v1_run2/bias_shrinkage/mito_log21_l2_1e-4_shift_4_-5/biological_replicates_full/depth_inputs.tsv'
    add('full_registry',full,'K562/HepG2 full-depth input registry')
    for row in pd.read_csv(full,sep='\t').itertuples():
        for suffix in ['bam','bai','dwm_expected']:
            add(row.sample+'_'+suffix,DATA/getattr(row,suffix),row.sample+' full development input '+suffix,getattr(row,suffix+'_sha256'))
    base=DATA/'data/public/processed/naked_dna_gse164997/samples/NakedDNA_rep1/alignment/NakedDNA_rep1.filtered.bam'
    add('NakedDNA_rep1_bam',base,'GSE164997 SRR13450940 naked DNA replicate1','a04cb6abc7bd06f00f0c813e31ebba774e24481d2d7075cba036e6e2638afd66')
    add('NakedDNA_rep1_bai',str(base)+'.bai','GSE164997 SRR13450940 naked DNA replicate1 index')
    dwm=DATA/'benchmarks/results/footprint_functional_v1/naked_dna_10m_s2026/ctcf_exact_geometry_safety/dwm'
    for suffix in ['AtacBias.pickle','expected.bw','uncorrected.bw']:
        add('control_old_'+suffix.replace('.','_'),dwm/('NakedDNA_rep1_'+suffix),'naked DNA replicate1 historical limited DWM support')
    def check(item):
        key,entry=item
        actual=digest(entry['path'])
        if entry['expected_sha256'] and actual != entry['expected_sha256']:
            raise ValueError('source identity changed '+key)
        entry['sha256']=actual;entry['bytes']=Path(entry['path']).stat().st_size
        print('Verified input',key,flush=True)
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(check,paths.items()))
    for directory in ['arrays','tables','figures','models','provenance','logs','experiments','historical_context','code']:
        (out/directory).mkdir(parents=True,exist_ok=True)
    dump(out/'provenance/input_access.local.json',dict(inputs=paths,protected_identity_denylist=DENY,
         denied='all nonallowlisted scientific inputs, all protected dataset aliases and resolved paths',
         replicate2_contents_allowed=False,final_controls_allowed=False))
    code={str(p.relative_to(ROOT)):digest(p) for p in sorted(Path(__file__).parent.glob('*.py'))}
    for name in ['conditional_protection_score.py','caliper_matching.py','fixed_shape_normalization.py','build_strand_functional_profiles.py']:
        p=ROOT/'benchmarks/scripts'/name;code[str(p.relative_to(ROOT))]=digest(p)
    freeze=dict(policy=POLICY,source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                created_utc=datetime.now(timezone.utc).isoformat(),initial_code=code,
                input_identities={k:{x:v[x] for x in ['identity','sha256','bytes']} for k,v in paths.items()},
                execution_commands=['python prepare.py arrays --outdir OUTPUT','python fragments.py --outdir OUTPUT --library LIBRARY',
                                    'python run.py --outdir OUTPUT','python verify.py OUTPUT','python package.py OUTPUT'],
                implementation_freezes='additional implementation modules must be hashed and fixtures passed before their first scientific scoring; no hypothesis, threshold or cohort changes')
    dump(out/'EXPERIMENT_PLAN.freeze.json',freeze)
    (out/'ORIGINAL_PROMOTION_REQUIREMENTS.md').write_text(PROMOTION)
    (out/'BASELINE_STATE.md').write_text('''# Baseline at the start of this round

Original research source: 0569e2f59def36d04c990e0494873ec396c3d3c4. Main is read-only.
The 22 September review packet and its separate original-source identities were checked.
The full-replicate-1 region calibration and count-screen stages are already complete;
their extraction is reused after checksums, not rerun. Previous full-depth controls
are not a controlled depth comparison with the old first-10M preparation.
The original training candidate pool exists (35,006 rows), and the frozen replay
has 8,106 matched rows. The replay remains retrospective and nonnested.
The current corrected-extraction evidence supersedes older historical model comparisons.
Historical ChromBPNet and 96-run results are not independently requalified here.
All DWM comparisons mean fp-tools expected-profile/residual comparators, not TOBIAS runs.

This round adds signed zero-centred shrinkage, a controlled overdispersion comparison,
training-chromosome held-out diagnostics, fragment accounting and power challenges.
All 21 tasks remain in ledgers. Full control DWM support, complete-fragment depth,
full biological replication, functional uncertainty and output/resource contracts
must be qualified separately. No locked controls or new labels will be opened.
''')
    print('Frozen before new scientific scoring:',out,flush=True)


def arrays(out):
    access=Access(out/'provenance/input_access.local.json',out/'logs/access.jsonl')
    sites=pd.read_csv(access.path('pool'),sep='\t')
    if not sites.TFBS_chr.isin(POLICY['training_chromosomes']).all(): raise ValueError('nontraining sites')
    replay=pd.read_csv(access.path('replay_sites'),sep='\t')
    tasks=pd.read_csv(access.path('tasks'),sep='\t')
    if len(tasks)!=21 or len(replay)!=8106: raise ValueError('denominator changed')
    # Stable identity includes task, chromosome, coordinates and motif strand;
    # locus identity omits task and strand, exposing duplicated windows.
    from hashlib import sha256
    def identify(table):
        table=table.copy()
        table['site_id']=[sha256('|'.join(map(str,r)).encode()).hexdigest()[:24] for r in
            table[['cell','tf','TFBS_chr','TFBS_start','TFBS_end','TFBS_strand']].itertuples(index=False,name=None)]
        table['center']=(table.TFBS_start+table.TFBS_end)//2
        table['parent_locus_id']=table.TFBS_chr+':'+table.center.astype(str)
        table['fold']=table.TFBS_chr
        table['overlap_group']=''
        for chrom,group in table.groupby('TFBS_chr'):
            right=-1;group_id=-1
            for index,row in group.sort_values(['center','site_id']).iterrows():
                if row.center-100>right: group_id+=1
                right=max(right,row.center+100)
                table.loc[index,'overlap_group']=f'{chrom}:{group_id}'
        return table
    sites=identify(sites);replay=identify(replay)
    sites.to_csv(out/'tables/pool.tsv.gz',sep='\t',index=False)
    replay.to_csv(out/'tables/replay_sites.tsv.gz',sep='\t',index=False)
    tasks.to_csv(out/'tables/task_ledger.tsv',sep='\t',index=False)
    for cell in ['HepG2','K562']:
        subset=sites[sites.cell.eq(cell)]
        with np.load(access.path(cell+'_wide'),allow_pickle=False) as z:
            indices=subset.wide_index.to_numpy(int)
            payload={key:z[key][indices,500:701].copy() for key in
                     ['plus_observed','minus_observed','plus_log_bias','minus_log_bias']}
            payload['valid']=z['bam_valid'][indices]
            payload['site_hash']=z['site_hash'][indices]
        reverse=subset.TFBS_strand.eq('-').to_numpy()
        for suffix in ['observed','log_bias']:
            p,m=payload['plus_'+suffix],payload['minus_'+suffix]
            saved=p[reverse].copy();p[reverse]=m[reverse,::-1];m[reverse]=saved[:,::-1]
        payload['site_id']=subset.site_id.to_numpy(str)
        np.savez_compressed(out/f'arrays/{cell}.replay_pool.npz',**payload)
        print('Packaged complete permitted training pool',cell,len(subset),flush=True)
    with np.load(access.path('replay_panels'),allow_pickle=False) as z:
        np.savez_compressed(out/'arrays/archived_replay.npz',raw=z['raw'],DWM_residual=z['DWM_expected'],LOG21_residual=z['parametric_direct'])
    controls=pd.read_csv(access.path('control_regions'),sep='\t')
    controls.to_csv(out/'tables/archived_control_regions.tsv.gz',sep='\t',index=False)
    keys=['cell','tf','TFBS_chr','TFBS_start','TFBS_end','TFBS_strand']
    control_sites=identify(controls[keys].drop_duplicates().reset_index(drop=True))
    control_sites.to_csv(out/'tables/control_sites.tsv.gz',sep='\t',index=False)
    if len(control_sites)!=4200: raise ValueError('control denominator changed')
    for cell in ['HepG2','K562']:
        with np.load(access.path(cell+'_control'),allow_pickle=False) as z:
            np.savez_compressed(out/f'arrays/{cell}.control_full.npz',**{k:z[k] for k in z.files})
    shutil.copyfile(access.path('LOG21'),out/'models/frozen_LOG21.npz')
    shutil.copyfile(access.path('resource_protocol'),out/'provenance/end_to_end_resource_protocol_v1.json')
    dump(out/'provenance/arrays_prepared.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),
        primary_scoring=False,protected_content_read=False,pool_sites=len(sites),replay_sites=len(replay),control_sites=4200))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','arrays']);p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args();globals()[a.stage](a.outdir.resolve())
