#!/usr/bin/env python3
"""Freeze and construct an all-site replicate-2 DWM reference; no rescreening."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
import audit_control_track_scope as audit
from audit_direct_bias_safety import record,write_immutable
from audit_chrombpnet_comparison import canonical_id,checked_path
from fp_tools.tools.atacorrect_functions import AtacBias
import fp_tools.cli

ROOT=audit.ROOT
DATA=audit.screen.DATA_ROOT
SOURCE=ROOT/'benchmarks/results/control_track_scope_20260922_v1'
PROFILE_CONFIG=ROOT/'benchmarks/results/corrected_safety_profiles_20260922_v1/configuration.freeze.json'
MODEL=DATA/'benchmarks/results/footprint_functional_v1/naked_dna_10m_rep2_s2026/ctcf_exact_geometry_safety/dwm/NakedDNA_rep2_AtacBias.pickle'
GENOME=DATA/'data/public/raw/genome/hg38.fa'
PROTOCOL=dict(control='previously_examined_NakedDNA_rep2',training_refit=False,
    output_regions='unique union of all label-free control motif intervals; native extension 100',
    matrix='existing DWM length 25',read_shift=[4,-5],window=100,k_flank=12,
    normalization='off; overlap qualification compares expected proportions and raw observed cuts',
    cores=4,split=100,blacklist=None,drop_chromosomes=['chrM','chrMT','M','MT','Mito'],
    candidate_scoring=False,thresholds_changed=False,promotion_eligible=False,
    overlap_required='expected proportions rtol 1e-5 atol 1e-7; report every drift before any rescreen')


def site_regions(tables):
    for table in tables:
        if any('label' in c.lower() or 'chip' in c.lower() for c in table.columns):
            raise ValueError('control regions must be label-free')
    result=pd.concat([t[['TFBS_chr','TFBS_start','TFBS_end']] for t in tables],ignore_index=True).drop_duplicates()
    if (result.TFBS_start<0).any() or (result.TFBS_end<=result.TFBS_start).any():raise ValueError('invalid motif region')
    return result.sort_values(['TFBS_chr','TFBS_start','TFBS_end']).reset_index(drop=True)


def freeze(out):
    if out.exists():raise ValueError('fresh directory required')
    audit.execute(SOURCE,verify=True)
    old=json.loads((SOURCE/'configuration.freeze.json').read_text());profile=json.loads(PROFILE_CONFIG.read_text())
    items=[i for i in profile['items'] if i['role']=='naked']
    if len(items)!=2 or {i['cell'] for i in items}!={'HepG2','K562'} or len({i['bam']['sha256'] for i in items})!=1:
        raise ValueError('control source identities changed')
    tables=[pd.read_csv(checked_path(i['sites'],ROOT),sep='\t') for i in items]
    if sum(map(len,tables))!=4200:raise ValueError('site denominator changed')
    regions=site_regions(tables)
    inputs=old['inputs']+[record(SOURCE/'configuration.freeze.json'),record(SOURCE/'coverage_manifest.json'),record(PROFILE_CONFIG),
        record(MODEL),record(GENOME),record(Path(str(GENOME)+'.fai'))]
    inputs += [i[k] for i in items for k in ['bam','bai','sites']]
    for entry in inputs:checked_path(entry,ROOT)
    model=AtacBias().from_pickle(str(MODEL))
    if model.stype!='DWM' or any(m.length!=25 for m in model.bias.values()):raise ValueError('bias model operator changed')
    paths={e['path'] for e in old['code']};paths.add(str(Path(__file__).resolve()))
    # Freeze the complete source/native implementation invoked by the CLI.
    paths.update(str(p.resolve()) for p in (ROOT/'src/fp_tools').rglob('*') if p.suffix in {'.py','.pyx','.so'})
    out.mkdir(parents=True);regions.to_csv(out/'all_control_sites.bed',sep='\t',index=False,header=False)
    command=[sys.executable,'-m','fp_tools.cli','--bams',items[0]['bam']['path'],
        '--genome',str(GENOME),'--peaks',old['bed']['path'],'--regions-out',str(out/'all_control_sites.bed'),
        '--bias-pkl',str(MODEL),'--norm-off','--write-tracks','uncorrected','expected','--skip-qc',
        '--scale-corrected','none','--extend','100','--window','100','--k_flank','12','--read_shift','4','-5',
        '--score_mat','DWM','--cores','4','--split','100','--outdir',str(out/'tracks'),'--prefix','NakedDNA_rep2_all_sites']
    value=dict(protocol=PROTOCOL,inputs=inputs+[record(out/'all_control_sites.bed')],code=[record(Path(p)) for p in sorted(paths)],
        command=command,items=items,original_model=record(MODEL),original_correction_factor=float(model.correction_factor),
        site_rows=4200,unique_regions=len(regions),runtime=audit.screen.prior.builder.extension.runtime())
    value['configuration_id']=canonical_id(value);write_immutable(out/'configuration.freeze.json',value)
    print('Frozen all-site reference command and',len(regions),'unique regions',flush=True)


def execute(out,verify=False):
    ext=audit.screen.prior.builder.extension
    config=ext.read_document(out/'configuration.freeze.json','configuration_id')
    if config['protocol']!=PROTOCOL or config['runtime']!=ext.runtime():raise ValueError('protocol or runtime changed')
    for r in config['inputs']+config['code']:checked_path(r,ROOT)
    target=out/'tracks_manifest.json'
    if target.exists():
        audit.screen.prior.builder.base.verify_stage(target,config['configuration_id'],ROOT)
        print('Verified constructed tracks; scientific overlap qualification is separate',flush=True);return
    if verify or (out/'tracks').exists():raise ValueError('missing or partial construction')
    subprocess.run(config['command'],check=True,cwd=ROOT)
    original=AtacBias().from_pickle(str(checked_path(config['original_model'],ROOT)))
    prefix=out/'tracks/NakedDNA_rep2_all_sites'
    produced=AtacBias().from_pickle(str(prefix)+'_AtacBias.pickle')
    if produced.stype!=original.stype or produced.correction_factor!=1.:raise ValueError('bias/normalization changed')
    for strand in original.bias:
        for attr in ['counts','bg_counts']:
            np.testing.assert_array_equal(getattr(original.bias[strand],attr),getattr(produced.bias[strand],attr))
    outputs={name:record(Path(str(prefix)+suffix)) for name,suffix in [('observed','_uncorrected.bw'),('expected','_expected.bw'),('bias_model','_AtacBias.pickle')]}
    result=dict(configuration_id=config['configuration_id'],outputs=outputs,model_refitted=False,
        overlap_qualified=False,candidates_scored=False,promotion_eligible=False)
    result['stage_id']=canonical_id(result);write_immutable(target,result)
    print('Constructed reference; overlap validation still required',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--outdir',type=Path,required=True)
    m=p.add_mutually_exclusive_group();m.add_argument('--freeze-only',action='store_true');m.add_argument('--verify-only',action='store_true')
    a=p.parse_args();out=a.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'):p.error('use research results')
    if a.freeze_only:freeze(out)
    else:
        with audit.screen.prior.builder.extension.exclusive_run(out):execute(out,a.verify_only)
