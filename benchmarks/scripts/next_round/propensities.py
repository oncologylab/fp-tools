"""Frozen sequence-only factorial controls, distinct from production expectations."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pysam
from access import Access,deny,digest
from prepare import ROOT,DATA,dump

sys.path.insert(0,str(ROOT/'src'))
from fp_tools.tools.parametric_bias import BiasFeatureSpec,ConditionalSequenceBiasModel,strand_log_bias


def freeze(out):
    prior=json.loads((out/'provenance/input_access.local.json').read_text())
    for cell in ['HepG2','K562']:
        for rep in [1,2,3]:
            sample=f'{cell}_rep{rep}'
            path=DATA/f'data/public/processed/encode_k562_hepg2_all_site_aggregate_20260817/work/atac_correct/{sample}/{sample}_AtacBias.pickle'
            deny(path,sample+' DWM model')
            prior['inputs'][sample+'_DWM_model']=dict(path=str(path),resolved_path=str(path.resolve()),
                identity=sample+' archived full-library DWM model; sequence propensity comparison only',sha256=digest(path),bytes=path.stat().st_size)
    prior['amendment']='metadata-only addition of explicitly recorded native DWM model paths; no candidate outcomes examined'
    dump(out/'provenance/propensity_access.local.json',prior)
    dump(out/'provenance/propensity.freeze.json',dict(code_sha256=digest(__file__),
        operator='DWM raw sequence log2 likelihood ratio times ln2; same count/shrinkage operator as LOG21',
        production_equivalence=False,archive_expected_equivalence=False,
        model_inputs={k:v['sha256'] for k,v in prior['inputs'].items() if 'DWM_model' in k},
        contract='201 motif-oriented positions plus12-base DWM sequence context; LOG21 context21; no occupancy refitting'))


def calculate(args):
    out,library=args;out=Path(out)
    access=Access(out/'provenance/propensity_access.local.json',out/'logs/access.jsonl')
    control=library=='NakedDNA_rep1'
    sites=pd.read_csv(out/'tables'/('control_sites.tsv.gz' if control else 'pool.tsv.gz'),sep='\t')
    if not control:sites=sites[sites.cell.eq(library.split('_')[0])].reset_index(drop=True)
    genome=access.path('genome')
    model_path=access.path('control_old_AtacBias_pickle' if control else library+'_DWM_model')
    from fp_tools.tools.atacorrect_functions import AtacBias
    from fp_tools.utils.sequences import nuc_to_num
    native=AtacBias().from_pickle(str(model_path))
    if native.stype!='DWM' or any(v.length!=25 for v in native.bias.values()):raise ValueError('DWM model contract changed')
    archive={}
    for strand in ['forward','reverse']:
        for field in ['bias_pwm_log','bg_pwm_log','bias_dwm_log','bg_dwm_log']:
            archive[strand+'_'+field]=np.asarray(getattr(native.bias[strand],field))
    np.savez_compressed(out/f'models/{library}.DWM_parameters.npz',**archive)
    p=np.full((len(sites),201),np.nan);m=p.copy();gc=np.zeros(len(sites))
    lp=p.copy();lm=p.copy()
    fixed=ConditionalSequenceBiasModel(BiasFeatureSpec.loglinear21())
    with np.load(access.path('LOG21'),allow_pickle=False) as z:
        fixed.main=z['main']
        for distance in fixed.feature_spec.pair_distances:fixed.pairs[distance]=z[f'pair_{distance}']
    with pysam.FastaFile(str(genome)) as fasta:
        sizes=dict(zip(fasta.references,fasta.lengths))
        for i,row in enumerate(sites.itertuples()):
            center=int(row.center);start=center-141;end=center+142
            if row.TFBS_chr not in sizes or start<0 or end>sizes[row.TFBS_chr]:continue
            sequence=fasta.fetch(row.TFBS_chr,start,end).upper()
            gc[i]=(sequence.count('G')+sequence.count('C'))/len(sequence)
            encoded=nuc_to_num(sequence)
            complement=np.asarray([1,0,3,2,4],np.int64)[encoded[::-1]]
            first=native.bias['forward'].score_sequence(encoded)[41:242]*np.log(2.)
            second=native.bias['reverse'].score_sequence(complement)[::-1][41:242]*np.log(2.)
            lfirst,lsecond,valid=strand_log_bias(fixed,sequence,np.arange(41,242))
            if row.TFBS_strand=='-':
                first,second=second[::-1],first[::-1]
                lfirst,lsecond=lsecond[::-1],lfirst[::-1]
            p[i],m[i]=first,second;lp[i],lm[i]=lfirst,lsecond
    # Training LOG21 sequence parity must hold without reading labels outside
    # the permitted pool or changing the five-seed bias coefficients.
    if not control:
        with np.load(out/f'arrays/{library.split("_")[0]}.replay_pool.npz',allow_pickle=False) as old:
            np.testing.assert_allclose(lp,old['plus_log_bias'],atol=1e-10,rtol=1e-10)
            np.testing.assert_allclose(lm,old['minus_log_bias'],atol=1e-10,rtol=1e-10)
    np.savez_compressed(out/f'arrays/{library}.propensities.npz',DWM_plus=p,DWM_minus=m,
                        LOG21_plus=lp,LOG21_minus=lm,site_id=sites.site_id.to_numpy(str),gc=gc)
    dump(out/f'provenance/{library}.propensities.json',dict(library=library,
        code_sha256=digest(__file__),finite_sites=int((np.isfinite(p).all(1)&np.isfinite(m).all(1)&np.isfinite(lp).all(1)&np.isfinite(lm).all(1)).sum()),
        sites=len(sites),production_expectation=False,LOG21_parity=not control))
    print('Sequence propensities complete',library,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);p.add_argument('--freeze',action='store_true');p.add_argument('--library')
    a=p.parse_args();out=a.outdir.resolve()
    if a.freeze:freeze(out)
    elif a.library:calculate((out,a.library))
    else:
        libraries=['NakedDNA_rep1']+[f'{c}_rep{r}' for c in ['HepG2','K562'] for r in [1,2,3]]
        with ProcessPoolExecutor(max_workers=3) as pool:list(pool.map(calculate,[(str(out),lib) for lib in libraries]))
