#!/usr/bin/env python3
"""Extract wider label-free training context without changing frozen studies."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
import build_global_parametric_training_profiles as old
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import checked_path, canonical_id

ROOT=old.ROOT
SOURCE=ROOT/'benchmarks/results/global_parametric_training_profiles_20260922_v1'
PROTOCOL=ROOT/'benchmarks/manifests/wide_background_training_v1.json'


def background_diagnostic(counts,log_bias,flank,outer):
    counts,log_bias=np.asarray(counts,float),np.asarray(log_bias,float)
    if counts.shape!=log_bias.shape or counts.ndim!=2 or counts.shape[1]!=2*flank+1 or not 50<outer<=flank:
        raise ValueError('invalid background diagnostic shape or radius')
    if not np.isfinite(counts).all() or (counts<0).any() or np.isinf(log_bias).any():
        raise ValueError('invalid counts or bias')
    left=np.arange(flank-outer,flank-49);right=np.arange(flank+50,flank+outer+1)
    both=np.r_[left,right];available=np.isfinite(log_bias[:,both]).all(axis=1)
    left_count=counts[:,left].sum(axis=1);right_count=counts[:,right].sum(axis=1)
    disagreement=np.full(len(counts),np.nan)
    bias=log_bias[available];offset=bias[:,both].max(axis=1)
    left_exposure=np.exp(bias[:,left]-offset[:,None]).sum(axis=1)
    right_exposure=np.exp(bias[:,right]-offset[:,None]).sum(axis=1)
    if (left_exposure<=0).any() or (right_exposure<=0).any(): raise ValueError('exposure underflow')
    disagreement[available]=np.abs(np.log2((left_count[available]+.5)/(right_count[available]+.5))
                                  +np.log2(right_exposure/left_exposure))
    return dict(available=available,left_count=left_count,right_count=right_count,disagreement=disagreement)


def infer(sites,model,genome,flank):
    if not isinstance(flank,int) or flank<200:
        raise ValueError('integer flank at least 200 required')
    width=2*flank+1
    margin=max(41,model.feature_spec.context_length//2+1)
    plus=np.full((len(sites),width),np.nan);minus=np.full_like(plus,np.nan)
    with old.strands.open_fasta(genome) as fasta:
        lengths=dict(zip(fasta.references,fasta.lengths))
        for index,row in enumerate(sites.itertuples(index=False)):
            center=(int(row.TFBS_start)+int(row.TFBS_end))//2
            start,end=center-flank-margin,center+flank+margin+1
            if row.TFBS_chr not in lengths or start<0 or end>lengths[row.TFBS_chr]: continue
            first,second,valid=old.strands.strand_log_bias(model,
                fasta.fetch(row.TFBS_chr,start,end).upper(),margin+np.arange(width))
            plus[index,valid],minus[index,valid]=first[valid],second[valid]
    return plus,minus


def freeze(out):
    if out.exists(): raise ValueError('fresh output directory required')
    original=SOURCE/'configuration.freeze.json'
    checked_path(dict(path=str(original),sha256='d903eec175130993afe35e6f3c64701e25399d1e29f6d2b98f8750825849b328'),ROOT)
    config=json.loads(original.read_text())
    inputs=config['inputs']+[record(original),record(PROTOCOL)]
    cells={}
    for cell,item in config['cells'].items():
        stage=SOURCE/(cell+'.json');source=json.loads(stage.read_text())
        inputs += [record(stage)]+list(source['outputs'].values())
        cells[cell]=dict(item,training_sites=source['outputs']['sites'],previous_profiles=source['outputs']['profiles'])
    value=dict(protocol=json.loads(PROTOCOL.read_text()),cells=cells,inputs=inputs,
        code=config['code']+[record(Path(__file__))],runtime=old.parent.extension.runtime(),
        labels_read=False,profiles_scored=False,model_refitted=False,promotion_eligible=False)
    for entry in value['inputs']+value['code']: checked_path(entry,ROOT)
    value['configuration_id']=canonical_id(value)
    out.mkdir(parents=True);write_immutable(out/'configuration.freeze.json',value)
    print('Frozen wider-context training extraction',flush=True)


def execute(out,verify=False):
    config=old.parent.extension.read_document(out/'configuration.freeze.json','configuration_id')
    if config['runtime']!=old.parent.extension.runtime() or config['protocol']!=json.loads(PROTOCOL.read_text()):
        raise ValueError('runtime or protocol changed')
    for entry in config['inputs']+config['code']: checked_path(entry,ROOT)
    stages=[];protocol=config['protocol'];flank=protocol['extraction_flank']
    for cell,item in config['cells'].items():
        receipt=out/(cell+'.json');target=out/(cell+'.npz')
        if receipt.exists(): old.parent.base.verify_stage(receipt,config['configuration_id'],ROOT)
        elif verify or target.exists() or (out/(cell+'.diagnostics.tsv')).exists():
            raise ValueError('missing result or partial output preserved')
        else:
            sites=pd.read_csv(item['training_sites']['path'],sep='\t')
            if list(sites.columns)!=old.SITE_KEYS or not sites.TFBS_chr.isin(protocol['chromosomes']).all():
                raise ValueError('label-free training support required')
            hashes=old.strands.site_hashes(sites)
            if hashes.tolist()!=item['training_site_hashes']: raise ValueError('site identity changed')
            plus,minus,valid=old.strands.extract_strand_cut_profiles(sites,item['bam']['path'],
                flank=flank,read_shift=tuple(protocol['read_shift']),minimum_mapq=protocol['minimum_mapq'],keep_duplicates=False)
            model=old.strands.ConditionalSequenceBiasModel.load(item['model']['path'])
            first,second=infer(sites,model,item['genome']['path'],flank)
            with np.load(item['previous_profiles']['path'],allow_pickle=False) as previous:
                np.testing.assert_array_equal(hashes,previous['site_hash'])
                shared=valid & previous['bam_valid']
                for key,array in [('plus_observed',plus),('minus_observed',minus)]:
                    np.testing.assert_array_equal(array[shared,flank-200:flank+201],previous[key][shared])
                for key,array in [('plus_log_bias',first),('minus_log_bias',second)]:
                    before=previous[key];after=array[:,flank-200:flank+201]
                    common=np.isfinite(before)&np.isfinite(after)
                    np.testing.assert_allclose(after[common],before[common],rtol=1e-10,atol=1e-10)
            rows=[]
            for strand,counts,bias in [('plus',plus,first),('minus',minus,second)]:
                for outer in protocol['background_outer_grid']:
                    diagnostic=background_diagnostic(counts,bias,flank,outer)
                    for tf in sorted(sites.tf.unique()):
                        task=sites.tf.eq(tf).to_numpy();supported=task&valid&diagnostic['available']
                        left=diagnostic['left_count'][supported];right=diagnostic['right_count'][supported]
                        disagreement=diagnostic['disagreement'][supported]
                        rows.append(dict(cell=cell,tf=tf,strand=strand,outer=outer,sites=int(task.sum()),
                            supported=int(supported.sum()),zero_background=int(((left+right)==0).sum()),
                            zero_left=int((left==0).sum()),zero_right=int((right==0).sum()),
                            median_background_count=float(np.median(left+right)) if len(left) else None,
                            median_absolute_log2_rate_disagreement=float(np.median(disagreement)) if len(left) else None,
                            p90_absolute_log2_rate_disagreement=float(np.quantile(disagreement,.9)) if len(left) else None))
            diagnostic_path=out/(cell+'.diagnostics.tsv')
            pd.DataFrame(rows).to_csv(diagnostic_path,sep='\t',index=False)
            np.savez_compressed(target,plus_observed=plus,minus_observed=minus,plus_log_bias=first,
                minus_log_bias=second,bam_valid=valid,site_hash=hashes)
            result=dict(configuration_id=config['configuration_id'],cell=cell,
                outputs=dict(profiles=record(target),sites=item['training_sites'],diagnostics=record(diagnostic_path)),
                audit=dict(sites=len(sites),bam_valid=int(valid.sum()),central_count_parity=True,shared_bias_parity=True),
                labels_read=False,profiles_scored=False,promotion_eligible=False)
            result['stage_id']=canonical_id(result);write_immutable(receipt,result)
        stages.append(record(receipt));print('Verified wider context:',cell,flush=True)
    result=dict(configuration_id=config['configuration_id'],stages=stages,outputs={},
        labels_read=False,profiles_scored=False,promotion_eligible=False)
    result['stage_id']=canonical_id(result);path=out/'profiles_manifest.json'
    if path.exists():
        if old.parent.base.verify_stage(path,config['configuration_id'],ROOT)!=result: raise ValueError('manifest differs')
    elif verify: raise ValueError('missing manifest')
    else: write_immutable(path,result)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--freeze-only',action='store_true');mode.add_argument('--verify-only',action='store_true')
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'): parser.error('use research results')
    if args.freeze_only: freeze(out)
    else:
        with old.parent.extension.exclusive_run(out): execute(out,args.verify_only)
