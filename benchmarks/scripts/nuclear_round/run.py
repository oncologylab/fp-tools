"""Frozen research stages. Existing outputs are never overwritten."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from hashlib import blake2b
import heapq
from itertools import groupby
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import pysam
from scipy.special import logsumexp
from scipy.stats import norm

from access import Access, deny, digest, dump
from protocol import CHROMS, POLICY, SEEDS

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'benchmarks/scripts'))
from fp_tools.tools.parametric_bias import (BiasFeatureSpec, ConditionalSequenceBiasModel,
    ensemble_sequence_bias_models, strand_log_bias, cut_position_from_alignment)
from evaluate_parametric_bias import ControlWindowDataset
from conditional_protection_score import depletion_score


def freeze(out):
    out.mkdir(parents=True,exist_ok=False)
    for name in ['logs','tables','arrays','models','scratch']: (out/name).mkdir()
    parent=ROOT.parent
    prior=parent/'fp-tools-next-round/benchmarks/results/next_round_20260922_v1'
    old=parent/'fp-tools-calibration/benchmarks/results/calibration_round_20261019_v1'
    registry=json.loads((prior/'provenance/input_access.local.json').read_text())['inputs']
    records={}
    for key in ['NakedDNA_rep1_bam','NakedDNA_rep1_bai','genome','genome_index']:
        r=registry[key]; p=Path(r['path']);deny(p,r['identity'])
        if str(p.resolve())!=r['resolved_path'] or digest(p)!=r['sha256']: raise ValueError(key)
        records[key]=dict(path=str(p),resolved_path=str(p.resolve()),identity=r['identity'],sha256=r['sha256'])
    blacklist=json.loads((ROOT/'benchmarks/manifests/direct_bias_replicate_depths_v1.json').read_text())['blacklist']
    extra=dict(blacklist=(parent/'fp-tools'/blacklist['path'],'hg38 blacklist v2',blacklist['sha256']),
        controls=(prior/'tables/control_counts.full.tsv.gz','naked DNA replicate1 control regional counts',None),
        control_arrays=(prior/'arrays/NakedDNA_rep1.full.npz','naked DNA replicate1 paired observed arrays',None),
        old_scores=(old/'tables/control_site_scores.tsv.gz','naked DNA replicate1 historical reference scores',None),
        ledger=(prior/'tables/task_ledger.tsv','all21 task metadata',None))
    for key,(p,identity,expected) in extra.items():
        deny(p,identity); h=digest(p)
        if expected is not None and h!=expected:raise ValueError(key)
        records[key]=dict(path=str(p),resolved_path=str(p.resolve()),identity=identity,sha256=h)
    code=list(Path(__file__).parent.glob('*.py'))+[ROOT/'src/fp_tools/tools/parametric_bias.py',
        ROOT/'benchmarks/scripts/evaluate_parametric_bias.py',ROOT/'benchmarks/scripts/conditional_protection_score.py']
    dump(out/'inputs.local.json',records)
    dump(out/'EXPERIMENT_PLAN.freeze.json',dict(policy=POLICY,
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        inputs={k:{'identity':v['identity'],'sha256':v['sha256']} for k,v in records.items()},
        code={str(p.relative_to(ROOT)):digest(p) for p in code},
        commands=['run.py freeze --outdir ATTEMPT','run.py extract --outdir ATTEMPT',
                  'run.py fit --outdir ATTEMPT --workers 16','run.py control --outdir ATTEMPT'],
        created_unix=time.time(),training_window_selection='genome-coordinate hash, independent of counts/labels'))


def check(out):
    frozen=json.loads((out/'EXPERIMENT_PLAN.freeze.json').read_text())
    if frozen['policy']!=POLICY:raise ValueError('policy changed; preserve attempt')
    for rel,h in frozen['code'].items():
        if digest(ROOT/rel)!=h:raise ValueError('code changed: '+rel)


def pairs(handle, inventory):
    for name, reads in groupby(handle.fetch(until_eof=True),key=lambda r:r.query_name):
        groups=defaultdict(list)
        for r in reads:
            inventory['records']+=1
            inventory['mapped_records']+=int(not r.is_unmapped)
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail or r.is_duplicate or r.mapping_quality<30:
                inventory['excluded_records']+=1
                continue
            groups[r.get_tag('RG') if r.has_tag('RG') else ''].append(r)
        for group in groups.values():
            if len(group)!=2:
                inventory['singleton_or_ambiguous_records']+=len(group);continue
            a,b=group
            if not (a.is_paired and b.is_paired and a.is_read1!=b.is_read1 and a.is_read2!=b.is_read2
                    and a.is_read1!=a.is_read2 and b.is_read1!=b.is_read2
                    and not a.mate_is_unmapped and not b.mate_is_unmapped
                    and a.next_reference_id==b.reference_id and b.next_reference_id==a.reference_id
                    and a.next_reference_start==b.reference_start and b.next_reference_start==a.reference_start):
                inventory['invalid_pair_records']+=2;continue
            inventory['complete_fragments']+=1
            yield a,b


def oracle(r):
    cigar=r.cigartuples or []
    edge=list(reversed(cigar)) if r.is_reverse else cigar
    soft=0
    for op,n in edge:
        if op==5:continue
        soft=n if op==4 else 0;break
    consumed=sum(n for op,n in cigar if op in (0,2,3,7,8))
    return r.reference_start+consumed+soft-5 if r.is_reverse else r.reference_start-soft+4


def extract(out):
    check(out);access=Access(out)
    bam=access.path('NakedDNA_rep1_bam');genome=access.path('genome')
    access.path('NakedDNA_rep1_bai');access.path('genome_index')
    bed=pd.read_csv(access.path('blacklist'),sep='\t',header=None,usecols=[0,1,2])
    if (out/'arrays/training_windows.npz').exists():raise FileExistsError('extraction complete; verify before reuse')
    heaps={c:[] for c in CHROMS}
    # Keep 3000 smallest hash windows per chromosome; this contains the global
    # 3000 for every leave-one-chromosome-out fold without storing the genome.
    with pysam.FastaFile(str(genome)) as fasta:
        for chrom in CHROMS:
            excluded=bed[bed[0].eq(chrom)].sort_values(1)
            starts=excluded[1].to_numpy();ends=excluded[2].to_numpy();j=0
            for start in range(200,fasta.get_reference_length(chrom)-241,200):
                while j<len(ends) and ends[j]<=start-41:j+=1
                if j<len(starts) and starts[j]<start+241:continue
                rank=int.from_bytes(blake2b(f'2026\0{chrom}\0{start}'.encode(),digest_size=8).digest(),'little')
                heap=heaps[chrom]
                if len(heap)==3000 and rank>=-heap[0][0]:continue
                sequence=fasta.fetch(chrom,start-41,start+241).upper()
                if set(sequence)-set('ACGT'):continue
                item=(-rank,start,sequence)
                if len(heap)<3000:heapq.heappush(heap,item)
                else:heapq.heapreplace(heap,item)
    candidates=sorted([(int(-rank),chrom,start,seq) for chrom,h in heaps.items() for rank,start,seq in h])
    chosen=set();folds={}
    for held in CHROMS:
        indices=[i for i,r in enumerate(candidates) if r[1]!=held][:3000]
        if len(indices)!=3000:raise ValueError('training windows unavailable')
        folds[held]=indices;chosen.update(indices)
    selected=sorted(chosen);mapping={old:new for new,old in enumerate(selected)}
    rows=[candidates[i] for i in selected]
    folds={c:[mapping[i] for i in idx] for c,idx in folds.items()}
    lookup={(r[1],r[2]):i for i,r in enumerate(rows)}
    counts=np.zeros((len(rows),2,200),np.uint32);inventory=Counter()
    collated=out/'scratch/name_collated.bam'
    if collated.exists():raise FileExistsError('preserve interrupted collate')
    print('Collating verified replicate1 BAM',flush=True)
    pysam.collate('-@','4','-o',str(collated),str(bam),str(out/'scratch/collate'))
    with pysam.AlignmentFile(str(collated),'rb') as handle:
        for pair in pairs(handle,inventory):
            for r in pair:
                cut=oracle(r)
                if cut!=cut_position_from_alignment(r,(4,-5)):raise AssertionError('coordinate oracle mismatch')
                inventory['oracle_comparisons']+=1
                idx=lookup.get((r.reference_name,(cut//200)*200))
                if idx is not None:counts[idx,int(r.is_reverse),cut%200]+=1
            if inventory['complete_fragments']%5000000==0:
                print('Complete fragments',inventory['complete_fragments'],flush=True)
    if inventory['complete_fragments']!=39179095 or inventory['mapped_records']!=78358765:
        raise ValueError('complete-fragment inventory differs from qualified source')
    np.savez_compressed(out/'arrays/training_windows.npz',chromosome=np.array([r[1] for r in rows]),
        start=np.array([r[2] for r in rows]),sequence=np.array([r[3] for r in rows]),counts=counts)
    dump(out/'models/folds.json',folds)
    dump(out/'tables/fragment_inventory.json',dict(inventory))
    dump(out/'extraction.complete.json',dict(arrays_sha256=digest(out/'arrays/training_windows.npz'),
        folds_sha256=digest(out/'models/folds.json'),inventory_sha256=digest(out/'tables/fragment_inventory.json')))
    print('Extraction complete',len(rows),'windows',int(counts.sum()),'assigned cuts',flush=True)


def fit_one(args):
    out,fold,l2,seed=args;out=Path(out)
    destination=out/f'models/{fold}.l2_{l2}.seed_{seed}.npz'
    if destination.exists():raise FileExistsError(destination)
    selected=json.loads((out/'models/folds.json').read_text())[fold]
    with np.load(out/'arrays/training_windows.npz',allow_pickle=False) as z:
        chrom=z['chromosome'][selected];seq=z['sequence'][selected];counts=z['counts'][selected]
        if np.any(chrom==fold):raise AssertionError('held chromosome in fit')
        ds=ControlWindowDataset('NakedDNA_rep1','fit_excluding_'+fold,'naked_dna',(4,-5),200,41,
            chrom,z['start'][selected],seq,counts[:,0],counts[:,1],np.zeros(len(selected)))
    model=ConditionalSequenceBiasModel(BiasFeatureSpec.loglinear21())
    x,y=ds.model_arrays(model.feature_spec)
    if not len(y):raise ValueError('no informative fitting windows')
    start=time.monotonic();model.fit(x,y,l2=l2,seed=seed,**POLICY['optimizer'])
    model.save(destination,dict(held_chromosome=fold,l2=l2,seed=seed,window_count=len(selected),
        informative_strand_windows=len(y),training_cuts=int(y.sum()),seconds=time.monotonic()-start,
        input_array_sha256=digest(out/'arrays/training_windows.npz')))
    return dict(fold=fold,l2=l2,seed=seed,status='COMPLETED',model=destination.name,sha256=digest(destination))


def fit(out,workers):
    check(out)
    if not 1<=workers<=32:raise ValueError('worker budget')
    receipt=json.loads((out/'extraction.complete.json').read_text())
    if digest(out/'arrays/training_windows.npz')!=receipt['arrays_sha256'] or digest(out/'models/folds.json')!=receipt['folds_sha256']:
        raise ValueError('extraction identity')
    jobs=[(str(out),fold,l2,seed) for l2 in [.01,.0001] for fold in CHROMS for seed in SEEDS]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for record in pool.map(fit_one,jobs):
            with (out/'logs/fits.jsonl').open('a') as log:log.write(json.dumps(record)+'\n')
            print(record,flush=True)
    for l2 in [.01,.0001]:
        for fold in CHROMS:
            members=[ConditionalSequenceBiasModel.load(out/f'models/{fold}.l2_{l2}.seed_{s}.npz') for s in SEEDS]
            ensemble_sequence_bias_models(members).save(out/f'models/{fold}.l2_{l2}.ensemble.npz',
                dict(held_chromosome=fold,l2=l2,seeds=SEEDS,weights=[.2]*5))
    dump(out/'fit.complete.json',dict(fits=150,ensembles=30))


def control(out):
    check(out);access=Access(out)
    if not (out/'fit.complete.json').exists():raise ValueError('fit prerequisite incomplete')
    frame=pd.read_csv(access.path('controls'),sep='\t')
    if len(frame)!=4200 or not set(frame.TFBS_chr)<=set(CHROMS):raise ValueError('control cohort')
    with np.load(access.path('control_arrays'),allow_pickle=False) as z:
        assert np.array_equal(z['site_id'],frame.site_id.to_numpy(str))
        plus=z['plus_observed'];minus=z['minus_observed']
    old=pd.read_csv(access.path('old_scores'),sep='\t')
    grid=np.arange(-100,101);center=abs(grid)<=7;region=center|((abs(grid)>=12)&(abs(grid)<=40))
    np.testing.assert_array_equal(plus[:,center].sum(1),frame.plus_center)
    np.testing.assert_array_equal(minus[:,center].sum(1),frame.minus_center)
    np.testing.assert_array_equal(plus[:,region].sum(1),frame.plus_total)
    np.testing.assert_array_equal(minus[:,region].sum(1),frame.minus_total)
    scores=[]
    with pysam.FastaFile(str(access.path('genome'))) as fasta:
        for l2 in [.01,.0001]:
            lp=np.full((4200,201),np.nan);lm=lp.copy()
            for chrom in CHROMS:
                model=ConditionalSequenceBiasModel.load(out/f'models/{chrom}.l2_{l2}.ensemble.npz')
                for i,row in frame[frame.TFBS_chr.eq(chrom)].iterrows():
                    sequence=fasta.fetch(chrom,int(row.center)-141,int(row.center)+142).upper()
                    p,m,valid=strand_log_bias(model,sequence,np.arange(41,242))
                    if row.TFBS_strand=='-':p,m=m[::-1],p[::-1]
                    lp[i],lm[i]=p,m
            np.savez_compressed(out/f'arrays/control.l2_{l2}.npz',site_id=frame.site_id.to_numpy(str),
                plus_observed=plus,minus_observed=minus,plus_log_bias=lp,minus_log_bias=lm)
            qp=np.exp(logsumexp(lp[:,center],axis=1)-logsumexp(lp[:,region],axis=1))
            qm=np.exp(logsumexp(lm[:,center],axis=1)-logsumexp(lm[:,region],axis=1))
            result=frame.copy();result['arm']='nuclear_primary' if l2==.01 else 'nuclear_ablation'
            result['q_plus']=qp;result['q_minus']=qm
            values=[]
            for row,p,m in zip(frame.itertuples(),qp,qm):
                if not np.isfinite([p,m]).all(): values.append((np.nan,np.nan));continue
                values.append(tuple(round(depletion_score(int(row.plus_total),int(row.minus_total),int(c),p,m),10)
                    for c in [row.plus_center+row.minus_center,0]))
            result['score']=[v[0] for v in values];result['ceiling']=[v[1] for v in values]
            result['status']=np.select([~np.isfinite(result.score),result.total.eq(0),result.ceiling.le(1.3010299957),
                result.score.gt(1.3010299957)],['UNAVAILABLE','NO_INFORMATION','ABSTAIN_UNREACHABLE','CALL'],default='NO_CALL')
            result['call']=result.status.eq('CALL');scores.append(result)
    scored=pd.concat(scores,ignore_index=True)
    scored.to_csv(out/'tables/control_site_scores.tsv.gz',sep='\t',index=False)
    evaluate_control(out,frame,scored,old)


def evaluate_control(out,frame,scored,old):
    weights=np.random.default_rng(2026).multinomial(15,np.repeat(1/15,15),2000)
    means=[];likelihood=[];calls=[]
    for arm in scored.arm.unique():
        g=scored[scored.arm.eq(arm)].set_index('site_id').loc[frame.site_id].reset_index()
        support=g.eligible & np.isfinite(g.q_plus) & np.isfinite(g.q_minus) & g.DWM_plus_valid & g.DWM_minus_valid & g.LOG21_plus_valid & g.LOG21_minus_valid
        g['support']=support.to_numpy()
        g['expected_center']=g.plus_total*g.q_plus+g.minus_total*g.q_minus
        for (cell,tf),task in g.groupby(['cell','tf']):
            valid=task[task.support];central=valid.plus_center+valid.minus_center
            sums=pd.DataFrame(dict(chrom=valid.TFBS_chr,observed=central,expected=valid.expected_center)).groupby('chrom').sum().reindex(CHROMS,fill_value=0)
            denom=weights@sums.expected.to_numpy();num=weights@sums.observed.to_numpy()
            ratios=np.divide(num,denom,out=np.full(len(num),np.nan),where=denom>0)
            finite=ratios[np.isfinite(ratios)];lo,hi=np.quantile(finite,[.025,.975]) if len(finite)>=1900 else (np.nan,np.nan)
            means.append(dict(arm=arm,cell=cell,tf=tf,all_sites=len(task),support_sites=len(valid),
                observed=int(central.sum()),expected=valid.expected_center.sum(),ratio_low=lo,ratio_high=hi,
                valid_draws=len(finite),discrepant=bool(lo>1 or hi<1),status='AVAILABLE' if len(finite)>=1900 else 'INCOMPLETE'))
            informative=valid[valid.total.gt(0)];n=len(informative);k=int(informative.call.sum());z=norm.ppf(.975)
            upper=((k/n+z*z/(2*n))+z*np.sqrt(k/n*(1-k/n)/n+z*z/(4*n*n)))/(1+z*z/n) if n else np.nan
            ref=old[old.arm.eq('DWM_count_tail')].set_index('site_id').loc[informative.site_id]
            calls.append(dict(arm=arm,cell=cell,tf=tf,all_sites=len(task),support_sites=len(valid),informative_sites=n,
                reachable_sites=int(informative.ceiling.gt(1.3010299957).sum()),calls=k,wilson_upper=upper,
                paired_DWM_increase=(k-int(ref.call.sum()))/n if n else np.nan,
                scope='examined replicate1; DWM propensity factorial, not qualified native DWM'))
        valid=g[g.support & g.deduplicated_fit_representative]
        def ll(table):
            return sum(table[s+'_center']*np.log(table['q_'+s])+(table[s+'_total']-table[s+'_center'])*np.log1p(-table['q_'+s]) for s in ['plus','minus'])
        for refarm in ['LOG21_count_tail','LOG21NC1_count_tail']:
            ref=old[old.arm.eq(refarm)].set_index('site_id').loc[valid.site_id].reset_index()
            delta=ll(valid).to_numpy()-ll(ref).to_numpy()
            sums=pd.Series(delta,index=valid.TFBS_chr).groupby(level=0).sum().reindex(CHROMS,fill_value=0)
            boot=weights@sums.to_numpy();lo,hi=np.quantile(boot,[.025,.975])
            likelihood.append(dict(arm=arm,reference=refarm,unique_windows=len(valid),gain=delta.sum(),low=lo,high=hi,
                nonempty_blocks=int(valid.TFBS_chr.nunique())))
    for name,records in [('control_task_means',means),('control_calls',calls),('paired_likelihood',likelihood)]:
        pd.DataFrame(records).to_csv(out/f'tables/{name}.tsv',sep='\t',index=False)
    m=pd.DataFrame(means).query("arm == 'nuclear_primary'");c=pd.DataFrame(calls).query("arm == 'nuclear_primary'")
    l=pd.DataFrame(likelihood).query("arm == 'nuclear_primary'")
    gates=dict(likelihood=bool(len(l)==2 and l.low.gt(0).all()),mean=bool(len(m)==21 and m.status.eq('AVAILABLE').all() and m.discrepant.sum()<=3),
        factorial_controls=bool(len(c)==21 and c.wilson_upper.le(.05).all() and c.paired_DWM_increase.le(.01).all()),
        native_DWM_support='INCOMPLETE',locked_controls='NOT_TESTED_LOCKED')
    failed=not all(gates[k] for k in ['likelihood','mean','factorial_controls'])
    dump(out/'control.decision.json',dict(gates=gates,decision='REJECTED_DEVELOPMENT' if failed else 'INCOMPLETE_NATIVE_DWM_QUALIFICATION',
        chromatin_authorized=False,reason='No method substitution or tuning after the frozen primary screen.'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['freeze','extract','fit','control'])
    parser.add_argument('--outdir',type=Path,required=True);parser.add_argument('--workers',type=int,default=16)
    a=parser.parse_args();out=a.outdir.resolve()
    if a.stage=='freeze':freeze(out);return
    receipt=dict(stage=a.stage,started=time.time(),command=sys.argv)
    import os
    receipt.update(pid=os.getpid(),boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
        process_start_ticks=Path('/proc/self/stat').read_text().split()[21])
    dump(out/f'logs/{a.stage}.launch.json',receipt)
    try:
        if a.stage=='fit':fit(out,a.workers)
        else:globals()[a.stage](out)
    except BaseException as e:
        dump(out/f'logs/{a.stage}.exit.json',dict(exit_code=1,error=repr(e),ended=time.time()));raise
    dump(out/f'logs/{a.stage}.exit.json',dict(exit_code=0,ended=time.time(),scientific_qualification=False))


if __name__=='__main__':main()
