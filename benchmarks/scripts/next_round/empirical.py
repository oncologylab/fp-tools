"""Whole-fragment synthetic protection on examined naked-DNA replicate1."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logit
from shrinkage import posterior
from access import digest
from prepare import dump


def perturb(args):
    out,block,row,endpoints,calibration=args
    out=Path(out);k=len(endpoints);rng=np.random.default_rng(np.random.SeedSequence([2026,block,741]))
    # columns endpoint1_position, strand, endpoint2_position, strand.
    region=np.zeros((k,2,3),np.int64)
    for end in [0,2]:
        pos=endpoints[:,end];strand=endpoints[:,end+1]
        for s in [0,1]:
            region[:,s,0]+=((abs(pos)<=7)&(strand==s)).astype(int)
            region[:,s,1]+=((pos>=-40)&(pos<=-12)&(strand==s)).astype(int)
            region[:,s,2]+=((pos>=12)&(pos<=40)&(strand==s)).astype(int)
    m=region[:,:,0].sum(1)
    original_center=float(m.sum());original_total=float(region.sum())
    all_counts=[];mean_outputs=[];call_outputs=[];rows=[]
    methods=['LOG21_binomial','LOG21_beta_binomial','DWM_binomial','DWM_beta_binomial']
    for r in [1.,.75,.5,.25,0.]:
        weights=np.power(r,m)
        if not k or weights.sum()==0:
            all_counts.append(np.full((2000,2,3),-1,np.int64));mean_outputs.append(np.full((2000,4),np.nan));call_outputs.append(np.zeros((2000,4),bool))
            for method in methods:rows.append(dict(block=block,site_id=row['site_id'],method=method,r=r,status='UNAVAILABLE_NO_SOURCE_OR_WEIGHT',source_fragments=k,repetitions=0))
            continue
        weights=weights/weights.sum()
        multiplicity=rng.multinomial(k,weights,size=2000)
        counts=(multiplicity@region.reshape(k,6)).reshape(2000,2,3)
        expected_region=(weights@region.reshape(k,6)).reshape(2,3)*k
        true_q=np.divide(expected_region[:,0],expected_region.sum(1),out=np.full(2,np.nan),where=expected_region.sum(1)>0)
        realised_ratio=expected_region[:,0].sum()/original_center if original_center else np.nan
        effect=counts.mean(0)
        means=np.full((2000,4),np.nan);calls=np.zeros((2000,4),bool)
        for j,method in enumerate(methods):
            arm=method.split('_')[0];record=calibration[(method,row['TFBS_chr'])]
            p=[row[arm+'_plus_p'],row[arm+'_minus_p']]
            if record['status']=='UNIDENTIFIED':
                rows.append(dict(block=block,site_id=row['site_id'],method=method,r=r,status='UNIDENTIFIED',source_fragments=k,repetitions=2000));continue
            cutoff=record['threshold'];rho=record['rho'];info=0
            for draw in range(2000):
                n=counts[draw].sum(1);c=counts[draw,:,0]
                result=posterior(*c,*n,*p,rho)
                means[draw,j]=result['mean'];info+=int(not result['prior_only'])
                calls[draw,j]=cutoff is not None and not result['prior_only'] and result['mean']>cutoff
            call_rate=calls[:,j].mean()
            rows.append(dict(block=block,site_id=row['site_id'],cell=row['cell'],tf=row['tf'],chromosome=row['TFBS_chr'],method=method,r=r,
                status='AVAILABLE' if cutoff is not None else 'CALLING_UNDERPOWERED',source_fragments=k,repetitions=2000,
                source_central_cuts=original_center,source_regional_cuts=original_total,
                realised_marginal_center_ratio=realised_ratio,expected_center=expected_region[:,0].sum(),
                expected_left=expected_region[:,1].sum(),expected_right=expected_region[:,2].sum(),
                expected_total=expected_region.sum(),sample_mean_center=effect[:,0].sum(),
                frozen_expected_center=row['LOG21_expected_center'],mean_score=float(np.mean(means[:,j])),
                calls=int(calls[:,j].sum()),call_rate=call_rate,informative_draws=info,
                monte_carlo_se=float(np.sqrt(call_rate*(1-call_rate)/2000)),
                true_q_plus=true_q[0],true_q_minus=true_q[1],coverage_status='UNDEFINED_NO_COMMON_SCALAR_TRUTH',
                independence='one empirical source window; pseudo-repetitions are not independent biological libraries'))
        all_counts.append(counts);mean_outputs.append(means);call_outputs.append(calls)
    np.savez_compressed(out/f'arrays/empirical_simulations/block_{block:03d}.npz',
        counts=np.stack(all_counts),mean=np.stack(mean_outputs),call=np.stack(call_outputs),
        r=np.asarray([1.,.75,.5,.25,0.]),methods=np.asarray(methods),source_endpoints=endpoints,site_id=np.asarray(row['site_id']))
    return pd.DataFrame(rows)


def execute(out,workers):
    controls=pd.read_csv(out/'tables/control_counts.full.tsv.gz',sep='\t')
    calibration=json.loads((out/'models/calibration.full.json').read_text())
    lookup={(r['method'],r['fold']):r for r in calibration['records']}
    with np.load(out/'arrays/control_complete_fragments.npz',allow_pickle=False) as z:
        indices=z['site_index'];endpoints=z['endpoints'];hashes=z['fragment_hash']
    eligible=controls[controls.eligible & controls.deduplicated_fit_representative].copy()
    eligible['selection_hash']=[sha256(str(x).encode()).hexdigest() for x in eligible.physical_id]
    selected=[];used={};stratum_count={}
    for index,row in eligible.sort_values('selection_hash').iterrows():
        stratum=(row.count_stratum,row.expected_stratum)
        if stratum_count.get(stratum,0)>=10:continue
        if any(abs(int(row.center)-center)<=200 for center in used.get(row.TFBS_chr,[])):continue
        selected.append(index);used.setdefault(row.TFBS_chr,[]).append(int(row.center));stratum_count[stratum]=stratum_count.get(stratum,0)+1
    panel=controls.loc[selected].copy();panel['block']=np.arange(len(panel))
    parent=list(range(len(panel)));owner={}
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    for index,row in panel.iterrows():
        block=int(row.block)
        for identity in hashes[indices==index]:
            identity=int(identity)
            if identity in owner:parent[root(block)]=root(owner[identity])
            else:owner[identity]=block
    panel['shared_fragment_component']=[root(i) for i in range(len(panel))]
    panel.to_csv(out/'tables/empirical_source_blocks.tsv',sep='\t',index=False)
    dump(out/'provenance/empirical_simulation.freeze.json',dict(code_sha256=digest(__file__),
        blocks_sha256=digest(out/'tables/empirical_source_blocks.tsv'),source_fragments_sha256=digest(out/'arrays/control_complete_fragments.npz'),
        calibration_sha256=digest(out/'models/calibration.full.json'),repetitions=2000,
        r=[1.,.75,.5,.25,0.],unit='whole fragments; fixed K, within-window normalized r**central_endpoint_count',
        protection_is_synthetic=True,source_block_count=len(panel),shared_fragment_components=int(panel.shared_fragment_component.nunique()),independent_libraries=1,
        selected_before_perturbation=True,no_common_scalar_truth_assumed=True))
    (out/'arrays/empirical_simulations').mkdir()
    jobs=[]
    for index,row in panel.iterrows():
        mask=indices==index
        if len(np.unique(hashes[mask]))!=int(mask.sum()):raise ValueError('duplicate whole fragment within a source window')
        jobs.append((str(out),int(row.block),row.to_dict(),endpoints[mask],lookup))
    summaries=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i,summary in enumerate(pool.map(perturb,jobs,chunksize=1)):
            summaries.append(summary);print('Empirical source block complete',i+1,'/',len(jobs),flush=True)
    pd.concat(summaries,ignore_index=True).to_csv(out/'tables/empirical_power.tsv',sep='\t',index=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);p.add_argument('--workers',type=int,default=8)
    a=p.parse_args();execute(a.outdir.resolve(),a.workers)
