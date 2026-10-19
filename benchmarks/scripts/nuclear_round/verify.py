"""Independent read-only verification; never fits, extracts or repairs outputs."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import binom, norm


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def verify(root,return_tables=False):
    hashes=root/'SHA256SUMS'
    files=0
    if hashes.exists():
        for line in hashes.read_text().splitlines():
            expected,name=line.split('  ',1)
            rel=Path(name)
            assert not rel.is_absolute() and '..' not in rel.parts
            p=root/rel
            assert p.is_file() and not p.is_symlink() and sha(p)==expected,name
            files+=1
    scientific=root/'nuclear' if (root/'nuclear').exists() else root
    if not (scientific/'control.decision.json').exists():
        raise ValueError('scientific control stage incomplete; no qualification claimed')
    scores=pd.read_csv(scientific/'tables/control_site_scores.tsv.gz',sep='\t')
    controls=pd.read_csv(root/'inputs/control_counts.full.tsv.gz',sep='\t')
    assert len(controls)==4200 and controls.groupby(['cell','tf']).size().eq(200).all()
    checked=0;score_replay=[];mean_replay=[];call_replay=[];likelihood_replay=[]
    for arm,l2 in [('nuclear_primary',.01),('nuclear_ablation',.0001)]:
        tab=scores[scores.arm.eq(arm)].set_index('site_id').loc[controls.site_id].reset_index()
        with np.load(scientific/f'arrays/control.l2_{l2}.npz',allow_pickle=False) as arrays:
            assert np.array_equal(arrays['site_id'],controls.site_id.to_numpy(str))
            if (root/'inputs/control_observed.full.npz').exists():
                with np.load(root/'inputs/control_observed.full.npz',allow_pickle=False) as original:
                    for name in ['site_id','plus_observed','minus_observed']:
                        np.testing.assert_array_equal(arrays[name],original[name])
            offset=np.arange(201)-100
            center=np.abs(offset)<=7
            region=center|((np.abs(offset)>=12)&(np.abs(offset)<=40))
            for strand in ['plus','minus']:
                observed=arrays[strand+'_observed'];logits=arrays[strand+'_log_bias']
                np.testing.assert_array_equal(observed[:,center].sum(1),tab[strand+'_center'])
                np.testing.assert_array_equal(observed[:,region].sum(1),tab[strand+'_total'])
                scaled=np.exp(logits-logits[:,region].max(1)[:,None])
                probability=scaled[:,center].sum(1)/scaled[:,region].sum(1)
                np.testing.assert_allclose(probability,tab['q_'+strand],atol=2e-14,rtol=2e-14)
                if (scientific/'tables/control_profile_plot_source.tsv.gz').exists():
                    plotted=pd.read_csv(scientific/'tables/control_profile_plot_source.tsv.gz',sep='\t')
                    expected=np.exp(logits-logsumexp(logits[:,region],axis=1)[:,None])*tab[strand+'_total'].to_numpy()[:,None]
                    for (cell,tf),group in tab.groupby(['cell','tf']):
                        rows=plotted[(plotted.arm==arm)&(plotted.strand==strand)&(plotted.cell==cell)&(plotted.tf==tf)].sort_values('offset')
                        assert len(rows)==201 and rows.sites.eq(len(group)).all()
                        np.testing.assert_allclose(rows.observed,observed[group.index].mean(0),atol=1e-12,rtol=1e-12)
                        np.testing.assert_allclose(rows.expected,expected[group.index].mean(0),atol=1e-12,rtol=1e-12)
            with np.load(scientific/'arrays/control_sequence_contexts.npz',allow_pickle=False) as contexts:
                assert np.array_equal(contexts['site_id'],arrays['site_id'])
                for chrom in [f'chr{i}' for i in range(1,16)]:
                    indices=np.flatnonzero(contexts['chromosome']==chrom)
                    if not len(indices):continue
                    with np.load(scientific/f'models/{chrom}.l2_{l2}.ensemble.npz',allow_pickle=False) as model:
                        members=[]
                        for seed in range(2026,2031):
                            path=scientific/f'models/{chrom}.l2_{l2}.seed_{seed}.npz'
                            meta=json.loads(path.with_suffix('.json').read_text())
                            assert meta['npz_sha256']==sha(path)
                            assert meta['metadata']['held_chromosome']==chrom and meta['metadata']['seed']==seed
                            with np.load(path,allow_pickle=False) as member:members.append({k:member[k] for k in member.files})
                        for key in model.files:
                            np.testing.assert_allclose(model[key],np.mean([member[key] for member in members],axis=0),atol=1e-14,rtol=1e-14)
                        encoded=np.array([[{'A':0,'C':1,'G':2,'T':3}[base] for base in s]
                                          for s in contexts['sequence'][indices]],dtype=np.int64)
                        x=np.lib.stride_tricks.sliding_window_view(encoded,21,axis=1)[:,31:232,:]
                        predictions=[]
                        for strand in [x,3-x[:,:,::-1]]:
                            predicted=sum(model['main'][j,strand[:,:,j]] for j in range(21))
                            for distance in [1,5,9]:
                                predicted+=sum(model[f'pair_{distance}'][j,strand[:,:,j],strand[:,:,j+distance]]
                                               for j in range(21-distance))
                            predictions.append(predicted)
                        p,m=predictions;reverse=contexts['motif_strand'][indices]=='-'
                        oriented_p=p.copy();oriented_m=m.copy()
                        oriented_p[reverse]=m[reverse,::-1];oriented_m[reverse]=p[reverse,::-1]
                        np.testing.assert_allclose(oriented_p,arrays['plus_log_bias'][indices],atol=1e-12,rtol=1e-12)
                        np.testing.assert_allclose(oriented_m,arrays['minus_log_bias'][indices],atol=1e-12,rtol=1e-12)
        for r in tab.itertuples():
            n,m,c=int(r.plus_total),int(r.minus_total),int(r.plus_center+r.minus_center)
            if not np.isfinite([r.q_plus,r.q_minus]).all():
                assert r.status=='UNAVAILABLE';continue
            # Independent explicit joint support summation, not scorer CDF code.
            x=np.arange(n+1);y=np.arange(m+1)
            joint=binom.logpmf(x,n,r.q_plus)[:,None]+binom.logpmf(y,m,r.q_minus)[None,:]
            value=min(300.,max(0.,-logsumexp(joint[(x[:,None]+y[None,:])<=c])/np.log(10)))
            ceiling=min(300.,max(0.,-joint[0,0]/np.log(10)))
            value,ceiling=round(float(value),10),round(float(ceiling),10)
            np.testing.assert_allclose([value,ceiling],[r.score,r.ceiling],atol=1.1e-10,rtol=0)
            status=('NO_INFORMATION' if n+m==0 else 'ABSTAIN_UNREACHABLE' if ceiling<=1.3010299957
                    else 'CALL' if value>1.3010299957 else 'NO_CALL')
            if not (r.eligible and r.DWM_plus_valid and r.DWM_minus_valid and r.LOG21_plus_valid and r.LOG21_minus_valid):status='EXCLUDED'
            assert status==r.status and (status=='CALL')==r.call
            score_replay.append(dict(site_id=r.site_id,cell=r.cell,tf=r.tf,chromosome=r.TFBS_chr,
                arm=arm,score=value,ceiling=ceiling,status=status,call=status=='CALL'))
            checked+=1
    folds=json.loads((scientific/'models/folds.json').read_text())
    with np.load(scientific/'arrays/training_windows.npz',allow_pickle=False) as z:
        for held,indices in folds.items():
            assert len(indices)==3000 and len(set(indices))==3000
            assert not np.any(z['chromosome'][indices]==held)
    decision=json.loads((scientific/'control.decision.json').read_text())
    assert not decision['chromatin_authorized']
    old=pd.read_csv(root/'inputs/control_reference_scores.tsv.gz',sep='\t')
    task_means=pd.read_csv(scientific/'tables/control_task_means.tsv',sep='\t')
    task_calls=pd.read_csv(scientific/'tables/control_calls.tsv',sep='\t')
    likelihood=pd.read_csv(scientific/'tables/paired_likelihood.tsv',sep='\t')
    chromosomes=[f'chr{i}' for i in range(1,16)]
    weights=np.random.default_rng(2026).multinomial(15,[1/15]*15,2000)
    recomputed={}
    for arm in scores.arm.unique():
        g=scores[scores.arm.eq(arm)].copy()
        valid=g.eligible & np.isfinite(g.q_plus) & np.isfinite(g.q_minus) & g.DWM_plus_valid & g.DWM_minus_valid & g.LOG21_plus_valid & g.LOG21_minus_valid
        g=g[valid];discrepant=0;controls_ok=True
        for (cell,tf),group in g.groupby(['cell','tf']):
            o=group.plus_center+group.minus_center
            e=group.plus_total*group.q_plus+group.minus_total*group.q_minus
            ob=np.array([o[group.TFBS_chr.eq(c)].sum() for c in chromosomes])
            ex=np.array([e[group.TFBS_chr.eq(c)].sum() for c in chromosomes])
            numerator=weights@ob;denominator=weights@ex;valid_draw=denominator>0
            assert valid_draw.sum()>=1900
            low,high=np.quantile(numerator[valid_draw]/denominator[valid_draw],[.025,.975])
            recorded=task_means[(task_means.arm==arm)&(task_means.cell==cell)&(task_means.tf==tf)].iloc[0]
            np.testing.assert_allclose([o.sum(),e.sum(),low,high],
                [recorded.observed,recorded.expected,recorded.ratio_low,recorded.ratio_high],rtol=1e-12,atol=1e-12)
            mean_replay.append(dict(arm=arm,cell=cell,tf=tf,observed=o.sum(),expected=e.sum(),
                ratio_low=low,ratio_high=high,discrepant=bool(low>1 or high<1),support_sites=len(group)))
            discrepant+=int(low>1 or high<1)
            info=group[group.total.gt(0)];n=len(info);k=int(info.call.sum());z=norm.ppf(.975)
            upper=((k/n+z*z/(2*n))+z*np.sqrt(k/n*(1-k/n)/n+z*z/(4*n*n)))/(1+z*z/n) if n else np.nan
            ref=old[old.arm.eq('DWM_count_tail')].set_index('site_id').loc[info.site_id]
            delta=(k-int(ref.call.sum()))/n if n else np.nan
            saved=task_calls[(task_calls.arm==arm)&(task_calls.cell==cell)&(task_calls.tf==tf)].iloc[0]
            np.testing.assert_allclose([n,k,upper,delta],[saved.informative_sites,saved.calls,saved.wilson_upper,saved.paired_DWM_increase],rtol=1e-12,atol=1e-12)
            call_replay.append(dict(arm=arm,cell=cell,tf=tf,informative_sites=n,calls=k,wilson_upper=upper,paired_DWM_increase=delta))
            controls_ok &= upper<=.05 and delta<=.01
        dedup=g[g.deduplicated_fit_representative].copy();positive=True
        for reference in ['LOG21_count_tail','LOG21NC1_count_tail']:
            ref=old[old.arm.eq(reference)].set_index('site_id').loc[dedup.site_id]
            delta=np.zeros(len(dedup))
            for strand in ['plus','minus']:
                n=dedup[strand+'_total'].to_numpy(int);c=dedup[strand+'_center'].to_numpy(int)
                delta+=binom.logpmf(c,n,dedup['q_'+strand])-binom.logpmf(c,n,ref['q_'+strand].to_numpy())
            blocks=np.array([delta[dedup.TFBS_chr.eq(c)].sum() for c in chromosomes])
            low,high=np.quantile(weights@blocks,[.025,.975])
            saved=likelihood[(likelihood.arm==arm)&(likelihood.reference==reference)].iloc[0]
            np.testing.assert_allclose([delta.sum(),low,high],[saved.gain,saved.low,saved.high],atol=1e-9,rtol=1e-11)
            likelihood_replay.append(dict(arm=arm,reference=reference,gain=delta.sum(),low=low,high=high,unique_windows=len(dedup)))
            positive &= low>0
        recomputed[arm]=dict(mean=bool(discrepant<=3),likelihood=bool(positive),factorial_controls=bool(controls_ok))
    for key,value in recomputed['nuclear_primary'].items():assert decision['gates'][key]==value
    power_checks=0
    if (scientific/'diagnostics.complete.json').exists():
        power=pd.read_csv(scientific/'tables/power.tsv',sep='\t')
        simroot=root/'historical_calibration/data/arrays/empirical_simulations'
        for block in range(138):
            with np.load(simroot/f'block_{block:03d}.npz',allow_pickle=False) as z:
                site=str(z['site_id']);draws=z['counts'];rs=z['r'];endpoints=z['source_endpoints']
            m=(np.abs(endpoints[:,0])<=7).astype(int)+(np.abs(endpoints[:,2])<=7).astype(int)
            original=m.sum();k=len(m)
            for arm in ['nuclear_primary','nuclear_ablation']:
                row=scores[(scores.arm==arm)&(scores.site_id==site)].iloc[0];cache={}
                for ix,r in enumerate(rs):
                    record=power[(power.block==block)&(power.arm==arm)&(power.r==r)].iloc[0]
                    if (draws[ix]<0).any():
                        assert record.status=='UNAVAILABLE_NO_SOURCE_OR_WEIGHT';power_checks+=1;continue
                    w=np.power(r,m);w=w/w.sum()
                    depletion=1-k*np.dot(w,m)/original if original else np.nan
                    np.testing.assert_allclose(depletion,record.realised_depletion,atol=1e-12,rtol=1e-12,equal_nan=True)
                    counts=draws[ix];n=counts[:,0].sum(1);minus=counts[:,1].sum(1);central=counts[:,:,0].sum(1)
                    triples,frequency=np.unique(np.column_stack([n,minus,central]),axis=0,return_counts=True)
                    calls=0;reachable=0
                    for (a,b,c),freq in zip(triples,frequency):
                        key=int(a),int(b),int(c)
                        if key not in cache:
                            distribution=np.convolve(binom.pmf(np.arange(a+1),a,row.q_plus),binom.pmf(np.arange(b+1),b,row.q_minus))
                            p=distribution[:c+1].sum();p0=distribution[0]
                            score=round(min(300.,max(0.,-np.log10(max(p,1e-300)))),10)
                            ceiling=round(min(300.,max(0.,-np.log10(max(p0,1e-300)))),10)
                            cache[key]=(score>1.3010299957,ceiling>1.3010299957)
                        call,reach=cache[key];calls+=int(freq)*call;reachable+=int(freq)*reach
                    assert calls==record.calls and reachable==record.reachable_draws
                    assert len(counts)==record.draws and int(((n+minus)>0).sum())==record.informative_draws
                    np.testing.assert_allclose(calls/len(counts),record.call_rate,rtol=1e-12,atol=1e-12)
                    power_checks+=1
    result=dict(checksummed_files=files,independent_score_and_call_rows=checked,
        fold_exclusion_checks=len(folds),status='VERIFIED_RECORDED_CONTROL_OUTPUTS',
        independent_power_rows=power_checks,scientific_qualification=False,decision=decision['decision'])
    if return_tables:
        result['tables']={name:pd.DataFrame(records) for name,records in
            [('per_site_scores',score_replay),('control_task_means',mean_replay),('control_calls',call_replay),('paired_likelihood',likelihood_replay)]}
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('package',type=Path)
    print(json.dumps(verify(p.parse_args().package.resolve()),indent=2))
