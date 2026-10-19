"""Prespecified control diagnostics; does not select, fit or advance candidates."""
import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import logsumexp

from access import deny,digest,dump
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'benchmarks/scripts'))
from conditional_protection_score import depletion_score


def freeze(out):
    prior=ROOT.parent/'fp-tools-calibration/benchmarks/results/calibration_round_20261019_v1/inputs.json'
    old=json.loads(prior.read_text())['permitted']
    keys=['tables/empirical_power.tsv']+[f'arrays/empirical_simulations/block_{i:03d}.npz' for i in range(138)]
    records={}
    for key in keys:
        r=old[key];p=Path(r['path']);deny(p,r['identity'])
        if str(p.resolve())!=r['resolved_path'] or digest(p)!=r['sha256']:raise ValueError(key)
        records[key]=r
    dump(out/'diagnostics.inputs.local.json',records)
    dump(out/'diagnostics.freeze.json',dict(code_sha256=digest(__file__),
        inputs={k:v['sha256'] for k,v in records.items()},
        empirical='reuse all 138 frozen whole-fragment blocks, 2000 draws and five r values; no independent endpoint resampling',
        effects='realized central depletion; undefined at zero source central count; r retained separately',
        endpoints=['central and both shoulder O/E by strand/task','fixed count and GC-stratum means',
                   'all-task full201 observed and expected curves','synthetic call rates and no-information/abstention'],
        gc_bins=[0,.3,.4,.5,.6,.7,1],status='frozen before nuclear control scoring'))


def run(out):
    frozen=json.loads((out/'diagnostics.freeze.json').read_text())
    expected=frozen['code_sha256']
    amendment=out/'DIAGNOSTICS_IMPLEMENTATION_CORRECTION.freeze.json'
    if amendment.exists():
        record=json.loads(amendment.read_text())
        assert record['original_sha256']==expected
        expected=record['corrected_sha256']
    assert digest(__file__)==expected
    inputs=json.loads((out/'diagnostics.inputs.local.json').read_text())
    def read_input(key):
        r=inputs[key];p=Path(r['path']);deny(p,r['identity'])
        assert str(p.resolve())==r['resolved_path'] and digest(p)==r['sha256']
        with (out/'logs/access.jsonl').open('a') as log:log.write(json.dumps(dict(input_id=key,sha256=r['sha256'],identity=r['identity']))+'\n')
        return p
    scores=pd.read_csv(out/'tables/control_site_scores.tsv.gz',sep='\t')
    registry=json.loads((out/'inputs.local.json').read_text());r=registry['old_scores'];p=Path(r['path']);deny(p,r['identity']);assert digest(p)==r['sha256']
    old=pd.read_csv(p,sep='\t');gc=old.drop_duplicates('site_id').set_index('site_id').gc
    offsets=np.arange(-100,101);regions=dict(center=abs(offsets)<=7,left=(offsets>=-40)&(offsets<=-12),right=(offsets>=12)&(offsets<=40))
    domain=regions['center']|regions['left']|regions['right'];regional=[];profiles=[]
    for arm,l2 in [('nuclear_primary',.01),('nuclear_ablation',.0001)]:
        s=scores[scores.arm.eq(arm)].reset_index(drop=True)
        with np.load(out/f'arrays/control.l2_{l2}.npz',allow_pickle=False) as z:
            assert np.array_equal(z['site_id'],s.site_id.to_numpy(str))
            for strand in ['plus','minus']:
                observed=z[strand+'_observed'];logits=z[strand+'_log_bias']
                expected=np.exp(logits-logsumexp(logits[:,domain],axis=1)[:,None])*s[strand+'_total'].to_numpy()[:,None]
                for name,mask in regions.items():
                    for i,row in s.iterrows():
                        regional.append(dict(site_id=row.site_id,cell=row.cell,tf=row.tf,chromosome=row.TFBS_chr,arm=arm,strand=strand,
                            region=name,observed=int(observed[i,mask].sum()),expected=expected[i,mask].sum(),
                            count_stratum=row.count_stratum,gc=float(gc.loc[row.site_id]),
                            total=int(row.total),status=row.status,eligible=row.eligible))
                for (cell,tf),group in s.groupby(['cell','tf']):
                    ix=group.index.to_numpy();o=observed[ix].mean(0);e=expected[ix].mean(0)
                    profiles.extend(dict(cell=cell,tf=tf,arm=arm,strand=strand,offset=x,observed=a,expected=b,
                        residual=a-b,sites=len(ix),units='cuts per site; conditional expected scaled by strand regional total')
                        for x,a,b in zip(offsets,o,e))
    regional=pd.DataFrame(regional);regional['gc_stratum']=pd.cut(regional.gc,[0,.3,.4,.5,.6,.7,1],include_lowest=True).astype(str)
    regional.to_csv(out/'tables/control_regional_diagnostics.tsv.gz',sep='\t',index=False)
    for grouping in ['count_stratum','gc_stratum']:
        summary=regional.groupby(['arm','strand','region',grouping],observed=False)[['observed','expected']].sum().reset_index()
        summary.to_csv(out/f'tables/calibration_by_{grouping}.tsv',sep='\t',index=False)
    profiles=pd.DataFrame(profiles);profiles.to_csv(out/'tables/control_profile_plot_source.tsv.gz',sep='\t',index=False)
    (out/'figures').mkdir(exist_ok=True)
    fig,axes=plt.subplots(3,7,figsize=(17,7),sharex=True)
    primary=profiles[profiles.arm.eq('nuclear_primary')].groupby(['cell','tf','offset'])[['observed','expected']].sum().reset_index()
    for ax,((cell,tf),g) in zip(axes.flat,primary.groupby(['cell','tf'])):
        ax.plot(g.offset,g.observed,label='observed',lw=.8);ax.plot(g.offset,g.expected,label='expected',lw=.8)
        ax.axvspan(-7,7,alpha=.1,color='grey');ax.set_title(f'{cell} {tf}',fontsize=8)
    axes.flat[0].legend(fontsize=7);fig.suptitle('Examined naked-DNA replicate 1: held-chromosome nuclear model; cuts per site')
    fig.tight_layout();fig.savefig(out/'figures/control_profiles_all_tasks.png',dpi=180);plt.close(fig)
    effects=pd.read_csv(read_input('tables/empirical_power.tsv'),sep='\t')
    power=[]
    for block in range(138):
        with np.load(read_input(f'arrays/empirical_simulations/block_{block:03d}.npz'),allow_pickle=False) as z:
            site=str(z['site_id']);counts=z['counts'];rs=z['r']
        for arm in ['nuclear_primary','nuclear_ablation']:
            site_score=scores[(scores.site_id==site)&(scores.arm==arm)]
            if len(site_score)!=1:raise ValueError('simulation site identity')
            row=site_score.iloc[0]
            for i,r in enumerate(rs):
                if not simulation_available(counts[i]):
                    power.append(dict(block=block,site_id=site,cell=row.cell,tf=row.tf,chromosome=row.TFBS_chr,
                        arm=arm,r=r,status='UNAVAILABLE_NO_SOURCE_OR_WEIGHT',draws=0))
                    continue
                c=counts[i];nplus=c[:,0].sum(1);nminus=c[:,1].sum(1);central=c[:,:,0].sum(1)
                triples=np.stack([nplus,nminus,central],axis=1);unique,inverse=np.unique(triples,axis=0,return_inverse=True)
                values=[];ceilings=[]
                for n,m,k in unique:
                    values.append(round(depletion_score(int(n),int(m),int(k),row.q_plus,row.q_minus),10))
                    ceilings.append(round(depletion_score(int(n),int(m),0,row.q_plus,row.q_minus),10))
                values=np.array(values)[inverse];ceilings=np.array(ceilings)[inverse]
                original=effects[(effects.block==block)&(effects.r==r)].iloc[0]
                power.append(dict(block=block,site_id=site,cell=row.cell,tf=row.tf,chromosome=row.TFBS_chr,arm=arm,r=r,
                    status='AVAILABLE',
                    realised_depletion=1-original.realised_marginal_center_ratio,draws=len(values),
                    informative_draws=int(((nplus+nminus)>0).sum()),reachable_draws=int((ceilings>1.3010299957).sum()),
                    calls=int((values>1.3010299957).sum()),call_rate=float((values>1.3010299957).mean()),
                    source_fragments=original.source_fragments,scope='one examined source window; 2000 draws not independent libraries'))
    power=pd.DataFrame(power);power.to_csv(out/'tables/power.tsv',sep='\t',index=False)
    fig,ax=plt.subplots(figsize=(6,4))
    for arm,g in power[np.isfinite(power.realised_depletion)].groupby('arm'):
        avg=g.groupby('r')[['realised_depletion','call_rate']].mean()
        ax.plot(avg.realised_depletion,avg.call_rate,'o-',label=arm)
    ax.set(xlabel='Realized central depletion (defined source windows)',ylabel='Call rate (same windows)',title='Whole-fragment synthetic challenge')
    ax.legend();fig.tight_layout();fig.savefig(out/'figures/power.png',dpi=180);plt.close(fig)
    dump(out/'diagnostics.complete.json',dict(regional_rows=len(regional),profile_rows=len(profiles),power_rows=len(power),independent_source_blocks=138))


def simulation_available(counts):
    counts=np.asarray(counts)
    if counts.ndim!=3 or counts.shape[1:]!=(2,3):raise ValueError('simulation count layout')
    return bool(np.all(counts>=0))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','run']);p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args();globals()[a.stage](a.outdir.resolve())
