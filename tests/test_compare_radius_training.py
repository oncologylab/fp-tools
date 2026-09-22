from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'benchmarks/scripts'))
import compare_radius_training as runner


def test_complete_synthetic_comparison_preserves_all_tasks_and_zero_ties(tmp_path,monkeypatch):
    protocol=json.loads(runner.PROTOCOL.read_text())
    sites=[];ledger=[]
    for index in range(21):
        cell='HepG2' if index<11 else 'K562';tf=f'TF{index}'
        ledger.append(dict(cell=cell,tf=tf,matched_pairs=1))
        for label in (1,0):
            sites.append(dict(cell=cell,tf=tf,TFBS_chr='chr1',TFBS_start=10000+index*2000+label*1000,
                TFBS_end=10010+index*2000+label*1000,TFBS_strand='-' if label else '+',chip_label=label,motif_score=1.))
    sites=pd.DataFrame(sites);panels={arm:np.zeros((42,201)) for arm in protocol['comparators']}
    monkeypatch.setattr(runner.previous.old,'load_split',lambda *args:(sites,panels.copy(),pd.DataFrame(ledger)))
    stages={};cohort={'outputs':{}}
    for cell in ('HepG2','K562'):
        part=sites.loc[sites.cell.eq(cell)]
        path=tmp_path/(cell+'.sites.tsv');part[runner.previous.KEYS].to_csv(path,sep='\t',index=False)
        archive=tmp_path/(cell+'.npz');zero=np.zeros((len(part),1201))
        np.savez(archive,plus_observed=zero,minus_observed=zero,plus_log_bias=zero,minus_log_bias=zero,
                 bam_valid=np.ones(len(part),bool),site_hash=runner.wide.old.strands.site_hashes(part))
        stages[cell]={'outputs':dict(sites=runner.record(path),profiles=runner.record(archive))}
        pairs=part[['cell','tf','chip_label']].copy();pairs['pair_id']=1;pairs['split']='train'
        pair_path=tmp_path/(cell+'.pairs.tsv');pairs.to_csv(pair_path,sep='\t',index=False)
        cohort['outputs'][cell+'_pairs']=runner.record(pair_path)
    cohort_path=tmp_path/'cohort.json';cohort_path.write_text(json.dumps(cohort))
    old_path=tmp_path/'old.json';old_path.write_text(json.dumps(dict(cohort_manifest=runner.record(cohort_path),dispersion=protocol['fixed_dispersion'])))
    monkeypatch.setattr(runner.previous,'OLD',old_path)
    scores=sites[runner.previous.KEYS+['chip_label']].copy();scores['arm']='gamma_fixed_dispersion';scores['score']=0.
    score_path=tmp_path/'prior.tsv';scores.to_csv(score_path,sep='\t',index=False)
    config=dict(configuration_id='synthetic',runtime={},protocol=protocol,inputs=[],code=[],stages=stages,prior_scores=runner.record(score_path))
    monkeypatch.setattr(runner.wide.old.parent.extension,'runtime',lambda:{})
    monkeypatch.setattr(runner.wide.old.parent.extension,'read_document',lambda *args:config)
    out=tmp_path/'result';out.mkdir();runner.execute(out)
    metrics=pd.read_csv(out/'metrics.tsv',sep='\t');parity=pd.read_csv(out/'radius100_parity.tsv',sep='\t')
    assert len(metrics)==126 and metrics.retained_pairs.eq(1).all()
    assert metrics.auroc.eq(.5).all() and metrics.auprc.eq(.5).all()
    assert metrics.status.eq('underpowered').all() and metrics.motif_smd.eq(0).all()
    assert not metrics.canonical_protection.any()
    assert len(parity)==21 and parity.max_absolute_score_difference.eq(0).all()
    result=json.loads((out/'comparison_manifest.json').read_text())
    assert not result['promotion_eligible'] and not result['candidate_selected']
