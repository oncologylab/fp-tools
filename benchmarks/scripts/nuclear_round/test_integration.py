"""Synthetic end-to-end output/verifier contract; not scientific evidence."""
import json
import numpy as np
import pandas as pd
from access import digest
from run import evaluate_control
from verify import verify
from conditional_protection_score import depletion_score


def test_control_tables_and_independent_verifier(tmp_path):
    out=tmp_path/'nuclear'
    for d in ['tables','arrays','models']:(out/d).mkdir(parents=True)
    (tmp_path/'inputs').mkdir()
    n=4200;position=np.arange(-100,101);center=abs(position)<=7
    q=15/73
    frame=pd.DataFrame(dict(site_id=[f'fixture_{i}' for i in range(n)],cell='fixture',
        tf=np.repeat([f'TF{i}' for i in range(21)],200),TFBS_chr=[f'chr{i%15+1}' for i in range(n)],
        eligible=True,LOG21_plus_valid=True,LOG21_minus_valid=True,DWM_plus_valid=True,DWM_minus_valid=True,
        deduplicated_fit_representative=True,plus_center=2,minus_center=2,plus_total=10,minus_total=10,total=20))
    frame.loc[:34,'eligible']=False
    frame.to_csv(tmp_path/'inputs/control_counts.full.tsv.gz',sep='\t',index=False)
    value=round(depletion_score(10,10,4,q,q),10);ceiling=round(depletion_score(10,10,0,q,q),10)
    score=pd.concat([frame.assign(arm=arm,q_plus=q,q_minus=q,score=value,ceiling=ceiling,status='NO_CALL',call=False)
                     for arm in ['nuclear_primary','nuclear_ablation']],ignore_index=True)
    score.loc[~score.eligible,'status']='EXCLUDED'
    score.to_csv(out/'tables/control_site_scores.tsv.gz',sep='\t',index=False)
    old=pd.concat([frame.assign(arm=a,q_plus=.2,q_minus=.2,call=False)
                   for a in ['LOG21_count_tail','LOG21NC1_count_tail','DWM_count_tail']],ignore_index=True)
    old.to_csv(tmp_path/'inputs/control_reference_scores.tsv.gz',sep='\t',index=False)
    evaluate_control(out,frame,score,old)
    observed=np.zeros((n,201),np.uint32);observed[:,100]=2;observed[:,70]=observed[:,130]=4
    for l2 in [.01,.0001]:
        np.savez_compressed(out/f'arrays/control.l2_{l2}.npz',site_id=frame.site_id.to_numpy(str),
            plus_observed=observed,minus_observed=observed,plus_log_bias=np.zeros((n,201)),minus_log_bias=np.zeros((n,201)))
        coefficients=dict(main=np.zeros((21,4)),**{f'pair_{d}':np.zeros((21-d,4,4)) for d in [1,5,9]})
        for c in range(1,16):
            for seed in range(2026,2031):
                p=out/f'models/chr{c}.l2_{l2}.seed_{seed}.npz';np.savez_compressed(p,**coefficients)
                p.with_suffix('.json').write_text(json.dumps(dict(npz_sha256=digest(p),metadata=dict(held_chromosome=f'chr{c}',seed=seed))))
            np.savez_compressed(out/f'models/chr{c}.l2_{l2}.ensemble.npz',**coefficients)
    np.savez_compressed(out/'arrays/control_sequence_contexts.npz',site_id=frame.site_id.to_numpy(str),
        chromosome=frame.TFBS_chr.to_numpy(str),motif_strand=np.array(['+','-']*(n//2)),sequence=np.array(['A'*283]*n))
    np.savez_compressed(out/'arrays/training_windows.npz',chromosome=np.array(['chr1']*3000+['chr2']*3000))
    (out/'models/folds.json').write_text(json.dumps({f'chr{i}':list(range(3000,6000)) if i==1 else list(range(3000)) for i in range(1,16)}))
    result=verify(tmp_path)
    assert result['independent_score_and_call_rows']==8400
    assert result['decision']=='REJECTED_DEVELOPMENT'
    assert not result['scientific_qualification']
