"""Write conclusions and missing-stage ledgers from completed recorded results."""
import argparse
import json
from pathlib import Path
import pandas as pd


def finalise(out):
    decision=json.loads((out/'control.decision.json').read_text())
    means=pd.read_csv(out/'tables/control_task_means.tsv',sep='\t')
    calls=pd.read_csv(out/'tables/control_calls.tsv',sep='\t')
    likelihood=pd.read_csv(out/'tables/paired_likelihood.tsv',sep='\t')
    inventory=json.loads((out/'tables/fragment_inventory.json').read_text())
    primary=means[means.arm.eq('nuclear_primary')]
    ablation=means[means.arm.eq('nuclear_ablation')]
    task=means[['cell','tf']].drop_duplicates()
    if len(task)!=21:raise ValueError('task denominator')
    task['control_status']='COMPLETED_EXAMINED_DEVELOPMENT'
    task['new_chromatin_status']='INTENTIONALLY_NOT_RUN_CONTROL_PREREQUISITE'
    task.loc[task.cell.eq('HepG2')&task.tf.eq('HMGA1'),'new_chromatin_status']='UNAVAILABLE_HISTORICAL_COHORT'
    task['depth_seed_replication_status']='INTENTIONALLY_NOT_RUN'
    task['resource_status']='INTENTIONALLY_NOT_RUN'
    task.to_csv(out/'tables/task_ledger.tsv',sep='\t',index=False)
    scores=pd.read_csv(out/'tables/control_site_scores.tsv.gz',sep='\t')
    scores['expected_center']=scores.plus_total*scores.q_plus+scores.minus_total*scores.q_minus
    scores['geometry']='7,12,40'
    scores['fitting_provenance']='naked DNA replicate1; exclude row TFBS_chr; equal coefficient seeds2026-2030'
    scores.to_csv(out/'tables/counts_and_expectations.tsv.gz',sep='\t',index=False)
    scores[['site_id','cell','tf','TFBS_chr','arm','geometry','fitting_provenance','score','ceiling','status','call','eligible','exclusion_reason']].to_csv(
        out/'tables/per_site_scores.tsv.gz',sep='\t',index=False)
    ledger_columns=['site_id','cell','tf','TFBS_chr','TFBS_start','TFBS_end','TFBS_strand','center',
        'parent_locus_id','fold','overlap_group','physical_id','eligible','library','depth','seed']
    sites=scores.drop_duplicates('site_id')[ledger_columns].copy()
    sites['cohort']='examined naked DNA replicate1 development controls'
    sites.to_csv(out/'tables/site_ledger.tsv.gz',sep='\t',index=False)
    means.merge(calls,on=['arm','cell','tf'],validate='one_to_one',suffixes=('_mean','_calls')).to_csv(out/'tables/calibration.tsv',sep='\t',index=False)
    unavailable=[]
    for row in task.itertuples():
        for c in range(1,16):
            unavailable.append(dict(cell=row.cell,tf=row.tf,fold=f'chr{c}',status=row.new_chromatin_status,
                auroc=None,average_precision=None,auroc_gain=None,relative_ap_gain=None,canonical_protection=None))
    for name in ['fold_metrics','paired_deltas','shape_metrics']:
        pd.DataFrame(unavailable).to_csv(out/f'tables/{name}.tsv',sep='\t',index=False)
    accounting=dict(mapped_records=inventory['mapped_records'],complete_fragments=inventory['complete_fragments'],
        excluded_records=inventory.get('excluded_records',0),singleton_or_ambiguous_retained_records=inventory.get('singleton_or_ambiguous_records',0),
        invalid_pair_records=inventory.get('invalid_pair_records',0),coordinate_oracle_comparisons=inventory['oracle_comparisons'])
    pd.DataFrame([accounting]).to_csv(out/'tables/fragment_accounting.tsv',sep='\t',index=False)
    gates=[]
    for key,value in decision['gates'].items():
        status=('DEVELOPMENT_CRITERION_MET' if value else 'FAIL') if isinstance(value,bool) else value
        gates.append(dict(gate=key,status=status,scope='development control screen; not promotion'))
    for gate in ['AUROC_gain','relative_AP_gain','canonical_difficult_family_protection','two_contexts',
                 'positive_control_noninferiority','depth_seed_biological_replication','end_to_end_CPU',
                 'absolute_peak_memory','complete_deployed_artifact','functional_uncertainty','occupancy_leakage']:
        gates.append(dict(gate=gate,status='INCOMPLETE',scope='not qualified by this control-only round'))
    pd.DataFrame(gates).to_csv(out/'tables/promotion_gates.tsv',sep='\t',index=False)
    runs=[]
    for l2 in [.01,.0001]:
        for c in range(1,16):
            for seed in range(2026,2031):
                p=out/f'models/chr{c}.l2_{l2}.seed_{seed}.npz'
                if not p.exists():raise ValueError('missing expected fit')
                runs.append(dict(stage='fit',fold=f'chr{c}',l2=l2,seed=seed,status='COMPLETED'))
    runs.extend(dict(stage=s,fold='',l2='',seed='',status='INTENTIONALLY_NOT_RUN_PREREQUISITE')
        for s in ['chromatin','robustness_panel','full_resource_matrix','locked_controls'])
    pd.DataFrame(runs).to_csv(out/'tables/run_ledger.tsv',sep='\t',index=False)
    panel=[];resources=[]
    for cell in ['K562','HepG2']:
        for replicate in [1,2,3]:
            library=f'{cell}_rep{replicate}'
            for depth in ['10M','25M','50M','full']:
                for seed in range(2026,2031):
                    panel.append(dict(library=library,depth=depth,seed=seed,
                        status='INTENTIONALLY_NOT_RUN_CONTROL_PREREQUISITE',evidence='no new chromatin scoring'))
            for depth in ['10M','full']:
                for repeat in range(1,7):
                    resources.append(dict(library=library,depth=depth,seed=2026,repeat=repeat,cores=4,
                        cpu_seconds=None,absolute_simultaneous_peak_bytes=None,
                        status='INTENTIONALLY_NOT_RUN_CONTROL_PREREQUISITE'))
    pd.DataFrame(panel).to_csv(out/'tables/depth_replication.tsv',sep='\t',index=False)
    pd.DataFrame(resources).to_csv(out/'tables/resource_results.tsv',sep='\t',index=False)
    ll=likelihood[likelihood.arm.eq('nuclear_primary')]
    likelihood_text='\n'.join(f"- Versus {r.reference}: likelihood gain {r.gain:.3f}, 95% chromosome-block interval [{r.low:.3f}, {r.high:.3f}]." for r in ll.itertuples())
    text=f'''# Nuclear sequence-bias follow-up: {decision['decision']}

This is an examined-development control experiment, not a promotion study.
All 21 control tasks remain in the denominator. No new chromatin ranking
results, held test/transfer labels, replicate-2/3 controls or new cell labels
were used. Main and the previously rejected calibration artifacts are unchanged.

## What ran

{inventory['complete_fragments']:,} complete naked-DNA replicate-1 fragments,
with {inventory['oracle_comparisons']:,} independent coordinate comparisons.
150 LOG21 fits: 15 held chromosomes × 5 seeds × 2 prespecified regularizations;
30 equal-coefficient ensembles. Primary L2=0.01; ablation L2=0.0001. The frozen
mitochondrial and intercept-only references were not refitted or replaced.
The detector, geometry (7,12,40), score rounding and strict threshold remain fixed.

## Control results

Primary: **{int(primary.discrepant.sum())}/21** task mean intervals exclude the
expected mean (required <=3). Ablation: {int(ablation.discrepant.sum())}/21; it
cannot replace a failed primary. Bootstrap intervals condition on fitted
predictions and only 15 chromosome blocks, not full refitting uncertainty.

{likelihood_text}

The independent verifier must recompute these results before this package is
called verified. An exit code alone does not establish scientific qualification.

## What this does and does not establish

The experiment tests transfer of a nuclear naked-DNA-trained sequence model to
held nuclear naked-DNA loci. It does not demonstrate improved chromatin
footprinting. Old ChIP-negative call rates are not known false-positive rates.
Nuclear training also changes source preparation/depth and window selection
relative to the historical mitochondrial model; differences cannot be attributed
uniquely to DNA source. The model has one-base prediction resolution; 21 bases
is its sequence context, not a footprint resolution or a TF-specific model.

No new transcription-factor-specific or chromatin-label fitting was performed.
Raw measured cuts do not change when their expectation or detector changes.
The historical count/confounding and lack-of-protection findings remain scoped
to tested models, observations, geometry and cohorts, not all possible algorithms.

## Stop decision and remaining limits

Decision: **{decision['decision']}**. No chromatin or robustness/resource matrix
was launched. Native production DWM support and model-training isolation remain
unqualified; raw DWM propensities are a distinct factorial comparator. Untouched
safety, full-workflow resources and functional uncertainty remain untested.
Control Wilson intervals are descriptive under dependent loci and shared fragments.
The extra nuclear model is not a release candidate.
'''
    (out/'EXECUTIVE_FINDINGS.md').write_text(text)
    (out/'NEXT_DEEP_RESEARCH_PROMPT.md').write_text('''Challenge this specific evidence and runnable verifier, not a generic literature review.
Check whether each held chromosome is excluded from all five training seeds;
whether genome-grid sampling, complete fragments and matched blacklist preserve
the declared observation contract; whether source/preparation/depth differences
confound attribution; whether central and shoulder means transfer; whether
likelihood gains coexist with motif-specific failures; and whether conservative
calls obscure weak power. Recompute counts, probabilities, exact tails, all-task
denominators and block intervals. Distinguish conditional fixed-prediction
uncertainty from full refitting uncertainty, and examined naked DNA from untouched
safety. Do not treat raw DWM propensities as native TOBIAS. Is stopping justified,
or does an identifiable, falsifiable limitation remain? Do not propose reopening
protected controls, relaxing gates or searching more models for a positive result.
''')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True)
    finalise(p.parse_args().outdir.resolve())
