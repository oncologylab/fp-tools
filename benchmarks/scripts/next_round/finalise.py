"""Consolidate complete denominators and decisions; never fit a new method."""
import argparse
import itertools
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from verify import read_table
from analyse import wilson
from prepare import dump
from access import digest


LIBRARIES=[f'{cell}_rep{r}' for cell in ['HepG2','K562'] for r in [1,2,3]]


def evidence_tables(out):
    tables=out/'tables'
    required=[f'panel_{library}_{depth}.attempt1' for library in LIBRARIES for depth in ['10M','full']]
    required += [f'links_{library}.attempt1' for library in LIBRARIES]
    for name in required:
        receipt=out/'logs'/name/'exit.json'
        if not receipt.exists() or json.loads(receipt.read_text())['returncode']!=0:
            raise RuntimeError('required initial evidence still incomplete: '+name)
    if not (tables/'task_ledger.cohort_freeze.tsv').exists():
        shutil.copyfile(tables/'task_ledger.tsv',tables/'task_ledger.cohort_freeze.tsv')
    task=read_table(tables/'task_ledger.cohort_freeze.tsv')
    task['family']=task.tf.map({'CTCF':'CTCF','MEF2A':'MEF2','MEF2D':'MEF2','MYC':'MYC/MAX','MAX':'MYC/MAX'}).fillna('NOT_ADJUDICATED_FOR_FAMILY_GATE')
    task['source_samples']=task.cell.map(lambda c:','.join(c+'_rep'+str(r) for r in [1,2,3]))
    task['expected_initial_runs']=6
    task['completed_initial_extractions']=6
    task['completed_initial_scoring_runs']=np.where(task.nested_pairs>0,6,0)
    task['stage_A']='INCOMPLETE_NATIVE_DWM_REFERENCE'
    task['stage_B']='DEVELOPMENT_DIAGNOSTICS'
    task['stage_C']=np.where(task.nested_pairs>0,'DEVELOPMENT_DIAGNOSTICS','UNAVAILABLE_NO_BALANCED_COHORT')
    task['stage_D']='INTENTIONALLY_NOT_RUN_SCREEN_NOT_MET'
    task['stage_E']='RESOURCE_GATE_INCOMPLETE'
    task.to_csv(tables/'task_ledger.tsv',sep='\t',index=False)
    all_counts=[];ceilings=[];ledger=[]
    nested=read_table(tables/'nested_sites.tsv.gz');replay=read_table(tables/'replay_sites.tsv.gz')
    for path in sorted(tables.glob('counts.*.tsv.gz'))+sorted(tables.glob('control_counts.*.tsv.gz')):
        data=read_table(path);library=data.library.iloc[0];depth=data.depth.iloc[0]
        links=tables/f'fragment_links.{library}.tsv.gz'
        if links.exists() and depth in ['10M','full']:
            source=read_table(links);columns=['site_id']+[c for c in source if c.endswith('_'+depth) or c=='shared_fragment_component']
            data=data.merge(source[columns],on='site_id',how='left',validate='one_to_one')
        data['observed_units']='motif-oriented cuts; two endpoints retained per complete fragment'
        data['sequence_expectation_scope']='frozen raw propensities, not native DWM corrected expectation'
        data['original_source_table']=path.name
        for arm in ['LOG21','DWM','uniform']:
            for s in ['plus','minus']:data[f'{arm}_{s}_expected_center']=data[s+'_total']*data[f'{arm}_{s}_p']
            value=np.zeros(len(data))
            for s in ['plus','minus']:
                n=data[s+'_total'].to_numpy();p=data[f'{arm}_{s}_p'].to_numpy()
                with np.errstate(divide='ignore',invalid='ignore'):
                    value+=np.where(n==0,0.,-n*np.log1p(-p)/np.log(10.))
            ceiling=data[['site_id','library','depth','seed','cell','tf','TFBS_chr','count_stratum','expected_stratum','total']].copy()
            ceiling['propensity']=arm;ceiling['uncapped_zero_center_ceiling']=value
            ceiling['historical_capped_rounded_ceiling']=np.round(np.minimum(300.,value),10)
            ceiling['strict_threshold']=1.3010299957
            ceiling['reachable']=ceiling.historical_capped_rounded_ceiling>1.3010299957
            ceiling['interpretation']='conditional attainable-score bound; not biological detection power'
            ceilings.append(ceiling)
        all_counts.append(data)
        columns=[c for c in ['site_id','parent_locus_id','TFBS_chr','TFBS_start','TFBS_end','TFBS_strand','center','motif_score','chip_label','central_201','outer_101_500','overlap_group','shared_fragment_component','fold','library','depth','supported','eligible'] if c in data]
        entry=data[columns].copy()
        entry['cohort_replay']=entry.site_id.isin(replay.site_id);entry['cohort_nested']=entry.site_id.isin(nested.site_id)
        entry=entry.merge(nested[['site_id','pair_id']],how='left',on='site_id',validate='one_to_one')
        entry['label_provenance']='no occupancy labels; examined naked-DNA replicate1' if library=='NakedDNA_rep1' else 'original permitted training-chromosome candidate pool'
        entry['blacklist_status']='archived permitted support/eligibility field; no inferred new blacklist'
        ledger.append(entry)
    pd.concat(all_counts,ignore_index=True).to_csv(tables/'counts_and_expectations.tsv.gz',sep='\t',index=False)
    pd.concat(ceilings,ignore_index=True).to_csv(tables/'attainable_ceiling.tsv.gz',sep='\t',index=False)
    pd.concat(ledger,ignore_index=True).to_csv(tables/'site_ledger.tsv.gz',sep='\t',index=False)
    inventories=[]
    for library in LIBRARIES+['NakedDNA_rep1']:
        value=json.loads((out/f'provenance/{library}.fragment_inventory.json').read_text())
        inventories.append(dict(library=library,**value))
    pd.DataFrame(inventories).to_csv(tables/'fragment_inventory.tsv',sep='\t',index=False)
    attempts=[]
    for directory in sorted((out/'logs').iterdir()):
        if not directory.is_dir():continue
        receipt=directory/'exit.json'
        value=json.loads(receipt.read_text()) if receipt.exists() else {}
        code=value.get('returncode')
        attempts.append(dict(attempt=directory.name,returncode=code,elapsed_seconds=value.get('elapsed_seconds'),
            status='EXIT_ZERO_NOT_SCIENTIFIC_QUALIFICATION' if code==0 else 'INTERRUPTED' if code is not None and code<0 else 'FAILED' if code is not None else 'INCOMPLETE',
            receipt_sha256=digest(receipt) if receipt.exists() else '',
            resource_gate='not a qualified complete-workflow resource measurement'))
    pd.DataFrame(attempts).to_csv(tables/'attempt_status.tsv',sep='\t',index=False)


def screen(out):
    folds=read_table(out/'tables/fold_metrics.tsv');folds=folds[folds.cohort.eq('nested')&folds.depth.eq('full')]
    tasks=read_table(out/'tables/task_ledger.tsv');rows=[];summary=[];loco=[]
    weights=np.random.default_rng(2026).multinomial(15,np.repeat(1/15,15),size=2000)
    chromosomes=[f'chr{i}' for i in range(1,16)]
    for method in ['LOG21_binomial','LOG21_beta_binomial']:
        for reference in ['archived_DWM_residual','LOG21_count_tail','LOG21_binomial','DWM_binomial']:
            if method==reference:continue
            boot=[];points=[];omissions={chrom:[] for chrom in chromosomes}
            for task in tasks.itertuples():
                subset=folds[folds.cell.eq(task.cell)&folds.tf.eq(task.tf)]
                a=subset[subset.method.eq(method)];b=subset[subset.method.eq(reference)]
                matched=a.merge(b,on=['library','fold'],suffixes=('','_ref'))
                gains=[];aps=[];bases=[];draws=[];blocks=[]
                for library in [task.cell+'_rep'+str(r) for r in [1,2,3]]:
                    g=matched[matched.library.eq(library)].set_index('fold').reindex(chromosomes)
                    delta=(g.auroc-g.auroc_ref).to_numpy();ap=(g.auprc-g.auprc_ref).to_numpy()
                    valid=np.isfinite(delta)&np.isfinite(ap);blocks.append(int(valid.sum()))
                    gains.append(np.mean(delta[valid]) if valid.any() else np.nan)
                    aps.append(np.mean(ap[valid]) if valid.any() else np.nan);bases.append(g.auprc_ref[valid].mean())
                    denom=weights@valid.astype(float)
                    draws.append(np.divide(weights@np.where(valid,delta,0.),denom,out=np.full(2000,np.nan),where=denom>0))
                complete=all(np.isfinite(gains))
                value=float(np.mean(gains)) if complete else np.nan
                apvalue=float(np.mean(aps)) if complete else np.nan
                base=float(np.mean(bases)) if complete else np.nan
                rows.append(dict(cell=task.cell,tf=task.tf,method=method,reference=reference,libraries_required=3,
                    libraries_available=int(np.isfinite(gains).sum()),minimum_nonempty_chromosomes=min(blocks),
                    auroc_gain=value,auprc_gain=apvalue,relative_auprc_gain=apvalue/base if base>0 else np.nan,
                    nested_gate_eligible_before_scoring=task.nested_gate_eligible_before_scoring,
                    status='DEVELOPMENT_DIAGNOSTIC' if complete else 'INCOMPLETE'))
                if complete:boot.append(np.mean(np.stack(draws),axis=0));points.append((value,apvalue/base))
                for omitted in chromosomes:
                    keep=matched[matched.fold.ne(omitted)].copy()
                    keep['delta']=keep.auroc-keep.auroc_ref
                    by_library=keep.groupby('library').delta.mean()
                    if len(by_library)==3 and by_library.notna().all():omissions[omitted].append(float(by_library.mean()))
            matrix=np.stack(boot) if boot else np.empty((0,2000));draw=matrix.mean(0) if len(matrix) else np.full(2000,np.nan)
            valid=np.isfinite(draw);lo,hi=np.quantile(draw[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan]
            summary.append(dict(method=method,reference=reference,available_tasks=len(points),original_tasks=21,
                mean_auroc_gain=float(np.mean([v[0] for v in points])) if points else np.nan,
                mean_relative_auprc_gain=float(np.mean([v[1] for v in points])) if points else np.nan,
                fixed_prediction_bootstrap_low=lo,fixed_prediction_bootstrap_high=hi,valid_draws=int(valid.sum()),
                principal_qualified=False,advancement=False,
                reason='missing qualified native DWM full-depth baseline, unavailable tasks, unchanged detector-only profiles and unmet calibration/uncertainty evidence'))
            for omitted,values in omissions.items():
                loco.append(dict(method=method,reference=reference,omitted_chromosome=omitted,available_tasks=len(values),original_tasks=21,
                                 mean_auroc_gain=float(np.mean(values)) if values else np.nan,scope='fixed-prediction block omission, not refitting'))
    pd.DataFrame(rows).to_csv(out/'tables/full_screen_task_deltas.tsv',sep='\t',index=False)
    pd.DataFrame(summary).to_csv(out/'tables/full_screen_summary.tsv',sep='\t',index=False)
    pd.DataFrame(loco).to_csv(out/'tables/leave_one_chromosome_sensitivity.tsv',sep='\t',index=False)


def control_pairs(out):
    rows=[]
    for depth in ['10M','full','archived_full']:
        frame=read_table(out/f'tables/control_scores.{depth}.tsv.gz')
        all_sites=read_table(out/f'tables/control_counts.{depth}.tsv.gz')
        for (cell,tf),task in all_sites.groupby(['cell','tf']):
            scores=frame[frame.cell.eq(cell)&frame.tf.eq(tf)]
            for method in ['LOG21_binomial','LOG21_beta_binomial']:
                # Same-operator control; explicitly not the archived DWM policy.
                reference='DWM_'+method.split('LOG21_',1)[1]
                a=scores[scores.method.eq(method)];b=scores[scores.method.eq(reference)]
                pair=a.merge(b,on='site_id',suffixes=('','_ref'),validate='one_to_one')
                pair=pair[pair.call_eligible&pair.call_eligible_ref]
                n=len(pair);k=int(pair.call.sum());kr=int(pair.call_ref.sum());lo,hi=wilson(k,n)
                rows.append(dict(depth=depth,cell=cell,tf=tf,method=method,reference=reference,
                    original_sites=len(task),eligible_original_sites=int(task.eligible.sum()),informative_original_sites=int((task.total>0).sum()),
                    paired_informative_eligible_sites=n,calls=k,reference_calls=kr,rate=k/n if n else np.nan,
                    wilson_high=hi,paired_rate_difference=(k-kr)/n if n else np.nan,
                    original_DWM_safety_comparison='INCOMPLETE_NATIVE_REFERENCE',
                    development_check='DEVELOPMENT_CRITERION_MET' if n and k/n<=.05 and hi<=.05 and (k-kr)/n<=.01 else 'FAIL',
                    dependence='same examined control library across tasks; Wilson descriptive, not independent-site guarantee'))
    pd.DataFrame(rows).to_csv(out/'tables/paired_control_rates.tsv',sep='\t',index=False)


def skipped(out):
    depth=[];resource=[]
    for library,d,seed in itertools.product(LIBRARIES,['10M','25M','50M','full'],range(2026,2031)):
        ran=d in ['10M','full'] and seed==2026
        depth.append(dict(library=library,depth=d,seed=seed,status='COMPLETED_INITIAL_DEVELOPMENT_PANEL' if ran else 'INTENTIONALLY_NOT_RUN',
            reason='initial frozen comparison/screen' if ran else 'advancement screen not met; do not expand search'))
    for library,d,r,arm in itertools.product(LIBRARIES,['10M','full'],range(1,7),['production_DWM','candidate_fit_and_apply']):
        resource.append(dict(library=library,depth=d,repetition=r,arm=arm,cores=4,numerical_threads=1,
            cpu_seconds=np.nan,absolute_simultaneous_peak_bytes=np.nan,status='INCOMPLETE',
            execution='INTENTIONALLY_NOT_RUN_SCREEN_NOT_MET',reason='equivalent production output contract and instrumentation not qualified'))
    pd.DataFrame(depth).to_csv(out/'tables/depth_replication.tsv',sep='\t',index=False)
    pd.DataFrame(resource).to_csv(out/'tables/resource_results.tsv',sep='\t',index=False)
    gates=[('mean_AUROC_gain_0.03','INCOMPLETE','qualified full-depth archived/native DWM comparator missing'),
        ('relative_AUPRC_gain_10_percent','INCOMPLETE','same limitation; existing LOG21 count-tail comparison is diagnostic'),
        ('canonical_shape_separation_gain_25_percent_two_difficult_families','INCOMPLETE','full-depth native DWM shape comparator unqualified; detector-only LOG21 profile change is exactly zero, which alone does not determine a correction-versus-DWM shape gain'),
        ('two_families_two_contexts','INCOMPLETE','no jointly qualified ranking/shape/safety result'),
        ('positive_control_loss_at_most_0.02','INCOMPLETE','original comparator incomplete'),
        ('principal_chromosome_interval_excludes_zero','INCOMPLETE','available bootstrap conditional on fixed predictions; original full comparator missing'),
        ('untouched_control_FPR_and_Wilson_at_most_0.05','NOT_TESTED_LOCKED','final replicate3 deliberately closed'),
        ('untouched_paired_DWM_FPR_increase_at_most_0.01','NOT_TESTED_LOCKED','final replicate3 deliberately closed'),
        ('new_context_holdouts','NOT_TESTED_LOCKED','SK-N-SH/GM23338 label contents deliberately closed'),
        ('bias_did_not_learn_occupancy','INCOMPLETE','mitochondrial coefficients unchanged, but this is not a new comprehensive transfer qualification'),
        ('depth_seed_replicate_stability','INCOMPLETE','six libraries at10M/full seed2026 only; expansion stopped'),
        ('complete_workflow_CPU_at_most_2x','INCOMPLETE','formal matrix intentionally not run'),
        ('absolute_process_tree_memory_at_most_1.5x','INCOMPLETE','instrumentation/output contract not qualified'),
        ('deployed_model_at_most_25MB','INCOMPLETE','score prototype is not a complete deployed workflow'),
        ('functional_uncertainty_coverage_calibrated','INCOMPLETE','site-effect working-model intervals are not functional bands'),
        ('GP_hybrid_beats_simpler_spline','INCOMPLETE','no GP/hybrid in this bounded round')]
    pd.DataFrame(gates,columns=['gate','status','reason']).to_csv(out/'tables/promotion_gates.tsv',sep='\t',index=False)
    record=dict(protocol_unchanged=True,execution='INTENTIONALLY_NOT_RUN',
        candidate_commands=None,baseline_commands=None,reason='No complete production-equivalent candidate output exists and scientific screen is not met. Supplying invented commands would not qualify the missing execution configuration.',
        required_instrumentation='isolated process-tree CPU and absolute simultaneous peak, including reaped descendants',
        core_budget=4,repeats=6,model_training='all five bias components and permitted control I/O must be included',
        timing_claims_made=False)
    path=out/'provenance/resource_execution_blockers.json'
    if path.exists():
        if json.loads(path.read_text())!=record:raise ValueError('resource blocker contract changed')
    else:dump(path,record)


def plots(out):
    model=read_table(out/'tables/model_power.tsv');empirical=read_table(out/'tables/empirical_power.tsv')
    uncertainty=[]
    for kind,table in [('model',model),('empirical_fragment',empirical)]:
        for r in table.itertuples():
            n=int(r.repetitions);k=int(r.calls) if np.isfinite(r.calls) else np.nan
            low,high=wilson(k,n) if n and np.isfinite(k) else (np.nan,np.nan)
            uncertainty.append(dict(layer=kind,source_id=getattr(r,'scenario',getattr(r,'block',None)),method=r.method,
                perturbation_r=getattr(r,'r',np.nan),status=r.status,calls=k,repetitions=n,conditional_call_rate=k/n if n and np.isfinite(k) else np.nan,
                conditional_MC_wilson_low=low,conditional_MC_wilson_high=high,
                uncertainty_scope='Monte Carlo draws conditional on generating model/source window; not uncertainty over biological libraries'))
    pd.DataFrame(uncertainty).to_csv(out/'tables/power_uncertainty.tsv',sep='\t',index=False)
    tasks=read_table(out/'tables/task_ledger.tsv')
    with PdfPages(out/'figures/all_task_power.pdf') as pdf:
        fig,axes=plt.subplots(7,3,figsize=(12,21))
        for ax,task in zip(axes.flat,tasks.itertuples()):
            sub=empirical[empirical.cell.eq(task.cell)&empirical.tf.eq(task.tf)]
            for method,g in sub.groupby('method'):
                points=g.groupby('r').call_rate.mean();ax.plot(points.index,points.values,marker='.',label=method)
            ax.set_title(task.cell+' '+task.tf);ax.set_xlabel('fragment weight r');ax.set_ylabel('mean window call fraction');ax.set_ylim(0,1)
            if not len(sub):ax.text(.5,.5,'NO SELECTED SOURCE BLOCK',transform=ax.transAxes,ha='center',fontsize=7)
        axes.flat[0].legend(fontsize=5);fig.suptitle('Synthetic fragment challenge: stratified source windows, not population FPR');fig.tight_layout(rect=[0,0,1,.98]);pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(1,2,figsize=(11,4))
        selected=model[model.p_plus.eq(.2)&model.p_minus.eq(.2)&model.shift_plus.eq(0)&model.shift_minus.eq(0)&model.generating_rho.eq(0)]
        for method,g in selected.groupby('method'):
            for theta in [0.,1.,2.]:
                v=g[g.theta.eq(theta)].sort_values('mean_expected_central')
                axes[0].plot(v.mean_expected_central,v.call_rate,marker='.',label=f'{method}, theta={theta}')
                axes[1].plot(v.mean_expected_central,v.coverage,marker='.',label=f'{method}, theta={theta}')
        for ax in axes:ax.set_xlabel('expected central count');ax.set_ylim(0,1.05)
        axes[0].set_ylabel('call fraction');axes[1].set_ylabel('conditional interval coverage');axes[0].legend(fontsize=6)
        fig.suptitle('Independent-binomial fixture slice; all210 scenarios retained in source tables');fig.tight_layout();pdf.savefig(fig);plt.close(fig)
    empirical.to_csv(out/'tables/power.tsv',sep='\t',index=False)
    sites=read_table(out/'tables/nested_sites.tsv.gz')
    sites=sites[sites.cell.eq('K562')&sites.tf.eq('CTCF')]
    with np.load(out/'arrays/K562_rep1.full.npz',allow_pickle=False) as z:
        ix=pd.Index(z['site_id']).get_indexer(sites.site_id)
        data=z['plus_observed'][ix]+z['minus_observed'][ix]
    figure_rows=[];fig,ax=plt.subplots(figsize=(7,4))
    for label in [0,1]:
        curve=data[sites.chip_label.to_numpy()==label].mean(0)
        ax.plot(np.arange(-100,101),curve,label=f'label {label}')
        figure_rows.extend(dict(label=label,position=x,value=y,sites=int(sites.chip_label.eq(label).sum())) for x,y in zip(range(-100,101),curve))
    pd.DataFrame(figure_rows).to_csv(out/'tables/representative_shape.tsv',sep='\t',index=False)
    ax.set(xlabel='base pairs from motif centre',ylabel='observed cuts per site per base',title='K562 CTCF, replicate1 full: complete held cohort')
    ax.legend();fig.tight_layout();fig.savefig(out/'figures/representative_shape.png',dpi=150);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,required=True);out=p.parse_args().outdir.resolve()
    for stage in [evidence_tables,screen,control_pairs,skipped,plots]:stage(out);print(stage.__name__,'complete',flush=True)
