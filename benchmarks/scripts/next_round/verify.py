"""Independent read-only verifier for the finished derived-data evidence package.

Never fits, extracts, repairs or writes missing scientific outputs. Optional
--figure-check writes only a temporary independent figure outside the package.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path,PurePosixPath
import tempfile
import numpy as np
import pandas as pd
from scipy.special import logsumexp, xlogy
from scipy.stats import rankdata, norm


def read_table(path):
    # Strict calls require exact round-tripping of the recorded float64 values.
    return pd.read_csv(path,sep='\t',float_precision='round_trip',low_memory=False)


def digest(path):
    h=sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def load_propensities(root,library):
    """Resolve losslessly shared LOG21 arrays in compact evidence packages."""
    with np.load(root/f'arrays/{library}.propensities.npz',allow_pickle=False) as z:
        values={k:z[k] for k in z.files}
    if 'LOG21_reference' in values:
        name=str(values['LOG21_reference']);path=PurePosixPath(name)
        if path.is_absolute() or '..' in path.parts or path.parts[0]!='arrays':raise AssertionError('unsafe shared expectation reference')
        with np.load(root/name,allow_pickle=False) as z:
            np.testing.assert_array_equal(values['site_id'],z['site_id'])
            for key in ['LOG21_plus','LOG21_minus']:values[key]=z[key]
    return values


def metrics(y,score):
    y=np.asarray(y,int);score=np.asarray(score,float)
    positive=y.sum();negative=len(y)-positive
    if not positive or not negative:return np.nan,np.nan
    ranks=rankdata(score,method='average')
    auc=(ranks[y==1].sum()-positive*(positive+1)/2)/(positive*negative)
    order=np.argsort(-score,kind='stable');ys=y[order];ss=score[order]
    ends=np.r_[np.flatnonzero(ss[:-1]!=ss[1:]),len(y)-1]
    tp=np.cumsum(ys)[ends];increment=np.diff(np.r_[0,tp])
    ap=np.sum(increment/positive*(tp/(ends+1)))
    return auc,ap


def verify_checksums(root):
    entries=[]
    for line in (root/'SHA256SUMS').read_text().splitlines():
        checksum,name=line.split('  ',1);p=PurePosixPath(name)
        if p.is_absolute() or '..' in p.parts:raise AssertionError('unsafe manifest path')
        target=root/name
        if target.is_symlink() or not target.is_file():raise AssertionError('missing or unsafe payload '+name)
        if digest(target)!=checksum:raise AssertionError('checksum mismatch '+name)
        entries.append(name)
    actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name!='SHA256SUMS'}
    if actual!=set(entries):raise AssertionError('unmanifested/missing payload')
    index=json.loads((root/'ORIGIN_INDEX.json').read_text())
    for item in index['files']:
        if item['packaged_path'] not in actual:raise AssertionError('unresolved packaged origin reference')
        if digest(root/item['packaged_path'])!=item['packaged_sha256']:raise AssertionError('origin packaged hash mismatch')
    return len(entries)


def verify_completion_records(root):
    """Check explicit missingness and stop decisions, without inferring success."""
    tasks=read_table(root/'tables/task_ledger.tsv')
    if len(tasks)!=21 or tasks[['cell','tf']].duplicated().any():
        raise AssertionError('original task denominator changed')
    missing=tasks[tasks.cell.eq('HepG2')&tasks.tf.eq('HMGA1')]
    if len(missing)!=1 or missing.nested_pairs.iloc[0]!=0:
        raise AssertionError('unavailable task silently replaced')
    folds=read_table(root/'tables/fold_metrics.tsv')
    if len(folds)!=30870 or folds[['cohort','library','depth','cell','tf','method','fold']].duplicated().any():
        raise AssertionError('incomplete or duplicate method/run/fold ledger')
    if folds.status.isna().any():raise AssertionError('unrecorded method status')
    gates=read_table(root/'tables/promotion_gates.tsv')
    if not set(gates.status)<=set(['NOT_TESTED_LOCKED','INCOMPLETE','FAIL','DEVELOPMENT_CRITERION_MET']):
        raise AssertionError('unqualified gate status')
    for field in ['depth_replication.tsv','resource_results.tsv','counts_and_expectations.tsv.gz',
                  'site_ledger.tsv.gz','principal_regional_shapes.tsv']:
        if not (root/'tables'/field).is_file():raise AssertionError('required evidence missing: '+field)
    links=list((root/'provenance').glob('*.fragment_links.json'))
    if len(links)!=6:raise AssertionError('fragment-sharing inventory incomplete')
    audit=json.loads((root/'provenance/access_scope_audit.json').read_text())
    if audit['prohibited_recorded_reads'] or audit['protected_content_packaged']:
        raise AssertionError('protected input policy violated')
    return dict(tasks=21,method_run_fold_records=len(folds),fragment_link_libraries=6,
                access_scope='recorded allowlists/traces only; not whole-host surveillance')


def verify_counts(root):
    maximum=0.;rows=0
    for path in sorted((root/'tables').glob('counts.*.tsv.gz'))+sorted((root/'tables').glob('control_counts.*.tsv.gz')):
        table=read_table(path);library=table.library.iloc[0];depth=table.depth.iloc[0];cell=library.split('_')[0]
        array_depth='full' if depth=='archived_full' else depth
        array_path=root/(f'arrays/{cell}.replay_pool.npz' if depth=='archived_10M' else f'arrays/{library}.{array_depth}.npz')
        with np.load(array_path,allow_pickle=False) as z:
            index=pd.Index(z['site_id']).get_indexer(table.site_id)
            if (index<0).any():raise AssertionError('missing observation ID')
            arrays={s:z[s+'_observed'][index] for s in ['plus','minus']}
        z=load_propensities(root,library)
        index=pd.Index(z['site_id']).get_indexer(table.site_id)
        if (index<0).any():raise AssertionError('missing expectation ID')
        for strand in ['plus','minus']:
            counts=arrays[strand]
            # Independent direct slices, not the implementation masks.
            center=counts[:,93:108].sum(1);left=counts[:,60:89].sum(1);right=counts[:,112:141].sum(1)
            for name,values in [('center',center),('left',left),('right',right),('total',center+left+right)]:
                np.testing.assert_array_equal(table[strand+'_'+name],values)
            for arm in ['LOG21','DWM']:
                logits=z[arm+'_'+strand][index]
                selected=np.concatenate([logits[:,60:89],logits[:,93:108],logits[:,112:141]],axis=1)
                values=np.exp(logsumexp(logits[:,93:108],axis=1)-logsumexp(selected,axis=1))
                actual=table[arm+'_'+strand+'_p'].to_numpy()
                np.testing.assert_allclose(actual,values,rtol=1e-10,atol=1e-12,equal_nan=True)
                maximum=max(maximum,float(np.nanmax(abs(actual-values))))
        if not table.TFBS_chr.isin([f'chr{i}' for i in range(1,16)]).all():raise AssertionError('nontraining labels')
        if depth=='archived_full':
            for s in ['plus','minus']:
                previous=[]
                for context in ['HepG2','K562']:
                    with np.load(root/f'arrays/{context}.control_full.npz',allow_pickle=False) as z:previous.append(z[s+'_observed'])
                np.testing.assert_array_equal(np.concatenate(previous),arrays[s])
        rows+=len(table)
    return dict(regional_rows=rows,maximum_probability_error=maximum)


def verify_metrics(root):
    scores=read_table(root/'tables/per_site_scores.tsv.gz')
    folds=read_table(root/'tables/fold_metrics.tsv')
    keys=['cohort','library','depth','cell','tf','method','fold']
    groups={key:group for key,group in scores.groupby(keys,dropna=False)}
    errors=[];checked=0
    for row in folds.itertuples():
        key=tuple(getattr(row,k) for k in keys);group=groups.get(key)
        if row.status!='AVAILABLE':continue
        if group is None:raise AssertionError('available fold lacks per-site evidence')
        finite=group[np.isfinite(group.score)]
        if len(finite)!=row.original_sites or len(finite)!=row.finite_sites:raise AssertionError('fold denominator changed')
        auc,ap=metrics(finite.chip_label,finite.score)
        errors.extend([abs(auc-row.auroc),abs(ap-row.auprc)]);checked+=1
    maximum=max(errors,default=0.)
    if maximum>1e-12:raise AssertionError('independent metric mismatch')
    valid=scores.call_eligible.eq(True)&~scores.prior_only.eq(True)&np.isfinite(scores.score)&np.isfinite(scores.threshold)
    expected=valid&(scores.score>scores.threshold)
    if not np.array_equal(expected.to_numpy(),scores.call.fillna(False).to_numpy()):raise AssertionError('strict call counts mismatch')
    total_calls=int(expected.sum())
    old=scores[scores.method.isin(['uniform_count_tail','LOG21_count_tail'])]
    if not np.allclose(old.score,np.round(old.score,10),rtol=0,atol=1e-13,equal_nan=True):raise AssertionError('old rounding changed')
    if not (old.threshold==1.3010299957).all():raise AssertionError('old threshold changed')
    task=read_table(root/'tables/task_ledger.tsv')
    if len(task)!=21 or len(task[task.cell.eq('HepG2')&task.tf.eq('HMGA1')])!=1:raise AssertionError('task denominator changed')
    deltas=read_table(root/'tables/paired_deltas.tsv')
    fold_groups={key:group for key,group in folds.groupby(['cohort','library','depth','cell','tf'])}
    for row in deltas.itertuples():
        group=fold_groups[(row.cohort,row.library,row.depth,row.cell,row.tf)]
        a=group[group.method.eq(row.method)].set_index('fold')
        b=group[group.method.eq(row.reference)].set_index('fold')
        valid=np.isfinite(a.auroc)&np.isfinite(b.auroc)&np.isfinite(a.auprc)&np.isfinite(b.auprc)
        expected=(a.auroc[valid]-b.auroc[valid]).mean()
        np.testing.assert_allclose(expected,row.auroc_gain,atol=1e-12,rtol=0,equal_nan=True)
        difference=(a.auprc[valid]-b.auprc[valid]).mean();relative=difference/b.auprc[valid].mean() if valid.any() else np.nan
        np.testing.assert_allclose([difference,relative],[row.auprc_gain,row.relative_auprc_gain],atol=1e-12,rtol=0,equal_nan=True)
    return dict(folds_recomputed=checked,maximum_metric_error=maximum,score_rows=len(scores),calls=total_calls)


def verify_simulation(root):
    checked=0
    if (root/'tables/model_power.tsv').exists():
        summary=read_table(root/'tables/model_power.tsv')
        for path in sorted((root/'tables/model_simulations').glob('*.tsv.gz')):
            data=read_table(path)
            for method,group in data.groupby('method'):
                calls=(group['mean']>group.threshold)&((group.plus_total+group.minus_total)>0)
                coverage=(group.low<=group.true_theta)&(group.high>=group.true_theta)
                row=summary[summary.scenario.eq(group.scenario.iloc[0])&summary.method.eq(method)].iloc[0]
                if calls.sum()!=row.calls:raise AssertionError('model-simulation call mismatch')
                np.testing.assert_allclose(coverage.mean(),row.coverage,rtol=0,atol=1e-12)
                np.testing.assert_allclose((group['mean']-group.true_theta).mean(),row.mean_error,rtol=0,atol=1e-12)
                checked+=1
    if (root/'tables/empirical_power.tsv').exists():
        table=read_table(root/'tables/empirical_power.tsv')
        for path in sorted((root/'arrays/empirical_simulations').glob('block_*.npz')):
            block=int(path.stem.split('_')[1])
            with np.load(path,allow_pickle=False) as z:
                for i,r in enumerate(z['r']):
                    for j,method in enumerate(z['methods']):
                        row=table[table.block.eq(block)&table.r.eq(r)&table.method.eq(method)].iloc[0]
                        if row.status=='AVAILABLE':
                            if z['call'][i,:,j].sum()!=row.calls:raise AssertionError('empirical calls mismatch')
                            np.testing.assert_allclose(z['mean'][i,:,j].mean(),row.mean_score,atol=1e-12,rtol=0)
                            np.testing.assert_allclose(z['counts'][i,:,:,0].sum(1).mean(),row.sample_mean_center,atol=1e-10,rtol=0)
                        checked+=1
    return checked


def verify_control_rates(root):
    if not (root/'tables/calibration.tsv').exists():return 'not yet summarized'
    summary=read_table(root/'tables/calibration.tsv');checked=0
    for depth,rows in summary.groupby('depth'):
        source=read_table(root/f'tables/control_scores.{depth}.tsv.gz')
        groups={key:group for key,group in source.groupby(['cell','tf','method'])}
        for r in rows.itertuples():
            g=groups[(r.cell,r.tf,r.method)]
            if r.stratum_type!='all':g=g[g[r.stratum_type].eq(r.stratum)]
            calls=g.call_eligible.eq(True)&~g.prior_only.eq(True)&(g['mean']>g.threshold)
            if len(g)!=r.supported_control_sites or int(calls.sum())!=r.calls:
                raise AssertionError('control denominator/calls mismatch')
            n=len(g);k=int(calls.sum())
            if n:
                z=norm.ppf(.975);p=k/n
                hi=(p+z*z/(2*n)+z*np.sqrt(p*(1-p)/n+z*z/(4*n*n)))/(1+z*z/n)
                np.testing.assert_allclose([p,hi],[r.rate,r.wilson_high],atol=1e-12,rtol=0)
            checked+=1
    return checked


def verify_calibration_order(root):
    path=root/'tables/calibration_fit_scores.full.tsv.gz'
    if not path.exists():return 'explicit calibration replay not included'
    fitting=read_table(path);controls=read_table(root/'tables/control_counts.full.tsv.gz').set_index('site_id')
    records=json.loads((root/'models/calibration.full.json').read_text())['records'];checked=0
    for r in records:
        frame=fitting[fitting.method.eq(r['method'])&fitting.fold.eq(r['fold'])]
        arm=r['method'].split('_')[0]
        expected=controls[controls.eligible&controls.deduplicated_fit_representative&controls.total.gt(0)&controls.TFBS_chr.ne(r['fold'])&controls[arm+'_plus_valid']&controls[arm+'_minus_valid']]
        if set(frame.site_id)!=set(expected.index) or frame.site_id.duplicated().any():raise AssertionError('calibration fitting membership changed')
        if len(frame)<200:raise AssertionError('insufficient calibration support')
        k=int(np.ceil(.95*(len(frame)+1)));cutoff=np.sort(frame.score.to_numpy())[k-1]
        np.testing.assert_allclose(cutoff,r['threshold'],rtol=0,atol=0);checked+=1
    return checked


def verify_null_means(root):
    table=read_table(root/'tables/control_mean_variance.tsv.gz')
    table=table[table.stratum_type.eq('all')&table.strand.eq('combined')&table.region.eq('center')]
    weights=np.random.default_rng(2026).multinomial(15,np.repeat(1/15,15),size=2000);checked=0
    for depth,summary in table.groupby('depth'):
        controls=read_table(root/f'tables/control_counts.{depth}.tsv.gz')
        for r in summary.itertuples():
            g=controls[controls.cell.eq(r.cell)&controls.tf.eq(r.tf)&controls.eligible]
            observed=g.plus_center+g.minus_center;total=g.plus_total+g.minus_total
            expected=g.plus_total*g[r.arm+'_plus_p']+g.minus_total*g[r.arm+'_minus_p']
            variance=sum(g[s+'_total']*g[r.arm+'_'+s+'_p']*(1-g[r.arm+'_'+s+'_p']) for s in ['plus','minus'])
            chrom=g.TFBS_chr.str[3:].astype(int).to_numpy()-1
            numerator=np.bincount(chrom,weights=(observed-expected).to_numpy(),minlength=15)
            denominator=np.bincount(chrom,weights=total.to_numpy(),minlength=15)
            den=weights@denominator;draw=np.divide(weights@numerator,den,out=np.full(2000,np.nan),where=den>0)
            valid=np.isfinite(draw);bounds=np.quantile(draw[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan]
            actual=[r.observed,r.expected,r.mean_error_fraction,r.independent_variance_ratio,r.mean_error_low,r.mean_error_high]
            want=[observed.sum(),expected.sum(),(observed-expected).sum()/total.sum(),((observed-expected)**2).sum()/variance.sum(),*bounds]
            np.testing.assert_allclose(actual,want,atol=1e-10,rtol=1e-12,equal_nan=True)
            if len(g)!=r.eligible_sites or int(valid.sum())!=r.valid_block_draws:raise AssertionError('null denominator changed')
            checked+=1
    return checked


def verify_confound_headlines(root):
    predictions=read_table(root/'tables/covariate_baseline_all_libraries.tsv.gz')
    reference=read_table(root/'tables/covariate_baseline_metrics.tsv').set_index(['library','depth','cell','tf','fold'])
    rows=[];maximum=0.
    for key,g in predictions.groupby(['library','depth','cell','tf','fold']):
        auc,ap=metrics(g.chip_label,g.score);r=reference.loc[key]
        np.testing.assert_allclose([auc,ap],[r.auroc,r.auprc],atol=1e-12,rtol=0,equal_nan=True)
        rows.append(dict(zip(['library','depth','cell','tf','fold'],key),auroc=auc,auprc=ap))
        if np.isfinite(auc):maximum=max(maximum,abs(auc-r.auroc),abs(ap-r.auprc))
    recomputed=pd.DataFrame(rows).groupby(['library','depth','cell','tf'])[['auroc','auprc']].mean().groupby(['depth','cell','tf']).mean()
    expected=read_table(root/'tables/covariate_baseline_task_metrics.tsv').set_index(['depth','cell','tf'])
    np.testing.assert_allclose(recomputed.loc[expected.index,['auroc','auprc']],expected[['auroc','auprc']],atol=1e-12,rtol=0,equal_nan=True)
    pool=read_table(root/'tables/pool.tsv.gz').set_index('site_id');cov=read_table(root/'tables/pool_library_covariates.tsv.gz')
    balance=read_table(root/'tables/actual_library_covariate_balance.tsv')
    balance=balance[balance.depth.eq('full')&balance.phase.eq('after')&balance.covariate.eq('actual_whole201')&balance.status.eq('AVAILABLE')]
    checked=0
    for library,frame in cov[cov.depth.eq('full')&cov.cohort_nested].groupby('library'):
        frame=frame.join(pool[['tf','chip_label']],on='site_id')
        with np.load(root/f'arrays/{library}.full.npz',allow_pickle=False) as z:
            ix=pd.Index(z['site_id']).get_indexer(frame.site_id)
            np.testing.assert_array_equal(frame.actual_whole201,z['plus_observed'][ix].sum(1)+z['minus_observed'][ix].sum(1))
        for tf,g in frame.groupby('tf'):
            a=np.log1p(g[g.chip_label.eq(1)].actual_whole201.to_numpy());b=np.log1p(g[g.chip_label.eq(0)].actual_whole201.to_numpy())
            variance=((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2)
            value=(a.mean()-b.mean())/np.sqrt(variance)
            row=balance[balance.library.eq(library)&balance.tf.eq(tf)].iloc[0]
            np.testing.assert_allclose(value,row.standardized_mean_difference,atol=1e-12,rtol=0);checked+=1
    return dict(baseline_folds=len(rows),maximum_metric_error=maximum,actual_full_depth_balance_rows=checked)


def verify_shapes(root):
    """Recalculate shape summaries without importing the production operators."""
    path=root/'tables/shape_metrics.tsv'
    if not path.exists():return 'not yet summarized'
    table=read_table(path);sites=read_table(root/'tables/nested_sites.tsv.gz')
    x=np.arange(-100,101);outer=abs(x)>=60;focus=abs(x)<=50
    checked=0;maximum=0.;regional_checked=0
    regional_path=root/'tables/principal_regional_shapes.tsv'
    regional=read_table(regional_path) if regional_path.exists() else None
    for (library,depth),rows in table.groupby(['library','depth']):
        cell=library.split('_')[0];held=sites[sites.cell.eq(cell)]
        p=root/(f'arrays/{cell}.replay_pool.npz' if depth=='archived_10M' else f'arrays/{library}.{depth}.npz')
        with np.load(p,allow_pickle=False) as z:
            ix=pd.Index(z['site_id']).get_indexer(held.site_id)
            if (ix<0).any():raise AssertionError('missing primary shape arrays')
            observed=np.stack([z[s+'_observed'][ix] for s in ['plus','minus']],axis=1)
        z=load_propensities(root,library)
        ix=pd.Index(z['site_id']).get_indexer(held.site_id)
        logits=np.stack([z['LOG21_'+s][ix] for s in ['plus','minus']],axis=1)
        strand_mu=np.exp(logits-logsumexp(logits,axis=2,keepdims=True))*observed.sum(2,keepdims=True)
        mu=strand_mu.sum(1)
        y=observed.sum(1);safe=np.maximum(mu,np.finfo(float).tiny)
        first=xlogy(y,np.maximum(y,np.finfo(float).tiny)/safe)
        dev=2*(first-(y-safe)) if cell=='HepG2' else 2*(first-(y+1/1.6144318781609228)*np.log((y+1/1.6144318781609228)/(safe+1/1.6144318781609228)))
        residual=np.sign(y-mu)*np.sqrt(np.maximum(np.nan_to_num(dev),0))
        values={'raw':y,'LOG21_expected':mu,'LOG21_residual':residual}
        for r in rows.itertuples():
            if r.status!='AVAILABLE':continue
            mask=held.tf.eq(r.tf).to_numpy();data=values[r.representation][mask];labels=held.loc[mask,'chip_label'].to_numpy(bool)
            # Symmetric outer coordinates make the least-squares intercept
            # the mean and the slope the dot product divided by x'x.
            intercept=data[:,outer].mean(1);slope=(data[:,outer]@x[outer])/(x[outer]@x[outer])
            detrended=data-intercept[:,None]-slope[:,None]*x
            scale=np.sqrt((detrended[:,outer]**2).mean(1));scale=np.where(scale>1e-8,scale,1.)
            normal=np.clip(detrended/scale[:,None],-12,12)
            delta=normal[labels].mean(0)-normal[~labels].mean(0)
            center=delta[95:106].mean();left=delta[60:86].mean();right=delta[115:141].mean()
            positive=normal[labels][:,focus];negative=normal[~labels][:,focus]
            pooled=((len(positive)-1)*positive.var(0,ddof=1)+(len(negative)-1)*negative.var(0,ddof=1))/(len(positive)+len(negative)-2)
            nonzero=pooled[pooled>0];floor=.01*np.median(nonzero) if len(nonzero) else 1.
            separation=np.sqrt(np.mean((positive.mean(0)-negative.mean(0))**2/np.maximum(pooled,floor)))
            expected=np.array([center,left,right,left-center,right-center,separation])
            actual=np.array([r.center_contrast,r.left_shoulder,r.right_shoulder,r.left_minus_center,r.right_minus_center,r.standardized_functional_separation])
            np.testing.assert_allclose(expected,actual,atol=1e-9,rtol=1e-9,equal_nan=True)
            maximum=max(maximum,float(np.nanmax(abs(expected-actual))));checked+=1
        if regional is not None:
            selected=regional[regional.library.eq(library)&regional.depth.eq(depth)]
            for r in selected.itertuples():
                mask=(held.tf.eq(r.tf)&held.chip_label.eq(r.label)).to_numpy()
                if r.strand=='combined':yy=y[mask];mm=mu[mask]
                else:
                    si=0 if r.strand=='plus' else 1
                    yy=observed[mask,si];mm=strand_mu[mask,si]
                if r.representation=='raw':data=yy
                elif r.representation=='LOG21_expected':data=mm
                else:
                    mm=np.maximum(mm,np.finfo(float).tiny)
                    first=xlogy(yy,np.maximum(yy,np.finfo(float).tiny)/mm)
                    dev=2*(first-(yy-mm)) if cell=='HepG2' else 2*(first-(yy+1/1.6144318781609228)*np.log((yy+1/1.6144318781609228)/(mm+1/1.6144318781609228)))
                    data=np.sign(yy-mm)*np.sqrt(np.maximum(np.nan_to_num(dev),0))
                weights=np.ones(mask.sum()) if r.aggregation=='site_equal' else y[mask].sum(1)
                curve=np.average(data,axis=0,weights=weights) if weights.sum()>0 else np.full(201,np.nan)
                c,inner=(5,15) if r.geometry=='canonical' else (7,12)
                center=curve[abs(x)<=c].mean();left=curve[(x>=-40)&(x<=-inner)].mean();right=curve[(x>=inner)&(x<=40)].mean()
                np.testing.assert_allclose([center,left,right,left-center,right-center,right-left],
                    [r.center_mean,r.left_mean,r.right_mean,r.left_minus_center,r.right_minus_center,r.right_minus_left],atol=1e-9,rtol=1e-9,equal_nan=True)
                regional_checked+=1
    return dict(rows=checked,maximum_absolute_error=maximum,regional_geometry_rows=regional_checked)


def verify_screen(root):
    path=root/'tables/full_screen_task_deltas.tsv'
    if not path.exists():return 'not yet summarized'
    folds=read_table(root/'tables/fold_metrics.tsv');folds=folds[folds.cohort.eq('nested')&folds.depth.eq('full')]
    task=read_table(path);checked=0;boot={}
    weights=np.random.default_rng(2026).multinomial(15,np.repeat(1/15,15),size=2000)
    for r in task.itertuples():
        g=folds[folds.cell.eq(r.cell)&folds.tf.eq(r.tf)]
        a=g[g.method.eq(r.method)].set_index(['library','fold']);b=g[g.method.eq(r.reference)].set_index(['library','fold'])
        paired=pd.DataFrame({'roc':a.auroc-b.auroc,'ap':a.auprc-b.auprc,'base':b.auprc}).dropna()
        means=paired.groupby(level=0).mean()
        if len(means)!=r.libraries_available:raise AssertionError('replicate denominator changed')
        expected=[means.roc.mean(),means.ap.mean(),means.ap.mean()/means.base.mean()] if len(means)==3 else [np.nan]*3
        np.testing.assert_allclose(expected,[r.auroc_gain,r.auprc_gain,r.relative_auprc_gain],atol=1e-12,rtol=0,equal_nan=True);checked+=1
        if len(means)==3:
            draws=[]
            for library in means.index:
                delta=paired.loc[library].roc.reindex([f'chr{i}' for i in range(1,16)]).to_numpy()
                valid=np.isfinite(delta);denom=weights@valid.astype(float)
                draws.append(np.divide(weights@np.nan_to_num(delta),denom,out=np.full(2000,np.nan),where=denom>0))
            boot.setdefault((r.method,r.reference),[]).append(np.mean(draws,axis=0))
    summary=read_table(root/'tables/full_screen_summary.tsv')
    for r in summary.itertuples():
        g=task[task.method.eq(r.method)&task.reference.eq(r.reference)&task.status.eq('DEVELOPMENT_DIAGNOSTIC')]
        if len(g)!=r.available_tasks or r.original_tasks!=21:raise AssertionError('screen task denominator changed')
        np.testing.assert_allclose([g.auroc_gain.mean(),g.relative_auprc_gain.mean()],
                                  [r.mean_auroc_gain,r.mean_relative_auprc_gain],atol=1e-12,rtol=0,equal_nan=True)
        draws=boot.get((r.method,r.reference),[])
        average=np.mean(draws,axis=0) if draws else np.full(2000,np.nan)
        valid=np.isfinite(average)
        interval=np.quantile(average[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan]
        if int(valid.sum())!=r.valid_draws:raise AssertionError('undefined bootstrap draw policy changed')
        np.testing.assert_allclose(interval,[r.fixed_prediction_bootstrap_low,r.fixed_prediction_bootstrap_high],atol=1e-12,rtol=0,equal_nan=True)
    return checked


def verify_effect_oracle(root):
    # Independent SciPy PMFs with fixed Gauss-Legendre quadrature, not the
    # adaptive likelihood/integration implementation used to make scores.
    from test_shrinkage import oracle
    scores=read_table(root/'tables/per_site_scores.tsv.gz')
    scores=scores[scores.method.isin(['LOG21_binomial','LOG21_beta_binomial','DWM_binomial','DWM_beta_binomial'])&np.isfinite(scores.score)]
    chosen=scores.sort_values('site_id').groupby(['library','depth','method']).head(1)
    checked=0;maximum=0.
    for (library,depth),group in chosen.groupby(['library','depth']):
        counts=read_table(root/f'tables/counts.{library}.{depth}.tsv.gz').set_index('site_id')
        for r in group.itertuples():
            c=counts.loc[r.site_id];arm=r.method.split('_')[0]
            expected=oracle([int(c.plus_center),int(c.minus_center)],[int(c.plus_total),int(c.minus_total)],
                            [c[arm+'_plus_p'],c[arm+'_minus_p']],r.rho)
            actual=[r.score,r.probability_positive,r.low,r.median,r.high]
            np.testing.assert_allclose(expected,actual,atol=1e-7,rtol=1e-7)
            maximum=max(maximum,float(np.max(np.abs(np.array(expected)-actual))));checked+=1
    return dict(independent_site_effects=checked,maximum_absolute_error=maximum)


def verify_cohorts(root):
    pool=read_table(root/'tables/pool.tsv.gz');held=read_table(root/'tables/nested_sites.tsv.gz')
    folds=read_table(root/'tables/cohort_folds.tsv')
    if pool.site_id.duplicated().any() or held.site_id.duplicated().any():raise AssertionError('duplicate site identity')
    if len(pool)!=35006 or len(held)!=6020:raise AssertionError('frozen cohort size changed')
    source=pool.set_index('site_id').loc[held.site_id]
    for key in ['TFBS_chr','TFBS_start','TFBS_end','TFBS_strand','chip_label','motif_score']:
        np.testing.assert_array_equal(source[key].to_numpy(),held[key].to_numpy())
    scales={}
    for r in folds.itertuples():
        fitting=pool[pool.cell.eq(r.cell)&pool.tf.eq(r.tf)&pool.supported&pool.TFBS_chr.ne(r.fold)]
        x=np.column_stack([fitting.motif_score,np.log1p(fitting.central_201),np.log1p(fitting.outer_101_500)])
        scale=x.std(0) if len(x) else np.ones(3);scale=np.where(scale>0,scale,1.)
        np.testing.assert_allclose(scale,[r.scale_motif,r.scale_central,r.scale_outer],atol=1e-12,rtol=1e-12)
        scales[(r.cell,r.tf,r.fold)]=scale
    for pair,group in held.groupby('pair_id'):
        if len(group)!=2 or group.chip_label.sum()!=1 or group.TFBS_chr.nunique()!=1:raise AssertionError('broken matched pair')
        row=group.iloc[0];scale=scales[(row.cell,row.tf,row.TFBS_chr)]
        x=np.column_stack([group.motif_score,np.log1p(group.central_201),np.log1p(group.outer_101_500)])
        distance=np.linalg.norm((x[0]-x[1])/scale)
        np.testing.assert_allclose(group.match_distance,distance,atol=1e-12,rtol=1e-10)
        if distance>.5+1e-12:raise AssertionError('matching caliper changed')
    receipts=[]
    for path in (root/'provenance').glob('*.fragment_links.json'):
        doc=json.loads(path.read_text());receipts.append(doc)
        if not doc['cross_chromosome_fit_leakage_qualified']:raise AssertionError('shared fragments cross fitting/evaluation chromosome')
    return dict(pool_sites=len(pool),held_sites=len(held),matched_pairs=held.pair_id.nunique(),fitting_scale_records=len(folds),fragment_link_libraries=len(receipts),
                scope='independent fitting-only scales, pair integrity and caliper; no rematching or refitting in verifier')


def main():
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--working',action='store_true',help='pre-package scientific verification only, no archive checksum assertion')
    p.add_argument('--figure-check',action='store_true');a=p.parse_args();root=a.root.resolve()
    result={'checksums':verify_checksums(root) if not a.working else 'not packaged yet',
            'counts':verify_counts(root),'metrics':verify_metrics(root),'simulation_rows':verify_simulation(root),
            'control_rate_rows':verify_control_rates(root),'shape':verify_shapes(root),'screen_task_rows':verify_screen(root),
            'effect_oracle':verify_effect_oracle(root),'cohorts':verify_cohorts(root),'calibration_order_records':verify_calibration_order(root),
            'null_mean_rows':verify_null_means(root),'confounding':verify_confound_headlines(root)}
    if not a.working:result['completion_records']=verify_completion_records(root)
    if a.figure_check:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        # Independent tiny figure from complete observation arrays, not a
        # repair of the submitted figure or table.
        sites=read_table(root/'tables/nested_sites.tsv.gz')
        sites=sites[sites.cell.eq('K562')&sites.tf.eq('CTCF')]
        with np.load(root/'arrays/K562_rep1.full.npz',allow_pickle=False) as z:
            index=pd.Index(z['site_id']).get_indexer(sites.site_id)
            data=z['plus_observed'][index]+z['minus_observed'][index]
        source=read_table(root/'tables/representative_shape.tsv')
        with tempfile.TemporaryDirectory() as directory:
            fig,ax=plt.subplots()
            for label in [0,1]:
                curve=data[sites.chip_label.to_numpy()==label].mean(0)
                np.testing.assert_array_equal(curve,source[source.label.eq(label)].sort_values('position').value)
                ax.plot(np.arange(-100,101),curve,label=str(label))
            ax.legend();ax.set_ylabel('cuts per site per base');ax.set_xlabel('base pairs from motif centre')
            target=Path(directory)/'independent_pool_profile.png';fig.savefig(target);plt.close(fig)
            if target.stat().st_size<1000:raise AssertionError('independent figure missing')
            result['independent_figure']='recomputed complete held K562 CTCF full-depth curves, matched all402 source values, independently rendered in a temporary directory'
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':main()
