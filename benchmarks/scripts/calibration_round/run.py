"""Execute the frozen calibration-round plan from next-round derived data.

Stages: naked-DNA recalibration and control endpoints, chromatin scoring on the
original nested cohort, discrimination metrics, regional shapes, synthetic
challenge re-scoring and the single advancement screen. Outputs are written
once; an existing output file stops execution instead of being overwritten.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.metrics import average_precision_score, roc_auc_score

from core import apply_recalibration, binomial_loglik, call_status, count_tail, fit_recalibration
from inputs import Access, check_frozen_code, digest, dump
from protocol import (ABLATION, ARMS, CHROMOSOMES, COMPARATOR, DEPTHS, LIBRARIES, POLICY, PRIMARY,
                      REFERENCE, THRESHOLD)

STRANDS = ['plus', 'minus']
COUNT_STRATA = ([-1, 0, 1, 4, 9, 19, 49, 99, np.inf], ['0', '1', '2-4', '5-9', '10-19', '20-49', '50-99', '100+'])
P_BINS = [0, .025, .05, .1, .15, .2, .25, .3, .4, 1.000001]
SEED = POLICY['bootstrap_seed']
DRAWS = POLICY['bootstrap_draws']


def read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip', low_memory=False)


def write(frame, path):
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    frame.to_csv(path, sep='\t', index=False)


def block_weights():
    return np.random.default_rng(SEED).multinomial(15, np.repeat(1 / 15, 15), size=DRAWS)


def wilson_upper(k, n):
    if not n:
        return np.nan
    z = norm.ppf(.975); p = k / n; d = 1 + z * z / n
    return (p + z * z / (2 * n)) / d + z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d


# ---------------------------------------------------------------- recalibration
def strand_rows(frame):
    central = np.concatenate([frame.plus_center.to_numpy(), frame.minus_center.to_numpy()])
    total = np.concatenate([frame.plus_total.to_numpy(), frame.minus_total.to_numpy()])
    p = np.concatenate([frame.LOG21_plus_p.to_numpy(), frame.LOG21_minus_p.to_numpy()])
    return central, total, p


def fit_folds(controls):
    training = controls[controls.eligible & controls.LOG21_plus_valid & controls.LOG21_minus_valid
                        & controls.deduplicated_fit_representative]
    records = []
    for fold in CHROMOSOMES + ['deployable_all_chr1_15']:
        rows = training if fold.startswith('deployable') else training[training.TFBS_chr.ne(fold)]
        central, total, p = strand_rows(rows)
        for arm, slope in [(PRIMARY, True), (ABLATION, False)]:
            fit = fit_recalibration(central, total, p, slope=slope)
            records.append(dict(fold=fold, arm=arm, fitting_windows=len(rows),
                                fitting_chromosomes=sorted(rows.TFBS_chr.unique()), **fit,
                                role='deployable record, never evaluated' if fold.startswith('deployable') else 'held-fold evaluation'))
    return records


def coefficient_lookup(records):
    return {(r['arm'], r['fold']): (r['a'], r['b']) for r in records}


def arm_probabilities(frame, arm, coefficients):
    if arm == 'uniform_count_tail':
        value = np.full(len(frame), 15 / 73)
        return value, value.copy()
    if arm == 'DWM_count_tail':
        return frame.DWM_plus_p.to_numpy(float), frame.DWM_minus_p.to_numpy(float)
    plus, minus = frame.LOG21_plus_p.to_numpy(float), frame.LOG21_minus_p.to_numpy(float)
    if arm == REFERENCE:
        return plus, minus
    a = frame.TFBS_chr.map(lambda c: coefficients[(arm, c)][0]).to_numpy(float)
    b = frame.TFBS_chr.map(lambda c: coefficients[(arm, c)][1]).to_numpy(float)
    return apply_recalibration(plus, a, b), apply_recalibration(minus, a, b)


def arm_valid(frame, arm):
    if arm == 'uniform_count_tail':
        return np.ones(len(frame), bool)
    if arm == 'DWM_count_tail':
        return (frame.DWM_plus_valid & frame.DWM_minus_valid).to_numpy(bool)
    return (frame.LOG21_plus_valid & frame.LOG21_minus_valid).to_numpy(bool)


def score_chunk(args):
    n_plus, n_minus, central, q_plus, q_minus = args
    return [count_tail(int(a), int(b), int(c), float(d), float(e)) if np.isfinite([d, e]).all() else (np.nan, np.nan)
            for a, b, c, d, e in zip(n_plus, n_minus, central, q_plus, q_minus)]


def score_frame(frame, coefficients, pool, base_columns):
    out = []
    central = (frame.plus_center + frame.minus_center).to_numpy(int)
    for arm in ARMS:
        q_plus, q_minus = arm_probabilities(frame, arm, coefficients)
        valid = arm_valid(frame, arm)
        q_plus = np.where(valid, q_plus, np.nan); q_minus = np.where(valid, q_minus, np.nan)
        chunks = np.array_split(np.arange(len(frame)), max(1, len(frame) // 250))
        jobs = [(frame.plus_total.to_numpy(int)[c], frame.minus_total.to_numpy(int)[c], central[c], q_plus[c], q_minus[c])
                for c in chunks]
        values = [v for part in pool.map(score_chunk, jobs) for v in part]
        result = frame[base_columns].copy()
        result['arm'] = arm
        result['q_plus'] = q_plus; result['q_minus'] = q_minus
        result['expected_center'] = frame.plus_total.to_numpy() * q_plus + frame.minus_total.to_numpy() * q_minus
        result['score'] = [v[0] for v in values]
        result['ceiling'] = [v[1] for v in values]
        result['status'] = [call_status(t, s, c) for t, s, c in zip(frame.total, result.score, result.ceiling)]
        result['call'] = result.status.eq('CALL')
        result['exclusion_reason'] = np.where(valid, '', 'sequence propensity unsupported for this arm')
        if arm in (PRIMARY, ABLATION):
            result['calibration_provenance'] = [
                f'naked_DNA_rep1_full:exclude_{c}:a={coefficients[(arm, c)][0]:.12g}:b={coefficients[(arm, c)][1]:.12g}'
                for c in frame.TFBS_chr]
        else:
            result['calibration_provenance'] = 'none (fixed propensities)'
        out.append(result)
    return pd.concat(out, ignore_index=True)


# ---------------------------------------------------------------- control endpoints
def control_stage(access, out, pool):
    controls = read(access.path('tables/control_counts.full.tsv.gz'))
    with np.load(access.path('arrays/NakedDNA_rep1.propensities.npz'), allow_pickle=False) as z:
        if not np.array_equal(z['site_id'], controls.site_id.to_numpy(str)):
            raise ValueError('control propensity order changed')
        controls['gc'] = z['gc']
    if len(controls) != 4200 or controls.groupby(['cell', 'tf']).size().ne(200).any():
        raise ValueError('control population changed')
    records = fit_folds(controls)
    dump(out / 'models/recalibration.json', dict(records=records, training_source='naked-DNA replicate 1 full depth',
                                                    labels_used=False, chromatin_used=False))
    coefficients = coefficient_lookup(records)
    support = (controls.eligible & controls.LOG21_plus_valid & controls.LOG21_minus_valid
               & controls.DWM_plus_valid & controls.DWM_minus_valid).to_numpy(bool)
    reason = np.where(~controls.eligible, 'ineligible: blacklist or unavailable model',
                      np.where(~(controls.LOG21_plus_valid & controls.LOG21_minus_valid), 'LOG21 unsupported',
                               np.where(~(controls.DWM_plus_valid & controls.DWM_minus_valid), 'DWM unsupported', '')))
    controls['common_support'] = support
    controls['support_exclusion'] = reason
    base = ['site_id', 'cell', 'tf', 'TFBS_chr', 'physical_id', 'deduplicated_fit_representative', 'eligible',
            'common_support', 'support_exclusion', 'gc', 'LOG21_plus_p', 'LOG21_minus_p',
            'plus_center', 'minus_center', 'plus_total', 'minus_total', 'total']
    scored = score_frame(controls, coefficients, pool, base)
    scored.loc[~scored.common_support, 'status'] = 'EXCLUDED'
    scored.loc[~scored.common_support, 'call'] = False
    write(scored, out / 'tables/control_site_scores.tsv.gz')

    weights = block_weights()
    loglik, variance, curves, means, calls = [], [], [], [], []
    dedup = scored[scored.common_support & scored.deduplicated_fit_representative]
    gc_edges = np.quantile(controls.gc, [1 / 3, 2 / 3])
    for arm in ARMS:
        rows = dedup[dedup.arm.eq(arm)]
        c = np.concatenate([rows.plus_center, rows.minus_center]).astype(float)
        n = np.concatenate([rows.plus_total, rows.minus_total]).astype(float)
        q = np.concatenate([rows.q_plus, rows.q_minus])
        chrom = np.concatenate([rows.TFBS_chr, rows.TFBS_chr])
        keep = n > 0
        loglik.append(dict(arm=arm, strand_observations=int(keep.sum()), central=int(c[keep].sum()),
                           expected=float((n * q)[keep].sum()), heldout_loglik=binomial_loglik(c[keep], n[keep], q[keep]),
                           **{f'loglik_{f}': binomial_loglik(c[keep & (chrom == f)], n[keep & (chrom == f)], q[keep & (chrom == f)])
                              for f in CHROMOSOMES}))
        pearson = (c - n * q) ** 2 / (n * q * (1 - q))
        for label, mask in [('all', keep)] + [(f'n={lo}-{hi}', keep & (n >= lo) & (n <= hi)) for lo, hi in [(1, 2), (3, 5), (6, 10), (11, 10**9)]]:
            variance.append(dict(arm=arm, stratum=label, observations=int(mask.sum()),
                                 pearson_dispersion=float(pearson[mask].sum() / mask.sum()) if mask.any() else np.nan))
        p21 = np.concatenate([rows.LOG21_plus_p, rows.LOG21_minus_p])
        gc = np.concatenate([rows.gc, rows.gc])
        for kind, values, edges in [('LOG21_p_bin', p21, P_BINS), ('gc_tertile', gc, [-np.inf, *gc_edges, np.inf])]:
            index = np.digitize(values, edges) - 1
            for k in range(len(edges) - 1):
                mask = keep & (index == k)
                curves.append(dict(arm=arm, kind=kind, bin=k, low=edges[k], high=edges[k + 1], observations=int(mask.sum()),
                                   cuts=float(n[mask].sum()), observed_central=float(c[mask].sum()),
                                   expected_central=float((n * q)[mask].sum()),
                                   observed_share=float(c[mask].sum() / n[mask].sum()) if mask.any() else np.nan,
                                   expected_share=float((n * q)[mask].sum() / n[mask].sum()) if mask.any() else np.nan))
        for (cell, tf), task in scored[scored.arm.eq(arm)].groupby(['cell', 'tf']):
            sup = task[task.common_support]
            o_chr = sup.groupby('TFBS_chr').apply(lambda g: (g.plus_center + g.minus_center).sum(), include_groups=False).reindex(CHROMOSOMES, fill_value=0).to_numpy(float)
            e_chr = sup.groupby('TFBS_chr').expected_center.sum().reindex(CHROMOSOMES, fill_value=0).to_numpy(float)
            O, E = o_chr.sum(), e_chr.sum()
            ow, ew = weights @ o_chr, weights @ e_chr
            ratio = np.divide(ow, ew, out=np.full(DRAWS, np.nan), where=ew > 0)
            valid = np.isfinite(ratio)
            low, high = np.quantile(ratio[valid], [.025, .975]) if valid.sum() >= 1900 else (np.nan, np.nan)
            means.append(dict(arm=arm, cell=cell, tf=tf, sites=len(task), supported_sites=len(sup),
                              excluded_sites=int((~task.common_support).sum()), observed_central=O, expected_central=E,
                              ratio=O / E if E > 0 else np.nan, ratio_low=low, ratio_high=high, valid_draws=int(valid.sum()),
                              discrepant=bool(np.isfinite(low) and (low > 1 or high < 1))))
            informative = sup[sup.total > 0]
            k = int(informative.call.sum())
            calls.append(dict(arm=arm, cell=cell, tf=tf, sites=len(task), eligible=len(sup), informative=len(informative),
                              reachable=int((informative.ceiling > THRESHOLD).sum()), calls=k,
                              eligible_rate=k / len(sup) if len(sup) else np.nan,
                              informative_rate=k / len(informative) if len(informative) else np.nan,
                              eligible_wilson_upper=wilson_upper(k, len(sup)),
                              informative_wilson_upper=wilson_upper(k, len(informative))))
    calls = pd.DataFrame(calls)
    dwm = calls[calls.arm.eq(COMPARATOR)].set_index(['cell', 'tf'])
    calls['paired_increase_over_DWM'] = [
        (r.calls - dwm.loc[(r.cell, r.tf), 'calls']) / r.informative if r.informative else np.nan for r in calls.itertuples()]
    write(pd.DataFrame(loglik), out / 'tables/control_heldout_loglik.tsv')
    write(pd.DataFrame(variance), out / 'tables/control_variance.tsv')
    write(pd.DataFrame(curves), out / 'tables/control_calibration_curves.tsv')
    write(pd.DataFrame(means), out / 'tables/control_task_means.tsv')
    write(calls, out / 'tables/control_calls.tsv')
    return controls, coefficients


# ---------------------------------------------------------------- chromatin scoring
def chromatin_stage(access, out, pool, coefficients):
    nested = read(access.path('tables/nested_sites.tsv.gz'))
    old = read(access.path('tables/per_site_scores.tsv.gz'))
    old = old[old.cohort.eq('nested') & old.method.isin(['LOG21_count_tail', 'uniform_count_tail'])]
    base = ['site_id', 'cell', 'tf', 'TFBS_chr', 'chip_label', 'pair_id', 'motif_score', 'central_201', 'library', 'depth',
            'plus_center', 'minus_center', 'plus_total', 'minus_total', 'total', 'count_stratum']
    frames, parity = [], []
    for library in LIBRARIES:
        cell = library.split('_')[0]
        expected = nested[nested.cell.eq(cell)]
        for depth in DEPTHS:
            counts = read(access.path(f'tables/counts.{library}.{depth}.tsv.gz'))
            sites = counts[counts.cohort_nested].merge(expected[['site_id', 'pair_id']], on='site_id', validate='one_to_one')
            if len(sites) != len(expected) or set(sites.site_id) != set(expected.site_id):
                raise ValueError(f'nested membership changed {library} {depth}')
            if not (sites.total == sites.plus_total + sites.minus_total).all():
                raise ValueError('regional totals inconsistent')
            sites = sites.sort_values(['tf', 'TFBS_chr', 'site_id']).reset_index(drop=True)
            result = score_frame(sites, coefficients, pool, base)
            for arm in ['LOG21_count_tail', 'uniform_count_tail']:
                mine = result[result.arm.eq(arm)].set_index('site_id')
                theirs = old[old.library.eq(library) & old.depth.eq(depth) & old.method.eq(arm)].set_index('site_id')
                theirs = theirs.loc[mine.index]
                np.testing.assert_array_equal(mine.score.to_numpy(), theirs.score.to_numpy())
                np.testing.assert_array_equal(mine.ceiling.to_numpy(), theirs.attainable_old_score_ceiling.to_numpy())
                np.testing.assert_array_equal(mine.call.to_numpy(), theirs.call.to_numpy())
                parity.append(dict(library=library, depth=depth, arm=arm, sites=len(mine), score_mismatches=0,
                                   ceiling_mismatches=0, call_mismatches=0))
            frames.append(result)
            print('scored', library, depth, flush=True)
    scores = pd.concat(frames, ignore_index=True)
    write(scores, out / 'tables/chromatin_site_scores.tsv.gz')
    write(pd.DataFrame(parity), out / 'tables/reference_parity.tsv')
    return scores, nested


# ---------------------------------------------------------------- metrics
def stratified_auc(labels, scores, strata):
    concordant = pairs = 0.
    for level in np.unique(strata):
        m = strata == level
        pos, neg = scores[m & (labels == 1)], scores[m & (labels == 0)]
        if len(pos) and len(neg):
            diff = pos[:, None] - neg[None, :]
            concordant += (diff > 0).sum() + .5 * (diff == 0).sum()
            pairs += len(pos) * len(neg)
    return concordant / pairs if pairs else np.nan, pairs


def fold_metrics(scores, covariates, task_ledger):
    rows = []
    covariates = covariates.rename(columns={'score': 'cov_score'})
    for (library, depth), group in scores.groupby(['library', 'depth']):
        cell = library.split('_')[0]
        cov = covariates[covariates.library.eq(library) & covariates.depth.eq(depth)]
        for t in task_ledger[task_ledger.cell.eq(cell)].itertuples():
            task = group[group.tf.eq(t.tf)]
            arms = list(ARMS) + ['covariates_only']
            for arm in arms:
                for fold in CHROMOSOMES:
                    if arm == 'covariates_only':
                        base = task[task.arm.eq(REFERENCE) & task.TFBS_chr.eq(fold)][['site_id', 'chip_label', 'count_stratum']]
                        sub = base.merge(cov[['site_id', 'cov_score']], on='site_id', how='left').rename(columns={'cov_score': 'score'})
                    else:
                        sub = task[task.arm.eq(arm) & task.TFBS_chr.eq(fold)]
                    finite = sub[np.isfinite(sub.score)]
                    ok = len(sub) > 0 and len(finite) == len(sub) and finite.chip_label.nunique() == 2
                    strat, pairs = stratified_auc(finite.chip_label.to_numpy(int), finite.score.to_numpy(float),
                                                  finite.count_stratum.to_numpy(str)) if ok else (np.nan, 0)
                    rows.append(dict(library=library, depth=depth, cell=cell, tf=t.tf, arm=arm, fold=fold, sites=len(sub),
                                     finite_sites=len(finite), positives=int(finite.chip_label.sum()),
                                     auroc=roc_auc_score(finite.chip_label, finite.score) if ok else np.nan,
                                     auprc=average_precision_score(finite.chip_label, finite.score) if ok else np.nan,
                                     count_stratified_auroc=strat, stratified_pairs=pairs,
                                     status='AVAILABLE' if ok else ('NO_SITES' if not len(sub) else 'INCOMPLETE')))
    return pd.DataFrame(rows)


def screen_pairs(folds, task_ledger, depth, method, reference):
    weights = block_weights()
    folds = folds[folds.depth.eq(depth)]
    rows, boot, points = [], [], []
    for t in task_ledger.itertuples():
        sub = folds[folds.cell.eq(t.cell) & folds.tf.eq(t.tf)]
        a = sub[sub.arm.eq(method)]; b = sub[sub.arm.eq(reference)]
        matched = a.merge(b, on=['library', 'fold'], suffixes=('', '_ref'))
        gains, aps, bases, draws, libs = [], [], [], [], []
        for library in [f'{t.cell}_rep{r}' for r in (1, 2, 3)]:
            g = matched[matched.library.eq(library)].set_index('fold').reindex(CHROMOSOMES)
            delta = (g.auroc - g.auroc_ref).to_numpy(); ap = (g.auprc - g.auprc_ref).to_numpy()
            valid = np.isfinite(delta) & np.isfinite(ap)
            gains.append(np.mean(delta[valid]) if valid.any() else np.nan)
            aps.append(np.mean(ap[valid]) if valid.any() else np.nan)
            bases.append(g.auprc_ref[valid].mean() if valid.any() else np.nan)
            libs.append(np.mean(delta[valid]) if valid.any() else np.nan)
            denom = weights @ valid.astype(float)
            draws.append(np.divide(weights @ np.where(valid, delta, 0.), denom, out=np.full(DRAWS, np.nan), where=denom > 0))
        complete = bool(np.isfinite(gains).all())
        value = float(np.mean(gains)) if complete else np.nan
        apvalue = float(np.mean(aps)) if complete else np.nan
        base = float(np.mean(bases)) if complete else np.nan
        rows.append(dict(depth=depth, cell=t.cell, tf=t.tf, method=method, reference=reference,
                         libraries_available=int(np.isfinite(gains).sum()), auroc_gain=value, auprc_gain=apvalue,
                         relative_auprc_gain=apvalue / base if complete and base > 0 else np.nan,
                         library_gains=json.dumps([None if not np.isfinite(x) else float(x) for x in libs]),
                         status='DEVELOPMENT_DIAGNOSTIC' if complete else 'UNAVAILABLE'))
        if complete:
            boot.append(np.mean(np.stack(draws), axis=0)); points.append((value, apvalue / base))
    matrix = np.stack(boot) if boot else np.full((1, DRAWS), np.nan)
    draw = matrix.mean(0); valid = np.isfinite(draw)
    low, high = np.quantile(draw[valid], [.025, .975]) if valid.sum() >= 1900 else (np.nan, np.nan)
    summary = dict(depth=depth, method=method, reference=reference, available_tasks=len(points), original_tasks=21,
                   mean_auroc_gain=float(np.mean([p[0] for p in points])) if points else np.nan,
                   mean_relative_auprc_gain=float(np.mean([p[1] for p in points])) if points else np.nan,
                   bootstrap_low=low, bootstrap_high=high, valid_draws=int(valid.sum()),
                   scope='fixed predictions; synchronized chromosome blocks; not refitting uncertainty')
    return rows, summary


def arm_means(folds, task_ledger):
    rows = []
    for (depth, arm), group in folds.groupby(['depth', 'arm']):
        per_task = []
        for t in task_ledger.itertuples():
            g = group[group.cell.eq(t.cell) & group.tf.eq(t.tf) & group.status.eq('AVAILABLE')]
            lib = g.groupby('library')[['auroc', 'auprc', 'count_stratified_auroc']].mean()
            ok = len(lib) == 3
            per_task.append(dict(depth=depth, arm=arm, cell=t.cell, tf=t.tf, libraries=len(lib),
                                 auroc=lib.auroc.mean() if ok else np.nan, auprc=lib.auprc.mean() if ok else np.nan,
                                 count_stratified_auroc=lib.count_stratified_auroc.mean() if ok else np.nan,
                                 status='AVAILABLE' if ok else 'UNAVAILABLE'))
        rows.extend(per_task)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- shape and information
def shape_stage(scores):
    weights = block_weights(); rows = []
    for (cell, tf, depth, arm), g in scores[scores.depth.isin(DEPTHS)].groupby(['cell', 'tf', 'depth', 'arm']):
        g = g[np.isfinite(g.expected_center)]
        stats = {}
        for label in (1, 0):
            h = g[g.chip_label.eq(label)]
            o = h.groupby('TFBS_chr').apply(lambda x: (x.plus_center + x.minus_center).sum(), include_groups=False).reindex(CHROMOSOMES, fill_value=0).to_numpy(float)
            e = h.groupby('TFBS_chr').expected_center.sum().reindex(CHROMOSOMES, fill_value=0).to_numpy(float)
            stats[label] = (o, e)
        (op, ep), (on, en) = stats[1], stats[0]
        rp = np.divide(weights @ op, weights @ ep, out=np.full(DRAWS, np.nan), where=weights @ ep > 0)
        rn = np.divide(weights @ on, weights @ en, out=np.full(DRAWS, np.nan), where=weights @ en > 0)
        rr = np.divide(rp, rn, out=np.full(DRAWS, np.nan), where=rn > 0)
        def interval(x):
            v = np.isfinite(x)
            return np.quantile(x[v], [.025, .975]) if v.sum() >= 1900 else (np.nan, np.nan)
        lp, hp = interval(rp); lr, hr = interval(rr)
        pos = op.sum() / ep.sum() if ep.sum() > 0 else np.nan
        neg = on.sum() / en.sum() if en.sum() > 0 else np.nan
        rows.append(dict(cell=cell, tf=tf, depth=depth, arm=arm, libraries_pooled=g.library.nunique(),
                         positive_sites=int(g[g.chip_label.eq(1)].site_id.nunique()), negative_sites=int(g[g.chip_label.eq(0)].site_id.nunique()),
                         positive_observed=op.sum(), positive_expected=ep.sum(), negative_observed=on.sum(), negative_expected=en.sum(),
                         positive_ratio=pos, positive_low=lp, positive_high=hp, negative_ratio=neg,
                         positive_over_negative=pos / neg if neg > 0 else np.nan, pon_low=lr, pon_high=hr,
                         canonical_protection=bool(np.isfinite(hp) and hp < 1 and np.isfinite(hr) and hr < 1)))
    return pd.DataFrame(rows)


def information_stage(scores):
    rows = []
    for (cell, tf, depth, arm, label), g in scores.groupby(['cell', 'tf', 'depth', 'arm', 'chip_label']):
        g = g.drop_duplicates(['site_id', 'library'])
        counts = g.status.value_counts()
        rows.append(dict(cell=cell, tf=tf, depth=depth, arm=arm, chip_label=label, site_library_rows=len(g),
                         median_total=float(g.total.median()), median_expected_center=float(g.expected_center.median()),
                         **{f'status_{k}': int(counts.get(k, 0)) for k in ['CALL', 'NO_CALL', 'ABSTAIN_UNREACHABLE', 'NO_INFORMATION', 'UNAVAILABLE']},
                         abstain_fraction=float(g.status.isin(['ABSTAIN_UNREACHABLE', 'NO_INFORMATION']).mean())))
    return pd.DataFrame(rows)


def power_stage(access, controls, coefficients):
    sim = read(access.path('tables/empirical_power.tsv'))
    sim = sim[['block', 'site_id']].drop_duplicates().set_index('block').site_id
    lookup = controls.set_index('site_id')
    rows = []
    for block in range(138):
        with np.load(access.path(f'arrays/empirical_simulations/block_{block:03d}.npz'), allow_pickle=False) as z:
            counts = z['counts']; rs = z['r']; site = str(z['site_id'])
        if sim.loc[block] != site:
            raise ValueError('synthetic block identity changed')
        if site not in lookup.index:
            raise ValueError('synthetic window missing from controls')
        row = lookup.loc[[site]]
        for arm in ARMS:
            qp, qm = arm_probabilities(row.reset_index(), arm, coefficients)
            qp, qm = float(qp[0]), float(qm[0])
            cache = {}
            for i, r in enumerate(rs):
                if (counts[i] < 0).any():
                    rows.append(dict(block=block, site_id=site, arm=arm, r=r, status='UNAVAILABLE_NO_SOURCE_OR_WEIGHT'))
                    continue
                calls = reach = info = 0
                for draw in counts[i]:
                    n = draw.sum(1); c = draw[:, 0].sum(); key = (int(n[0]), int(n[1]), int(c))
                    if key not in cache:
                        cache[key] = count_tail(*key, qp, qm)
                    s, ceil = cache[key]
                    info += int(n.sum() > 0); calls += int(n.sum() > 0 and ceil > THRESHOLD and s > THRESHOLD)
                    reach += int(n.sum() > 0 and ceil > THRESHOLD)
                rows.append(dict(block=block, site_id=site, cell=row.cell.iloc[0], tf=row.tf.iloc[0], chromosome=row.TFBS_chr.iloc[0],
                                 arm=arm, r=float(r), status='AVAILABLE', draws=len(counts[i]), informative_draws=info,
                                 reachable_draws=reach, calls=calls, call_rate=calls / len(counts[i]),
                                 scope='one examined source window; resamples are not independent libraries'))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- screen
def advancement(out, means, calls, summaries, deltas, shape, task_ledger):
    s = {(r['depth'], r['method'], r['reference']): r for r in summaries}
    ref_disc = int(means[means.arm.eq(REFERENCE)].discrepant.sum())
    pri_disc = int(means[means.arm.eq(PRIMARY)].discrepant.sum())
    ll = read(out / 'tables/control_heldout_loglik.tsv').set_index('arm').heldout_loglik
    gates = []
    def gate(key, passed, evidence):
        gates.append(dict(gate=key, status='PASS' if passed else 'FAIL', evidence=evidence))
    gate('A_calibration', ll[PRIMARY] > ll[REFERENCE] and pri_disc <= 3 and pri_disc < ref_disc,
         f'held-out loglik primary {ll[PRIMARY]:.3f} vs reference {ll[REFERENCE]:.3f}; discrepant tasks {pri_disc}/21 vs {ref_disc}/21')
    main = s[('full', PRIMARY, COMPARATOR)]
    gate('B_accuracy', main['mean_auroc_gain'] >= .03 and main['mean_relative_auprc_gain'] >= .10 and main['bootstrap_low'] > 0,
         f"AUROC gain {main['mean_auroc_gain']:.5f} [{main['bootstrap_low']:.5f}, {main['bootstrap_high']:.5f}], "
         f"relative AP {main['mean_relative_auprc_gain']:.4%}, {main['available_tasks']}/21 tasks")
    vs_ref = s[('full', PRIMARY, REFERENCE)]
    gate('C_reference', vs_ref['mean_auroc_gain'] >= -.01, f"AUROC change versus LOG21 count-tail {vs_ref['mean_auroc_gain']:.5f}")
    full = deltas[deltas.depth.eq('full') & deltas.method.eq(PRIMARY) & deltas.reference.eq(COMPARATOR)]
    ctcf = full[full.tf.eq('CTCF')]
    gate('D_positive_control', len(ctcf) == 2 and (ctcf.auroc_gain >= -.02).all(),
         '; '.join(f'{r.cell} CTCF {r.auroc_gain:.5f}' for r in ctcf.itertuples()))
    fam = {tf: f for f, members in POLICY['difficult_families'].items() for tf in members}
    fam['CTCF'] = 'CTCF'
    improved = full[full.auroc_gain > 0]
    families = improved.tf.map(fam).nunique(); cells = improved.cell.nunique()
    gate('E_contexts', families >= 2 and cells >= 2, f'{families} families and {cells} cells with positive task gain')
    sh = shape[shape.depth.eq('full') & shape.arm.eq(PRIMARY) & shape.canonical_protection & ~shape.tf.eq('CTCF')]
    shfam = sh.tf.map(fam).nunique()
    gate('F_shape', shfam >= 2, f'canonical regional protection in {shfam} difficult families: ' + ', '.join(f'{r.cell} {r.tf}' for r in sh.itertuples()))
    # The fixed count-tail detector needs no chromatin training, so all 21
    # control tasks are evaluable here, including HepG2 HMGA1.
    available = calls[calls.arm.eq(PRIMARY)]
    ok = (available.informative_wilson_upper <= .05).all() and (available.paired_increase_over_DWM <= .01).all()
    gate('G_controls', bool(ok), f'max informative Wilson upper {available.informative_wilson_upper.max():.4f}; '
         f'max paired increase {available.paired_increase_over_DWM.max():.4f}; examined replicate-1, not untouched safety')
    direction = []
    for depth in DEPTHS:
        d = deltas[deltas.depth.eq(depth) & deltas.method.eq(PRIMARY) & deltas.reference.eq(COMPARATOR) & deltas.status.eq('DEVELOPMENT_DIAGNOSTIC')]
        libs = np.array([json.loads(x) for x in d.library_gains], float)
        cell = d.cell.to_numpy()
        per_library = [np.nanmean(libs[cell == c, r]) for c in ['HepG2', 'K562'] for r in range(3)]
        direction.append(sum(v > 0 for v in per_library))
    gate('H_direction', all(v >= 5 for v in direction), f'positive libraries full {direction[0]}/6, 10M {direction[1]}/6')
    ablation_gap = ll[PRIMARY] - ll[ABLATION]
    gates.append(dict(gate='ablation_preference', status='PREFER_ABLATION' if ablation_gap < 1.92 else 'PRIMARY_JUSTIFIED',
                      evidence=f'held-out loglik primary minus ablation {ablation_gap:.3f}'))
    g = pd.DataFrame(gates)
    core = g[g.gate.str.match('[A-H]_')]
    if (core.status == 'PASS').all():
        decision = POLICY['advancement_screen']['outcomes']['all']
    elif g.set_index('gate').status['A_calibration'] == 'PASS':
        decision = POLICY['advancement_screen']['outcomes']['A_only']
    else:
        decision = POLICY['advancement_screen']['outcomes']['A_fails']
    g = pd.concat([g, pd.DataFrame([dict(gate='DECISION', status='', evidence=decision)])], ignore_index=True)
    write(g, out / 'tables/advancement_screen.tsv')
    return decision


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=16)
    args = parser.parse_args()
    out = args.outdir.resolve()
    plan = check_frozen_code(out)
    for sub in ['tables', 'models', 'logs']:
        (out / sub).mkdir(exist_ok=False)
    dump(out / 'logs/launch.json', dict(command=sys.argv, started_utc=datetime.now(timezone.utc).isoformat(),
                                        plan_sha256=digest(out / 'PLAN.freeze.json')))
    access = Access(out / 'inputs.json', out / 'logs/access.jsonl')
    task_ledger = read(access.path('tables/task_ledger.tsv'))
    if len(task_ledger) != 21:
        raise ValueError('task denominator changed')
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        controls, coefficients = control_stage(access, out, pool)
        print('control stage complete', flush=True)
        scores, nested = chromatin_stage(access, out, pool, coefficients)
    covariates = read(access.path('tables/covariate_baseline_all_libraries.tsv.gz'))
    folds = fold_metrics(scores, covariates, task_ledger)
    write(folds, out / 'tables/fold_metrics.tsv')
    write(arm_means(folds, task_ledger), out / 'tables/arm_task_metrics.tsv')
    deltas, summaries = [], []
    for depth in DEPTHS:
        for method, reference in [(PRIMARY, COMPARATOR), (ABLATION, COMPARATOR), (REFERENCE, COMPARATOR),
                                  (PRIMARY, REFERENCE), (ABLATION, REFERENCE), (PRIMARY, ABLATION),
                                  ('uniform_count_tail', COMPARATOR), (PRIMARY, 'covariates_only')]:
            rows, summary = screen_pairs(folds, task_ledger, depth, method, reference)
            deltas.extend(rows); summaries.append(summary)
    deltas = pd.DataFrame(deltas)
    write(deltas, out / 'tables/task_deltas.tsv')
    write(pd.DataFrame(summaries), out / 'tables/screen_summary.tsv')
    shape = shape_stage(scores)
    write(shape, out / 'tables/regional_shapes.tsv')
    write(information_stage(scores), out / 'tables/information.tsv')
    write(power_stage(access, controls, coefficients), out / 'tables/synthetic_power.tsv')
    means = read(out / 'tables/control_task_means.tsv'); calls = read(out / 'tables/control_calls.tsv')
    decision = advancement(out, means, calls, summaries, deltas, shape, task_ledger)
    outputs = {str(p.relative_to(out)): digest(p) for p in sorted(out.rglob('*')) if p.is_file()
               and p.name not in ('result.json',) and not p.name.endswith('.jsonl') and p.parent.name != 'logs'}
    dump(out / 'result.json', dict(decision=decision, plan_sha256=digest(out / 'PLAN.freeze.json'), outputs=outputs,
                                   completed_utc=datetime.now(timezone.utc).isoformat(), evidence='development only'))
    print('DECISION:', decision, flush=True)


if __name__ == '__main__':
    main()
