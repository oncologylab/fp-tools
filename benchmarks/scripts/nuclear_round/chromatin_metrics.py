"""Metric aggregation copied unchanged from the frozen calibration comparison."""
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
CHROMOSOMES=[f"chr{i}" for i in range(1,16)]
SEED=2026
DRAWS=2000
REFERENCE="LOG21_count_tail"
ARMS=["uniform_count_tail","LOG21_count_tail","DWM_count_tail","LOG21NC1_count_tail","nuclear_primary","nuclear_ablation"]

def block_weights():
    return np.random.default_rng(SEED).multinomial(15, np.repeat(1 / 15, 15), size=DRAWS)

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

