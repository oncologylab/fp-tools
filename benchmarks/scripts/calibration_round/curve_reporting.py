"""Descriptive curves with the same fold/library/task weights as scalar metrics.

Average precision is computed from original predictions, never by integrating
the displayed average curve. Tied scores are grouped by sklearn. PR lookup uses
the first achieved recall at or above the grid value, without a precision envelope.
"""
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve


def fold_curve(labels, scores, grid):
    labels, scores = np.asarray(labels), np.asarray(scores, float)
    if len(labels) == 0 or not np.isfinite(scores).all() or len(np.unique(labels)) != 2:
        return None
    fpr, tpr, _ = roc_curve(labels, scores)
    precision, recall, _ = precision_recall_curve(labels, scores)
    # Reverse sklearn's decreasing-recall arrays. searchsorted selects the first
    # threshold reaching the requested recall, rather than a favorable precision.
    r, p = recall[::-1], precision[::-1]
    pr = p[np.searchsorted(r, grid, side='left').clip(max=len(r)-1)]
    return dict(roc=np.interp(grid, fpr, tpr), pr=pr,
                auroc=roc_auc_score(labels, scores), ap=average_precision_score(labels, scores),
                prevalence=labels.mean(), sites=len(labels))


def curves(scores, covariates, ledger, arms, reference):
    grid = np.linspace(0, 1, 201)
    full = scores[scores.depth.eq('full')].copy()
    cov = covariates[covariates.depth.eq('full')][['site_id', 'library', 'score']]
    c = full[full.arm.eq(reference)].drop(columns='score').merge(
        cov, on=['site_id', 'library'], how='left', validate='one_to_one')
    c['arm'] = 'covariates_only'
    full = pd.concat([full, c], ignore_index=True)
    source, metrics, rows = [], [], {}
    for task in ledger.itertuples():
        for arm in arms:
            g = full[full.cell.eq(task.cell) & full.tf.eq(task.tf) & full.arm.eq(arm)]
            libs = []
            for library, lib in g.groupby('library', sort=True):
                folds = []
                for chrom, fold in lib.groupby('TFBS_chr', sort=True):
                    v = fold_curve(fold.chip_label, fold.score, grid)
                    metrics.append(dict(cell=task.cell, tf=task.tf, arm=arm, library=library,
                        chromosome=chrom, status='AVAILABLE' if v is not None else 'UNDEFINED',
                        **({k: v[k] for k in ['auroc', 'ap', 'prevalence', 'sites']} if v else {})))
                    if v is None:
                        continue
                    folds.append(v)
                    source.extend(dict(cell=task.cell, tf=task.tf, arm=arm, library=library,
                        chromosome=chrom, grid=x, tpr=y, precision=z)
                        for x,y,z in zip(grid, v['roc'], v['pr']))
                if folds:
                    libs.append({k: np.mean([f[k] for f in folds], axis=0)
                                 for k in ['roc', 'pr', 'auroc', 'ap', 'prevalence']})
            if not libs:
                rows[task.cell, task.tf, arm] = None
                metrics.append(dict(cell=task.cell,tf=task.tf,arm=arm,library='',chromosome='',status='UNAVAILABLE'))
            else:
                v = {k: np.mean([lib[k] for lib in libs], axis=0) for k in libs[0]}
                rows[task.cell, task.tf, arm] = dict(fpr=grid, tpr=v['roc'], rec=grid,
                    prec=v['pr'], roc_grid=v['roc'], pr_grid=v['pr'], prevalence=v['prevalence'],
                    auroc=v['auroc'], ap=v['ap'])
    return grid, rows, pd.DataFrame(source), pd.DataFrame(metrics)


def realised_power(power, empirical):
    keys = ['block', 'site_id', 'r']
    effect = empirical[keys + ['realised_marginal_center_ratio']].drop_duplicates()
    if effect.duplicated(keys).any():
        raise ValueError('conflicting empirical effect definitions')
    result = power.merge(effect, on=keys, how='left', validate='many_to_one')
    result['realised_depletion'] = 1-result.realised_marginal_center_ratio
    # Undefined effects (zero source central counts) stay in the exported ledger.
    result['effect_status'] = np.where(np.isfinite(result.realised_depletion),
                                     'DEFINED', 'UNDEFINED_ZERO_SOURCE_CENTER')
    return result
