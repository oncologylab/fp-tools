"""Read-only figures from completed calibration-round tables (no scientific fitting)."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import precision_recall_curve, roc_curve  # noqa: E402

from protocol import COMPARATOR, PRIMARY, REFERENCE  # noqa: E402

COLORS = {PRIMARY: '#2a78d6', REFERENCE: '#eb6834', COMPARATOR: '#1baf7a', 'covariates_only': '#eda100'}
LABELS = {PRIMARY: 'LOG21 + naked-DNA recalibration (candidate)', REFERENCE: 'LOG21 count-tail (reference)',
          COMPARATOR: 'DWM propensities, same detector', 'covariates_only': 'Covariates only (no central shape)'}
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 8, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': .6,
                     'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': '#fcfcfb',
                     'axes.facecolor': '#fcfcfb', 'savefig.facecolor': '#fcfcfb'})


def read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip', low_memory=False)


def save(fig, out, name):
    for ext in ['pdf', 'png']:
        fig.savefig(out / f'figures/{name}.{ext}', dpi=200, bbox_inches='tight')
    plt.close(fig)


def curve_data(scores, covariates, ledger):
    """Per task, pool the three full-depth libraries' held sites (scores are unfitted except fold-held recalibration)."""
    full = scores[scores.depth.eq('full')]
    cov = covariates[covariates.depth.eq('full')][['site_id', 'library', 'score']]
    grid = np.linspace(0, 1, 201)
    rows = {}
    for t in ledger.itertuples():
        task = full[full.cell.eq(t.cell) & full.tf.eq(t.tf)]
        for arm in [PRIMARY, REFERENCE, COMPARATOR, 'covariates_only']:
            if arm == 'covariates_only':
                g = task[task.arm.eq(REFERENCE)][['site_id', 'library', 'chip_label']].merge(cov, on=['site_id', 'library'], how='left')
            else:
                g = task[task.arm.eq(arm)]
            g = g[np.isfinite(g.score)]
            if g.chip_label.nunique() != 2:
                rows[(t.cell, t.tf, arm)] = None
                continue
            fpr, tpr, _ = roc_curve(g.chip_label, g.score)
            prec, rec, _ = precision_recall_curve(g.chip_label, g.score)
            # Stepwise PR: precision at recall r is the best precision achieved at recall >= r.
            order = np.argsort(rec)
            rec_s, prec_s = rec[order], prec[order]
            step = np.array([prec_s[rec_s >= r].max() for r in grid])
            rows[(t.cell, t.tf, arm)] = dict(fpr=fpr, tpr=tpr, rec=rec, prec=prec, roc_grid=np.interp(grid, fpr, tpr),
                                             pr_grid=step, prevalence=g.chip_label.mean())
    return grid, rows


def roc_pr(out, grid, rows, ledger):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    for arm in [PRIMARY, REFERENCE, COMPARATOR, 'covariates_only']:
        avail = [v for (c, t, a), v in rows.items() if a == arm and v is not None]
        axes[0].plot(grid, np.mean([v['roc_grid'] for v in avail], 0), color=COLORS[arm], lw=2, label=LABELS[arm])
        axes[1].step(grid, np.mean([v['pr_grid'] for v in avail], 0), where='post', color=COLORS[arm], lw=2)
    axes[0].plot([0, 1], [0, 1], color=MUTED, lw=1, ls='--')
    prevalence = np.mean([v['prevalence'] for (c, t, a), v in rows.items() if a == REFERENCE and v is not None])
    axes[1].axhline(prevalence, color=MUTED, lw=1, ls='--')
    axes[0].set(xlabel='False-positive rate', ylabel='True-positive rate', title='ROC, equal-task average')
    axes[1].set(xlabel='Recall', ylabel='Precision (stepwise)', title='Precision–recall, equal-task average', ylim=(.4, 1))
    fig.legend(loc='lower center', ncol=2, frameon=False, bbox_to_anchor=(.5, -.12))
    fig.suptitle('Full depth, nested held cohort, 20 available tasks (HepG2 HMGA1 unavailable), three libraries pooled per task',
                 fontsize=7, color=MUTED)
    save(fig, out, 'roc_pr_average')
    for kind in ['roc', 'pr']:
        fig, axes = plt.subplots(3, 7, figsize=(14, 6.4), sharex=True, sharey=True)
        for ax, t in zip(axes.flat, ledger.itertuples()):
            ax.set_title(f'{t.cell} {t.tf}', fontsize=7)
            for arm in [PRIMARY, REFERENCE, COMPARATOR, 'covariates_only']:
                v = rows[(t.cell, t.tf, arm)]
                if v is None:
                    ax.text(.5, .5, 'unavailable:\nno balanced pairs', ha='center', va='center', color=MUTED, fontsize=7)
                    break
                if kind == 'roc':
                    ax.plot(v['fpr'], v['tpr'], color=COLORS[arm], lw=1.2)
                else:
                    ax.step(v['rec'], v['prec'], color=COLORS[arm], lw=1.2, where='post')
            if kind == 'roc':
                ax.plot([0, 1], [0, 1], color=MUTED, lw=.6, ls='--')
        handles = [plt.Line2D([], [], color=COLORS[a], lw=2) for a in COLORS]
        fig.legend(handles, [LABELS[a] for a in COLORS], loc='lower center', ncol=4, frameon=False)
        fig.suptitle(f'All 21 tasks, full depth, {"ROC" if kind == "roc" else "precision–recall"} (pooled libraries)', color=INK)
        fig.tight_layout(rect=[0, .05, 1, .96])
        save(fig, out, f'{kind}_all_tasks')


def control_plot(out, means):
    arms = [REFERENCE, PRIMARY, COMPARATOR]
    tasks = means[['cell', 'tf']].drop_duplicates().reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(6.4, 6))
    for j, arm in enumerate(arms):
        m = means[means.arm.eq(arm)].set_index(['cell', 'tf']).loc[list(map(tuple, tasks.values))]
        y = np.arange(len(tasks)) + (j - 1) * .25
        ax.errorbar(m.ratio, y, xerr=[m.ratio - m.ratio_low, m.ratio_high - m.ratio], fmt='o', ms=4, color=COLORS[arm],
                    lw=1.2, label=f'{LABELS[arm]}: {int(m.discrepant.sum())}/21 discrepant')
    ax.axvline(1, color=MUTED, lw=1, ls='--')
    ax.set_yticks(range(len(tasks)), [f'{c} {t}' for c, t in tasks.values]); ax.invert_yaxis()
    ax.set_xlabel('Naked-DNA central cuts, observed / expected (95% chromosome-block interval)')
    ax.set_title('Null-mean calibration on examined naked-DNA replicate 1 (held chromosomes)', fontsize=8)
    ax.legend(loc='lower center', bbox_to_anchor=(.5, -.22), frameon=False, fontsize=7)
    save(fig, out, 'control_null_mean')


def shape_plot(out, shape):
    full = shape[shape.depth.eq('full') & shape.arm.eq(PRIMARY)].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    y = np.arange(len(full))
    ax.errorbar(full.positive_over_negative, y, xerr=[full.positive_over_negative - full.pon_low, full.pon_high - full.positive_over_negative],
                fmt='o', ms=4, color=COLORS[PRIMARY], lw=1.2)
    ax.axvline(1, color=MUTED, lw=1, ls='--')
    ax.set_yticks(y, [f'{c} {t}' for c, t in zip(full.cell, full.tf)]); ax.invert_yaxis()
    ax.set_xlabel('Central O/E, ChIP-positive ÷ ChIP-negative (below 1 = relative protection)')
    ax.set_title('Measured regional protection, full depth, three libraries pooled (HepG2 HMGA1 unavailable)', fontsize=8)
    save(fig, out, 'regional_protection')


def power_plot(out, power):
    p = power[power.status.eq('AVAILABLE')]
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    for arm in [PRIMARY, REFERENCE, COMPARATOR]:
        g = p[p.arm.eq(arm)].groupby('r').call_rate.mean()
        ax.plot(1 - g.index, g.values, marker='o', ms=4, lw=2, color=COLORS[arm], label=LABELS[arm])
    ax.set(xlabel='Synthetic fraction of central-endpoint fragments removed', ylabel='Call rate (mean over windows)',
           title='Whole-fragment challenge, 138 examined windows')
    ax.legend(frameon=False, fontsize=6)
    save(fig, out, 'synthetic_power')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--covariates', type=Path, required=True)
    parser.add_argument('--ledger', type=Path, required=True)
    args = parser.parse_args()
    out = args.outdir.resolve()
    (out / 'figures').mkdir(exist_ok=True)
    ledger = read(args.ledger)[['cell', 'tf']]
    if len(ledger) != 21:
        raise ValueError('task denominator changed')
    scores = read(out / 'tables/chromatin_site_scores.tsv.gz')
    grid, rows = curve_data(scores, read(args.covariates), ledger)
    roc_pr(out, grid, rows, ledger)
    control_plot(out, read(out / 'tables/control_task_means.tsv'))
    shape_plot(out, read(out / 'tables/regional_shapes.tsv'))
    power_plot(out, read(out / 'tables/synthetic_power.tsv'))


if __name__ == '__main__':
    main()
