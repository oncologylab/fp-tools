"""All-task, all-run and promotion-gate ledgers from completed tables (read-only)."""
import argparse
from pathlib import Path

import pandas as pd

from protocol import ARMS, COMPARATOR, DEPTHS, LIBRARIES, PRIMARY

GATE_OUT = 'tables/promotion_gate_ledger.tsv'


def read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip', low_memory=False)


def task_ledger(out, ledger):
    metrics = read(out / 'tables/arm_task_metrics.tsv')
    deltas = read(out / 'tables/task_deltas.tsv')
    deltas = deltas[deltas.reference.eq(COMPARATOR)]
    means = read(out / 'tables/control_task_means.tsv')
    calls = read(out / 'tables/control_calls.tsv')
    shape = read(out / 'tables/regional_shapes.tsv')
    rows = []
    for t in ledger.itertuples():
        for depth in DEPTHS:
            for arm in list(ARMS) + ['covariates_only']:
                m = metrics[metrics.cell.eq(t.cell) & metrics.tf.eq(t.tf) & metrics.depth.eq(depth) & metrics.arm.eq(arm)]
                d = deltas[deltas.cell.eq(t.cell) & deltas.tf.eq(t.tf) & deltas.depth.eq(depth) & deltas.method.eq(arm)]
                c = means[means.cell.eq(t.cell) & means.tf.eq(t.tf) & means.arm.eq(arm)]
                k = calls[calls.cell.eq(t.cell) & calls.tf.eq(t.tf) & calls.arm.eq(arm)]
                s = shape[shape.cell.eq(t.cell) & shape.tf.eq(t.tf) & shape.depth.eq(depth) & shape.arm.eq(arm)]
                first = lambda f, col: f[col].iloc[0] if len(f) else None  # noqa: E731
                rows.append(dict(cell=t.cell, tf=t.tf, depth=depth, arm=arm,
                                 chromatin_status=first(m, 'status') or 'UNAVAILABLE',
                                 auroc=first(m, 'auroc'), auprc=first(m, 'auprc'),
                                 count_stratified_auroc=first(m, 'count_stratified_auroc'),
                                 auroc_gain_vs_DWM=first(d, 'auroc_gain'), relative_ap_gain_vs_DWM=first(d, 'relative_auprc_gain'),
                                 control_ratio=first(c, 'ratio'), control_discrepant=first(c, 'discrepant'),
                                 control_calls=first(k, 'calls'), control_informative=first(k, 'informative'),
                                 control_informative_wilson_upper=first(k, 'informative_wilson_upper'),
                                 positive_over_negative_central_oe=first(s, 'positive_over_negative'),
                                 canonical_protection=first(s, 'canonical_protection'),
                                 note='no balanced held pairs; never counted as zero or pass'
                                 if (t.cell, t.tf) == ('HepG2', 'HMGA1') else ''))
    frame = pd.DataFrame(rows)
    if frame[['cell', 'tf']].drop_duplicates().shape[0] != 21:
        raise ValueError('task denominator changed')
    frame.to_csv(out / 'tables/all_task_ledger.tsv', sep='\t', index=False)


def run_ledger(out):
    folds = read(out / 'tables/fold_metrics.tsv')
    rows = []
    for library in LIBRARIES:
        for depth in DEPTHS:
            for arm in list(ARMS) + ['covariates_only']:
                f = folds[folds.library.eq(library) & folds.depth.eq(depth) & folds.arm.eq(arm)]
                rows.append(dict(library=library, depth=depth, seed=2026, arm=arm, attempt=1, exit_code=0,
                                 fold_records=len(f), available_folds=int(f.status.eq('AVAILABLE').sum()),
                                 incomplete_folds=int(f.status.eq('INCOMPLETE').sum()), no_site_folds=int(f.status.eq('NO_SITES').sum()),
                                 sites=int(f.sites.sum())))
    for combo, reason in [('seeds other than 2026', 'NOT_RUN: candidate failed advancement screen'),
                          ('naked-DNA replicate 2 calibration transfer', 'NOT_RUN: examined data; screen failed'),
                          ('naked-DNA replicate 3', 'NOT_TESTED_LOCKED'),
                          ('SK-N-SH / GM23338 labels', 'NOT_TESTED_LOCKED'),
                          ('complete-workflow resource protocol', 'NOT_RUN: candidate failed advancement screen')]:
        rows.append(dict(library=combo, depth='', seed='', arm='', attempt='', exit_code='', fold_records=0,
                         available_folds=0, incomplete_folds=0, no_site_folds=0, sites=0, status=reason))
    pd.DataFrame(rows).to_csv(out / 'tables/run_ledger.tsv', sep='\t', index=False)


def gate_ledger(out):
    screen = read(out / 'tables/advancement_screen.tsv').set_index('gate')
    summary = read(out / 'tables/screen_summary.tsv')
    main = summary[summary.depth.eq('full') & summary.method.eq(PRIMARY) & summary.reference.eq(COMPARATOR)].iloc[0]
    ev = lambda g: screen.loc[g, 'evidence']  # noqa: E731
    arm = read(out / 'tables/arm_task_metrics.tsv')
    strat = arm[arm.depth.eq('full') & arm.status.eq('AVAILABLE')].groupby('arm').count_stratified_auroc.mean()
    rows = [
        ('Mean AUROC gain >= 0.03 over DWM', 'FAIL',
         f"+{main.mean_auroc_gain:.5f} versus raw-DWM-propensity factorial control is numerically above 0.03, but the "
         "qualified native DWM residual comparator is unavailable at full depth and the gain is coverage-confounded "
         f"(count-stratified AUROC {strat[PRIMARY]:.4f} versus covariate-only {strat['covariates_only']:.4f}); "
         "development screen failed elsewhere"),
        ('Relative AUPRC gain >= 10%', 'FAIL', f'{main.mean_relative_auprc_gain:.4%} versus DWM factorial control'),
        ('>=25% stronger standardized separation with canonical protection in two difficult families', 'FAIL', ev('F_shape')),
        ('Improvement in two families and two cell contexts', 'DEVELOPMENT_CRITERION_MET',
         ev('E_contexts') + '; ranking only, not protection'),
        ('No CTCF/positive-control AUROC loss > 0.02', 'DEVELOPMENT_CRITERION_MET', ev('D_positive_control')),
        ('Principal chromosome-block interval excludes zero', 'DEVELOPMENT_CRITERION_MET',
         f'[{main.bootstrap_low:.5f}, {main.bootstrap_high:.5f}] fixed predictions, factorial DWM control only'),
        ('Untouched naked-DNA FPR and Wilson upper <= 5%, <= 1 point over DWM', 'NOT_TESTED_LOCKED',
         'replicate 3 closed. Examined replicate 1: ' + ev('G_controls')),
        ('No evidence bias component learned occupancy', 'DEVELOPMENT_CRITERION_MET',
         'recalibration fitted only on protein-free naked DNA, no labels or chromatin; LOG21 unchanged'),
        ('Stable direction across depths, seeds, biological replicates', 'FAIL', ev('H_direction') + '; seeds not run'),
        ('Complete-workflow CPU <= 2x DWM', 'INCOMPLETE', 'not run: candidate failed the advancement screen'),
        ('Absolute process-tree peak memory <= 1.5x DWM', 'INCOMPLETE', 'not run: candidate failed the advancement screen'),
        ('Deployed model and sidecars <= 25 MB', 'DEVELOPMENT_CRITERION_MET',
         'LOG21 NPZ 7,519 bytes plus two recalibration scalars in JSON; whole deployment not packaged'),
        ('Calibrated functional uncertainty', 'INCOMPLETE', 'not addressed; regional ratios only'),
        ('GP/hybrid must outperform simpler model', 'NOT_APPLICABLE', 'no GP or hybrid model proposed'),
        ('Two new cell contexts (SK-N-SH, GM23338)', 'NOT_TESTED_LOCKED', 'labels closed'),
        ('Calibration-round screen A: null mean', screen.loc['A_calibration', 'status'], ev('A_calibration')),
        ('Overall decision', 'FAIL', screen.loc['DECISION', 'evidence'] + '; general replacement rejected'),
    ]
    pd.DataFrame(rows, columns=['requirement', 'status', 'evidence']).to_csv(out / GATE_OUT, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--ledger', type=Path, required=True)
    a = parser.parse_args()
    out = a.outdir.resolve()
    task_ledger(out, read(a.ledger))
    run_ledger(out)
    gate_ledger(out)


if __name__ == '__main__':
    main()
