#!/usr/bin/env python3
"""Independently check geometry selection, held scores, shapes and summaries."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crossfit_stable_geometry as experiment


def auc(y, score):
    positive = y == 1
    n = int(positive.sum())
    return (rankdata(score)[positive].sum() - n * (n + 1) / 2) / (n * (len(y) - n))


def weights(g):
    x = np.arange(-100, 101)
    c, start, end = g
    inner = abs(x) <= c
    outer = (abs(x) >= start) & (abs(x) <= end)
    return outer / outer.sum() - inner / inner.sum()


def contrasts(delta, g):
    x = np.arange(-100, 101)
    c, start, end = g
    a, b, d = delta[abs(x) <= c].mean(), delta[(x >= -end) & (x <= -start)].mean(), delta[(x >= start) & (x <= end)].mean()
    return [a, b, d], bool(a < -1e-10 and b - a > 1e-10 and d - a > 1e-10)


def verify(out):
    experiment.execute(out, verify=True)
    config = json.loads((out / 'configuration.freeze.json').read_text())
    paths = config['source_outputs']
    read = lambda path: pd.read_csv(path, sep='\t', float_precision='round_trip')
    sites = read(paths['sites']['path'])
    fits, folds, summary, scores = [read(out / (name + '.tsv')) for name in ['fits', 'folds', 'summary', 'scores']]
    assert len(fits) == 2205 and len(folds) == 6615 and len(summary) == 441
    key = ['cell', 'tf', 'arm', 'chromosome']
    lookup = folds.set_index(key + ['detector'])
    score_groups = {k: v for k, v in scores.groupby(['cell', 'tf', 'arm', 'TFBS_chr', 'detector'])}
    with np.load(paths['panels']['path'], allow_pickle=False) as archive:
        panels = {arm: archive[arm] for arm in experiment.source.ARMS}
    x = np.arange(-100, 101); outer = abs(x) >= 60
    design = np.column_stack([np.ones(201), x]); normal = {}
    for arm, p in panels.items():
        coefficients = np.linalg.lstsq(design[outer], p[:, outer].T, rcond=None)[0]
        residual = p - (design @ coefficients).T
        scale = np.sqrt(np.mean(residual[:, outer] ** 2, axis=1))
        normal[arm] = np.clip(residual / np.where(scale > 1e-8, scale, 1)[:, None], -12, 12)
    checked = 0
    sensitivity = []
    for r in fits.itertuples(index=False):
        ix = np.flatnonzero(sites.cell.eq(r.cell) & sites.tf.eq(r.tf))
        table = sites.iloc[ix]; y = table.chip_label.to_numpy(int)
        held = table.TFBS_chr.eq(r.chromosome).to_numpy()
        excluded = set(table.loc[held, 'pair_id'])
        fit = np.array([p not in excluded for p in table.pair_id], dtype=bool)
        assert r.fit_positive == (y[fit] == 1).sum() and r.fit_negative == (y[fit] == 0).sum()
        assert r.excluded_pair_mates == (~fit & ~held).sum()
        raw, alternative = panels[r.arm][ix], normal[r.arm][ix]
        # Reproduce the specified arithmetic for rank-sensitive selection;
        # also retain the independent least-squares solver as a sensitivity.
        coefficients = raw[:, outer] @ np.linalg.pinv(design[outer]).T
        residual = raw - (coefficients[:, [0]] + coefficients[:, [1]] * x[None, :])
        scale = np.sqrt(np.mean(residual[:, outer] ** 2, axis=1))
        norm = np.clip(residual / np.where(scale > 1e-8, scale, 1)[:, None], -12, 12)
        np.testing.assert_allclose(norm, alternative, atol=1e-11, rtol=1e-10)
        zero = panels['raw'][ix].sum(1) == 0
        geometry = tuple(json.loads(r.geometry))
        if min(r.fit_positive, r.fit_negative) < 200:
            assert r.reason == 'underpowered_fallback' and geometry == (7, 12, 40)
        else:
            delta = norm[fit][y[fit] == 1].mean(0) - norm[fit][y[fit] == 0].mean(0)
            eligible = {}
            for g in experiment.GEOMETRIES:
                if contrasts(delta, g)[1]:
                    s = np.round(norm[fit] @ weights(g), 10); s[zero[fit]] = 0
                    eligible[g] = auc(y[fit], s)
            if eligible:
                assert r.reason == 'selected' and geometry in eligible
                np.testing.assert_allclose(eligible[geometry], max(eligible.values()), atol=1e-12, rtol=0)
            else:
                assert r.reason == 'no_protection_fallback' and geometry == (7, 12, 40)
            alternate_delta = alternative[fit][y[fit] == 1].mean(0) - alternative[fit][y[fit] == 0].mean(0)
            alternate_eligible = {}
            for g in experiment.GEOMETRIES:
                if contrasts(alternate_delta, g)[1]:
                    s = np.round(alternative[fit] @ weights(g), 10); s[zero[fit]] = 0
                    alternate_eligible[g] = auc(y[fit], s)
            assert alternate_eligible == eligible, 'solver changed fitting metrics or eligibility'
            if eligible:
                assert max(alternate_eligible, key=alternate_eligible.get) == geometry, 'solver changed selected geometry'
        for detector in experiment.PROTOCOL['detectors']:
            g = geometry if detector == 'normalized_selected' else (7, 12, 40)
            row = lookup.loc[(r.cell, r.tf, r.arm, r.chromosome, detector)]
            v = raw if detector == 'raw_fixed' else norm
            expected = np.round(v[held] @ weights(g), 10); expected[zero[held]] = 0
            yy = y[held]
            assert row.positive == (yy == 1).sum() and row.negative == (yy == 0).sum()
            if held.any():
                frame = score_groups[(r.cell, r.tf, r.arm, r.chromosome, detector)]
                cols = experiment.source.cohorts.prior.previous.KEYS + ['chip_label', 'pair_id']
                pd.testing.assert_frame_equal(frame[cols].reset_index(drop=True), table.loc[held, cols].reset_index(drop=True))
                np.testing.assert_array_equal(frame.score, expected)
            if row.positive and row.negative:
                assert row.evaluation_status == 'available'
                archived = frame.score.to_numpy()
                np.testing.assert_allclose([row.auroc, row.auprc], [auc(yy, archived), average_precision_score(yy, archived)], atol=1e-12)
                alternate_scores = np.round((raw if detector == 'raw_fixed' else alternative)[held] @ weights(g), 10)
                alternate_scores[zero[held]] = 0
                np.testing.assert_array_equal(alternate_scores, archived)
                sensitivity.append(dict(cell=r.cell, tf=r.tf, arm=r.arm, chromosome=r.chromosome,
                    detector=detector, maximum_score_difference=float(np.max(abs(alternate_scores - archived))),
                    auroc_difference=auc(yy, alternate_scores) - row.auroc,
                    auprc_difference=average_precision_score(yy, alternate_scores) - row.auprc))
                delta = norm[held][yy == 1].mean(0) - norm[held][yy == 0].mean(0)
                values, canonical = contrasts(delta, g)
                np.testing.assert_allclose([row.center, row.left, row.right], values, atol=1e-11)
                assert row.canonical == canonical
                checked += 1
            else:
                assert row.evaluation_status == 'missing_class' and pd.isna(row.auroc)
    for r in summary.itertuples(index=False):
        f = folds.loc[folds.cell.eq(r.cell) & folds.tf.eq(r.tf) & folds.arm.eq(r.arm) & folds.detector.eq(r.detector)]
        good = f.loc[f.evaluation_status.eq('available')]
        assert len(f) == 15 and len(good) == r.evaluable_chromosomes
        vals = [good.auroc.mean(), good.auprc.mean(), good.canonical.mean()]
        np.testing.assert_allclose([r.mean_auroc, r.mean_auprc, r.canonical_fraction], vals, equal_nan=True)
    print(f'Independently verified 2205 fits, 6615 folds ({checked} evaluable), all saved scores and 441 summaries', flush=True)
    diagnostic = pd.DataFrame(sensitivity)
    changed = diagnostic.auroc_difference.abs().gt(1e-12) | diagnostic.auprc_difference.abs().gt(1e-12)
    print('Equivalent-solver sensitivity:', int(changed.sum()), 'changed folds;',
          diagnostic[['maximum_score_difference', 'auroc_difference', 'auprc_difference']].abs().max().to_dict(), flush=True)
    assert not changed.any(), 'rank sensitivity persists'
    return diagnostic


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--sensitivity-output', type=Path)
    args = parser.parse_args()
    diagnostic = verify(args.outdir.resolve())
    if args.sensitivity_output:
        target = args.sensitivity_output.resolve()
        if not target.is_relative_to(experiment.ROOT / 'benchmarks/results'):
            parser.error('use research results')
        with target.open('x') as handle:
            diagnostic.to_csv(handle, sep='\t', index=False)
