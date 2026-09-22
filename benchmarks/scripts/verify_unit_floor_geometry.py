#!/usr/bin/env python3
"""Check unit-floor fits and held scores with independent least-squares arithmetic."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crossfit_unit_floor_geometry as study
from verify_stable_geometry_crossfit import auc, weights, contrasts


def verify(out):
    study.execute(out, verify=True)
    config = json.loads((out / 'configuration.freeze.json').read_text())
    read = lambda p: pd.read_csv(p, sep='\t', float_precision='round_trip')
    sites = read(config['source_outputs']['sites']['path'])
    fits, folds, scores, summary, curves, paired = [read(out / (name + '.tsv')) for name in
                                                 ['fits', 'folds', 'scores', 'summary', 'curves', 'paired']]
    assert len(fits) == 2205 and len(folds) == 4410 and len(summary) == 294 and len(paired) == 4410
    with np.load(config['source_outputs']['panels']['path'], allow_pickle=False) as archive:
        panels = {arm: archive[arm] for arm in study.PROTOCOL['arms']}
    x = np.arange(-100, 101)
    outer = abs(x) >= 60
    design = np.column_stack([np.ones(201), x])
    normal = {}
    for arm, values in panels.items():
        coef = np.linalg.lstsq(design[outer], values[:, outer].T, rcond=None)[0]
        residual = values - (design @ coef).T
        scale = np.sqrt((residual[:, outer] ** 2).mean(1))
        normal[arm] = np.clip(residual / np.maximum(scale, 1)[:, None], -12, 12)
    keys = ['cell', 'tf', 'arm', 'chromosome']
    lookup = folds.set_index(keys + ['detector'])
    score_groups = {k: v for k, v in scores.groupby(['cell', 'tf', 'arm', 'TFBS_chr', 'detector'])}
    checked = 0
    score_differences = []
    for fit_row in fits.itertuples(index=False):
        ix = np.flatnonzero(sites.cell.eq(fit_row.cell) & sites.tf.eq(fit_row.tf))
        table = sites.iloc[ix]
        y = table.chip_label.to_numpy(int)
        held = table.TFBS_chr.eq(fit_row.chromosome).to_numpy()
        excluded = set(table.loc[held, 'pair_id'])
        fit = np.array([v not in excluded for v in table.pair_id], dtype=bool)
        assert fit_row.fit_positive == (y[fit] == 1).sum() and fit_row.fit_negative == (y[fit] == 0).sum()
        assert fit_row.excluded_pair_mates == (~fit & ~held).sum()
        values = normal[fit_row.arm][ix]
        zero = panels['raw'][ix].sum(1) == 0
        geometry = tuple(json.loads(fit_row.geometry))
        if min(fit_row.fit_positive, fit_row.fit_negative) < 200:
            assert fit_row.reason == 'underpowered_fallback' and geometry == (7, 12, 40)
        else:
            delta = values[fit][y[fit] == 1].mean(0) - values[fit][y[fit] == 0].mean(0)
            candidates = {}
            for g in study.prior.GEOMETRIES:
                if contrasts(delta, g)[1]:
                    s = np.round(values[fit] @ weights(g), 10)
                    s[zero[fit]] = 0
                    candidates[g] = auc(y[fit], s)
            if candidates:
                assert fit_row.reason == 'selected' and geometry == max(candidates, key=candidates.get)
            else:
                assert fit_row.reason == 'no_protection_fallback' and geometry == (7, 12, 40)
        for detector in study.PROTOCOL['detectors']:
            g = geometry if detector == 'unit_floor_selected' else (7, 12, 40)
            row = lookup.loc[(fit_row.cell, fit_row.tf, fit_row.arm, fit_row.chromosome, detector)]
            yy = y[held]
            assert row.positive == (yy == 1).sum() and row.negative == (yy == 0).sum()
            expected = np.round(values[held] @ weights(g), 10)
            expected[zero[held]] = 0
            if held.any():
                saved = score_groups[(fit_row.cell, fit_row.tf, fit_row.arm, fit_row.chromosome, detector)]
                cols = study.prior.source.cohorts.prior.previous.KEYS + ['chip_label', 'pair_id']
                pd.testing.assert_frame_equal(saved[cols].reset_index(drop=True), table.loc[held, cols].reset_index(drop=True))
                discrepancy = np.abs(expected - saved.score.to_numpy())
                np.testing.assert_allclose(expected, saved.score, rtol=0, atol=1.01e-10)
                if np.any(discrepancy):
                    score_differences.append(dict(cell=fit_row.cell, tf=fit_row.tf, arm=fit_row.arm,
                        chromosome=fit_row.chromosome, detector=detector, changed_scores=int((discrepancy > 0).sum()),
                        maximum_difference=float(discrepancy.max())))
            if row.positive and row.negative:
                np.testing.assert_allclose([row.auroc, row.auprc], [auc(yy, expected), average_precision_score(yy, expected)], atol=1e-12, rtol=0)
                delta = values[held][yy == 1].mean(0) - values[held][yy == 0].mean(0)
                for prefix, geom in [('', g), ('fixed_', (5, 15, 40))]:
                    v, canonical = contrasts(delta, geom)
                    np.testing.assert_allclose([getattr(row, prefix + k) for k in ['center', 'left', 'right']], v, atol=1e-11)
                    assert getattr(row, prefix + 'canonical') == canonical
                checked += 1
            else:
                assert row.evaluation_status == 'missing_class' and pd.isna(row.auroc)
    for row in summary.itertuples(index=False):
        chosen = folds[folds.cell.eq(row.cell) & folds.tf.eq(row.tf) & folds.arm.eq(row.arm) & folds.detector.eq(row.detector)]
        assert len(chosen) == 15 and chosen.auroc.notna().sum() == row.evaluable_chromosomes
        np.testing.assert_allclose([row.mean_auroc, row.mean_auprc, row.canonical_fraction, row.fixed_canonical_fraction],
                                  [chosen.auroc.mean(), chosen.auprc.mean(), chosen.canonical.mean(), chosen.fixed_canonical.mean()], equal_nan=True)
    for (cell, tf, arm), frame in curves.groupby(['cell', 'tf', 'arm']):
        ix = np.flatnonzero(sites.cell.eq(cell) & sites.tf.eq(tf))
        y = sites.iloc[ix].chip_label.to_numpy()
        delta = normal[arm][ix][y == 1].mean(0) - normal[arm][ix][y == 0].mean(0)
        np.testing.assert_allclose(frame.contrast, delta, atol=1e-11)
    baseline = read(config['baseline_outputs']['folds']['path']).set_index(keys + ['detector'])
    for row in paired.itertuples(index=False):
        k = (row.cell, row.tf, row.arm, row.chromosome)
        new, old = lookup.loc[k + (row.detector_new,)], baseline.loc[k + (row.detector_old,)]
        np.testing.assert_allclose([row.auroc_gain, row.auprc_gain], [new.auroc - old.auroc, new.auprc - old.auprc], equal_nan=True, rtol=0, atol=1e-12)
    print(f'Independent verification passed: 2205 fits, 4410 folds ({checked} evaluable), scores within 1.01e-10, 294 summaries, 140 curves and 4410 paired effects within 1e-12', flush=True)
    print('Exact rounded-score parity:', not score_differences,
          '; one-decimal-unit sensitivity:', json.dumps(score_differences), flush=True)
    return score_differences


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    verify(parser.parse_args().outdir.resolve())
