#!/usr/bin/env python3
"""Independently verify probability convolution, geometry choices and held ranks."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import binom
from sklearn.metrics import average_precision_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crossfit_conditional_protection as study
from verify_stable_geometry_crossfit import auc, contrasts


def verify(out):
    study.execute(out, verify=True)
    config = study.preflight(out)
    read = lambda p: pd.read_csv(p, sep='\t', float_precision='round_trip')
    sites = read(config['source_outputs']['sites']['path'])
    with np.load(out / 'score_grid.npz', allow_pickle=False) as archive:
        saved = {key: archive[key] for key in archive.files}
    independent = {arm: np.zeros_like(saved[arm]) for arm in study.ARMS}
    x = np.arange(-100, 101)
    changed, maximum = 0, 0.
    for item in config['items']:
        ix = np.flatnonzero(sites.cell.eq(item['cell']))
        original = read(item['sites']['path'])
        order = pd.Index(study.safety.site_hashes(original)).get_indexer(study.safety.site_hashes(sites.iloc[ix]))
        with np.load(item['profiles']['path'], allow_pickle=False) as z:
            counts = [z['plus_observed'][order], z['minus_observed'][order]]
        with np.load(item['bias_arrays']['path'], allow_pickle=False) as z:
            logits = [z['plus_log_bias'][order], z['minus_log_bias'][order]]
        for arm in study.ARMS:
            for j, (center, left, right) in enumerate(study.prior.GEOMETRIES):
                inner = abs(x) <= center
                region = inner | ((abs(x) >= left) & (abs(x) <= right))
                n = [c[:, region].sum(1).astype(int) for c in counts]
                k = sum(c[:, inner].sum(1).astype(int) for c in counts)
                ps = []
                for log in logits:
                    if arm == 'conditional_uniform':
                        ps.append(np.repeat(inner.sum() / region.sum(), len(ix)))
                    else:
                        propensity = np.exp(log - log[:, region].max(1)[:, None])
                        ps.append(propensity[:, inner].sum(1) / propensity[:, region].sum(1))
                for i, row in enumerate(ix):
                    # Direct convolution, rather than log-pmf plus log-CDF.
                    mass = np.convolve(binom.pmf(np.arange(n[0][i]+1), n[0][i], ps[0][i]),
                                       binom.pmf(np.arange(n[1][i]+1), n[1][i], ps[1][i]))
                    probability = min(1., max(1e-300, mass[:k[i]+1].sum()))
                    independent[arm][row, j] = np.round(-np.log10(probability), 10)
            difference = abs(independent[arm][ix] - saved[arm][ix])
            changed += int((difference > 0).sum())
            maximum = max(maximum, float(difference.max()))
            np.testing.assert_allclose(independent[arm][ix], saved[arm][ix], rtol=0, atol=1.01e-10)
        print('Independent probability grid:', item['cell'], flush=True)
    fits, folds, scores, summary = [read(out / (name + '.tsv')) for name in ['fits', 'folds', 'scores', 'summary']]
    assert len(fits) == 630 and len(folds) == 1260 and len(summary) == 84
    lookup = folds.set_index(['cell', 'tf', 'arm', 'chromosome', 'detector'])
    score_groups = {k: v for k, v in scores.groupby(['cell', 'tf', 'arm', 'TFBS_chr', 'detector'])}
    for r in fits.itertuples(index=False):
        ix = np.flatnonzero(sites.cell.eq(r.cell) & sites.tf.eq(r.tf))
        table = sites.iloc[ix]
        labels = table.chip_label.to_numpy(int)
        held = table.TFBS_chr.eq(r.chromosome).to_numpy()
        excluded = set(table.loc[held, 'pair_id'])
        fit = np.array([v not in excluded for v in table.pair_id], bool)
        assert r.fit_positive == (labels[fit] == 1).sum() and r.fit_negative == (labels[fit] == 0).sum()
        assert r.excluded_pair_mates == (~fit & ~held).sum()
        shape = saved[r.arm + '_shape'][ix]
        expected_index, reason = 0, 'underpowered_fallback'
        if min(r.fit_positive, r.fit_negative) >= 200:
            delta = shape[fit][labels[fit] == 1].mean(0) - shape[fit][labels[fit] == 0].mean(0)
            options = {j: auc(labels[fit], independent[r.arm][ix][fit, j]) for j, geometry in enumerate(study.prior.GEOMETRIES)
                       if contrasts(delta, geometry)[1]}
            expected_index, reason = (max(options, key=options.get), 'selected') if options else (0, 'no_protection_fallback')
        assert tuple(json.loads(r.geometry)) == study.prior.GEOMETRIES[expected_index] and r.reason == reason
        for detector in study.PROTOCOL['detectors']:
            j = expected_index if detector == 'conditional_selected' else 0
            row = lookup.loc[(r.cell, r.tf, r.arm, r.chromosome, detector)]
            yy = labels[held]
            expected = independent[r.arm][ix][held, j]
            assert row.positive == (yy == 1).sum() and row.negative == (yy == 0).sum()
            if held.any():
                group = score_groups[(r.cell, r.tf, r.arm, r.chromosome, detector)]
                cols = study.prior.source.cohorts.prior.previous.KEYS + ['chip_label', 'pair_id']
                pd.testing.assert_frame_equal(group[cols].reset_index(drop=True), table.loc[held, cols].reset_index(drop=True))
                np.testing.assert_allclose(group.score, expected, rtol=0, atol=1.01e-10)
            if row.positive and row.negative:
                np.testing.assert_allclose([row.auroc, row.auprc], [auc(yy, expected), average_precision_score(yy, expected)], rtol=0, atol=1e-12)
                delta = shape[held][yy == 1].mean(0) - shape[held][yy == 0].mean(0)
                values, canonical = contrasts(delta, study.prior.GEOMETRIES[j])
                np.testing.assert_allclose([row.center, row.left, row.right], values, rtol=0, atol=1e-12)
                assert row.canonical == canonical and row.fixed_canonical == contrasts(delta, (5, 15, 40))[1]
            else:
                assert row.evaluation_status == 'missing_class' and pd.isna(row.auroc)
    for r in summary.itertuples(index=False):
        f = folds[folds.cell.eq(r.cell) & folds.tf.eq(r.tf) & folds.arm.eq(r.arm) & folds.detector.eq(r.detector)]
        assert len(f) == 15 and f.auroc.notna().sum() == r.evaluable_chromosomes
        np.testing.assert_allclose([r.mean_auroc, r.mean_auprc, r.canonical_fraction, r.fixed_canonical_fraction],
                                  [f.auroc.mean(), f.auprc.mean(), f.canonical.mean(), f.fixed_canonical.mean()], equal_nan=True)
    print(f'Verified 178332 grid scores, 630 fits, 1260 folds and 84 summaries; numerical score differences {changed}, maximum {maximum}', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir', type=Path, required=True)
    verify(p.parse_args().outdir.resolve())
