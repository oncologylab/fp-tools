#!/usr/bin/env python3
"""Recheck regional calibration, resampling and predictive intervals independently."""
import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_rep1_region_calibration as study


def bernoulli_interval(counts, probabilities):
    mass = np.array([1.])
    for n, p in zip(counts.astype(int), probabilities):
        for _ in range(n):
            mass = np.r_[mass * (1-p), 0.] + np.r_[0., mass * p]
    cdf = (mass / mass.sum()).cumsum()
    return np.searchsorted(cdf, [.025, .975])


def verify(out):
    study.execute(out, verify=True)
    config = study.preflight(out)
    read = lambda p: pd.read_csv(p, sep='\t', float_precision='round_trip')
    rows, calibration, likelihood, blocks = [read(out / (name + '.tsv')) for name in ['sites', 'calibration', 'likelihood', 'blocks']]
    weights = np.load(out / 'bootstrap_weights.npy', allow_pickle=False)
    np.testing.assert_array_equal(weights, np.random.default_rng(2026).multinomial(15, np.repeat(1/15, 15), size=2000))
    x = np.arange(-100, 101)
    masks = [abs(x) <= 7, (x >= -40) & (x <= -12), (x >= 12) & (x <= 40)]
    union = np.logical_or.reduce(masks)
    for item in config['items']:
        sites = read(item['sites']['path'])
        with np.load(config['profiles'][item['cell']]['path'], allow_pickle=False) as z:
            observed = {strand: z[strand + '_observed'] for strand in ['plus', 'minus']}
        with np.load(item['arrays']['path'], allow_pickle=False) as z:
            logs = {strand: z[strand + '_log_bias'] for strand in ['plus', 'minus']}
        for strand in ['plus', 'minus']:
            proportions = np.exp(logs[strand] - logs[strand][:, union].max(1)[:, None])
            counts = np.column_stack([observed[strand][:, mask].sum(1) for mask in masks])
            for arm in ['parametric', 'uniform']:
                expected_prob = np.column_stack([proportions[:, mask].sum(1) / proportions[:, union].sum(1) for mask in masks]) if arm == 'parametric' else np.tile(np.array([15, 29, 29])/73, (len(sites), 1))
                for j, region in enumerate(['center', 'left', 'right']):
                    frame = rows[rows.cell.eq(item['cell']) & rows.strand.eq(strand) & rows.arm.eq(arm) & rows.region.eq(region)].reset_index(drop=True)
                    pd.testing.assert_frame_equal(frame[study.source.safety.SITE_KEYS], sites[study.source.safety.SITE_KEYS])
                    np.testing.assert_array_equal(frame.total, counts.sum(1))
                    np.testing.assert_array_equal(frame.observed, counts[:, j])
                    np.testing.assert_allclose(frame.probability, expected_prob[:, j], rtol=1e-12, atol=1e-12)
                    np.testing.assert_allclose(frame.expected, counts.sum(1)*expected_prob[:, j], rtol=1e-12, atol=1e-12)
    for r in calibration.itertuples(index=False):
        frame = rows[rows.cell.eq(r.cell) & rows.tf.eq(r.tf) & rows.arm.eq(r.arm) & rows.region.eq(r.region) & rows.eligible]
        if r.strand != 'pooled':
            frame = frame[frame.strand.eq(r.strand)]
        observed, expected, total = frame.observed.sum(), frame.expected.sum(), frame.total.sum()
        assert len(frame) == r.site_strand_rows and (frame.total > 0).sum() == r.informative_site_strands
        np.testing.assert_allclose([r.total_cuts, r.observed, r.expected], [total, observed, expected], rtol=0, atol=1e-10)
        np.testing.assert_array_equal([r.independent_null_count_low, r.independent_null_count_high], bernoulli_interval(frame.total.to_numpy(), frame.probability.to_numpy()))
        b = frame.groupby('TFBS_chr')[['total', 'observed', 'expected']].sum().reindex(study.CHROMS, fill_value=0)
        saved = blocks[blocks.cell.eq(r.cell) & blocks.tf.eq(r.tf) & blocks.strand.eq(r.strand) & blocks.arm.eq(r.arm) & blocks.region.eq(r.region)].set_index('chromosome').loc[study.CHROMS]
        np.testing.assert_allclose(saved[['total', 'observed', 'expected']], b, rtol=0, atol=1e-10)
        denominator = weights @ b.total.to_numpy()
        valid = denominator > 0
        residual = weights @ (b.observed - b.expected).to_numpy()
        interval = np.percentile(residual[valid]/denominator[valid], [2.5, 97.5]) if valid.any() else [np.nan, np.nan]
        assert valid.sum() == r.valid_bootstrap_draws
        np.testing.assert_allclose([r.block_error_low, r.block_error_high], interval, rtol=0, atol=1e-12, equal_nan=True)
    for r in likelihood.itertuples(index=False):
        frame = rows[rows.cell.eq(r.cell) & rows.tf.eq(r.tf) & rows.arm.eq('parametric') & rows.eligible].copy()
        if r.strand != 'pooled':
            frame = frame[frame.strand.eq(r.strand)]
        frame['gain'] = frame.observed * (np.log(frame.probability) - np.log(frame.region.map({'center': 15/73, 'left': 29/73, 'right': 29/73})))
        gain = frame.groupby('TFBS_chr').gain.sum().reindex(study.CHROMS, fill_value=0).to_numpy()
        counts = frame[frame.region.eq('center')].groupby('TFBS_chr').total.sum().reindex(study.CHROMS, fill_value=0).to_numpy()
        np.testing.assert_allclose([r.total_cuts, r.log_likelihood_gain, r.gain_per_cut], [counts.sum(), gain.sum(), gain.sum()/counts.sum()], rtol=0, atol=1e-10)
        denominator = weights @ counts
        valid = denominator > 0
        interval = np.percentile((weights @ gain)[valid]/denominator[valid], [2.5, 97.5])
        assert valid.sum() == r.valid_bootstrap_draws
        np.testing.assert_allclose([r.block_gain_low, r.block_gain_high], interval, rtol=0, atol=1e-12)
    assert len(rows) == 50400 and len(calibration) == 378 and len(likelihood) == 63 and len(blocks) == 5670
    print('Independent verification: 50400 site-region rows, 378 calibration/predictive intervals, 63 likelihood comparisons and all 5670 chromosome blocks', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir', type=Path, required=True)
    verify(p.parse_args().outdir.resolve())
