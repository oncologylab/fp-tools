#!/usr/bin/env python3
"""Independently recompute saved normalization components with least squares."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import diagnose_low_count_normalization as study


def verify(out):
    config = study.preflight(out)
    base = study.extraction.base
    result = base.verify_stage(out / 'diagnostic_manifest.json', config['configuration_id'], study.ROOT)
    profiles = base.verify_stage(out / 'profiles_manifest.json', config['configuration_id'], study.ROOT)
    rows = pd.read_csv(result['outputs']['rows']['path'], sep='\t', float_precision='round_trip')
    fitted = json.loads(Path(config['models']['path']).read_text())['models']
    fitted = {(m['cell'], m['tf'], m['arm']): m for m in fitted}
    thresholds = pd.read_csv(config['thresholds']['path'], sep='\t', float_precision='round_trip').set_index(['cell', 'tf', 'arm'])
    training = pd.read_csv(config['training']['sites']['path'], sep='\t')
    with np.load(config['training']['panels']['path'], allow_pickle=False) as z:
        training_panels = {arm: z[arm] for arm in study.ARMS}
    cohorts = [('ATAC_training', training, training_panels)]
    for item in config['items']:
        sites = pd.read_csv(item['sites']['path'], sep='\t')
        mask = sites.TFBS_chr.isin(study.PROTOCOL['chromosomes']).to_numpy()
        with np.load(profiles['outputs'][item['cell']]['path'], allow_pickle=False) as z:
            observed, expected = z['observed'], z['expected']
            np.testing.assert_array_equal(observed, z['plus_observed'] + z['minus_observed'])
            np.testing.assert_allclose(observed.sum(1), expected.sum(1), rtol=1e-12, atol=1e-12)
            reconstructed = np.zeros_like(expected)
            with np.load(item['arrays']['path'], allow_pickle=False) as cached:
                for strand in ['plus', 'minus']:
                    logits = cached[strand + '_log_bias'][mask]
                    finite = np.isfinite(logits)
                    valid = finite.any(1)
                    logits = np.where(finite[valid], logits[valid], -np.inf)
                    probabilities = np.exp(logits - logsumexp(logits, axis=1)[:, None])
                    reconstructed[valid] += probabilities * z[strand + '_observed'][valid].sum(1)[:, None]
            np.testing.assert_allclose(reconstructed, expected, rtol=1e-12, atol=1e-12)
        _, _, residual = study.safety.direct_scores(observed, expected, config['dispersion'][item['cell']])
        cohorts.append(('naked_rep1', sites.loc[mask].reset_index(drop=True), dict(raw=observed, parametric_direct=residual)))
    checked = 0
    x = np.arange(-100, 101)
    outer = abs(x) >= 60
    design = np.column_stack([np.ones(outer.sum()), x[outer]])
    for source, sites, panels in cohorts:
        for (cell, tf), task in sites.groupby(['cell', 'tf']):
            ix = task.index.to_numpy()
            for arm in study.ARMS:
                saved = rows.loc[rows.source.eq(source) & rows.cell.eq(cell) & rows.tf.eq(tf) & rows.arm.eq(arm)].reset_index(drop=True)
                pd.testing.assert_frame_equal(saved[study.safety.SITE_KEYS], task[study.safety.SITE_KEYS].reset_index(drop=True))
                values = panels[arm][ix]
                coefficients = np.linalg.lstsq(design, values[:, outer].T, rcond=None)[0]
                residual = values - coefficients[0, :, None] - coefficients[1, :, None] * x
                scale = np.sqrt((residual[:, outer] ** 2).mean(1))
                denominator = np.where(scale > 1e-8, scale, 1.)
                normalized = np.clip(residual / denominator[:, None], -12, 12)
                np.testing.assert_allclose(scale, saved.scale, rtol=1e-10, atol=1e-12)
                np.testing.assert_array_equal(scale <= 1e-8, saved.fallback)
                np.testing.assert_allclose(1. / denominator, saved.amplification, rtol=1e-10, atol=1e-10)
                np.testing.assert_allclose((abs(residual / denominator[:, None]) > 12).mean(1), saved.clipped_fraction, rtol=0, atol=1e-15)
                model = fitted[(cell, tf, arm)]
                if model['status'] != 'unavailable':
                    center, left, right = model['geometry']
                    weights = ((abs(x) >= left) & (abs(x) <= right)).astype(float)
                    weights /= weights.sum()
                    inner = (abs(x) <= center).astype(float)
                    weights -= inner / inner.sum()
                    score = np.round(normalized @ weights, 10)
                    score[panels['raw'][ix].sum(1) == 0] = 0
                    np.testing.assert_allclose(score, saved.score, rtol=0, atol=1.01e-10)
                    cutoff = thresholds.loc[(cell, tf, arm)].candidate_threshold
                    np.testing.assert_array_equal(saved.threshold_evaluated & (score > cutoff), saved.exceeds_frozen_threshold)
                else:
                    assert saved.score.isna().all() and not saved.threshold_evaluated.any()
                checked += len(saved)
    if checked != len(rows):
        raise ValueError('rows not independently checked')
    print(f'Independent least-squares verification: {checked} normalization rows, scores and threshold decisions', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    verify(parser.parse_args().outdir.resolve())
