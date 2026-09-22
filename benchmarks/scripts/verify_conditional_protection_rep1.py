#!/usr/bin/env python3
"""Verify frozen full-training geometries and replicate-1 count-tail decisions."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import binom, binomtest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import screen_conditional_protection_rep1 as screen
from verify_stable_geometry_crossfit import auc, contrasts


def verify(out):
    screen.execute(out, verify=True)
    config = screen.preflight(out)
    read = lambda p: pd.read_csv(p, sep='\t', float_precision='round_trip')
    models = json.loads((out / 'models.freeze.json').read_text())['models']
    rows, summary = read(out / 'sites.tsv'), read(out / 'summary.tsv')
    sites = read(config['training']['sites']['path'])
    with np.load(config['grid']['path'], allow_pickle=False) as grid:
        for m in models:
            ix = np.flatnonzero(sites.cell.eq(m['cell']) & sites.tf.eq(m['tf']))
            y = sites.iloc[ix].chip_label.to_numpy()
            if not len(ix):
                assert m['status'] == 'unavailable' and m['geometry'] is None
                continue
            best, reason = 0, 'underpowered_fallback'
            if min((y == 1).sum(), (y == 0).sum()) >= 200:
                shape = grid[m['arm'] + '_shape'][ix]
                delta = shape[y == 1].mean(0) - shape[y == 0].mean(0)
                options = {j: auc(y, grid[m['arm']][ix, j]) for j, g in enumerate(screen.study.prior.GEOMETRIES)
                           if contrasts(delta, g)[1]}
                best, reason = (max(options, key=options.get), 'selected') if options else (0, 'no_protection_fallback')
            assert tuple(m['geometry']) == screen.study.prior.GEOMETRIES[best] and m['status'] == reason
    blacklist = pd.read_csv(config['blacklist']['path'], sep='\t', header=None)
    x = np.arange(-100, 101)
    checked, score_difference = 0, 0.
    for item in config['items']:
        original = read(item['sites']['path'])
        with np.load(config['control_profiles'][item['cell']]['path'], allow_pickle=False) as z:
            counts = [z['plus_observed'], z['minus_observed']]
            profile_valid = z['valid']
        with np.load(item['arrays']['path'], allow_pickle=False) as z:
            logits = [z['plus_log_bias'], z['minus_log_bias']]
        for m in [m for m in models if m['cell'] == item['cell']]:
            ix = np.flatnonzero(original.tf.eq(m['tf']))
            saved = rows[rows.cell.eq(m['cell']) & rows.tf.eq(m['tf']) & rows.arm.eq(m['arm'])].reset_index(drop=True)
            pd.testing.assert_frame_equal(saved[screen.study.safety.SITE_KEYS], original.iloc[ix][screen.study.safety.SITE_KEYS].reset_index(drop=True))
            g = m['geometry'] or screen.study.prior.GEOMETRIES[0]
            center, left, right = g
            inner = abs(x) <= center
            region = inner | ((abs(x) >= left) & (abs(x) <= right))
            n = [v[ix][:, region].sum(1).astype(int) for v in counts]
            k = sum(v[ix][:, inner].sum(1).astype(int) for v in counts)
            for name, expected in [('n_plus', n[0]), ('n_minus', n[1]), ('central', k)]:
                np.testing.assert_array_equal(saved[name], expected)
            for j, site in enumerate(original.iloc[ix].itertuples(index=False)):
                position = (int(site.TFBS_start) + int(site.TFBS_end)) // 2
                black = ((blacklist[0] == site.TFBS_chr) & (blacklist[1] < position + 101) & (blacklist[2] > position - 100)).any()
                assert bool(saved.iloc[j].window_blacklisted) == bool(black)
                probabilities = []
                for log in logits:
                    if m['arm'] == 'conditional_uniform':
                        probabilities.append(inner.sum() / region.sum())
                    else:
                        prop = np.exp(log[ix[j]] - np.max(log[ix[j], region]))
                        probabilities.append(prop[inner].sum() / prop[region].sum())
                np.testing.assert_allclose([saved.iloc[j].p_plus, saved.iloc[j].p_minus], probabilities, rtol=1e-12, atol=1e-12)
                eligible = bool(profile_valid[ix[j]] and not black and m['status'] != 'unavailable')
                informative = eligible and n[0][j] + n[1][j] > 0
                assert bool(saved.iloc[j].eligible) == eligible and bool(saved.iloc[j].informative) == informative
                if m['status'] != 'unavailable':
                    mass = np.convolve(binom.pmf(np.arange(n[0][j]+1), n[0][j], probabilities[0]),
                                       binom.pmf(np.arange(n[1][j]+1), n[1][j], probabilities[1]))
                    score = np.round(-np.log10(min(1., max(1e-300, mass[:k[j]+1].sum()))), 10)
                    difference = abs(score - saved.iloc[j].score)
                    score_difference = max(score_difference, difference)
                    assert difference <= 1.01e-10
                    assert bool(saved.iloc[j].called) == bool(informative and score > screen.PROTOCOL['threshold'])
                else:
                    assert pd.isna(saved.iloc[j].score) and not saved.iloc[j].called
                checked += 1
            r = summary[summary.cell.eq(m['cell']) & summary.tf.eq(m['tf']) & summary.arm.eq(m['arm'])].iloc[0]
            calls = int(saved.called.sum())
            assert r.sites == len(saved) and r.calls == calls
            assert r.eligible == saved.eligible.sum() and r.informative == saved.informative.sum()
            assert r.blacklisted == saved.window_blacklisted.sum()
            for prefix, mask in [('all_', saved.eligible), ('informative_', saved.informative)]:
                n_sites = int(mask.sum())
                if n_sites:
                    np.testing.assert_allclose(r[prefix + 'rate'], calls / n_sites, rtol=0, atol=1e-12)
                    np.testing.assert_allclose(r[prefix + 'wilson_upper'], binomtest(calls, n_sites).proportion_ci(method='wilson').high, rtol=0, atol=1e-12)
                else:
                    assert pd.isna(r[prefix + 'rate']) and pd.isna(r[prefix + 'wilson_upper'])
            assert bool(r.nominal_fpr_bound_met) == bool(m['status'] != 'unavailable' and r.all_wilson_upper <= .05 and r.informative_wilson_upper <= .05)
    assert checked == 8400
    print(f'Independent verification: 42 models, 8400 site rows and 42 rate/Wilson summaries; maximum score difference {score_difference}', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir', type=Path, required=True)
    verify(p.parse_args().outdir.resolve())
