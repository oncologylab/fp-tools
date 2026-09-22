#!/usr/bin/env python3
"""Verify ceilings through probability products and observed scores by convolution."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_conditional_score_ceiling as study


def verify(out):
    study.execute(out, verify=True)
    config = study.preflight(out)
    result = json.loads((out / 'capacity_manifest.json').read_text())
    rows = pd.read_csv(result['outputs']['sites']['path'], sep='\t', float_precision='round_trip', low_memory=False)
    summary = pd.read_csv(result['outputs']['summary']['path'], sep='\t', float_precision='round_trip')
    n1, n2 = rows.n_plus.to_numpy(int), rows.n_minus.to_numpy(int)
    p1, p2 = rows.p_plus.to_numpy(), rows.p_minus.to_numpy()
    probability = (1-p1)**n1 * (1-p2)**n2
    ceiling = np.round(-np.log10(np.clip(probability, 1e-300, 1)), 10)
    available = rows.score.notna().to_numpy()
    ceiling[~available] = np.nan
    np.testing.assert_allclose(ceiling, rows.ceiling, rtol=0, atol=1.01e-10, equal_nan=True)
    eligible = rows.eligible.to_numpy()
    reachable = eligible & (ceiling > study.THRESHOLD)
    np.testing.assert_array_equal(reachable, rows.reachable)
    np.testing.assert_array_equal(eligible & (ceiling <= study.THRESHOLD), rows.count_limited)
    np.testing.assert_array_equal(eligible & ((n1+n2) > 0), rows.informative)
    max_difference = 0.
    for i in np.flatnonzero(available):
        mass = np.convolve(binom.pmf(np.arange(n1[i]+1), n1[i], p1[i]),
                           binom.pmf(np.arange(n2[i]+1), n2[i], p2[i]))
        tail = min(1., max(1e-300, mass[:int(rows.iloc[i].central)+1].sum()))
        score = np.round(-np.log10(tail), 10)
        max_difference = max(max_difference, abs(score - rows.iloc[i].score))
        assert abs(score - rows.iloc[i].score) <= 1.01e-10
        assert score <= ceiling[i] + 1.01e-10
        assert bool(rows.iloc[i].called) == bool(eligible[i] and score > study.THRESHOLD)
    np.testing.assert_array_equal(rows.reachable_uncalled, reachable & ~rows.called.to_numpy())
    assert not (rows.called & ~rows.reachable).any()
    keys = ['source', 'cell', 'tf', 'arm', 'detector', 'group']
    for r in summary.itertuples(index=False):
        mask = np.ones(len(rows), bool)
        for key in keys:
            mask &= rows[key].eq(getattr(r, key)).to_numpy()
        group = rows.loc[mask]
        assert len(group) == r.sites
        for key in ['eligible', 'informative', 'reachable', 'called', 'count_limited', 'reachable_uncalled']:
            assert int(group[key].sum()) == getattr(r, key)
        for key in ['ceiling', 'n_plus', 'n_minus']:
            expected = group[key].median() if len(group) and group[key].notna().any() else np.nan
            np.testing.assert_allclose(getattr(r, 'median_' + key), expected, equal_nan=True)
    assert len(rows) == 40824 and len(summary) == 210 and summary.sites.sum() == len(rows)
    print(f'Independent capacity verification: 40824 rows and 210 summaries; maximum observed-score difference {max_difference}', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir', type=Path, required=True)
    verify(p.parse_args().outdir.resolve())
