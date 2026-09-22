#!/usr/bin/env python3
"""Fixed-geometry cleavage calibration on examined naked-DNA replicate 1."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parent))
import diagnose_low_count_normalization as source
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = source.ROOT
SOURCE = ROOT / 'benchmarks/results/low_count_normalization_20260922_v1'
POSITIONS = np.arange(-100, 101)
REGIONS = {'center': abs(POSITIONS) <= 7, 'left': (POSITIONS >= -40) & (POSITIONS <= -12),
           'right': (POSITIONS >= 12) & (POSITIONS <= 40)}
CHROMS = [f'chr{i}' for i in range(1, 16)]
PROTOCOL = dict(evidence='examined_replicate1_fixed_geometry_calibration', geometry=[7, 12, 40],
                regions=list(REGIONS), arms=['parametric', 'uniform'], strands=['plus', 'minus', 'pooled'],
                chromosomes=CHROMS, bootstrap_draws=2000, bootstrap_seed=2026,
                resampling='shared multinomial chromosome-block weights; no redraw of zero-count draws',
                predictive_distribution='sum independent site-strand binomials; diagnostic assumption, not validated biological coverage',
                paired_metric='parametric minus uniform conditional log likelihood per cut, natural logarithm',
                model_fitted=False, thresholds_changed=False, promotion_eligible=False,
                final_control_read=False, new_holdout_labels_read=False)


def region_probabilities(logits):
    values = np.asarray(logits, float)
    if values.ndim != 2 or values.shape[1] != 201:
        raise ValueError('201-position sequence logits required')
    region = np.logical_or.reduce(list(REGIONS.values()))
    valid = np.isfinite(values[:, region]).all(1)
    result = np.full((len(values), 3), np.nan)
    denominator = logsumexp(values[valid][:, region], axis=1)
    for j, mask in enumerate(REGIONS.values()):
        result[valid, j] = np.exp(logsumexp(values[valid][:, mask], axis=1) - denominator)
    return result, valid


def predictive_interval(counts, probabilities):
    counts, probabilities = np.asarray(counts), np.asarray(probabilities, float)
    if (counts.shape != probabilities.shape or counts.ndim != 1 or not np.isfinite(counts).all()
            or (counts < 0).any() or (counts != np.floor(counts)).any() or not np.isfinite(probabilities).all()
            or (probabilities < 0).any() or (probabilities > 1).any()):
        raise ValueError('aligned integer counts and probabilities required')
    mass = np.array([1.])
    for n, p in zip(counts.astype(int), probabilities):
        if n:
            mass = np.convolve(mass, binom.pmf(np.arange(n+1), n, p))
    mass /= mass.sum()
    cumulative = mass.cumsum()
    return int(np.searchsorted(cumulative, .025)), int(np.searchsorted(cumulative, .975))


def freeze(out):
    if out.exists():
        raise ValueError('fresh output required')
    source.execute(SOURCE, verify=True)
    old = source.preflight(SOURCE)
    profiles = json.loads((SOURCE / 'profiles_manifest.json').read_text())
    config = dict(protocol=PROTOCOL, items=old['items'], profiles=profiles['outputs'], blacklist=old['blacklist'],
                  inputs=old['inputs'] + [record(SOURCE / 'configuration.freeze.json'), record(SOURCE / 'profiles_manifest.json')] + list(profiles['outputs'].values()),
                  code=old['code'] + [record(Path(__file__))], runtime=old['runtime'])
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / 'configuration.freeze.json', config)


def preflight(out):
    config = source.extraction.extension.read_document(out / 'configuration.freeze.json', 'configuration_id')
    if config['protocol'] != PROTOCOL or config['runtime'] != source.extraction.extension.runtime():
        raise ValueError('protocol/runtime changed')
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    return config


def execute(out, verify=False):
    config = preflight(out)
    target = out / 'calibration_manifest.json'
    if target.exists():
        source.extraction.base.verify_stage(target, config['configuration_id'], ROOT)
        print('Verified fixed-region replicate-1 calibration', flush=True)
        return
    if verify or list(out.glob('*.tsv')):
        raise ValueError('missing or partial result preserved')
    intervals = [(v[0], int(v[1]), int(v[2])) for v in
                 (line.split('\t') for line in Path(config['blacklist']['path']).read_text().splitlines())]
    rows = []
    for item in config['items']:
        sites = pd.read_csv(item['sites']['path'], sep='\t')
        if not sites.TFBS_chr.isin(CHROMS).all() or any('label' in c.lower() or 'chip' in c.lower() for c in sites.columns):
            raise ValueError('only label-free training-chromosome controls allowed')
        black = source.exclusions.overlap_flags(sites, intervals, 100)
        with np.load(config['profiles'][item['cell']]['path'], allow_pickle=False) as z:
            np.testing.assert_array_equal(z['source_index'], np.arange(len(sites)))
            np.testing.assert_array_equal(z['site_hash'], source.safety.site_hashes(sites))
            counts = {strand: z[strand + '_observed'] for strand in ['plus', 'minus']}
            profile_valid = z['valid']
        with np.load(item['arrays']['path'], allow_pickle=False) as z:
            logits = {strand: z[strand + '_log_bias'] for strand in ['plus', 'minus']}
        for strand in ['plus', 'minus']:
            probabilities, valid = region_probabilities(logits[strand])
            observed = np.column_stack([counts[strand][:, mask].sum(1) for mask in REGIONS.values()])
            total = observed.sum(1)
            eligible = profile_valid & valid & ~black
            for arm in ['parametric', 'uniform']:
                p = probabilities if arm == 'parametric' else np.tile([mask.sum()/73 for mask in REGIONS.values()], (len(sites), 1))
                for j, region in enumerate(REGIONS):
                    frame = sites[source.safety.SITE_KEYS].copy()
                    for name, value in dict(strand=strand, arm=arm, region=region, eligible=eligible,
                                           window_blacklisted=black, total=total, observed=observed[:, j],
                                           probability=p[:, j], expected=total*p[:, j]).items():
                        frame[name] = value
                    rows.append(frame)
    sites = pd.concat(rows, ignore_index=True)
    if len(sites) != 50400:
        raise ValueError('site denominator changed')
    weights = np.random.default_rng(PROTOCOL['bootstrap_seed']).multinomial(15, np.repeat(1/15, 15), size=PROTOCOL['bootstrap_draws'])
    summaries, likelihoods, blocks = [], [], []
    for (cell, tf), task in sites.groupby(['cell', 'tf'], sort=True):
        for strand in PROTOCOL['strands']:
            eligible = task[task.eligible & (True if strand == 'pooled' else task.strand.eq(strand))]
            for arm in PROTOCOL['arms']:
                for region in REGIONS:
                    subset = eligible[eligible.arm.eq(arm) & eligible.region.eq(region)]
                    block = subset.groupby('TFBS_chr')[['total', 'observed', 'expected']].sum().reindex(CHROMS, fill_value=0)
                    for chrom, r in block.iterrows():
                        blocks.append(dict(cell=cell, tf=tf, strand=strand, arm=arm, region=region, chromosome=chrom,
                                           total=r.total, observed=r.observed, expected=r.expected))
                    total, observed, expected = block.sum().to_numpy()
                    sampled_total = weights @ block.total.to_numpy()
                    sampled_error = weights @ (block.observed - block.expected).to_numpy()
                    usable = sampled_total > 0
                    error = sampled_error[usable] / sampled_total[usable]
                    lower, upper = np.quantile(error, [.025, .975]) if len(error) else (np.nan, np.nan)
                    pl, pu = predictive_interval(subset.total.to_numpy(), subset.probability.to_numpy())
                    summaries.append(dict(cell=cell, tf=tf, strand=strand, arm=arm, region=region,
                        site_strand_rows=len(subset), informative_site_strands=int((subset.total > 0).sum()),
                        total_cuts=total, observed=observed, expected=expected,
                        observed_fraction=observed/total if total else np.nan, expected_fraction=expected/total if total else np.nan,
                        error_fraction=(observed-expected)/total if total else np.nan,
                        block_error_low=lower, block_error_high=upper, valid_bootstrap_draws=int(usable.sum()),
                        independent_null_count_low=pl, independent_null_count_high=pu))
            # Both arms use identical observed categories and eligibility.
            parametric = eligible[eligible.arm.eq('parametric')].copy()
            uniform_p = parametric.region.map({key: value.sum()/73 for key, value in REGIONS.items()}).to_numpy()
            terms = np.zeros(len(parametric))
            nonzero = parametric.observed.to_numpy() > 0
            if (parametric.loc[nonzero, 'probability'] <= 0).any():
                raise ValueError('positive count has zero predicted probability')
            terms[nonzero] = parametric.observed.to_numpy()[nonzero] * np.log(parametric.probability.to_numpy()[nonzero] / uniform_p[nonzero])
            parametric['log_gain'] = terms
            gain = parametric.groupby('TFBS_chr').log_gain.sum().reindex(CHROMS, fill_value=0).to_numpy()
            # Each strand total occurs in three category rows; count it once.
            count = parametric[parametric.region.eq('center')].groupby('TFBS_chr').total.sum().reindex(CHROMS, fill_value=0).to_numpy()
            sampled_count = weights @ count
            usable = sampled_count > 0
            gains = (weights @ gain)[usable] / sampled_count[usable]
            lo, hi = np.quantile(gains, [.025, .975]) if len(gains) else (np.nan, np.nan)
            likelihoods.append(dict(cell=cell, tf=tf, strand=strand, total_cuts=int(count.sum()),
                log_likelihood_gain=float(gain.sum()), gain_per_cut=gain.sum()/count.sum() if count.sum() else np.nan,
                block_gain_low=lo, block_gain_high=hi, valid_bootstrap_draws=int(usable.sum())))
    if len(summaries) != 378 or len(likelihoods) != 63 or len(blocks) != 5670:
        raise ValueError('task/strand/region denominator changed')
    outputs = {}
    for name, frame in [('sites', sites), ('calibration', pd.DataFrame(summaries)),
                        ('likelihood', pd.DataFrame(likelihoods)), ('blocks', pd.DataFrame(blocks))]:
        path = out / (name + '.tsv')
        frame.to_csv(path, sep='\t', index=False)
        outputs[name] = record(path)
    path = out / 'bootstrap_weights.npy'
    np.save(path, weights)
    outputs['weights'] = record(path)
    result = dict(configuration_id=config['configuration_id'], outputs=outputs, site_rows=len(sites),
                  calibration_rows=len(summaries), likelihood_rows=len(likelihoods),
                  multiple_comparisons_adjusted=False, biological_coverage_calibrated=False,
                  model_refitted=False, promotion_eligible=False)
    result['stage_id'] = canonical_id(result)
    write_immutable(target, result)
    print(f'Completed calibration: {len(sites)} site-region rows, {len(summaries)} calibration rows, {len(likelihoods)} likelihood rows', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--stage', choices=['freeze', 'execute', 'verify'], required=True)
    args = p.parse_args()
    out = args.outdir.resolve()
    if not out.is_relative_to(ROOT / 'benchmarks/results'):
        p.error('use research results')
    if args.stage == 'freeze':
        freeze(out)
    else:
        with source.extraction.extension.exclusive_run(out):
            execute(out, args.stage == 'verify')
