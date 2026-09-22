#!/usr/bin/env python3
"""Frozen score-capacity diagnostic at observed strand totals; not power prediction."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import screen_conditional_protection_rep1 as source
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = source.ROOT
TRAINING = ROOT / 'benchmarks/results/conditional_protection_20260922_v1'
CONTROL = ROOT / 'benchmarks/results/conditional_protection_rep1_20260922_v1'
THRESHOLD = source.PROTOCOL['threshold']
PROTOCOL = dict(evidence='training_and_examined_replicate1_conditional_capacity',
                definition='maximum lower-tail score occurs at zero central cuts, fixed observed strand totals and propensities',
                threshold=THRESHOLD, threshold_changed=False, geometry_changed=False,
                depth_prediction=False, biological_power_estimated=False, promotion_eligible=False)


def score_ceiling(n_plus, n_minus, p_plus, p_minus):
    n1, n2, p1, p2 = np.broadcast_arrays(*[np.asarray(v, float) for v in [n_plus, n_minus, p_plus, p_minus]])
    for n in [n1, n2]:
        if not np.isfinite(n).all() or (n < 0).any() or (n != np.floor(n)).any():
            raise ValueError('nonnegative integer totals required')
    for p in [p1, p2]:
        if not np.isfinite(p).all() or (p < 0).any() or (p > 1).any():
            raise ValueError('finite probabilities in [0,1] required')
    total = np.zeros(n1.shape)
    for n, p in [(n1, p1), (n2, p2)]:
        with np.errstate(divide='ignore'):
            contribution = np.zeros(n.shape)
            np.multiply(n, -np.log1p(-p) / np.log(10.), out=contribution, where=n > 0)
        total += contribution
    return np.round(np.minimum(total, source.study.scoring.MAX_SCORE), 10)


def annotate(frame, components, scores, eligible):
    table = frame.copy()
    for key in ['n_plus', 'n_minus', 'central', 'p_plus', 'p_minus']:
        table[key] = components[key]
    ceiling = score_ceiling(*(components[k] for k in ['n_plus', 'n_minus', 'p_plus', 'p_minus']))
    score = np.asarray(scores, float)
    eligible = np.asarray(eligible, bool)
    available = np.isfinite(score)
    if ((score[available] - ceiling[available]) > 1.01e-10).any():
        raise ValueError('observed score exceeds analytic ceiling')
    ceiling[~available] = np.nan
    table['score'], table['ceiling'], table['eligible'] = score, ceiling, eligible
    table['informative'] = eligible & ((table.n_plus + table.n_minus) > 0)
    table['reachable'] = eligible & (ceiling > THRESHOLD)
    table['called'] = eligible & (score > THRESHOLD)
    table['count_limited'] = eligible & (ceiling <= THRESHOLD)
    table['reachable_uncalled'] = table.reachable & ~table.called
    if (table.called & ~table.reachable).any():
        raise ValueError('called site classified unreachable')
    return table


def freeze(out):
    if out.exists():
        raise ValueError('fresh output required')
    source.study.execute(TRAINING, verify=True)
    source.execute(CONTROL, verify=True)
    train = source.study.preflight(TRAINING)
    control = source.preflight(CONTROL)
    first = json.loads((TRAINING / 'crossfit_manifest.json').read_text())
    second = json.loads((CONTROL / 'screen_manifest.json').read_text())
    config = dict(protocol=PROTOCOL, training=train,
                  train_outputs=first['outputs'], control_outputs=second['outputs'],
                  inputs=train['inputs'] + control['inputs'] + list(first['outputs'].values()) + list(second['outputs'].values()) +
                  [record(TRAINING / 'crossfit_manifest.json'), record(CONTROL / 'screen_manifest.json'), second['model_artifact']],
                  code=train['code'] + control['code'] + [record(Path(__file__))], runtime=train['runtime'])
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / 'configuration.freeze.json', config)


def preflight(out):
    ext = source.study.prior.source.cohorts.prior.wide.old.parent.extension
    config = ext.read_document(out / 'configuration.freeze.json', 'configuration_id')
    if config['protocol'] != PROTOCOL or config['runtime'] != ext.runtime():
        raise ValueError('protocol/runtime changed')
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    return config


def execute(out, verify=False):
    config = preflight(out)
    base = source.study.prior.source.cohorts.prior.wide.old.parent.base
    target = out / 'capacity_manifest.json'
    if target.exists():
        base.verify_stage(target, config['configuration_id'], ROOT)
        print('Verified conditional score-capacity diagnostic', flush=True)
        return
    if verify or list(out.glob('*.tsv')):
        raise ValueError('missing or partial result preserved')
    read = lambda entry: pd.read_csv(entry['path'], sep='\t', float_precision='round_trip')
    sites = read(config['training']['source_outputs']['sites'])
    folds, scores = read(config['train_outputs']['folds']), read(config['train_outputs']['scores'])
    groups = {k: v for k, v in scores.groupby(['cell', 'tf', 'arm', 'TFBS_chr', 'detector'])}
    rows = []
    for item in config['training']['items']:
        subset = sites[sites.cell.eq(item['cell'])].reset_index(drop=True)
        all_sites = read(item['sites'])
        order = pd.Index(source.study.safety.site_hashes(all_sites)).get_indexer(source.study.safety.site_hashes(subset))
        if (order < 0).any():
            raise ValueError('site identity missing')
        with np.load(item['profiles']['path'], allow_pickle=False) as z:
            plus, minus = z['plus_observed'][order], z['minus_observed'][order]
        with np.load(item['bias_arrays']['path'], allow_pickle=False) as z:
            bias = [z['plus_log_bias'][order], z['minus_log_bias'][order]]
        cache = {}
        for row in folds[folds.cell.eq(item['cell'])].itertuples(index=False):
            ix = np.flatnonzero(subset.tf.eq(row.tf) & subset.TFBS_chr.eq(row.chromosome))
            if not len(ix):
                continue
            geometry = tuple(json.loads(row.geometry))
            key = (row.arm, geometry)
            if key not in cache:
                first, second = bias if row.arm == 'conditional_parametric' else [np.zeros_like(plus), np.zeros_like(minus)]
                cache[key] = source.study.scoring.region_components(plus, minus, first, second, geometry)
            comp = {k: v[ix] for k, v in cache[key].items()}
            saved = groups[(row.cell, row.tf, row.arm, row.chromosome, row.detector)]
            cols = source.study.prior.source.cohorts.prior.previous.KEYS + ['chip_label', 'pair_id']
            pd.testing.assert_frame_equal(saved[cols].reset_index(drop=True), subset.iloc[ix][cols].reset_index(drop=True))
            frame = saved[cols + ['arm', 'detector']].copy()
            frame['source'], frame['group'], frame['geometry'] = 'training_held', np.where(frame.chip_label == 1, 'positive', 'negative'), row.geometry
            rows.append(annotate(frame, comp, saved.score.to_numpy(), np.ones(len(ix), bool)))
    control = read(config['control_outputs']['sites'])
    frame = control[source.study.safety.SITE_KEYS + ['arm', 'model_status', 'window_blacklisted']].copy()
    frame['source'], frame['group'], frame['detector'] = 'naked_rep1', 'control', 'full_training_selected'
    rows.append(annotate(frame, {k: control[k].to_numpy() for k in ['n_plus', 'n_minus', 'central', 'p_plus', 'p_minus']},
                         control.score.to_numpy(), control.eligible.to_numpy()))
    data = pd.concat(rows, ignore_index=True)
    if len(data) != 40824:
        raise ValueError('row denominator changed')
    keys = ['source', 'cell', 'tf', 'arm', 'detector', 'group']
    summary = data.groupby(keys).agg(sites=('ceiling', 'size'), eligible=('eligible', 'sum'),
        informative=('informative', 'sum'), reachable=('reachable', 'sum'), called=('called', 'sum'),
        count_limited=('count_limited', 'sum'), reachable_uncalled=('reachable_uncalled', 'sum'),
        median_ceiling=('ceiling', 'median'), median_n_plus=('n_plus', 'median'), median_n_minus=('n_minus', 'median')).reset_index()
    # Retain unavailable training tasks explicitly despite their empty site sets.
    missing = folds[folds.matched_pairs.eq(0)][['cell', 'tf', 'arm', 'detector']].drop_duplicates()
    empty = []
    for row in missing.to_dict('records'):
        for group in ['positive', 'negative']:
            empty.append(dict(row, source='training_held', group=group, sites=0, eligible=0, informative=0,
                              reachable=0, called=0, count_limited=0, reachable_uncalled=0,
                              median_ceiling=np.nan, median_n_plus=np.nan, median_n_minus=np.nan))
    summary = pd.concat([summary, pd.DataFrame(empty)], ignore_index=True).sort_values(keys)
    if len(summary) != 210:
        raise ValueError('task/group denominator changed')
    outputs = {}
    for name, table in [('sites', data), ('summary', summary)]:
        path = out / (name + '.tsv')
        table.to_csv(path, sep='\t', index=False)
        outputs[name] = record(path)
    result = dict(configuration_id=config['configuration_id'], outputs=outputs, rows=len(data),
                  summary_rows=len(summary), thresholds_changed=False, depth_prediction=False, promotion_eligible=False)
    result['stage_id'] = canonical_id(result)
    write_immutable(target, result)
    print(f'Completed capacity audit: {len(data)} rows, {len(summary)} task/group summaries', flush=True)


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
        with source.study.prior.source.cohorts.prior.wide.old.parent.extension.exclusive_run(out):
            execute(out, args.stage == 'verify')
