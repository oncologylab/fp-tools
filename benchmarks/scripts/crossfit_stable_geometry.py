#!/usr/bin/env python3
"""Frozen successor with explicit numerical ties for training geometry ranks."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import evaluate_radius_full_pool as source
from corrected_shape_candidates import GEOMETRIES, POSITIONS, geometry_weights
from fixed_shape_normalization import normalize_fixed
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path
import stable_geometry_scores as stable

ROOT = source.ROOT
SOURCE = ROOT / 'benchmarks/results/radius_full_pool_comparison_20260922_v1'
PROTOCOL = dict(
    evidence='training_only_leave_one_chromosome_out_stable_geometry_diagnostic',
    score_precision_decimals=stable.DECIMALS, protection_margin=stable.SHAPE_MARGIN,
    numerical_policy='round scores to 10 decimal places before fitting ranks and evaluation; exact rank-sum fitting AUROC; shape inequalities require margin 1e-10',
    chromosomes=[f'chr{i}' for i in range(1, 16)], arms=source.ARMS,
    geometries=[list(g) for g in GEOMETRIES], minimum_fit_per_class=200,
    selection='largest fitting AUROC among central-negative contrasts below both shoulders; first grid entry wins exact ties',
    fallback='conventional geometry when underpowered or no protective fitting contrast',
    fit_exclusion='held chromosome and both members of any pair touching it',
    normalization='per-site linear outer detrending, RMS scale, fallback 1, clip 12',
    detectors=['raw_fixed', 'normalized_fixed', 'normalized_selected'],
    evaluation='within-chromosome AUROC and AUPRC; no pooling differently fitted scores',
    summary='equal-weight mean over shared evaluable chromosomes; retain missing-class rows',
    matching='fixed full-training matched cohorts, not nested; descriptive development evidence only',
    zero_read_rule='all detector scores zero for zero-read 201-position windows',
    candidate_selection=False, promotion_eligible=False, uncertainty_calibrated=False)


def fitting_mask(sites, chromosome):
    """Called within one cell/TF task; exclude held sites and their pair mates."""
    held = sites.TFBS_chr.eq(chromosome)
    return ~sites.pair_id.isin(sites.loc[held, 'pair_id'])


def geometry_contrast(delta, geometry):
    center, start, end = geometry
    inner = float(delta[abs(POSITIONS) <= center].mean())
    left = float(delta[(POSITIONS >= -end) & (POSITIONS <= -start)].mean())
    right = float(delta[(POSITIONS >= start) & (POSITIONS <= end)].mean())
    margin = stable.SHAPE_MARGIN
    return inner, left, right, bool(inner < -margin and left - inner > margin and right - inner > margin)


def fit_geometry(profiles, labels, zero_read, minimum=200):
    """Only fitting rows enter selection; normalized inputs are per-site fixed."""
    values, y = np.asarray(profiles, float), np.asarray(labels)
    zero = np.asarray(zero_read, bool)
    if (values.shape != (len(y), 201) or zero.shape != y.shape
            or not np.isfinite(values).all() or not np.isin(y, [0, 1]).all()):
        raise ValueError('finite aligned fitting profiles and binary labels required')
    if min(int((y == 0).sum()), int((y == 1).sum())) < minimum:
        return GEOMETRIES[0], 'underpowered_fallback'
    delta = values[y == 1].mean(0) - values[y == 0].mean(0)
    eligible = []
    for geometry in GEOMETRIES:
        if geometry_contrast(delta, geometry)[3]:
            score = stable.stable_scores(values, geometry_weights(geometry), zero)
            eligible.append((stable.rank_auc(y, score), geometry))
    if not eligible:
        return GEOMETRIES[0], 'no_protection_fallback'
    return max(eligible, key=lambda item: item[0])[1], 'selected'


def freeze(out):
    if out.exists():
        raise ValueError('fresh directory required')
    source.execute(SOURCE, verify=True)
    old = json.loads((SOURCE / 'configuration.freeze.json').read_text())
    result = json.loads((SOURCE / 'comparison_manifest.json').read_text())
    config = dict(protocol=PROTOCOL, inputs=old['inputs'] + [record(SOURCE / 'configuration.freeze.json'),
        record(SOURCE / 'comparison_manifest.json')] + list(result['outputs'].values()),
        source_outputs=result['outputs'], code=old['code'] + [record(Path(__file__)), record(Path(stable.__file__))],
        runtime=source.cohorts.prior.wide.old.parent.extension.runtime(), promotion_eligible=False)
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / 'configuration.freeze.json', config)
    print('Frozen chromosome-wise geometry diagnostic', flush=True)


def execute(out, verify=False):
    config = source.cohorts.prior.wide.old.parent.extension.read_document(out / 'configuration.freeze.json', 'configuration_id')
    if config['protocol'] != PROTOCOL or config['runtime'] != source.cohorts.prior.wide.old.parent.extension.runtime():
        raise ValueError('protocol or runtime changed')
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    target = out / 'crossfit_manifest.json'
    if target.exists():
        source.cohorts.prior.wide.old.parent.base.verify_stage(target, config['configuration_id'], ROOT)
        print('Verified chromosome-wise geometry diagnostic', flush=True)
        return
    if verify or list(out.glob('*.tsv')):
        raise ValueError('missing or partial output')
    paths = config['source_outputs']
    read = lambda name: pd.read_csv(paths[name]['path'], sep='\t', float_precision='round_trip')
    sites, metrics = read('sites'), read('metrics')
    if not sites.TFBS_chr.isin(PROTOCOL['chromosomes']).all():
        raise ValueError('nontraining coordinate')
    with np.load(paths['panels']['path'], allow_pickle=False) as archive:
        panels = {arm: archive[arm] for arm in source.ARMS}
    rows, score_rows, fits = [], [], []
    for task in metrics[['cell', 'tf', 'matched_pairs', 'status']].drop_duplicates().itertuples(index=False):
        ix = np.flatnonzero(sites.cell.eq(task.cell) & sites.tf.eq(task.tf))
        table = sites.iloc[ix].reset_index(drop=True)
        y = table.chip_label.to_numpy(int)
        zero = panels['raw'][ix].sum(1) == 0
        for arm in source.ARMS:
            raw = panels[arm][ix]
            norm = normalize_fixed(raw, POSITIONS, 1.) if len(ix) else raw.copy()
            for chrom in PROTOCOL['chromosomes']:
                fit = fitting_mask(table, chrom).to_numpy()
                held = table.TFBS_chr.eq(chrom).to_numpy()
                geometry, reason = fit_geometry(norm[fit], y[fit], zero[fit], PROTOCOL['minimum_fit_per_class'])
                common = dict(cell=task.cell, tf=task.tf, arm=arm, chromosome=chrom,
                              cohort_status=task.status, matched_pairs=task.matched_pairs)
                fits.append(dict(common, geometry=json.dumps(geometry), reason=reason,
                    fit_positive=int((y[fit] == 1).sum()), fit_negative=int((y[fit] == 0).sum()),
                    excluded_pair_mates=int((~fit & ~held).sum())))
                for detector in PROTOCOL['detectors']:
                    g = geometry if detector == 'normalized_selected' else GEOMETRIES[0]
                    v = raw if detector == 'raw_fixed' else norm
                    scores = stable.stable_scores(v[held], geometry_weights(g), zero[held])
                    yy = y[held]
                    available = bool((yy == 1).any() and (yy == 0).any())
                    row = dict(common, detector=detector, geometry=json.dumps(g),
                        positive=int((yy == 1).sum()), negative=int((yy == 0).sum()),
                        evaluation_status='available' if available else 'missing_class')
                    if available:
                        delta = norm[held][yy == 1].mean(0) - norm[held][yy == 0].mean(0)
                        center, left, right, canonical = geometry_contrast(delta, g)
                        row.update(auroc=roc_auc_score(yy, scores), auprc=average_precision_score(yy, scores),
                                   center=center, left=left, right=right, canonical=canonical)
                    rows.append(row)
                    frame = table.loc[held, source.cohorts.prior.previous.KEYS + ['chip_label', 'pair_id']].copy()
                    frame['arm'], frame['detector'], frame['score'] = arm, detector, scores
                    score_rows.append(frame)
        print('Completed cross-fitting:', task.cell, task.tf, flush=True)
    if len(rows) != 21 * 7 * 15 * 3:
        raise ValueError('task denominator changed')
    results = pd.DataFrame(rows)
    summary = results.groupby(['cell', 'tf', 'arm', 'detector'], sort=True).agg(
        evaluable_chromosomes=('auroc', 'count'), mean_auroc=('auroc', 'mean'),
        mean_auprc=('auprc', 'mean'), canonical_fraction=('canonical', 'mean')).reset_index()
    outputs = {}
    for name, table in [('folds', results), ('fits', pd.DataFrame(fits)), ('summary', summary),
                        ('scores', pd.concat(score_rows, ignore_index=True))]:
        path = out / (name + '.tsv')
        table.to_csv(path, sep='\t', index=False)
        outputs[name] = record(path)
    result = dict(configuration_id=config['configuration_id'], outputs=outputs,
        training_only=True, nested_matching=False, promotion_eligible=False, uncertainty_calibrated=False)
    result['stage_id'] = canonical_id(result)
    write_immutable(target, result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--freeze-only', action='store_true')
    mode.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    out = args.outdir.resolve()
    if not out.is_relative_to(ROOT / 'benchmarks/results'):
        parser.error('use research results')
    if args.freeze_only:
        freeze(out)
    else:
        with source.cohorts.prior.wide.old.parent.extension.exclusive_run(out):
            execute(out, args.verify_only)
