#!/usr/bin/env python3
"""Freeze a training-only strand-conditioned central-depletion comparison."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crossfit_stable_geometry as prior
import conditional_protection_score as scoring
import audit_direct_bias_safety as safety
from crossfit_unit_floor_geometry import normalize_unit_floor
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = prior.ROOT
SOURCE = ROOT / 'benchmarks/results/stable_geometry_crossfit_20260922_v1'
PROFILE_DIR = ROOT / 'benchmarks/results/corrected_safety_profiles_20260922_v1'
DIAGNOSTIC = ROOT / 'benchmarks/results/low_count_normalization_20260922_v1/configuration.freeze.json'
ARMS = ['conditional_parametric', 'conditional_uniform']
PROTOCOL = dict(evidence='training_only_conditional_protection_crossfit', arms=ARMS,
                chromosomes=prior.PROTOCOL['chromosomes'], geometries=[list(g) for g in prior.GEOMETRIES],
                score='capped negative log10 lower tail of two independent strand-conditioned binomials; round 10 decimals',
                conditioning='per-strand cuts in center plus both shoulders; gap and other flanks excluded',
                assumptions='fixed correct sequence propensities and independent cuts; biological calibration unproven',
                geometry_rule='maximum training AUROC among unit-floor residual contrasts protective in both shoulders; minimum 200/class',
                fallback='conventional geometry for insufficient labels or no protective contrast',
                detectors=['conditional_fixed', 'conditional_selected'],
                support='all frozen matched sites required to have finite geometry scores; fail rather than silently filter',
                nested_matching=False, promotion_eligible=False, safety_established=False,
                max_score=scoring.MAX_SCORE, thresholds_fitted=False)


def freeze(out):
    if out.exists():
        raise ValueError('fresh output required')
    prior.execute(SOURCE, verify=True)
    old = json.loads((SOURCE / 'configuration.freeze.json').read_text())
    baseline = json.loads((SOURCE / 'crossfit_manifest.json').read_text())
    profile = json.loads((PROFILE_DIR / 'configuration.freeze.json').read_text())
    dispersion = json.loads(DIAGNOSTIC.read_text())['dispersion']
    items = []
    inputs = old['inputs'] + [record(SOURCE / 'configuration.freeze.json'), record(SOURCE / 'crossfit_manifest.json'),
                             record(PROFILE_DIR / 'configuration.freeze.json'), record(DIAGNOSTIC)] + list(baseline['outputs'].values())
    for item in profile['items']:
        if item['role'] != 'validation':
            continue
        path = PROFILE_DIR / (item['cell'] + '.validation.json')
        stage = json.loads(path.read_text())
        copied = dict(cell=item['cell'], sites=item['sites'], bias_arrays=item['arrays'], profiles=stage['outputs']['profiles'])
        items.append(copied)
        inputs += [record(path), item['manifest'], copied['sites'], copied['bias_arrays'], copied['profiles']]
    code = {r['path'] for r in old['code'] + profile['code']}
    code.update([str(Path(__file__).resolve()), str(Path(scoring.__file__).resolve()),
                 str(ROOT / 'benchmarks/scripts/crossfit_unit_floor_geometry.py')])
    config = dict(protocol=PROTOCOL, inputs=inputs, items=items, source_outputs=old['source_outputs'],
                  baseline_outputs=baseline['outputs'], dispersion=dispersion,
                  code=[record(Path(p)) for p in sorted(code)], runtime=old['runtime'])
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / 'configuration.freeze.json', config)


def preflight(out):
    ext = prior.source.cohorts.prior.wide.old.parent.extension
    config = ext.read_document(out / 'configuration.freeze.json', 'configuration_id')
    if config['protocol'] != PROTOCOL or config['runtime'] != ext.runtime():
        raise ValueError('protocol/runtime changed')
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    return config


def prepare(config, sites, panels):
    count_scores = {arm: np.zeros((len(sites), len(prior.GEOMETRIES))) for arm in ARMS}
    shape = {arm: np.zeros((len(sites), 201)) for arm in ARMS}
    for item in config['items']:
        cell = item['cell']
        ix = np.flatnonzero(sites.cell.eq(cell))
        all_sites = pd.read_csv(item['sites']['path'], sep='\t')
        order = pd.Index(safety.site_hashes(all_sites)).get_indexer(safety.site_hashes(sites.iloc[ix]))
        if (order < 0).any():
            raise ValueError('matched site missing')
        np.testing.assert_array_equal(all_sites.iloc[order].chip_label, sites.iloc[ix].chip_label)
        with np.load(item['profiles']['path'], allow_pickle=False) as archive:
            np.testing.assert_array_equal(archive['site_hash'][order], safety.site_hashes(sites.iloc[ix]))
            plus, minus = archive['plus_observed'][order], archive['minus_observed'][order]
            if not archive['valid'][order].all():
                raise ValueError('invalid matched count profile')
        np.testing.assert_array_equal(plus + minus, panels['raw'][ix])
        with np.load(item['bias_arrays']['path'], allow_pickle=False) as archive:
            np.testing.assert_array_equal(archive['site_hash'][order], safety.site_hashes(sites.iloc[ix]))
            plus_bias, minus_bias = archive['plus_log_bias'][order], archive['minus_log_bias'][order]
        for arm in ARMS:
            first, second = (plus_bias, minus_bias) if arm == 'conditional_parametric' else (np.zeros_like(plus), np.zeros_like(minus))
            for j, geometry in enumerate(prior.GEOMETRIES):
                score, components = scoring.profile_scores(plus, minus, first, second, geometry)
                if not components['valid'].all() or not np.isfinite(score).all():
                    raise ValueError('incomplete matched probability support; no automatic exclusion')
                count_scores[arm][ix, j] = score
            if arm == 'conditional_parametric':
                residual = panels['parametric_direct'][ix]
            else:
                _, _, residual = safety.direct_scores(plus + minus, np.ones_like(plus), config['dispersion'][cell])
            shape[arm][ix] = normalize_unit_floor(residual)
        print('Prepared count-conditioned scores:', cell, flush=True)
    return count_scores, shape


def select_geometry(values, scores, labels):
    if min(int((labels == 0).sum()), int((labels == 1).sum())) < 200:
        return 0, 'underpowered_fallback'
    delta = values[labels == 1].mean(0) - values[labels == 0].mean(0)
    candidates = [(prior.stable.rank_auc(labels, scores[:, j]), j) for j, g in enumerate(prior.GEOMETRIES)
                  if prior.geometry_contrast(delta, g)[3]]
    return (max(candidates, key=lambda x: x[0])[1], 'selected') if candidates else (0, 'no_protection_fallback')


def execute(out, verify=False):
    config = preflight(out)
    base = prior.source.cohorts.prior.wide.old.parent.base
    target = out / 'crossfit_manifest.json'
    if target.exists():
        base.verify_stage(target, config['configuration_id'], ROOT)
        print('Verified conditional-count crossfit', flush=True)
        return
    if verify or (out / 'score_grid.npz').exists():
        raise ValueError('missing or partial result preserved')
    sites = pd.read_csv(config['source_outputs']['sites']['path'], sep='\t')
    tasks = pd.read_csv(config['source_outputs']['metrics']['path'], sep='\t')[['cell', 'tf', 'matched_pairs', 'status']].drop_duplicates()
    if len(sites) != 8106 or len(tasks) != 21 or not sites.TFBS_chr.isin(PROTOCOL['chromosomes']).all():
        raise ValueError('training membership changed')
    with np.load(config['source_outputs']['panels']['path'], allow_pickle=False) as archive:
        panels = {arm: archive[arm] for arm in ['raw', 'parametric_direct']}
    count_scores, shape = prepare(config, sites, panels)
    np.savez_compressed(out / 'score_grid.npz', site_hash=safety.site_hashes(sites), **count_scores,
                        **{arm + '_shape': values for arm, values in shape.items()})
    folds, fits, held_scores = [], [], []
    for task in tasks.itertuples(index=False):
        ix = np.flatnonzero(sites.cell.eq(task.cell) & sites.tf.eq(task.tf))
        table = sites.iloc[ix].reset_index(drop=True)
        y = table.chip_label.to_numpy(int)
        for arm in ARMS:
            scores, values = count_scores[arm][ix], shape[arm][ix]
            for chrom in PROTOCOL['chromosomes']:
                fit = prior.fitting_mask(table, chrom).to_numpy()
                held = table.TFBS_chr.eq(chrom).to_numpy()
                index, reason = select_geometry(values[fit], scores[fit], y[fit])
                common = dict(cell=task.cell, tf=task.tf, arm=arm, chromosome=chrom,
                              cohort_status=task.status, matched_pairs=task.matched_pairs)
                fits.append(dict(common, geometry=json.dumps(prior.GEOMETRIES[index]), reason=reason,
                                 fit_positive=int((y[fit] == 1).sum()), fit_negative=int((y[fit] == 0).sum()),
                                 excluded_pair_mates=int((~fit & ~held).sum())))
                for detector in PROTOCOL['detectors']:
                    j = index if detector == 'conditional_selected' else 0
                    score, yy = scores[held, j], y[held]
                    available = bool((yy == 0).any() and (yy == 1).any())
                    row = dict(common, detector=detector, geometry=json.dumps(prior.GEOMETRIES[j]),
                               positive=int((yy == 1).sum()), negative=int((yy == 0).sum()),
                               evaluation_status='available' if available else 'missing_class')
                    if available:
                        delta = values[held][yy == 1].mean(0) - values[held][yy == 0].mean(0)
                        c, l, r, canonical = prior.geometry_contrast(delta, prior.GEOMETRIES[j])
                        row.update(auroc=roc_auc_score(yy, score), auprc=average_precision_score(yy, score),
                                   center=c, left=l, right=r, canonical=canonical,
                                   fixed_canonical=prior.geometry_contrast(delta, (5, 15, 40))[3])
                    folds.append(row)
                    frame = table.loc[held, prior.source.cohorts.prior.previous.KEYS + ['chip_label', 'pair_id']].copy()
                    frame['arm'], frame['detector'], frame['score'] = arm, detector, score
                    held_scores.append(frame)
    folds = pd.DataFrame(folds)
    if len(folds) != 1260 or len(fits) != 630:
        raise ValueError('task denominator changed')
    summary = folds.groupby(['cell', 'tf', 'arm', 'detector']).agg(
        evaluable_chromosomes=('auroc', 'count'), mean_auroc=('auroc', 'mean'), mean_auprc=('auprc', 'mean'),
        canonical_fraction=('canonical', 'mean'), fixed_canonical_fraction=('fixed_canonical', 'mean')).reset_index()
    outputs = {'grid': record(out / 'score_grid.npz')}
    for name, frame in [('folds', folds), ('fits', pd.DataFrame(fits)), ('summary', summary), ('scores', pd.concat(held_scores, ignore_index=True))]:
        path = out / (name + '.tsv')
        frame.to_csv(path, sep='\t', index=False)
        outputs[name] = record(path)
    result = dict(configuration_id=config['configuration_id'], outputs=outputs, training_only=True,
                  nested_matching=False, promotion_eligible=False, biological_null_calibration=False)
    result['stage_id'] = canonical_id(result)
    write_immutable(target, result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--stage', choices=['freeze', 'execute', 'verify'], required=True)
    args = parser.parse_args()
    out = args.outdir.resolve()
    if not out.is_relative_to(ROOT / 'benchmarks/results'):
        parser.error('use research results')
    if args.stage == 'freeze':
        freeze(out)
    else:
        with prior.source.cohorts.prior.wide.old.parent.extension.exclusive_run(out):
            execute(out, args.stage == 'verify')
