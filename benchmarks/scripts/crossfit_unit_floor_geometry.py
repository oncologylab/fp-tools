#!/usr/bin/env python3
"""Training-only comparison of fixed unit-floor normalization and geometry."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crossfit_stable_geometry as prior
from fixed_shape_normalization import detrend
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = prior.ROOT
SOURCE = ROOT / 'benchmarks/results/stable_geometry_crossfit_20260922_v1'
PROTOCOL = dict(evidence='training_only_unit_floor_geometry_crossfit',
                chromosomes=prior.PROTOCOL['chromosomes'], arms=prior.source.ARMS,
                normalization='outer linear detrend; denominator=max(RMS,1); clip at 12',
                floor=1., floor_fitted=False, geometry_rule=prior.PROTOCOL,
                detectors=['unit_floor_fixed', 'unit_floor_selected'],
                comparison='paired archived current-normalization folds; raw conventional DWM retained',
                shape='selected-geometry and fixed center +/-5, shoulders 15..40; descriptive held-fold contrasts',
                promotion_eligible=False, uncertainty_calibrated=False, nested_matching=False)


def normalize_unit_floor(values):
    residual, scale = detrend(values, prior.POSITIONS)
    return np.clip(residual / np.maximum(scale, 1.)[:, None], -12., 12.)


def freeze(out):
    if out.exists():
        raise ValueError('fresh output required')
    prior.execute(SOURCE, verify=True)
    old = json.loads((SOURCE / 'configuration.freeze.json').read_text())
    manifest = json.loads((SOURCE / 'crossfit_manifest.json').read_text())
    config = dict(protocol=PROTOCOL, inputs=old['inputs'] + [record(SOURCE / 'configuration.freeze.json'),
                  record(SOURCE / 'crossfit_manifest.json')] + list(manifest['outputs'].values()),
                  source_outputs=old['source_outputs'], baseline_outputs=manifest['outputs'],
                  code=old['code'] + [record(Path(__file__))], runtime=old['runtime'])
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


def execute(out, verify=False):
    config = preflight(out)
    target = out / 'crossfit_manifest.json'
    if target.exists():
        prior.source.cohorts.prior.wide.old.parent.base.verify_stage(target, config['configuration_id'], ROOT)
        print('Verified unit-floor crossfit', flush=True)
        return
    if verify or list(out.glob('*.tsv')):
        raise ValueError('missing or partial result preserved')
    paths = config['source_outputs']
    sites = pd.read_csv(paths['sites']['path'], sep='\t')
    metrics = pd.read_csv(paths['metrics']['path'], sep='\t')
    if len(sites) != 8106 or not sites.TFBS_chr.isin(PROTOCOL['chromosomes']).all():
        raise ValueError('training membership changed')
    with np.load(paths['panels']['path'], allow_pickle=False) as archive:
        panels = {arm: archive[arm] for arm in PROTOCOL['arms']}
    folds, fits, scores, curves = [], [], [], []
    for task in metrics[['cell', 'tf', 'matched_pairs', 'status']].drop_duplicates().itertuples(index=False):
        ix = np.flatnonzero(sites.cell.eq(task.cell) & sites.tf.eq(task.tf))
        table = sites.iloc[ix].reset_index(drop=True)
        labels = table.chip_label.to_numpy(int)
        zero = panels['raw'][ix].sum(1) == 0
        for arm in PROTOCOL['arms']:
            values = normalize_unit_floor(panels[arm][ix]) if len(ix) else panels[arm][ix].copy()
            if len(ix):
                delta = values[labels == 1].mean(0) - values[labels == 0].mean(0)
                curves.append(pd.DataFrame(dict(cell=task.cell, tf=task.tf, arm=arm,
                                                position=prior.POSITIONS, contrast=delta)))
            for chromosome in PROTOCOL['chromosomes']:
                fit = prior.fitting_mask(table, chromosome).to_numpy()
                held = table.TFBS_chr.eq(chromosome).to_numpy()
                geometry, reason = prior.fit_geometry(values[fit], labels[fit], zero[fit], 200)
                common = dict(cell=task.cell, tf=task.tf, arm=arm, chromosome=chromosome,
                              cohort_status=task.status, matched_pairs=task.matched_pairs)
                fits.append(dict(common, geometry=json.dumps(geometry), reason=reason,
                                 fit_positive=int((labels[fit] == 1).sum()), fit_negative=int((labels[fit] == 0).sum()),
                                 excluded_pair_mates=int((~fit & ~held).sum())))
                y = labels[held]
                available = bool((y == 1).any() and (y == 0).any())
                for detector in PROTOCOL['detectors']:
                    g = geometry if detector == 'unit_floor_selected' else prior.GEOMETRIES[0]
                    score = prior.stable.stable_scores(values[held], prior.geometry_weights(g), zero[held])
                    row = dict(common, detector=detector, geometry=json.dumps(g),
                               positive=int((y == 1).sum()), negative=int((y == 0).sum()),
                               evaluation_status='available' if available else 'missing_class')
                    if available:
                        delta = values[held][y == 1].mean(0) - values[held][y == 0].mean(0)
                        center, left, right, canonical = prior.geometry_contrast(delta, g)
                        fc, fl, fr, fp = prior.geometry_contrast(delta, (5, 15, 40))
                        row.update(auroc=roc_auc_score(y, score), auprc=average_precision_score(y, score),
                                   center=center, left=left, right=right, canonical=canonical,
                                   fixed_center=fc, fixed_left=fl, fixed_right=fr, fixed_canonical=fp)
                    folds.append(row)
                    frame = table.loc[held, prior.source.cohorts.prior.previous.KEYS + ['chip_label', 'pair_id']].copy()
                    frame['arm'], frame['detector'], frame['score'] = arm, detector, score
                    scores.append(frame)
        print('Completed:', task.cell, task.tf, flush=True)
    if len(folds) != 21 * 7 * 15 * 2 or len(fits) != 21 * 7 * 15:
        raise ValueError('task denominator changed')
    folds = pd.DataFrame(folds)
    summary = folds.groupby(['cell', 'tf', 'arm', 'detector']).agg(
        evaluable_chromosomes=('auroc', 'count'), mean_auroc=('auroc', 'mean'), mean_auprc=('auprc', 'mean'),
        canonical_fraction=('canonical', 'mean'), fixed_canonical_fraction=('fixed_canonical', 'mean')).reset_index()
    baseline = pd.read_csv(config['baseline_outputs']['folds']['path'], sep='\t')
    comparisons = []
    for detector, old_detector in [('unit_floor_fixed', 'normalized_fixed'), ('unit_floor_selected', 'normalized_selected')]:
        new = folds[folds.detector.eq(detector)]
        old = baseline[baseline.detector.eq(old_detector)]
        paired = new.merge(old, on=['cell', 'tf', 'arm', 'chromosome'], suffixes=('_new', '_old'), validate='one_to_one')
        if len(paired) != 21 * 7 * 15 or not (paired.positive_new == paired.positive_old).all() or not (paired.negative_new == paired.negative_old).all():
            raise ValueError('paired held-site denominator changed')
        paired['auroc_gain'] = paired.auroc_new - paired.auroc_old
        paired['auprc_gain'] = paired.auprc_new - paired.auprc_old
        comparisons.append(paired)
    outputs = {}
    for name, frame in [('folds', folds), ('fits', pd.DataFrame(fits)), ('summary', summary),
                        ('scores', pd.concat(scores, ignore_index=True)), ('curves', pd.concat(curves, ignore_index=True)),
                        ('paired', pd.concat(comparisons, ignore_index=True))]:
        path = out / (name + '.tsv')
        frame.to_csv(path, sep='\t', index=False)
        outputs[name] = record(path)
    result = dict(configuration_id=config['configuration_id'], outputs=outputs, training_only=True,
                  nested_matching=False, promotion_eligible=False, uncertainty_calibrated=False)
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
