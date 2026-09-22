#!/usr/bin/env python3
"""Fit all fixed-rule research geometries on chr1-15; serialize exact scorers."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crossfit_stable_geometry as source
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = source.ROOT
SOURCE = ROOT / 'benchmarks/results/stable_geometry_crossfit_20260922_v1'
PROTOCOL = dict(training_chromosomes=[f'chr{i}' for i in range(1, 16)],
    rule=source.PROTOCOL, task_count=21, arms=source.source.ARMS,
    selection='unchanged full-training per-task fitting rule; no selection between corrections',
    missing_cohort='unavailable, not fitted', underpowered='fixed-geometry fallback, explicitly marked',
    calibration_performed=False, safety_established=False, promotion_eligible=False)


def score_model(model, profiles, zero_read):
    if model['status'] == 'unavailable':
        raise ValueError('unavailable model has no scored cohort')
    geometry = tuple(model['geometry'])
    if geometry not in source.GEOMETRIES or model.get('normalization_fallback') != 1. or model.get('decimals') != 10:
        raise ValueError('model operator identity changed')
    values = source.normalize_fixed(profiles, source.POSITIONS, model['normalization_fallback'])
    return source.stable.stable_scores(values, source.geometry_weights(geometry), zero_read)


def freeze(out):
    if out.exists(): raise ValueError('fresh directory required')
    source.execute(SOURCE, verify=True)
    old = json.loads((SOURCE / 'configuration.freeze.json').read_text())
    result = json.loads((SOURCE / 'crossfit_manifest.json').read_text())
    config = dict(protocol=PROTOCOL, inputs=old['inputs'] + [record(SOURCE / 'configuration.freeze.json'),
        record(SOURCE / 'crossfit_manifest.json')] + list(result['outputs'].values()),
        source_outputs=old['source_outputs'], code=old['code'] + [record(Path(__file__))],
        runtime=source.source.cohorts.prior.wide.old.parent.extension.runtime())
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True); write_immutable(out / 'configuration.freeze.json', config)


def execute(out, verify=False):
    ext = source.source.cohorts.prior.wide.old.parent.extension
    config = ext.read_document(out / 'configuration.freeze.json', 'configuration_id')
    if canonical_id(config['protocol']) != canonical_id(PROTOCOL) or config['runtime'] != ext.runtime():
        raise ValueError('protocol or runtime changed')
    for entry in config['inputs'] + config['code']: checked_path(entry, ROOT)
    target = out / 'models.freeze.json'
    if target.exists():
        source.source.cohorts.prior.wide.old.parent.base.verify_stage(target, config['configuration_id'], ROOT)
        print('Verified full-training stable geometry artifacts', flush=True); return
    if verify or (out / 'detectors.json').exists(): raise ValueError('missing or partial output')
    paths = config['source_outputs']; read = lambda key: pd.read_csv(paths[key]['path'], sep='\t', float_precision='round_trip')
    sites, metrics = read('sites'), read('metrics')
    if not sites.TFBS_chr.isin(PROTOCOL['training_chromosomes']).all(): raise ValueError('nontraining coordinate')
    with np.load(paths['panels']['path'], allow_pickle=False) as archive: panels = {arm:archive[arm] for arm in source.source.ARMS}
    models = []
    for task in metrics[['cell','tf','matched_pairs','status']].drop_duplicates().itertuples(index=False):
        ix = np.flatnonzero(sites.cell.eq(task.cell) & sites.tf.eq(task.tf)); y = sites.iloc[ix].chip_label.to_numpy(int)
        zero = panels['raw'][ix].sum(1) == 0
        for arm in source.source.ARMS:
            geometry, reason = None, 'unavailable'
            if len(ix):
                norm = source.normalize_fixed(panels[arm][ix], source.POSITIONS, 1.)
                geometry, reason = source.fit_geometry(norm, y, zero, minimum=200)
            model = dict(cell=task.cell, tf=task.tf, arm=arm, geometry=list(geometry) if geometry else None,
                status=reason, positive=int((y==1).sum()), negative=int((y==0).sum()),
                cohort_status=task.status, normalization_fallback=1., decimals=10)
            if len(ix):
                expected = source.stable.stable_scores(norm, source.geometry_weights(geometry), zero)
                loaded = json.loads(json.dumps(model))
                np.testing.assert_array_equal(score_model(loaded, panels[arm][ix], zero), expected)
                # Inference must not depend on batch composition.
                for j in (0, len(ix)//2, len(ix)-1):
                    np.testing.assert_array_equal(score_model(loaded, panels[arm][ix[j:j+1]], zero[j:j+1]), expected[j:j+1])
            models.append(model)
    if len(models) != 147: raise ValueError('task denominator changed')
    artifact = dict(schema='fp-tools-research-stable-geometry-models-v1', configuration_id=config['configuration_id'],
        models=models, probabilities=False, safety_established=False, promotion_eligible=False)
    artifact['artifact_id'] = canonical_id(artifact)
    write_immutable(out / 'detectors.json', artifact)
    table = pd.DataFrame(models); table.geometry = table.geometry.map(json.dumps); table.to_csv(out / 'training.tsv', sep='\t', index=False)
    result = dict(configuration_id=config['configuration_id'], outputs={key:record(out/name) for key,name in
        [('detectors','detectors.json'), ('training','training.tsv')]}, training_only=True,
        calibration_performed=False, safety_established=False, promotion_eligible=False)
    result['stage_id'] = canonical_id(result); write_immutable(target, result)
    print('Frozen all 147 full-training detector records', flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--outdir',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group(); mode.add_argument('--freeze-only',action='store_true'); mode.add_argument('--verify-only',action='store_true')
    args=parser.parse_args(); out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'): parser.error('use research results')
    if args.freeze_only: freeze(out)
    else:
        with source.source.cohorts.prior.wide.old.parent.extension.exclusive_run(out): execute(out,args.verify_only)
