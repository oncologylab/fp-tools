#!/usr/bin/env python3
"""Freeze full-training count geometries, then audit nominal tails on replicate 1."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crossfit_conditional_protection as study
import diagnose_low_count_normalization as controls
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = study.ROOT
SOURCE = ROOT / 'benchmarks/results/conditional_protection_20260922_v1'
CONTROL = ROOT / 'benchmarks/results/low_count_normalization_20260922_v1'
PROTOCOL = dict(evidence='previously_examined_replicate1_development_null_audit',
                nominal_alpha=.05, threshold=float(np.round(-np.log10(.05), 10)),
                calls='strict score greater than rounded threshold; require positive conditioned counts',
                geometry='unchanged full-training rule; freeze all 42 model records before control scoring',
                controls='all original 4200 label-free sites; retain blacklist and missing-model exclusions',
                safety_established=False, promotion_eligible=False, thresholds_fitted=False,
                final_control_read=False, new_holdout_labels_read=False)


def freeze(out):
    if out.exists():
        raise ValueError('fresh output required')
    study.execute(SOURCE, verify=True)
    controls.execute(CONTROL, verify=True)
    source = study.preflight(SOURCE)
    control = controls.preflight(CONTROL)
    result = json.loads((SOURCE / 'crossfit_manifest.json').read_text())
    profiles = json.loads((CONTROL / 'profiles_manifest.json').read_text())
    inputs = source['inputs'] + control['inputs'] + [record(SOURCE / 'configuration.freeze.json'),
              record(SOURCE / 'crossfit_manifest.json'), record(CONTROL / 'configuration.freeze.json'),
              record(CONTROL / 'profiles_manifest.json')] + list(result['outputs'].values()) + list(profiles['outputs'].values())
    config = dict(protocol=PROTOCOL, inputs=inputs,
                  code=source['code'] + control['code'] + [record(Path(__file__))],
                  runtime=source['runtime'], training=source['source_outputs'], grid=result['outputs']['grid'],
                  items=control['items'], control_profiles=profiles['outputs'], blacklist=control['blacklist'])
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / 'configuration.freeze.json', config)


def preflight(out):
    ext = study.prior.source.cohorts.prior.wide.old.parent.extension
    config = ext.read_document(out / 'configuration.freeze.json', 'configuration_id')
    if config['protocol'] != PROTOCOL or config['runtime'] != ext.runtime():
        raise ValueError('protocol/runtime changed')
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    return config


def fit(out, config):
    target = out / 'models.freeze.json'
    ext = study.prior.source.cohorts.prior.wide.old.parent.extension
    if target.exists():
        doc = ext.read_document(target, 'artifact_id')
        if doc['configuration_id'] != config['configuration_id']:
            raise ValueError('model configuration changed')
        return doc
    sites = pd.read_csv(config['training']['sites']['path'], sep='\t')
    if not sites.TFBS_chr.isin(study.PROTOCOL['chromosomes']).all():
        raise ValueError('nontraining labels')
    tasks = pd.read_csv(config['training']['metrics']['path'], sep='\t')[['cell', 'tf']].drop_duplicates()
    rows = []
    with np.load(config['grid']['path'], allow_pickle=False) as z:
        np.testing.assert_array_equal(z['site_hash'], study.safety.site_hashes(sites))
        for task in tasks.itertuples(index=False):
            ix = np.flatnonzero(sites.cell.eq(task.cell) & sites.tf.eq(task.tf))
            y = sites.iloc[ix].chip_label.to_numpy(int)
            for arm in study.ARMS:
                j, reason = study.select_geometry(z[arm + '_shape'][ix], z[arm][ix], y) if len(ix) else (0, 'unavailable')
                rows.append(dict(cell=task.cell, tf=task.tf, arm=arm, status=reason,
                                 geometry=list(study.prior.GEOMETRIES[j]) if len(ix) else None,
                                 positive=int((y == 1).sum()), negative=int((y == 0).sum())))
    if len(rows) != 42:
        raise ValueError('model denominator changed')
    doc = dict(configuration_id=config['configuration_id'], models=rows, promotion_eligible=False,
               threshold=PROTOCOL['threshold'], control_profiles_scored=False)
    doc['artifact_id'] = canonical_id(doc)
    write_immutable(target, doc)
    return doc


def execute(out, verify=False):
    config = preflight(out)
    base = study.prior.source.cohorts.prior.wide.old.parent.base
    target = out / 'screen_manifest.json'
    if target.exists():
        result = base.verify_stage(target, config['configuration_id'], ROOT)
        checked_path(result['model_artifact'], ROOT)
        print('Verified replicate-1 count-tail diagnostic', flush=True)
        return
    if verify or (out / 'sites.tsv').exists():
        raise ValueError('missing or partial screen preserved')
    model = fit(out, config)
    lock = dict(configuration_id=config['configuration_id'], models=record(out / 'models.freeze.json'),
                threshold=PROTOCOL['threshold'], promotion_eligible=False)
    lock['stage_id'] = canonical_id(lock)
    write_immutable(out / 'screen.freeze.json', lock)
    intervals = [(v[0], int(v[1]), int(v[2])) for v in
                 (line.split('\t') for line in Path(config['blacklist']['path']).read_text().splitlines())]
    rows, summaries = [], []
    for item in config['items']:
        cell = item['cell']
        sites = pd.read_csv(item['sites']['path'], sep='\t')
        if any('label' in c.lower() or 'chip' in c.lower() for c in sites.columns):
            raise ValueError('control labels prohibited')
        if not sites.TFBS_chr.isin(study.PROTOCOL['chromosomes']).all():
            raise ValueError('control scope changed')
        black = controls.exclusions.overlap_flags(sites, intervals, 100)
        with np.load(config['control_profiles'][cell]['path'], allow_pickle=False) as z:
            np.testing.assert_array_equal(z['source_index'], np.arange(len(sites)))
            np.testing.assert_array_equal(z['site_hash'], study.safety.site_hashes(sites))
            plus, minus, profile_valid = z['plus_observed'], z['minus_observed'], z['valid']
        with np.load(item['arrays']['path'], allow_pickle=False) as z:
            bias = [z['plus_log_bias'], z['minus_log_bias']]
        for m in [m for m in model['models'] if m['cell'] == cell]:
            ix = np.flatnonzero(sites.tf.eq(m['tf']))
            available = m['status'] != 'unavailable'
            first, second = [v[ix] for v in bias] if m['arm'] == 'conditional_parametric' else [np.zeros_like(plus[ix]), np.zeros_like(minus[ix])]
            geometry = tuple(m['geometry']) if available else study.prior.GEOMETRIES[0]
            score, comp = study.scoring.profile_scores(plus[ix], minus[ix], first, second, geometry)
            valid = profile_valid[ix] & comp['valid'] & ~black[ix] & available
            informative = valid & ((comp['n_plus'] + comp['n_minus']) > 0)
            calls = informative & (score > PROTOCOL['threshold'])
            frame = sites.iloc[ix][study.safety.SITE_KEYS].copy()
            for name, values in dict(arm=m['arm'], model_status=m['status'], score=score if available else np.nan,
                                     window_blacklisted=black[ix], eligible=valid, informative=informative,
                                     called=calls, n_plus=comp['n_plus'], n_minus=comp['n_minus'],
                                     central=comp['central'], p_plus=comp['p_plus'], p_minus=comp['p_minus']).items():
                frame[name] = values
            rows.append(frame)
            row = dict(cell=cell, tf=m['tf'], arm=m['arm'], model_status=m['status'], sites=len(ix),
                       blacklisted=int(black[ix].sum()), eligible=int(valid.sum()), informative=int(informative.sum()), calls=int(calls.sum()))
            for prefix, mask in [('all_', valid), ('informative_', informative)]:
                n = int(mask.sum())
                row[prefix + 'rate'] = calls.sum() / n if n else np.nan
                row[prefix + 'wilson_upper'] = study.safety.wilson_interval(int(calls.sum()), n)[1]
            row['nominal_fpr_bound_met'] = bool(available and row['all_wilson_upper'] <= .05 and row['informative_wilson_upper'] <= .05)
            summaries.append(row)
    table, summary = pd.concat(rows, ignore_index=True), pd.DataFrame(summaries)
    if len(table) != 8400 or len(summary) != 42:
        raise ValueError('control denominator changed')
    outputs = {}
    for name, frame in [('sites', table), ('summary', summary)]:
        path = out / (name + '.tsv')
        frame.to_csv(path, sep='\t', index=False)
        outputs[name] = record(path)
    result = dict(configuration_id=config['configuration_id'], model_artifact=record(out / 'models.freeze.json'),
                  screen_freeze=record(out / 'screen.freeze.json'), outputs=outputs,
                  promotion_eligible=False, paired_DWM_safety_established=False, final_control_read=False)
    result['stage_id'] = canonical_id(result)
    write_immutable(target, result)
    print(summary.to_string(index=False), flush=True)


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
        with study.prior.source.cohorts.prior.wide.old.parent.extension.exclusive_run(out):
            execute(out, args.stage == 'verify')
