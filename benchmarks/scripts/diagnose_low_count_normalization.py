#!/usr/bin/env python3
"""Frozen descriptive normalization audit on training and naked-DNA replicate 1."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rebuild_corrected_safety_profiles as extraction
import freeze_stable_geometry_models as models
import audit_direct_bias_safety as safety
import audit_control_blacklist as exclusions
from fixed_shape_normalization import detrend, normalize_fixed
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = models.ROOT
DATA = Path('/home/exouser/projects/fp-tools')
PROFILE_CONFIG = ROOT / 'benchmarks/results/corrected_safety_profiles_20260922_v1/configuration.freeze.json'
MODEL_DIR = ROOT / 'benchmarks/results/stable_geometry_models_20260922_v1'
SCREEN = ROOT / 'benchmarks/results/stable_geometry_screen_20260922_v2'
BAM = DATA / 'benchmarks/results/footprint_functional_v1/naked_dna_10m_s2026/NakedDNA_rep1.first10m.bam'
POSITIONS = np.arange(-100, 101)
ARMS = ['raw', 'parametric_direct']
PROTOCOL = dict(evidence='training_and_replicate1_descriptive_diagnostic',
                chromosomes=[f'chr{i}' for i in range(1, 16)], arms=ARMS,
                count_bins=[0, 1, 2, 5, 10, 20, 50, 100], scale_bins=[0, 1e-8, .01, .1, 1, 10],
                normalization='unchanged outer linear detrend, RMS, fallback 1, clip 12',
                score='unchanged serialized geometry and ten-decimal rounding',
                thresholds='existing frozen thresholds only; no calibration or fitting',
                final_control_read=False, new_holdout_labels_read=False, promotion_eligible=False)


def components(values):
    residual, scale = detrend(values, POSITIONS)
    denominator = np.where(scale > 1e-8, scale, 1.)
    before_clip = residual / denominator[:, None]
    return dict(scale=scale, denominator=denominator, fallback=scale <= 1e-8,
                amplification=1. / denominator, clipped_fraction=(np.abs(before_clip) > 12).mean(1),
                outer_variation=(np.abs(residual[:, abs(POSITIONS) >= 60]) > 1e-8).sum(1))


def count_stratum(counts):
    counts = np.asarray(counts, float)
    if not np.isfinite(counts).all() or (counts < 0).any():
        raise ValueError('finite nonnegative counts required')
    return np.searchsorted(PROTOCOL['count_bins'], counts, side='right') - 1


def freeze(out):
    if out.exists():
        raise ValueError('fresh output required')
    profile = extraction.extension.read_document(PROFILE_CONFIG, 'configuration_id')
    model_config = extraction.extension.read_document(MODEL_DIR / 'configuration.freeze.json', 'configuration_id')
    model_manifest = extraction.extension.read_document(MODEL_DIR / 'models.freeze.json', 'stage_id')
    thresholds = extraction.extension.read_document(SCREEN / 'thresholds.freeze.json', 'stage_id')
    screen = extraction.extension.read_document(SCREEN / 'configuration.freeze.json', 'configuration_id')
    direct = screen['prior']['direct']['source_configuration']['cells']
    items = [{k: item[k] for k in ['cell', 'manifest', 'arrays', 'sites']} for item in profile['items'] if item['role'] == 'naked']
    inputs = [record(p) for p in [PROFILE_CONFIG, MODEL_DIR / 'configuration.freeze.json', MODEL_DIR / 'models.freeze.json',
              SCREEN / 'thresholds.freeze.json', SCREEN / 'configuration.freeze.json', BAM, Path(str(BAM) + '.bai'), exclusions.BLACKLIST]]
    inputs += list(model_config['source_outputs'].values()) + list(model_manifest['outputs'].values()) + list(thresholds['outputs'].values())
    inputs += [item[k] for item in items for k in ['manifest', 'arrays', 'sites']]
    for entry in inputs:
        checked_path(entry, ROOT)
    if record(exclusions.BLACKLIST)['sha256'] != exclusions.BLACKLIST_SHA:
        raise ValueError('blacklist changed')
    for item in items:
        sites = pd.read_csv(item['sites']['path'], sep='\t')
        with np.load(item['arrays']['path'], allow_pickle=False) as archive:
            extraction.validate_sites(sites, archive, item['cell'], 'naked')
    paths = {entry['path'] for entry in profile['code'] + model_config['code']}
    paths.update([str(Path(__file__).resolve()), str(Path(exclusions.__file__).resolve())])
    config = dict(protocol=PROTOCOL, items=items, inputs=inputs, code=[record(Path(p)) for p in sorted(paths)],
                  runtime=extraction.extension.runtime(), bam=record(BAM), blacklist=record(exclusions.BLACKLIST),
                  training=model_config['source_outputs'], models=model_manifest['outputs']['detectors'],
                  thresholds=thresholds['outputs']['thresholds'],
                  dispersion={cell: direct[cell]['dispersion'] for cell in ['HepG2', 'K562']})
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / 'configuration.freeze.json', config)


def preflight(out):
    config = extraction.extension.read_document(out / 'configuration.freeze.json', 'configuration_id')
    if config['protocol'] != PROTOCOL or config['runtime'] != extraction.extension.runtime():
        raise ValueError('protocol/runtime changed')
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    return config


def extract(out, config):
    target = out / 'profiles_manifest.json'
    if target.exists():
        return extraction.base.verify_stage(target, config['configuration_id'], ROOT)
    if any(out.glob('*.rep1.npz')):
        raise ValueError('partial extraction preserved')
    outputs = {}
    for item in config['items']:
        sites = pd.read_csv(item['sites']['path'], sep='\t')
        mask = sites.TFBS_chr.isin(PROTOCOL['chromosomes']).to_numpy()
        with np.load(item['arrays']['path'], allow_pickle=False) as archive:
            cached = {k: archive[k][mask] for k in ['plus_log_bias', 'minus_log_bias', 'valid']}
        subset = sites.loc[mask].reset_index(drop=True)
        arrays = extraction.rebuild(subset, cached, checked_path(config['bam'], ROOT))
        arrays['source_index'] = np.flatnonzero(mask)
        path = out / (item['cell'] + '.rep1.npz')
        np.savez_compressed(path, **arrays)
        outputs[item['cell']] = record(path)
    result = dict(configuration_id=config['configuration_id'], outputs=outputs, replicate=1, candidate_fitted=False)
    result['stage_id'] = canonical_id(result)
    write_immutable(target, result)
    return result


def row_diagnostics(sites, panels, totals, valid, blacklist, fitted, thresholds, source, groups):
    rows = []
    for arm in ARMS:
        for (cell, tf), task in sites.groupby(['cell', 'tf'], sort=True):
            ix = task.index.to_numpy()
            model = fitted[(cell, tf, arm)]
            values = panels[arm][ix]
            parts = components(values)
            available = model['status'] != 'unavailable'
            score = models.score_model(model, values, totals[ix] == 0) if available else np.full(len(ix), np.nan)
            cutoff = thresholds.loc[(cell, tf, arm)]
            threshold = float(cutoff.candidate_threshold)
            eligible = valid[ix] & ~blacklist[ix] & (totals[ix] > 0) & available
            evaluated = eligible & bool(cutoff.threshold_eligible) & np.isfinite(threshold)
            table = task[safety.SITE_KEYS].copy()
            for name, array in dict(source=source, group=groups[ix], arm=arm, model_status=model['status'],
                                    total_cuts=totals[ix], count_stratum=count_stratum(totals[ix]),
                                    profile_valid=valid[ix], window_blacklisted=blacklist[ix],
                                    diagnostic_eligible=eligible, threshold_evaluated=evaluated,
                                    score=score, threshold=threshold,
                                    exceeds_frozen_threshold=evaluated & (score > threshold), **parts).items():
                table[name] = array
            rows.append(table)
    return pd.concat(rows, ignore_index=True)


def execute(out, verify=False):
    config = preflight(out)
    target = out / 'diagnostic_manifest.json'
    if target.exists():
        extraction.base.verify_stage(target, config['configuration_id'], ROOT)
        extraction.base.verify_stage(out / 'profiles_manifest.json', config['configuration_id'], ROOT)
        print('Verified frozen normalization diagnostic', flush=True)
        return
    if verify or (out / 'rows.tsv').exists():
        raise ValueError('missing or partial diagnostic preserved')
    profiles = extract(out, config)
    artifact = extraction.extension.read_document(Path(config['models']['path']), 'artifact_id')
    fitted = {(m['cell'], m['tf'], m['arm']): m for m in artifact['models']}
    thresholds = pd.read_csv(config['thresholds']['path'], sep='\t').set_index(['cell', 'tf', 'arm'])
    intervals = [(v[0], int(v[1]), int(v[2])) for v in
                 (line.split('\t') for line in Path(config['blacklist']['path']).read_text().splitlines())]
    sites = pd.read_csv(config['training']['sites']['path'], sep='\t')
    if not sites.TFBS_chr.isin(PROTOCOL['chromosomes']).all() or len(sites) != 8106:
        raise ValueError('training membership changed')
    with np.load(config['training']['panels']['path'], allow_pickle=False) as archive:
        panels = {arm: archive[arm] for arm in ARMS}
    totals = panels['raw'].sum(1)
    groups = np.where(sites.chip_label.to_numpy() == 1, 'training_positive', 'training_negative')
    frames = [row_diagnostics(sites, panels, totals, np.ones(len(sites), bool),
              exclusions.overlap_flags(sites, intervals, 100), fitted, thresholds, 'ATAC_training', groups)]
    membership = []
    for item in config['items']:
        original = pd.read_csv(item['sites']['path'], sep='\t')
        mask = original.TFBS_chr.isin(PROTOCOL['chromosomes']).to_numpy()
        subset = original.loc[mask].reset_index(drop=True)
        ledger = original[safety.SITE_KEYS].copy()
        ledger['training_chromosome'] = mask
        ledger['window_blacklisted'] = exclusions.overlap_flags(original, intervals, 100)
        membership.append(ledger)
        with np.load(profiles['outputs'][item['cell']]['path'], allow_pickle=False) as archive:
            np.testing.assert_array_equal(archive['source_index'], np.flatnonzero(mask))
            np.testing.assert_array_equal(archive['site_hash'], safety.site_hashes(subset))
            observed, expected, valid = archive['observed'], archive['expected'], archive['valid']
        _, _, residual = safety.direct_scores(observed, expected, config['dispersion'][item['cell']])
        panels = dict(raw=observed, parametric_direct=residual)
        frames.append(row_diagnostics(subset, panels, observed.sum(1), valid,
                      ledger.loc[mask, 'window_blacklisted'].to_numpy(), fitted, thresholds,
                      'naked_rep1', np.repeat('control', len(subset))))
    rows = pd.concat(frames, ignore_index=True)
    members = pd.concat(membership, ignore_index=True)
    if len(members) != 4200:
        raise ValueError('control denominator changed')
    keys = ['source', 'cell', 'tf', 'arm', 'group', 'count_stratum']
    summary = rows.groupby(keys).agg(sites=('score', 'size'), eligible=('diagnostic_eligible', 'sum'),
        threshold_evaluated=('threshold_evaluated', 'sum'), exceedances=('exceeds_frozen_threshold', 'sum'),
        fallback=('fallback', 'sum'), median_scale=('scale', 'median'), median_amplification=('amplification', 'median'),
        median_score=('score', 'median'), maximum_amplification=('amplification', 'max'),
        mean_clipped_fraction=('clipped_fraction', 'mean')).reset_index()
    outputs = {}
    for key, frame in [('rows', rows), ('summary', summary), ('control_membership', members)]:
        path = out / (key + '.tsv')
        frame.to_csv(path, sep='\t', index=False)
        outputs[key] = record(path)
    result = dict(configuration_id=config['configuration_id'], outputs=outputs, training_sites=8106,
                  control_sites=4200, control_training_sites=int(members.training_chromosome.sum()),
                  rows=len(rows), models_refitted=False, thresholds_changed=False, promotion_eligible=False)
    result['stage_id'] = canonical_id(result)
    write_immutable(target, result)
    print(json.dumps({k: v for k, v in result.items() if isinstance(v, (int, bool))}), flush=True)


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
        with extraction.extension.exclusive_run(out):
            execute(out, args.stage == 'verify')
