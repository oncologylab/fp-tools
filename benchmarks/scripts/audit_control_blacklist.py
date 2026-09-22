#!/usr/bin/env python3
"""Frozen blacklist intervention on the historical naked-DNA region layout."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import isolate_control_region_layout as prior
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = prior.ROOT
SOURCE = ROOT / 'benchmarks/results/control_region_layout_20260922_v1'
BLACKLIST = prior.builder.DATA / 'data/public/processed/encode_k562_hepg2_atac_replicates/reference/hg38-blacklist.v2.bed'
BLACKLIST_SHA = '31c69342df43bbc19dd8ef2886611a8150cb53e90f341de14d30e742c9251737'
PROTOCOL = dict(evidence='examined_replicate2_blacklist_intervention', candidates_scored=False,
                promotion_eligible=False, preserve_all_sites=True, flank=100)


def blacklist_command(command, out, blacklist):
    result = prior.historical_command(command, out)
    if '--blacklist' in result:
        raise ValueError('source already specifies a blacklist')
    return result + ['--blacklist', str(blacklist)]


def overlap_flags(sites, intervals, flank=0):
    result = []
    for site in sites.itertuples():
        start, end = int(site.TFBS_start), int(site.TFBS_end)
        if flank:
            center = (start + end) // 2
            start, end = center - flank, center + flank + 1
        result.append(any(chrom == site.TFBS_chr and left < end and right > start
                          for chrom, left, right in intervals))
    return np.asarray(result, dtype=bool)


def freeze(out):
    if out.exists():
        raise ValueError('fresh output required')
    prior.review(SOURCE, verify=True)
    old = prior.preflight(SOURCE)
    blacklist = record(BLACKLIST)
    if blacklist['sha256'] != BLACKLIST_SHA:
        raise ValueError('blacklist identity changed')
    config = dict(protocol=PROTOCOL, inputs=old['inputs'] + [blacklist,
                  record(SOURCE / 'configuration.freeze.json'), record(SOURCE / 'layout_manifest.json')],
                  code=old['code'] + [record(Path(__file__))], runtime=old['runtime'],
                  original_model=old['original_model'], blacklist=blacklist,
                  command=blacklist_command(old['command'], out, BLACKLIST))
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / 'configuration.freeze.json', config)


def preflight(out):
    ext = prior.builder.audit.screen.prior.builder.extension
    config = ext.read_document(out / 'configuration.freeze.json', 'configuration_id')
    if config['protocol'] != PROTOCOL or config['runtime'] != ext.runtime():
        raise ValueError('protocol/runtime changed')
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    return config


def execute(out, verify=False):
    config = preflight(out)
    base = prior.builder.audit.screen.prior.builder.base
    manifest = out / 'tracks_manifest.json'
    if manifest.exists():
        tracks = base.verify_stage(manifest, config['configuration_id'], ROOT)
    else:
        if verify or (out / 'tracks').exists():
            raise ValueError('missing or partial construction preserved')
        subprocess.run(config['command'], cwd=ROOT, check=True)
        prefix = out / 'tracks/NakedDNA_rep2_historical'
        old = prior.builder.AtacBias().from_pickle(str(checked_path(config['original_model'], ROOT)))
        new = prior.builder.AtacBias().from_pickle(str(prefix) + '_AtacBias.pickle')
        if new.stype != old.stype or new.correction_factor != 1.:
            raise ValueError('model settings changed')
        for strand in old.bias:
            for attr in ['counts', 'bg_counts']:
                np.testing.assert_array_equal(getattr(old.bias[strand], attr), getattr(new.bias[strand], attr))
        tracks = dict(configuration_id=config['configuration_id'], model_refitted=False,
                      outputs={k: record(Path(str(prefix) + suffix)) for k, suffix in
                               [('observed', '_uncorrected.bw'), ('expected', '_expected.bw'),
                                ('model', '_AtacBias.pickle')]})
        tracks['stage_id'] = canonical_id(tracks)
        write_immutable(manifest, tracks)
    target = out / 'review_manifest.json'
    if target.exists():
        base.verify_stage(target, config['configuration_id'], ROOT)
        print('Verified blacklist diagnostic', flush=True)
        return
    if verify or (out / 'sites.tsv').exists():
        raise ValueError('missing or partial review preserved')
    intervals = []
    for line in checked_path(config['blacklist'], ROOT).read_text().splitlines():
        chrom, start, end, *_ = line.split('\t')
        intervals.append((chrom, int(start), int(end)))
    original = json.loads((prior.builder.SOURCE / 'configuration.freeze.json').read_text())
    safety = prior.builder.audit.screen.prior.safety
    tables = []
    for cell in ['HepG2', 'K562']:
        entry = original['screen']['prior']['direct']['source_configuration']['cells'][cell]['naked_dwm']
        sites, old_observed, old_expected, _, _ = safety.load_artifact(entry, prior.builder.DATA, cell, control=True)
        observed, ov = safety.extract_prediction_profiles(sites, checked_path(tracks['outputs']['observed'], ROOT), 100, require_dense=False)
        expected, ev = safety.extract_prediction_profiles(sites, checked_path(tracks['outputs']['expected'], ROOT), 100, require_dense=False)
        positive, match, error = prior.compare_expectations(old_expected, expected)
        op, om, _ = prior.compare_expectations(old_observed, observed)
        om |= ~op & (observed.sum(1) == 0)
        table = sites[safety.SITE_KEYS].copy()
        for name, values in dict(old_positive=positive, new_positive=expected.sum(1) > 0,
                                 expected_shape_match=match, expected_max_error=error,
                                 observed_shape_match=om, coordinate_valid=ov & ev,
                                 motif_blacklisted=overlap_flags(sites, intervals),
                                 window_blacklisted=overlap_flags(sites, intervals, 100)).items():
            table[name] = values
        tables.append(table)
    data = pd.concat(tables, ignore_index=True)
    if len(data) != 4200:
        raise ValueError('site denominator changed')
    path = out / 'sites.tsv'
    data.to_csv(path, sep='\t', index=False)
    result = dict(configuration_id=config['configuration_id'], sites=len(data),
                  old_positive=int(data.old_positive.sum()), new_positive=int(data.new_positive.sum()),
                  support_changes=int((data.old_positive != data.new_positive).sum()),
                  expected_shape_failures=int((~data.loc[data.old_positive, 'expected_shape_match']).sum()),
                  observed_shape_failures=int((~data.observed_shape_match).sum()),
                  motif_blacklisted=int(data.motif_blacklisted.sum()), window_blacklisted=int(data.window_blacklisted.sum()),
                  candidates_scored=False, promotion_eligible=False, outputs={'sites': record(path)})
    result['stage_id'] = canonical_id(result)
    write_immutable(target, result)
    print(json.dumps(result), flush=True)


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
        with prior.builder.audit.screen.prior.builder.extension.exclusive_run(out):
            execute(out, args.stage == 'verify')
