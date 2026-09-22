#!/usr/bin/env python3
"""Separate blacklist exclusion from expanded control output-layout effects."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_control_blacklist as prior
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = prior.ROOT
EXPANDED = ROOT / 'benchmarks/results/all_site_control_dwm_20260922_v1'
HISTORICAL = ROOT / 'benchmarks/results/control_blacklist_20260922_v1'
PROTOCOL = dict(evidence='examined_replicate2_blacklisted_layout_diagnostic',
                candidates_scored=False, promotion_eligible=False,
                preserve_all_sites=True, expectation_rtol=1e-5, expectation_atol=1e-7)


def expanded_blacklist_command(command, out, blacklist):
    result = list(command)
    if '--blacklist' in result:
        raise ValueError('source already has blacklist')
    for option, value in [('--outdir', str(out / 'tracks')), ('--prefix', 'NakedDNA_rep2_all_sites')]:
        if result.count(option) != 1:
            raise ValueError('unique output option required')
        result[result.index(option) + 1] = value
    return result + ['--blacklist', str(blacklist)]


def freeze(out):
    if out.exists():
        raise ValueError('fresh output required')
    prior.execute(HISTORICAL, verify=True)
    prior.prior.builder.execute(EXPANDED, verify=True)
    historical = prior.preflight(HISTORICAL)
    expanded = json.loads((EXPANDED / 'configuration.freeze.json').read_text())
    tracks = json.loads((HISTORICAL / 'tracks_manifest.json').read_text())
    review = json.loads((HISTORICAL / 'review_manifest.json').read_text())
    refs = {'historical_observed': tracks['outputs']['observed'],
            'historical_expected': tracks['outputs']['expected'],
            'historical_sites': review['outputs']['sites']}
    expanded_tracks = json.loads((EXPANDED / 'tracks_manifest.json').read_text())
    refs.update({'expanded_' + k: expanded_tracks['outputs'][k] for k in ['observed', 'expected']})
    config = dict(protocol=PROTOCOL, inputs=historical['inputs'] + expanded['inputs'] + list(refs.values()) +
                  [record(EXPANDED / 'configuration.freeze.json'), record(EXPANDED / 'tracks_manifest.json'),
                   record(HISTORICAL / 'tracks_manifest.json'), record(HISTORICAL / 'review_manifest.json')],
                  code=historical['code'] + [record(Path(__file__))], runtime=historical['runtime'],
                  original_model=expanded['original_model'], references=refs,
                  command=expanded_blacklist_command(expanded['command'], out, prior.BLACKLIST))
    for entry in config['inputs']:
        checked_path(entry, ROOT)
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / 'configuration.freeze.json', config)


def execute(out, verify=False):
    base = prior.prior.builder.audit.screen.prior.builder.base
    ext = prior.prior.builder.audit.screen.prior.builder.extension
    config = ext.read_document(out / 'configuration.freeze.json', 'configuration_id')
    if config['protocol'] != PROTOCOL or config['runtime'] != ext.runtime():
        raise ValueError('protocol/runtime changed')
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    manifest = out / 'tracks_manifest.json'
    if manifest.exists():
        tracks = base.verify_stage(manifest, config['configuration_id'], ROOT)
    else:
        if verify or (out / 'tracks').exists():
            raise ValueError('missing or partial tracks preserved')
        subprocess.run(config['command'], cwd=ROOT, check=True)
        prefix = out / 'tracks/NakedDNA_rep2_all_sites'
        old = prior.prior.builder.AtacBias().from_pickle(str(checked_path(config['original_model'], ROOT)))
        new = prior.prior.builder.AtacBias().from_pickle(str(prefix) + '_AtacBias.pickle')
        if new.stype != old.stype or new.correction_factor != 1.:
            raise ValueError('model settings changed')
        for strand in old.bias:
            for attr in ['counts', 'bg_counts']:
                np.testing.assert_array_equal(getattr(old.bias[strand], attr), getattr(new.bias[strand], attr))
        tracks = dict(configuration_id=config['configuration_id'], model_refitted=False,
                      outputs={k: record(Path(str(prefix) + suffix)) for k, suffix in
                               [('observed', '_uncorrected.bw'), ('expected', '_expected.bw'), ('model', '_AtacBias.pickle')]})
        tracks['stage_id'] = canonical_id(tracks)
        write_immutable(manifest, tracks)
    target = out / 'review_manifest.json'
    if target.exists():
        base.verify_stage(target, config['configuration_id'], ROOT)
        print('Verified blacklisted layout comparison', flush=True)
        return
    if verify or (out / 'sites.tsv').exists():
        raise ValueError('missing or partial review preserved')
    sites = pd.read_csv(checked_path(config['references']['historical_sites'], ROOT), sep='\t')
    safety = prior.prior.builder.audit.screen.prior.safety
    arrays = {}
    for name, entry in {**config['references'], **tracks['outputs']}.items():
        if name.endswith(('observed', 'expected')):
            arrays[name], valid = safety.extract_prediction_profiles(sites, checked_path(entry, ROOT), 100, require_dense=False)
            if not valid.all():
                raise ValueError('coordinate coverage changed')
    table = sites[safety.SITE_KEYS + ['motif_blacklisted', 'window_blacklisted']].copy()
    table['new_positive'] = arrays['expected'].sum(1) > 0
    for label, source in [('historical', 'historical_expected'), ('unblacklisted', 'expanded_expected')]:
        positive, match, error = prior.prior.compare_expectations(arrays[source], arrays['expected'])
        table[label + '_positive'] = positive
        table[label + '_shape_match'] = match
        table[label + '_max_error'] = error
    table['unblacklisted_observed_exact'] = (arrays['observed'] == arrays['expanded_observed']).all(1)
    if len(table) != 4200:
        raise ValueError('site denominator changed')
    outputs = {}
    for cell in ['HepG2', 'K562']:
        mask = sites.cell.to_numpy() == cell
        path = out / (cell + '.profiles.npz')
        np.savez_compressed(path, observed=arrays['observed'][mask], expected=arrays['expected'][mask],
                            site_hash=safety.site_hashes(sites.loc[mask]),
                            window_blacklisted=table.loc[mask, 'window_blacklisted'].to_numpy())
        outputs[cell + '_profiles'] = record(path)
    path = out / 'sites.tsv'
    table.to_csv(path, sep='\t', index=False)
    outputs['sites'] = record(path)
    eligible = ~table.window_blacklisted
    result = dict(configuration_id=config['configuration_id'], sites=len(table),
                  blacklisted=int((~eligible).sum()), new_positive=int(table.new_positive.sum()),
                  blacklisted_positive=int(table.loc[~eligible, 'new_positive'].sum()),
                  historical_shape_failures=int((~table.loc[table.historical_positive, 'historical_shape_match']).sum()),
                  historical_positive_lost=int((table.historical_positive & ~table.new_positive).sum()),
                  eligible_blacklist_shape_changes=int((eligible & table.unblacklisted_positive & ~table.unblacklisted_shape_match).sum()),
                  eligible_observed_changes=int((eligible & ~table.unblacklisted_observed_exact).sum()),
                  blacklist_observed_changes=int((~eligible & ~table.unblacklisted_observed_exact).sum()),
                  candidates_scored=False, promotion_eligible=False, qualified_for_rescreen=False, outputs=outputs)
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
        with prior.prior.builder.audit.screen.prior.builder.extension.exclusive_run(out):
            execute(out, args.stage == 'verify')
