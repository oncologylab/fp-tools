#!/usr/bin/env python3
"""Frozen higher-count replicate-1 calibration with unchanged regional analysis."""
import argparse
import json
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_rep1_region_calibration as calibration
import verify_rep1_region_calibration as independent
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = calibration.ROOT
PRIOR = ROOT / 'benchmarks/results/rep1_region_calibration_20260922_v1'
AUDIT = ROOT / 'benchmarks/results/rep1_depth_input_audit_20260922_v1'
POLICY = dict(evidence='examined_replicate1_higher_count_development_calibration',
              depth_only_comparison=False, model_refitted=False, thresholds_changed=False,
              preprocessing='full prepared alignment: read-level blacklist and mitochondrial filtering; full-input deduplication',
              site_selection='unchanged 4200 training-chromosome controls, original row order',
              new_holdout_labels_read=False, final_control_read=False, promotion_eligible=False)


def freeze(out):
    if out.exists():
        raise ValueError('fresh output required')
    calibration.execute(PRIOR, verify=True)
    old = calibration.preflight(PRIOR)
    audit = json.loads((AUDIT / 'configuration.freeze.json').read_text())
    result = json.loads((AUDIT / 'audit.json').read_text())
    if result['configuration_id'] != audit['configuration_id']:
        raise ValueError('input audit identity mismatch')
    bam = next(r for r in audit['inputs'] if r['path'].endswith('/alignment/NakedDNA_rep1.filtered.bam'))
    if bam['sha256'] != 'a04cb6abc7bd06f00f0c813e31ebba774e24481d2d7075cba036e6e2638afd66':
        raise ValueError('audited replicate-1 alignment changed')
    total = 0
    for item in old['items']:
        sites = pd.read_csv(item['sites']['path'], sep='\t')
        if not sites.TFBS_chr.isin(calibration.CHROMS).all() or any('label' in c.lower() or 'chip' in c.lower() for c in sites.columns):
            raise ValueError('only label-free training controls allowed')
        total += len(sites)
    if total != 4200:
        raise ValueError('site denominator changed')
    config = dict(policy=POLICY, protocol=calibration.PROTOCOL,
                  items=old['items'], blacklist=old['blacklist'], bam=bam,
                  runtime=old['runtime'],
                  inputs=old['inputs'] + audit['inputs'] + [record(AUDIT / 'configuration.freeze.json'),
                    record(AUDIT / 'audit.json'), record(PRIOR / 'configuration.freeze.json')],
                  code=old['code'] + [record(Path(__file__)), record(Path(independent.__file__))])
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / 'configuration.freeze.json', config)


def execute(out, verify=False):
    config = calibration.source.extraction.extension.read_document(out / 'configuration.freeze.json', 'configuration_id')
    if config['policy'] != POLICY or config['protocol'] != calibration.PROTOCOL or config['runtime'] != calibration.source.extraction.extension.runtime():
        raise ValueError('frozen study changed')
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    if verify and not (out / 'profiles_manifest.json').exists():
        raise ValueError('missing profiles')
    # Reuse the previously verified extractor, not archived control counts.
    profiles = calibration.source.extract(out, config)
    destination = out / 'calibration'
    child_config = dict(protocol=config['protocol'], policy=config['policy'],
                        items=config['items'], profiles=profiles['outputs'],
                        blacklist=config['blacklist'], runtime=config['runtime'], code=config['code'],
                        inputs=config['inputs'] + [record(out / 'configuration.freeze.json'),
                          record(out / 'profiles_manifest.json')] + list(profiles['outputs'].values()))
    child_config['configuration_id'] = canonical_id(child_config)
    if destination.exists():
        actual = json.loads((destination / 'configuration.freeze.json').read_text())
        if actual != child_config:
            raise ValueError('derived calibration configuration changed')
    else:
        if verify:
            raise ValueError('missing calibration')
        destination.mkdir()
        write_immutable(destination / 'configuration.freeze.json', child_config)
    calibration.execute(destination, verify=verify)
    independent.verify(destination)


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
        with calibration.source.extraction.extension.exclusive_run(out):
            execute(out, args.stage == 'verify')
