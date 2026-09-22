#!/usr/bin/env python3
"""Apply the existing frozen count detector to higher-count development controls."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import screen_conditional_protection_rep1 as screen
import verify_conditional_protection_rep1 as independent
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = screen.ROOT
OLD = ROOT / 'benchmarks/results/conditional_protection_rep1_20260922_v1'
FULL = ROOT / 'benchmarks/results/full_rep1_region_calibration_20260922_v1'
POLICY = dict(evidence='examined_replicate1_higher_count_development_screen',
              models='exact existing 42 records, no fitting', threshold=screen.PROTOCOL['threshold'],
              depth_only_comparison=False, independent_safety=False,
              promotion_eligible=False, final_control_read=False, new_holdout_labels_read=False)


def freeze(out):
    if out.exists():
        raise ValueError('fresh output required')
    screen.execute(OLD, verify=True)
    old = screen.preflight(OLD)
    full = screen.controls.extraction.extension.read_document(FULL / 'configuration.freeze.json', 'configuration_id')
    models = record(OLD / 'models.freeze.json')
    if models['sha256'] != '4e1aeded40b4bf84629c1b9ce0382c682642554ddc2c5181d3fd0a73610a7592':
        raise ValueError('frozen models changed')
    if old['items'] != full['items'] or old['blacklist'] != full['blacklist']:
        raise ValueError('site or exclusion definitions changed')
    config = dict(policy=POLICY, source_config=old, profile_configuration_id=full['configuration_id'],
                  models=models, inputs=old['inputs'] + full['inputs'] + [models,
                    record(OLD / 'configuration.freeze.json'), record(FULL / 'configuration.freeze.json')],
                  code=old['code'] + full['code'] + [record(Path(__file__)), record(Path(independent.__file__))],
                  runtime=old['runtime'])
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out / 'configuration.freeze.json', config)


def execute(out, verify=False):
    ext = screen.controls.extraction.extension
    config = ext.read_document(out / 'configuration.freeze.json', 'configuration_id')
    if config['policy'] != POLICY or config['runtime'] != ext.runtime():
        raise ValueError('frozen screen changed')
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    profiles = screen.controls.extraction.base.verify_stage(
        FULL / 'profiles_manifest.json', config['profile_configuration_id'], ROOT)
    old = config['source_config']
    child = dict(old, control_profiles=profiles['outputs'], policy=POLICY,
                 inputs=config['inputs'] + [record(out / 'configuration.freeze.json'),
                    record(FULL / 'profiles_manifest.json')] + list(profiles['outputs'].values()),
                 code=config['code'])
    child.pop('configuration_id')
    child['configuration_id'] = canonical_id(child)
    target = out / 'screen'
    if not target.exists():
        if verify:
            raise ValueError('missing screen')
        target.mkdir()
        write_immutable(target / 'configuration.freeze.json', child)
        original = ext.read_document(checked_path(config['models'], ROOT), 'artifact_id')
        artifact = dict(configuration_id=child['configuration_id'], models=original['models'],
                        promotion_eligible=False, threshold=POLICY['threshold'],
                        control_profiles_scored=False, source_models=config['models'])
        artifact['artifact_id'] = canonical_id(artifact)
        write_immutable(target / 'models.freeze.json', artifact)
    else:
        if json.loads((target / 'configuration.freeze.json').read_text()) != child:
            raise ValueError('derived screen changed')
    frozen = ext.read_document(target / 'models.freeze.json', 'artifact_id')
    original = ext.read_document(checked_path(config['models'], ROOT), 'artifact_id')
    if frozen['models'] != original['models'] or frozen['threshold'] != POLICY['threshold']:
        raise ValueError('model records or threshold changed')
    screen.execute(target, verify=verify)
    independent.verify(target)


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
        with screen.controls.extraction.extension.exclusive_run(out):
            execute(out, args.stage == 'verify')
