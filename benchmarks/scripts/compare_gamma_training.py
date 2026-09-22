#!/usr/bin/env python3
"""Frozen training-only comparison; no candidate selection or promotion."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_corrected_shape_development as old
import gamma_poisson_global_expectation as gamma
from fixed_shape_normalization import normalize_fixed
from fp_tools.tools.parametric_bias import calibrated_residuals
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = old.ROOT
SOURCE = ROOT / 'benchmarks/results/global_parametric_training_profiles_20260922_v1'
OLD = ROOT / 'benchmarks/results/corrected_shape_development_20260922_v1/configuration.freeze.json'
KEYS = ['cell', 'tf', 'TFBS_chr', 'TFBS_start', 'TFBS_end', 'TFBS_strand']


def load_tracks(archive, order):
    """Decompress each archive member once, preserving requested site order."""
    return {key:archive[key][order] for key in
            ('plus_observed','minus_observed','plus_log_bias','minus_log_bias')}


def align_indices(target, source):
    left = pd.MultiIndex.from_frame(target[KEYS])
    right = pd.MultiIndex.from_frame(source[KEYS])
    if not right.is_unique or not left.is_unique:
        raise ValueError('duplicate site identity')
    index = right.get_indexer(left)
    if (index < 0).any():
        raise ValueError('missing training site')
    return index


def complete_pairs(pairs, supported):
    if len(pairs) != len(supported):
        raise ValueError('pair support length mismatch')
    frame = pairs[['cell', 'tf', 'pair_id', 'chip_label']].copy()
    frame['supported'] = np.asarray(supported, bool)
    grouped = frame.groupby(['cell', 'tf', 'pair_id'], sort=False)
    if not (grouped.size().eq(2).all() and grouped.chip_label.sum().eq(1).all()):
        raise ValueError('each pair needs one positive and one negative')
    return grouped.supported.transform('all').to_numpy()


def freeze(out):
    if out.exists():
        raise ValueError('fresh output directory required')
    previous = json.loads(OLD.read_text())
    inputs = previous['inputs'] + [record(OLD), previous['cohort_manifest'], previous['cohort_configuration']]
    for cell in ('HepG2', 'K562'):
        stage = SOURCE / (cell + '.json')
        inputs += [record(stage)] + list(json.loads(stage.read_text())['outputs'].values())
    paths = [Path(__file__), Path(gamma.__file__), ROOT/'benchmarks/scripts/fixed_shape_normalization.py',
             ROOT/'benchmarks/scripts/global_parametric_scoring.py', ROOT/'src/fp_tools/tools/parametric_bias.py']
    config = dict(inputs=inputs, code=previous['code']+[record(p) for p in paths],
        runtime=old.base.source.builder.extension.runtime(), training_chromosomes=[f'chr{i}' for i in range(1,16)],
        geometry=[7,12,40], normalization_fallback=1., model_selection=False,
        evidence='training_only_diagnostic', promotion_eligible=False,
        source_label_scope='existing chr1-18 source tables; only chr1-15 selected or evaluated',
        pairing='existing balanced v3 pairs; retain both members only when all arms have complete profile support',
        arms=['raw','DWM_expected','parametric_direct','gamma_fixed_dispersion','gamma_predictive'],
        uncertainty_calibrated=False, score_semantics='ranks, not binding probabilities')
    for item in inputs + config['code']:
        checked_path(item, ROOT)
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True)
    write_immutable(out/'configuration.freeze.json', config)


def execute(out, verify=False):
    config = old.base.source.builder.extension.read_document(out/'configuration.freeze.json','configuration_id')
    if config['runtime'] != old.base.source.builder.extension.runtime():
        raise ValueError('runtime changed')
    for item in config['inputs'] + config['code']:
        checked_path(item, ROOT)
    manifest = out/'comparison_manifest.json'
    if manifest.exists():
        old.base.source.builder.base.verify_stage(manifest,config['configuration_id'],ROOT)
        print('Verified completed training comparison', flush=True)
        return
    if verify or list(out.glob('*.tsv')):
        raise ValueError('missing result or partial output preserved')
    previous = json.loads(OLD.read_text())
    sites, panels, ledger = old.load_split(previous,ROOT,'train')
    if not sites.TFBS_chr.isin(config['training_chromosomes']).all():
        raise ValueError('nontraining chromosome')
    cohort = json.loads(checked_path(previous['cohort_manifest'],ROOT).read_text())
    pairs = []
    new = {key:np.full((len(sites),201),np.nan) for key in ('gamma_fixed_dispersion','gamma_predictive')}
    for cell in ('HepG2','K562'):
        pair = pd.read_csv(checked_path(cohort['outputs'][cell+'_pairs'],ROOT),sep='\t')
        pair = pair.loc[pair.split.eq('train')].copy()
        pairs.append(pair)
        ix = np.flatnonzero(sites.cell.eq(cell))
        np.testing.assert_array_equal(pair.chip_label,sites.iloc[ix].chip_label)
        source = pd.read_csv(SOURCE/(cell+'.sites.tsv.gz'),sep='\t')
        order = align_indices(sites.iloc[ix],source)
        with np.load(SOURCE/(cell+'.npz'),allow_pickle=False) as data:
            # NPZ indexing decompresses a whole member. Load each member once,
            # then select the fixed cohort before iterating over site windows.
            tracks = load_tracks(data, order)
            for row, dest in enumerate(ix):
                result = gamma.score_tracks(tracks['plus_observed'][row],tracks['minus_observed'][row],
                    tracks['plus_log_bias'][row],tracks['minus_log_bias'][row])
                observed = result['observed'][100:301]
                mean = result['expected'][100:301]
                predictive = result['corrected'][100:301]
                reverse = sites.iloc[dest].TFBS_strand == '-'
                np.testing.assert_array_equal(observed[::-1] if reverse else observed,panels['raw'][dest])
                if not np.isfinite(mean).all():
                    continue
                fixed = calibrated_residuals(observed,mean,'deviance',dispersion=previous['dispersion'][cell])
                new['gamma_fixed_dispersion'][dest] = fixed[::-1] if reverse else fixed
                new['gamma_predictive'][dest] = predictive[::-1] if reverse else predictive
        print(f'Computed {cell} training expectations',flush=True)
    panels.update(new)
    pair = pd.concat(pairs,ignore_index=True)
    support = np.logical_and.reduce([np.isfinite(p).all(axis=1) for p in panels.values()])
    keep = complete_pairs(pair,support)
    rows, curves, score_rows = [], [], []
    weights = old.candidates.geometry_weights(config['geometry'])
    for task in ledger.itertuples(index=False):
        ix = np.flatnonzero((sites.cell.eq(task.cell)&sites.tf.eq(task.tf)).to_numpy() & keep)
        y = sites.iloc[ix].chip_label.to_numpy(int)
        for arm, values in panels.items():
            row = dict(cell=task.cell,tf=task.tf,arm=arm,original_pairs=task.matched_pairs,
                retained_pairs=len(ix)//2,cohort_status=task.status,
                status='underpowered' if len(ix)<400 else 'training_diagnostic')
            if len(np.unique(y))==2:
                score = values[ix] @ weights
                score[panels['raw'][ix].sum(axis=1)==0] = 0.
                norm = normalize_fixed(values[ix],old.candidates.POSITIONS,config['normalization_fallback'])
                delta = norm[y==1].mean(axis=0)-norm[y==0].mean(axis=0)
                x = old.candidates.POSITIONS
                row.update(auroc=roc_auc_score(y,score),auprc=average_precision_score(y,score),
                    center_contrast=delta[abs(x)<=5].mean(),left_shoulder=delta[(x>=-40)&(x<=-15)].mean(),
                    right_shoulder=delta[(x>=15)&(x<=40)].mean(),canonical_protection=old.candidates.protection_direction(delta))
                curves.append(pd.DataFrame(dict(cell=task.cell,tf=task.tf,arm=arm,position=x,
                    positive_mean=norm[y==1].mean(axis=0),negative_mean=norm[y==0].mean(axis=0))))
                scores = sites.iloc[ix][KEYS+['chip_label']].copy()
                scores['arm'],scores['score']=arm,score
                score_rows.append(scores)
            rows.append(row)
    if len(rows)!=105:
        raise ValueError('full 21-task denominator changed')
    outputs = {}
    for name,table in [('metrics',pd.DataFrame(rows)),('curves',pd.concat(curves)),('scores',pd.concat(score_rows))]:
        path=out/(name+'.tsv');table.to_csv(path,sep='\t',index=False);outputs[name]=record(path)
    result=dict(configuration_id=config['configuration_id'],outputs=outputs,promotion_eligible=False,
        evidence='training_only_diagnostic',uncertainty_calibrated=False,candidate_selected=False)
    result['stage_id']=canonical_id(result)
    write_immutable(manifest,result)
    print('Completed 105 training comparison rows',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir',type=Path,required=True)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--freeze-only',action='store_true')
    mode.add_argument('--verify-only',action='store_true')
    args=parser.parse_args();out=args.outdir.resolve()
    if not out.is_relative_to(ROOT/'benchmarks/results'):
        parser.error('use research results')
    if args.freeze_only: freeze(out)
    else:
        with old.base.source.builder.extension.exclusive_run(out): execute(out,args.verify_only)
