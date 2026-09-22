#!/usr/bin/env python3
"""Frozen descriptive chromosome-direction ledger; not confidence intervals."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crossfit_stable_geometry as source
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id, checked_path

ROOT = source.ROOT
SOURCE = ROOT / 'benchmarks/results/stable_geometry_crossfit_20260922_v1'
EFFECTS = {
    'normalization': ('normalized_fixed', None, 'raw_fixed'),
    'detector': ('normalized_selected', None, 'normalized_fixed'),
    'correction': ('normalized_selected', 'DWM_expected', 'normalized_selected'),
    'combined': ('normalized_selected', 'DWM_expected', 'raw_fixed'),
}
PROTOCOL = dict(effects=EFFECTS, chromosomes=[f'chr{i}' for i in range(1, 16)],
    fixed_center=[-5, 5], left=[-40, -15], right=[15, 40], shape_margin=1e-10,
    gain_margin=1e-12, summaries='equal-chromosome means, signs, minimum and maximum, and leave-one-chromosome-out mean ranges',
    intervals=False, training_only=True, multiplicity_adjusted=False,
    nested_matching=False, promotion_eligible=False,
    plots=['DWM_expected', 'parametric_direct', 'radius_500'])


def direction_summary(gains, protection):
    g, p = np.asarray(gains, float), np.asarray(protection, bool)
    if g.shape != p.shape or g.ndim != 1:
        raise ValueError('aligned one-dimensional gains and protection required')
    good = np.isfinite(g); g, p = g[good], p[good]
    row = dict(evaluable=int(len(g)), positive=int((g > 1e-12).sum()),
        negative=int((g < -1e-12).sum()), tied=int((abs(g) <= 1e-12).sum()),
        fixed_protection=int(p.sum()), positive_and_protective=int(((g > 1e-12) & p).sum()))
    omitted = (g.sum() - g) / (len(g) - 1) if len(g) > 1 else np.array([])
    row.update(mean=float(g.mean()) if len(g) else np.nan,
        minimum=float(g.min()) if len(g) else np.nan, maximum=float(g.max()) if len(g) else np.nan,
        omitted_mean_min=float(omitted.min()) if len(omitted) else np.nan,
        omitted_mean_max=float(omitted.max()) if len(omitted) else np.nan)
    return row


def freeze(out):
    if out.exists():
        raise ValueError('fresh directory required')
    source.execute(SOURCE, verify=True)
    old = json.loads((SOURCE / 'configuration.freeze.json').read_text())
    result = json.loads((SOURCE / 'crossfit_manifest.json').read_text())
    config = dict(protocol=PROTOCOL, inputs=old['inputs'] + [record(SOURCE / 'configuration.freeze.json'),
        record(SOURCE / 'crossfit_manifest.json')] + list(result['outputs'].values()),
        source_outputs=old['source_outputs'], fitted_outputs=result['outputs'],
        code=old['code'] + [record(Path(__file__))], runtime=source.source.cohorts.prior.wide.old.parent.extension.runtime())
    config['configuration_id'] = canonical_id(config)
    out.mkdir(parents=True); write_immutable(out / 'configuration.freeze.json', config)


def execute(out, verify=False):
    extension = source.source.cohorts.prior.wide.old.parent.extension
    config = extension.read_document(out / 'configuration.freeze.json', 'configuration_id')
    if canonical_id(config['protocol']) != canonical_id(PROTOCOL) or config['runtime'] != extension.runtime():
        raise ValueError('protocol or runtime changed')
    for entry in config['inputs'] + config['code']:
        checked_path(entry, ROOT)
    target = out / 'direction_manifest.json'
    if target.exists():
        source.source.cohorts.prior.wide.old.parent.base.verify_stage(target, config['configuration_id'], ROOT)
        print('Verified chromosome-direction ledger', flush=True); return
    if verify or list(out.glob('*.tsv')):
        raise ValueError('missing or partial output')
    read = lambda entry: pd.read_csv(entry['path'], sep='\t', float_precision='round_trip')
    folds = read(config['fitted_outputs']['folds']); sites = read(config['source_outputs']['sites'])
    if not sites.TFBS_chr.isin(PROTOCOL['chromosomes']).all():
        raise ValueError('nontraining coordinate')
    lookup = folds.set_index(['cell', 'tf', 'arm', 'chromosome', 'detector'])
    tasks = list(folds[['cell', 'tf']].drop_duplicates().itertuples(index=False, name=None))
    with np.load(config['source_outputs']['panels']['path'], allow_pickle=False) as archive:
        panels = {arm: archive[arm] for arm in source.source.ARMS}
    rows = []; x = np.arange(-100, 101)
    for cell, tf in tasks:
        ix = np.flatnonzero(sites.cell.eq(cell) & sites.tf.eq(tf)); table = sites.iloc[ix]
        for arm in source.source.ARMS:
            values = source.normalize_fixed(panels[arm][ix], x, 1.) if len(ix) else np.empty((0, 201))
            for chrom in PROTOCOL['chromosomes']:
                held = table.TFBS_chr.eq(chrom).to_numpy(); y = table.loc[held, 'chip_label'].to_numpy(int)
                available = bool((y == 0).any() and (y == 1).any())
                center = left = right = np.nan; protective = False
                if available:
                    delta = values[held][y == 1].mean(0) - values[held][y == 0].mean(0)
                    center, left, right = delta[abs(x) <= 5].mean(), delta[(x >= -40) & (x <= -15)].mean(), delta[(x >= 15) & (x <= 40)].mean()
                    protective = bool(center < -1e-10 and left - center > 1e-10 and right - center > 1e-10)
                for effect, (candidate, reference_arm, reference) in EFFECTS.items():
                    a = lookup.loc[(cell, tf, arm, chrom, candidate)]
                    b = lookup.loc[(cell, tf, reference_arm or arm, chrom, reference)]
                    assert a.positive == b.positive == (y == 1).sum() and a.negative == b.negative == (y == 0).sum()
                    assert (a.evaluation_status == 'available') == (b.evaluation_status == 'available') == available
                    rows.append(dict(cell=cell, tf=tf, arm=arm, chromosome=chrom, effect=effect,
                        cohort_status=a.cohort_status, matched_pairs=a.matched_pairs, available=available,
                        auroc_gain=a.auroc-b.auroc, relative_auprc_gain=a.auprc/b.auprc-1 if available else np.nan,
                        fixed_center=center, fixed_left=left, fixed_right=right, fixed_protection=protective,
                        detector_protection=bool(a.canonical) if available else False))
    frame = pd.DataFrame(rows)
    if len(frame) != 8820:
        raise ValueError('task denominator changed')
    summaries = []
    for (cell, tf, arm, effect), group in frame.groupby(['cell', 'tf', 'arm', 'effect'], sort=True):
        row = dict(cell=cell, tf=tf, arm=arm, effect=effect, cohort_status=group.cohort_status.iloc[0], matched_pairs=int(group.matched_pairs.iloc[0]))
        for metric in ['auroc_gain', 'relative_auprc_gain']:
            row.update({metric+'_'+k:v for k,v in direction_summary(group[metric], group.fixed_protection).items()})
        summaries.append(row)
    outputs = {}
    for name, data in [('folds', frame), ('summary', pd.DataFrame(summaries))]:
        path = out / (name+'.tsv'); data.to_csv(path, sep='\t', index=False); outputs[name] = record(path)
    pdf = out / 'all_task_detector_direction.pdf'
    with PdfPages(pdf) as pages:
        for arm in PROTOCOL['plots']:
            fig, ax = plt.subplots(figsize=(11.7, 10), constrained_layout=True)
            matrix = np.full((21, 15), np.nan)
            for i, (cell, tf) in enumerate(tasks):
                group = frame.loc[frame.cell.eq(cell) & frame.tf.eq(tf) & frame.arm.eq(arm) & frame.effect.eq('detector')].set_index('chromosome').loc[PROTOCOL['chromosomes']]
                matrix[i] = group.auroc_gain
                for j, r in enumerate(group.itertuples()):
                    if r.available and r.fixed_protection:
                        ax.plot(j, i, 'o', color='black', markersize=3)
            bound = max(.05, float(np.nanmax(abs(matrix))))
            im = ax.imshow(np.ma.masked_invalid(matrix), cmap='RdBu', vmin=-bound, vmax=bound, aspect='auto')
            ax.set_xticks(range(15), range(1, 16)); ax.set_yticks(range(21), [f'{c} {t}' for c,t in tasks])
            ax.set_xlabel('Held training chromosome'); ax.set_title(f'{arm}: selected minus fixed normalized detector\nBlack dot: fixed-center protection; white: unavailable. Training diagnostic only.')
            fig.colorbar(im, ax=ax, label='Within-chromosome AUROC gain')
            pages.savefig(fig); path = out / (arm+'.png'); fig.savefig(path, dpi=110); plt.close(fig); outputs[arm+'_plot'] = record(path)
    outputs['pdf'] = record(pdf)
    result = dict(configuration_id=config['configuration_id'], outputs=outputs, training_only=True,
        confidence_intervals=False, promotion_eligible=False, visual_review_complete=False)
    result['stage_id'] = canonical_id(result); write_immutable(target, result)
    print('Completed all 8820 paired rows and 588 task/effect summaries', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--outdir', type=Path, required=True)
    mode = parser.add_mutually_exclusive_group(); mode.add_argument('--freeze-only', action='store_true'); mode.add_argument('--verify-only', action='store_true')
    args = parser.parse_args(); out = args.outdir.resolve()
    if not out.is_relative_to(ROOT / 'benchmarks/results'):
        parser.error('use research results')
    if args.freeze_only: freeze(out)
    else:
        with source.source.cohorts.prior.wide.old.parent.extension.exclusive_run(out): execute(out, args.verify_only)
