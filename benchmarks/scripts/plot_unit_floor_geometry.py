#!/usr/bin/env python3
"""Plot all unit-floor training contrasts without claiming uncertainty coverage."""
import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crossfit_unit_floor_geometry as study
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import canonical_id


def plot(source, out):
    if out.exists():
        raise ValueError('fresh plot directory required')
    study.execute(source, verify=True)
    manifest = json.loads((source / 'crossfit_manifest.json').read_text())
    curves = pd.read_csv(manifest['outputs']['curves']['path'], sep='\t')
    summary = pd.read_csv(manifest['outputs']['summary']['path'], sep='\t')
    tasks = summary[['cell', 'tf']].drop_duplicates()
    fig, axes = plt.subplots(7, 3, figsize=(15, 21), constrained_layout=True)
    arms = [('DWM_expected', 'Dinucleotide weight matrix', '#555555'),
            ('parametric_direct', 'Conditional parametric', '#1874bd'),
            ('radius_500', '500-base background', '#a83c83')]
    for ax, task in zip(axes.flat, tasks.itertuples(index=False)):
        sub = curves[curves.cell.eq(task.cell) & curves.tf.eq(task.tf)]
        ax.axhline(0, color='#999999', lw=.6)
        ax.axvspan(-5, 5, color='#aaaaaa', alpha=.15)
        for arm, label, color in arms:
            data = sub[sub.arm.eq(arm)]
            if len(data):
                ax.plot(data.position, data.contrast, label=label, color=color, lw=.8)
        if sub.empty:
            ax.text(.5, .5, 'No matched cohort', ha='center', transform=ax.transAxes)
        ax.set(title=f'{task.cell} — {task.tf}', xlim=(-100, 100), xlabel='Bases from motif center',
               ylabel='Positive − negative contrast')
        ax.tick_params(labelsize=8)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=3, frameon=False)
    fig.suptitle('Training contrasts after unit-floor normalization\nFull matched cohorts; unsmoothed profiles; no uncertainty bands', fontsize=13)
    out.mkdir(parents=True)
    outputs = {}
    for extension in ['pdf', 'png']:
        path = out / ('training_contrasts.' + extension)
        fig.savefig(path, dpi=130)
        outputs[extension] = record(path)
    plt.close(fig)
    result = dict(source_manifest=record(source / 'crossfit_manifest.json'), code=record(Path(__file__)),
                  outputs=outputs, training_only=True, confidence_bands=False, promotion_eligible=False)
    result['stage_id'] = canonical_id(result)
    write_immutable(out / 'plots_manifest.json', result)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    args = p.parse_args()
    if not args.outdir.resolve().is_relative_to(study.ROOT / 'benchmarks/results'):
        p.error('use research results')
    plot(args.source.resolve(), args.outdir.resolve())
