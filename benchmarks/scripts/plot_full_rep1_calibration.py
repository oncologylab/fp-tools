#!/usr/bin/env python3
"""Plot all tasks/regions without interpreting null intervals as biological coverage."""
from pathlib import Path
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import checked_path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'benchmarks/results/full_rep1_region_calibration_20260922_v1/calibration'
OUT = ROOT / 'benchmarks/results/full_rep1_calibration_plots_20260922_v1'


def main():
    manifest = checked_path(dict(path=str(SOURCE / 'calibration_manifest.json'),
        sha256='5767c81e3a6464a3cc48d08e2e7d6b2a7075514af566f0cd0f81a9779d82e7c8'), ROOT)
    entry = json.loads(manifest.read_text())['outputs']['calibration']
    frame = pd.read_csv(checked_path(entry, ROOT), sep='\t')
    frame = frame[frame.strand.eq('pooled') & frame.arm.eq('parametric')].copy()
    tasks = list(frame[['cell', 'tf']].drop_duplicates().itertuples(index=False, name=None))
    assert len(tasks) == 21 and len(frame) == 63
    OUT.mkdir(parents=True, exist_ok=False)
    write_immutable(OUT / 'inputs.json', dict(inputs=[record(manifest), entry], code=record(Path(__file__)),
        interpretation='examined replicate-1 diagnostic; predictive intervals assume independent cuts, not validated biological coverage'))
    fig, axes = plt.subplots(1, 3, figsize=(15, 10), sharey=True)
    for ax, region in zip(axes, ['center', 'left', 'right']):
        rows = frame[frame.region.eq(region)].set_index(['cell', 'tf']).loc[tasks]
        y = list(range(len(tasks)))
        ax.hlines(y, rows.independent_null_count_low/rows.total_cuts,
                  rows.independent_null_count_high/rows.total_cuts,
                  color='#80a5c7', linewidth=4, label='95% model-predictive count interval')
        ax.scatter(rows.expected_fraction, y, marker='|', color='#164d7d', s=100, label='Parametric expectation', zorder=3)
        ax.scatter(rows.observed_fraction, y, color='#b74225', s=27, label='Observed fraction', zorder=4)
        ax.axvline((15 if region == 'center' else 29)/73, color='#555555', ls='--', label='Uniform expectation')
        ax.set_title(region.capitalize()); ax.set_xlabel('Fraction of cuts across center and shoulders')
        ax.set_xlim(0, .7); ax.grid(axis='x', alpha=.2)
    axes[0].set_yticks(range(len(tasks)), [f'{cell} {tf}' for cell, tf in tasks])
    axes[0].invert_yaxis()
    fig.suptitle('Higher-count naked-DNA replicate 1: fixed-region calibration', fontsize=16)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, bbox_to_anchor=(.5, .035))
    fig.text(.5, .012, 'Examined development control; unadjusted independent-cut intervals; preprocessing differs from sparse-count input.', ha='center', fontsize=10)
    fig.tight_layout(rect=(0, .10, 1, .95))
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(OUT / ('calibration.' + ext), dpi=150)
    plt.close(fig)
    write_immutable(OUT / 'outputs.json', {ext: record(OUT / ('calibration.' + ext)) for ext in ['png', 'pdf', 'svg']})


if __name__ == '__main__':
    main()
