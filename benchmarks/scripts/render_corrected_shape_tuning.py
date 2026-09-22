#!/usr/bin/env python3
"""Render every corrected tuning task without implying confidence or promotion."""

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
import run_corrected_shape_development as experiment
from audit_direct_bias_safety import record, write_immutable
from audit_chrombpnet_comparison import checked_path, canonical_id

CORRECTIONS=("raw","DWM_expected","parametric_direct")
LABELS=("Raw counts","DWM expected-signal residual","Parametric residual")


def validate_tables(metrics,curves):
    keys=["cell","tf","correction"]
    if metrics.duplicated(keys+["detector"]).any() or curves.duplicated(keys+["position"]).any():
        raise ValueError("duplicate task/detector or curve position")
    tasks=sorted(set(map(tuple,metrics[["cell","tf"]].to_numpy())))
    if set(map(tuple,curves[["cell","tf"]].to_numpy()))!=set(tasks):
        raise ValueError("curve and metric task denominators differ")
    for cell,tf in tasks:
        rows=metrics.loc[metrics.cell.eq(cell)&metrics.tf.eq(tf)]
        if (len(rows)!=9 or set(rows.correction)!=set(CORRECTIONS)
                or set(rows.detector)!={"conventional","trained_geometry","hierarchical_spline"}
                or rows.matched_pairs.nunique()!=1 or rows.cohort_status.nunique()!=1):
            raise ValueError("incomplete or inconsistent task factorial")
        for correction in CORRECTIONS:
            curve=curves.loc[curves.cell.eq(cell)&curves.tf.eq(tf)&curves.correction.eq(correction)].sort_values("position")
            if not np.array_equal(curve.position.to_numpy(),np.arange(-100,101)):
                raise ValueError("incomplete aggregate curve")
            values=curve[["positive_mean","negative_mean","positive_minus_negative"]].to_numpy(float)
            if not np.isfinite(values).all() or not np.allclose(values[:,0]-values[:,1],values[:,2],rtol=1e-8,atol=1e-10):
                raise ValueError("aggregate contrast differs from means")
    return tasks


def render(metrics,curves,out):
    tasks=validate_tables(metrics,curves)
    if out.exists():
        raise ValueError("render into a new directory")
    out.mkdir(parents=True)
    pdf=out/"all_task_aggregates.pdf"
    pages=[]
    with PdfPages(pdf) as document:
        for start in range(0,len(tasks),3):
            fig,axes=plt.subplots(3,3,figsize=(13,9),sharex=True,squeeze=False)
            fig.suptitle("Corrected training/tuning study: aggregate shape by correction",fontsize=14)
            for row in range(3):
                if start+row>=len(tasks):
                    for ax in axes[row]:ax.set_visible(False)
                    continue
                cell,tf=tasks[start+row]
                selected=metrics.loc[metrics.cell.eq(cell)&metrics.tf.eq(tf)].iloc[0]
                data=curves.loc[curves.cell.eq(cell)&curves.tf.eq(tf)]
                limit=max(.03,float(data.loc[data.position.abs().le(60),"positive_minus_negative"].abs().max())*1.12)
                for column,correction in enumerate(CORRECTIONS):
                    ax=axes[row,column]
                    curve=data.loc[data.correction.eq(correction)].sort_values("position")
                    ax.axhline(0,color="0.5",lw=.7)
                    ax.axvspan(-5,5,color="0.85",alpha=.5)
                    ax.axvspan(-40,-15,color="#e9f2f7",alpha=.6)
                    ax.axvspan(15,40,color="#e9f2f7",alpha=.6)
                    ax.plot(curve.position,curve.positive_minus_negative,color=("#777777","#2369a1","#a34032")[column],lw=1.3)
                    ax.set(xlim=(-60,60),ylim=(-limit,limit),xlabel="Position relative to motif center (base pairs)")
                    ax.set_title(f"{cell} {tf} | {LABELS[column]}\n{int(selected.matched_pairs)} pairs; {selected.cohort_status}",fontsize=9)
                    if column==0:ax.set_ylabel("Positive − negative\nnormalized profile")
                    ax.tick_params(labelsize=8)
            fig.text(.5,.018,"Tuning chromosomes only. No confidence bands. Detector choice changes ranks, not these input profiles.\nDWM: dinucleotide weight matrix. Underpowered tasks are retained; no promotion claim.",ha="center",fontsize=9)
            fig.tight_layout(rect=(.02,.065,.99,.95))
            document.savefig(fig)
            path=out/f"page_{len(pages)+1:02d}.png"
            fig.savefig(path,dpi=120)
            pages.append(path)
            plt.close(fig)
    return [pdf]+pages


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root",type=Path,required=True)
    parser.add_argument("--result-dir",type=Path,required=True)
    parser.add_argument("--outdir",type=Path,required=True)
    args=parser.parse_args()
    out=args.outdir.resolve()
    if not out.is_relative_to(experiment.ROOT/"benchmarks/results"):
        parser.error("use research results")
    experiment.evaluate(args.data_root.resolve(),args.result_dir.resolve(),verify_only=True)
    source=args.result_dir.resolve()/"tuning_manifest.json"
    doc=json.loads(source.read_text())
    metrics=pd.read_csv(checked_path(doc["outputs"]["metrics"],args.data_root),sep="\t")
    curves=pd.read_csv(checked_path(doc["outputs"]["curves"],args.data_root),sep="\t")
    paths=render(metrics,curves,out)
    manifest=dict(source=record(source),code=record(Path(__file__)),tasks=len(metrics[["cell","tf"]].drop_duplicates()),
                  outputs=[record(p) for p in paths],confidence_bands=False,promotion_eligible=False)
    manifest["report_id"]=canonical_id(manifest)
    write_immutable(out/"render_manifest.json",manifest)
    print(out/"all_task_aggregates.pdf")
