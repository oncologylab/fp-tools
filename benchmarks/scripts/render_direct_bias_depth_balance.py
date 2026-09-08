#!/usr/bin/env python3
"""Render all-TF primary and current-depth sensitivity profiles on shared axes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_direct_bias_depth_balance as audit
from audit_chrombpnet_comparison import canonical_id, checked_path
from audit_direct_bias_safety import record, write_immutable


def render(metrics, profiles, retention, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    primary = "motif_accessibility_observed_total"
    methods = {audit.base.DWM: ("#777777", "Fixed DWM"), audit.base.SAMPLE_DWM: ("#7B3294", "Sample DWM"),
               audit.base.PARAMETRIC: ("#D55E00", "Parametric")}
    pdf_path = outdir / "current_depth_balance_all_tf.pdf"
    outputs, pages = {}, []
    with PdfPages(pdf_path, metadata={"CreationDate": None, "ModDate": None}) as pdf:
        for (cell, tf), retained in retention.groupby(["cell", "tf"], sort=True):
            task = profiles.loc[profiles.cell.eq(cell) & profiles.tf.eq(tf)]
            fig, axes = plt.subplots(2, 4, sharey=True, figsize=(11.69, 8.27))
            fig.suptitle(f"{cell} {tf}: sensitivity to current-depth cut balance\n"
                         "Original partners; fixed 0.2-SD training-scale caliper; unchanged candidate", fontsize=12)
            for row, (cohort, size_field, label) in enumerate(((primary, "original_pairs", "Original pairs"),
                                                           (audit.SECONDARY, "retained_pairs", "Current-cut caliper"))):
                for column, depth in enumerate(("10m", "25m", "50m", "full")):
                    ax = axes[row, column]
                    counts = retained.loc[retained.depth.eq(depth), size_field]
                    panel = task.loc[task.cohort.eq(cohort) & task.depth.eq(depth)]
                    for method, (color, name) in methods.items():
                        data = panel.loc[panel.method.eq(method)]
                        for _, group in data.groupby(["sample", "seed"], sort=True):
                            group = group.sort_values("position")
                            ax.plot(group.position, group.positive_minus_negative, color=color, alpha=.15, linewidth=.6)
                        if not data.empty:
                            mean = data.groupby("position").positive_minus_negative.mean()
                            ax.plot(mean.index, mean.values, color=color, linewidth=1.1, label=name)
                    candidate = metrics.loc[metrics.cell.eq(cell) & metrics.tf.eq(tf) & metrics.depth.eq(depth)
                                            & metrics.cohort.eq(cohort) & metrics.method.eq(audit.base.PARAMETRIC)]
                    gain = (f"ΔAUROC vs raw {candidate.auroc_gain_over_raw.min():+.3f} to {candidate.auroc_gain_over_raw.max():+.3f}"
                            if len(candidate) else "No estimable comparison")
                    power = f"{int((counts < 200).sum())}/{len(counts)} runs underpowered"
                    ax.set(title=f"{depth}: {label}; {int(counts.min())}–{int(counts.max())} pairs\n{power}\n{gain}",
                           xlabel="Motif offset (bp)", xlim=(-100, 100))
                    if panel.empty:
                        ax.text(.5, .5, "Fewer than two pairs\nProfile not estimated", ha="center", va="center",
                                transform=ax.transAxes, fontsize=8)
                    ax.title.set_fontsize(7.5)
                    ax.tick_params(labelsize=7)
                    ax.axhline(0, color="#888888", linewidth=.5)
                    ax.axvspan(-5, 5, alpha=.15, color="#aaaaaa")
                    ax.spines[["top", "right"]].set_visible(False)
                    if column == 0:
                        ax.set_ylabel("Positive − negative\n(normalized deviance)", fontsize=8)
            if axes[0, 0].get_legend_handles_labels()[0]:
                axes[0, 0].legend(fontsize=7, frameon=False)
            fig.text(.065, .045,
                     "Shared vertical scale. Thin curves: individual runs; thick curves: their mean. Each within-run method comparison uses identical sites.\n"
                     "Five nested seeds of replicate 1 at 10M/25M/50M; three biological replicates at full. Per-run pointwise bands and metric intervals are in tables.\n"
                     "Post-hoc examined chr19–22/X. Different subsets can differ in task difficulty; conditioning on ATAC counts does not establish a causal correction benefit.",
                     fontsize=7.3)
            fig.subplots_adjust(left=.07, right=.98, top=.8, bottom=.16, wspace=.28, hspace=.58)
            pdf.savefig(fig)
            pages.append({"cell": cell, "tf": tf})
            if tf == "CTCF":
                path = outdir / f"{cell}_CTCF_current_depth_balance.png"
                fig.savefig(path, dpi=160)
                outputs[path.name] = record(path)
            plt.close(fig)
    outputs[pdf_path.name] = record(pdf_path)
    return outputs, pages


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--calculation-commit", required=True)
    args = parser.parse_args(argv)
    manifest = args.outdir / "figure_manifest.json"
    if manifest.exists():
        done = audit.read_document(manifest, "figure_id")
        if done["source_result"] != record(args.result) or done["code"] != record(Path(__file__)):
            raise ValueError("figure source or renderer differs from the completed report")
        for entry in done["outputs"].values():
            checked_path(entry, args.data_root)
        print("Verified completed depth-balance figures; no rendering repeated.")
        return 0
    result = audit.verify_result(args.result.parent, args.data_root, args.calculation_commit)
    if result != json.loads(args.result.read_text()):
        raise ValueError("requested result is not the verified sensitivity manifest")
    tables = {name: pd.read_csv(checked_path(result["outputs"][name], args.data_root), sep="\t")
              for name in ("metrics", "profiles", "retention")}
    args.outdir.mkdir(parents=True, exist_ok=False)
    outputs, pages = render(tables["metrics"], tables["profiles"], tables["retention"], args.outdir)
    done = {"source_result": record(args.result), "calculation_commit": args.calculation_commit,
            "code": record(Path(__file__)), "outputs": outputs, "pages": pages,
            "evidence": result["evidence"], "promotion_eligible": False, "causal_effect_established": False}
    done["figure_id"] = canonical_id(done)
    write_immutable(manifest, done)
    print(json.dumps({"pages": len(pages), "figure_id": done["figure_id"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
