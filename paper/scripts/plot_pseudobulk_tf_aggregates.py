#!/usr/bin/env python
"""Plot real-data pseudobulk TF aggregate cut-site profiles."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig

DEFAULT_GROUPS = ["B_cell", "CD4_T", "NK_T_cytotoxic", "CD14_Monocyte", "FCGR3A_Monocyte", "Dendritic_cell"]
DEFAULT_TFS = ["PAX5", "TCF7", "CEBPB", "CTCF"]
GROUP_COLORS = {
    "B_cell": "#1f77b4",
    "CD4_T": "#2ca02c",
    "NK_T_cytotoxic": "#9467bd",
    "CD14_Monocyte": "#d62728",
    "FCGR3A_Monocyte": "#ff7f0e",
    "Dendritic_cell": "#8c564b",
}


def load_sites(path: Path) -> list[tuple[str, int, int]]:
    sites = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            center = (int(start) + int(end)) // 2
            sites.append((chrom, center, center + 1))
    return sites


def mean_profile(bigwig: Path, sites: list[tuple[str, int, int]], flank: int) -> list[float]:
    bw = pyBigWig.open(str(bigwig))
    values = [0.0] * (2 * flank)
    used = 0
    try:
        chroms = bw.chroms()
        for chrom, center, _end in sites:
            start = center - flank
            end = center + flank
            if chrom not in chroms or start < 0 or end > chroms[chrom]:
                continue
            row = bw.values(chrom, start, end)
            if len(row) != 2 * flank:
                continue
            for index, value in enumerate(row):
                if value == value:  # skip NaN positions in sparse bigWigs
                    values[index] += float(value)
            used += 1
    finally:
        bw.close()
    if used:
        values = [value / used for value in values]
    return values


def smooth_profile(values: list[float], window: int) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if window <= 1:
        return arr
    if window % 2 == 0:
        window += 1
    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(arr, kernel, mode="same")


def protection_profile(values: list[float], center_half_width: int, flank_inner: int, flank_outer: int) -> np.ndarray:
    arr = smooth_profile(values, center_half_width * 2 + 1)
    n = arr.shape[0]
    protected = np.zeros(n, dtype=float)
    for index in range(n):
        left_start = max(0, index - flank_outer)
        left_end = max(0, index - flank_inner)
        right_start = min(n, index + flank_inner + 1)
        right_end = min(n, index + flank_outer + 1)
        flanks = []
        if left_end > left_start:
            flanks.append(arr[left_start:left_end])
        if right_end > right_start:
            flanks.append(arr[right_start:right_end])
        if not flanks:
            continue
        protected[index] = float(np.nanmean(np.concatenate(flanks)) - arr[index])
    return protected


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--tf-site-dir", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--groups", default=",".join(DEFAULT_GROUPS))
    parser.add_argument("--tfs", default=",".join(DEFAULT_TFS))
    parser.add_argument("--flank", type=int, default=500)
    parser.add_argument("--footprint-like-output", default=None, help="Optional PNG/PDF prefix for flank-minus-center protection-score aggregate plots.")
    parser.add_argument("--protection-center-half-width", type=int, default=10)
    parser.add_argument("--protection-flank-inner", type=int, default=25)
    parser.add_argument("--protection-flank-outer", type=int, default=100)
    args = parser.parse_args(argv)

    manifest = pd.read_csv(args.manifest, sep="\t")
    groups = [group.strip() for group in args.groups.split(",") if group.strip()]
    tfs = [tf.strip() for tf in args.tfs.split(",") if tf.strip()]
    records = []
    profiles = {}

    for tf in tfs:
        sites = load_sites(Path(args.tf_site_dir) / f"{tf}.motif_peaks.bed")
        for group in groups:
            rows = manifest[(manifest["group"] == group) & (manifest["passes_filters"].astype(bool))]
            if rows.empty:
                continue
            bigwig = Path(str(rows.iloc[0]["cutsite_bigwig"]))
            if not bigwig.exists():
                continue
            profile = mean_profile(bigwig, sites, args.flank)
            profiles[(tf, group)] = profile
            for offset, value in zip(range(-args.flank, args.flank), profile):
                records.append({"tf": tf, "group": group, "offset_bp": offset, "cutsite_cpm": value, "n_sites": len(sites)})

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    table_path = out_prefix.with_suffix(".tsv")
    pd.DataFrame(records).to_csv(table_path, sep="\t", index=False)

    ncols = 2
    nrows = math.ceil(len(tfs) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.2, 2.7 * nrows), sharex=True)
    axes = list(axes.flat if hasattr(axes, "flat") else [axes])
    xvals = list(range(-args.flank, args.flank))
    for ax, tf in zip(axes, tfs):
        for group in groups:
            profile = profiles.get((tf, group))
            if profile is None:
                continue
            ax.plot(xvals, profile, label=group.replace("_", " "), linewidth=1.2, color=GROUP_COLORS.get(group))
        ax.axvline(0, color="black", linewidth=0.7, alpha=0.5)
        ax.set_title(tf)
        ax.set_ylabel("Cut sites (CPM)")
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[len(tfs):]:
        ax.axis("off")
    axes[min(len(tfs), len(axes)) - 1].legend(frameon=False, fontsize=7, loc="upper right")
    for ax in axes[-ncols:]:
        ax.set_xlabel("Distance from motif-associated peak center (bp)")
    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(out_prefix.with_suffix(".pdf"))

    if args.footprint_like_output:
        protection_prefix = Path(args.footprint_like_output)
        protection_prefix.parent.mkdir(parents=True, exist_ok=True)
        protection_records = []
        fig2, axes2 = plt.subplots(nrows, ncols, figsize=(7.2, 2.7 * nrows), sharex=True)
        axes2 = list(axes2.flat if hasattr(axes2, "flat") else [axes2])
        for ax, tf in zip(axes2, tfs):
            for group in groups:
                profile = profiles.get((tf, group))
                if profile is None:
                    continue
                protected = protection_profile(
                    profile,
                    center_half_width=args.protection_center_half_width,
                    flank_inner=args.protection_flank_inner,
                    flank_outer=args.protection_flank_outer,
                )
                ax.plot(xvals, protected, label=group.replace("_", " "), linewidth=1.2, color=GROUP_COLORS.get(group))
                for offset, value in zip(xvals, protected):
                    protection_records.append({"tf": tf, "group": group, "offset_bp": offset, "protection_score": value})
            ax.axhline(0, color="0.65", linewidth=0.7)
            ax.axvline(0, color="black", linewidth=0.7, alpha=0.5)
            ax.set_title(tf)
            ax.set_ylabel("Protection score\n(flank - center CPM)")
            ax.spines[["top", "right"]].set_visible(False)
        for ax in axes2[len(tfs):]:
            ax.axis("off")
        axes2[min(len(tfs), len(axes2)) - 1].legend(frameon=False, fontsize=7, loc="upper right")
        for ax in axes2[-ncols:]:
            ax.set_xlabel("Distance from motif-associated peak center (bp)")
        fig2.suptitle("Footprint-like pseudobulk protection score", y=1.01)
        fig2.tight_layout()
        fig2.savefig(protection_prefix.with_suffix(".png"), dpi=300)
        fig2.savefig(protection_prefix.with_suffix(".pdf"))
        pd.DataFrame(protection_records).to_csv(protection_prefix.with_suffix(".tsv"), sep="\t", index=False)

    print(f"Wrote {out_prefix.with_suffix('.png')} and {table_path}")
    if args.footprint_like_output:
        print(f"Wrote {Path(args.footprint_like_output).with_suffix('.png')} and {Path(args.footprint_like_output).with_suffix('.tsv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
