#!/usr/bin/env python
"""Plot real-data pseudobulk TF aggregate cut-site profiles."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--tf-site-dir", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--groups", default=",".join(DEFAULT_GROUPS))
    parser.add_argument("--tfs", default=",".join(DEFAULT_TFS))
    parser.add_argument("--flank", type=int, default=500)
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
    print(f"Wrote {out_prefix.with_suffix('.png')} and {table_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
