#!/usr/bin/env python
"""Plot motif-centered pseudobulk TF aggregate cut-site profiles."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig

DEFAULT_GROUPS = ["B_cell", "CD4_T", "NK_T_cytotoxic", "CD14_Monocyte", "FCGR3A_Monocyte", "Dendritic_cell", "Mixed_myeloid"]
DEFAULT_TFS = ["BACH2", "RORA", "CEBPB", "CTCF"]
LINEAGE_BY_TF = {
    "PAX5": "B_cell", "EBF1": "B_cell", "POU2F2": "B_cell", "SPIB": "B_cell", "MEF2C": "B_cell", "BHLHE41": "B_cell", "BACH2": "B_cell",
    "TCF7": "T_NK", "TCF7L2": "T_NK", "LEF1": "T_NK", "Gata3": "T_NK", "GATA3": "T_NK", "RUNX3": "T_NK", "RORA": "T_NK", "TBX21": "T_NK", "EOMES": "T_NK", "IRF4": "T_NK",
    "Spi1": "Myeloid", "SPI1": "Myeloid", "CEBPB": "Myeloid", "CEBPA": "Myeloid", "CEBPD": "Myeloid", "JUNB": "Myeloid", "FOS": "Myeloid", "BATF": "Myeloid", "MAF": "Myeloid", "REL": "Myeloid",
    "CTCF": "Control",
}
LINEAGE_LABELS = {
    "B_cell": "B-cell-associated",
    "T_NK": "T/NK-associated",
    "Myeloid": "myeloid-associated",
    "Control": "ubiquitous control",
}
EXPECTED_GROUPS = {
    "B_cell": {"B_cell"},
    "T_NK": {"CD4_T", "NK_T_cytotoxic"},
    "Myeloid": {"CD14_Monocyte", "FCGR3A_Monocyte", "Dendritic_cell", "Mixed_myeloid"},
}
GROUP_COLORS = {
    "B_cell": "#1f77b4",
    "CD4_T": "#2ca02c",
    "NK_T_cytotoxic": "#9467bd",
    "CD14_Monocyte": "#d62728",
    "FCGR3A_Monocyte": "#ff7f0e",
    "Dendritic_cell": "#8c564b",
    "Mixed_myeloid": "#17becf",
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


def site_path(site_dir: Path, tf: str) -> Path:
    for suffix in ("motif_hits", "motif_peaks"):
        path = site_dir / f"{tf}.{suffix}.bed"
        if path.exists():
            return path
    raise FileNotFoundError(f"No motif-centered BED found for {tf} in {site_dir}")


def discover_tfs(site_dir: Path) -> list[str]:
    names = []
    for path in sorted(site_dir.glob("*.motif_hits.bed")) + sorted(site_dir.glob("*.motif_peaks.bed")):
        names.append(path.name.split(".motif_", 1)[0])
    return sorted(dict.fromkeys(names))


def read_site_summary(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None or not path.exists():
        return {}
    table = pd.read_csv(path, sep="\t")
    return {str(row["tf"]): {key: str(row[key]) for key in table.columns} for _, row in table.iterrows()}


def truthy(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


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
                if value == value:
                    values[index] += float(value)
            used += 1
    finally:
        bw.close()
    if used:
        values = [value / used for value in values]
    return values


def smooth_profile(values: list[float] | np.ndarray, window: int) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if window <= 1:
        return arr
    if window % 2 == 0:
        window += 1
    kernel = np.ones(window, dtype=float) / float(window)
    return np.convolve(arr, kernel, mode="same")


def protection_profile(values: list[float] | np.ndarray, center_half_width: int, flank_inner: int, flank_outer: int) -> np.ndarray:
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
        if flanks:
            protected[index] = float(np.nanmean(np.concatenate(flanks)) - arr[index])
    return protected


def lineage_for_tf(tf: str, summary: dict[str, dict[str, str]]) -> str:
    return summary.get(tf, {}).get("lineage") or LINEAGE_BY_TF.get(tf, "Control")


def line_style(lineage: str, group: str) -> dict[str, float]:
    expected = EXPECTED_GROUPS.get(lineage)
    if expected is None:
        return {"linewidth": 1.25, "alpha": 0.95}
    if group in expected:
        return {"linewidth": 1.9, "alpha": 1.0}
    return {"linewidth": 0.85, "alpha": 0.38}


def score_tf(tf: str, profiles: dict[tuple[str, str], list[float]], groups: list[str], lineage: str, center_half_width: int, flank_inner: int, flank_outer: int) -> dict[str, float | str]:
    group_scores = {}
    mid = len(next(iter(profiles.values()))) // 2 if profiles else 0
    for group in groups:
        profile = profiles.get((tf, group))
        if profile is None:
            group_scores[group] = np.nan
            continue
        protected = protection_profile(profile, center_half_width, flank_inner, flank_outer)
        group_scores[group] = float(np.nanmean(protected[max(0, mid - center_half_width): mid + center_half_width + 1]))
    expected = EXPECTED_GROUPS.get(lineage, set())
    if expected:
        expected_values = [group_scores[group] for group in expected if group in group_scores]
        other_values = [value for group, value in group_scores.items() if group not in expected]
        expected_mean = float(np.nanmean(expected_values)) if expected_values else float("nan")
        other_mean = float(np.nanmean(other_values)) if other_values else float("nan")
    else:
        expected_mean = float(np.nanmean(list(group_scores.values())))
        other_mean = expected_mean
    row: dict[str, float | str] = {"tf": tf, "lineage": lineage, "expected_mean_protection": expected_mean, "other_mean_protection": other_mean, "expected_contrast": expected_mean - other_mean}
    row.update({f"{group}_center_protection": value for group, value in group_scores.items()})
    return row


def numeric_column(frame: pd.DataFrame, column: str, default: float = 0.0) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(default, index=frame.index)
    return pd.to_numeric(frame[column], errors="coerce").fillna(default)


def select_tfs(screen: pd.DataFrame, min_total_sites: int, min_plotted_sites: int, per_lineage: int, controls: list[str], min_center_protection: float) -> list[str]:
    selected: list[str] = []
    total_sites = numeric_column(screen, "total_motif_hits")
    plotted_sites = numeric_column(screen, "plotted_sites", default=np.nan).fillna(numeric_column(screen, "n_sites"))
    center = numeric_column(screen, "expected_mean_protection")
    contrast = numeric_column(screen, "expected_contrast")
    base = (total_sites >= min_total_sites) & (plotted_sites >= min_plotted_sites) & (center > min_center_protection)
    for lineage in ["B_cell", "T_NK", "Myeloid"]:
        subset = screen[(screen["lineage"] == lineage) & base & (contrast > 0)].copy()
        subset = subset.sort_values(["expected_contrast", "expected_mean_protection", "plotted_sites"], ascending=[False, False, False]).head(per_lineage)
        selected.extend(subset["tf"].astype(str).tolist())
    for control in controls:
        row = screen[screen["tf"].astype(str) == control]
        if not row.empty and bool((base.loc[row.index] & (center.loc[row.index] > min_center_protection)).iloc[0]) and control not in selected:
            selected.append(control)
    if controls and not any(control in selected for control in controls):
        fill = screen[base & (contrast > 0) & ~screen["tf"].isin(selected)].copy()
        fill = fill.sort_values(["expected_contrast", "expected_mean_protection", "plotted_sites"], ascending=[False, False, False]).head(1)
        selected.extend(fill["tf"].astype(str).tolist())
    return selected


def _summary_int(summary: dict[str, dict[str, str]], tf: str, key: str, fallback: int) -> int:
    try:
        value = summary.get(tf, {}).get(key, "")
        return int(float(value)) if value not in ("", "nan") else fallback
    except (TypeError, ValueError):
        return fallback


def tf_label(tf: str, lineage: str, n_sites: int, summary: dict[str, dict[str, str]]) -> str:
    role = LINEAGE_LABELS.get(lineage, lineage)
    prefixes = summary.get(tf, {}).get("motif_prefixes", "")
    motif_text = prefixes.replace(",", ", ") if prefixes else tf
    total = _summary_int(summary, tf, "total_motif_hits", n_sites)
    plotted = _summary_int(summary, tf, "plotted_sites", n_sites)
    if total != plotted:
        count_text = f"{total:,} hits; top {plotted:,} plotted"
    else:
        count_text = f"{plotted:,} motif sites"
    return f"{tf} ({role})\n{motif_text}; {count_text}"


def plot_profiles(out_prefix: Path, tfs: list[str], groups: list[str], profiles: dict[tuple[str, str], list[float]], counts: dict[str, int], summary: dict[str, dict[str, str]], flank: int, ylabel: str, protection: bool, center_half_width: int, flank_inner: int, flank_outer: int) -> None:
    ncols = 2
    nrows = math.ceil(len(tfs) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.6, 2.85 * nrows), sharex=True)
    axes = list(axes.flat if hasattr(axes, "flat") else [axes])
    xvals = list(range(-flank, flank))
    for ax, tf in zip(axes, tfs):
        lineage = lineage_for_tf(tf, summary)
        for group in groups:
            profile = profiles.get((tf, group))
            if profile is None:
                continue
            values = protection_profile(profile, center_half_width, flank_inner, flank_outer) if protection else smooth_profile(profile, 5)
            ax.plot(xvals, values, label=group.replace("_", " "), color=GROUP_COLORS.get(group), **line_style(lineage, group))
        ax.axvline(0, color="black", linewidth=0.7, alpha=0.5)
        if protection:
            ax.axhline(0, color="0.65", linewidth=0.7)
        ax.set_title(tf_label(tf, lineage, counts.get(tf, 0), summary), fontsize=8.5, fontweight="bold")
        ax.set_ylabel(ylabel, fontsize=8.5)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[len(tfs):]:
        ax.axis("off")
    axes[min(len(tfs), len(axes)) - 1].legend(frameon=False, fontsize=7, loc="upper right")
    for ax in axes[-ncols:]:
        ax.set_xlabel("Distance from motif center (bp)", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(out_prefix.with_suffix(".pdf"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--tf-site-dir", required=True)
    parser.add_argument("--site-summary", default=None)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--groups", default=",".join(DEFAULT_GROUPS))
    parser.add_argument("--tfs", default=",".join(DEFAULT_TFS), help="Comma-separated TF names, or 'auto'.")
    parser.add_argument("--flank", type=int, default=250)
    parser.add_argument("--screen-output", default=None)
    parser.add_argument("--signal-column", default="cutsite_bigwig", help="Manifest column containing the bigWig signal to aggregate (default: cutsite_bigwig).")
    parser.add_argument("--footprint-like-output", default=None, help="Optional PNG/PDF prefix for flank-minus-center protection-score aggregate plots.")
    parser.add_argument("--protection-center-half-width", type=int, default=10)
    parser.add_argument("--protection-flank-inner", type=int, default=25)
    parser.add_argument("--protection-flank-outer", type=int, default=100)
    parser.add_argument("--auto-min-sites", type=int, default=50, help="Backward-compatible minimum plotted sites for auto-selection.")
    parser.add_argument("--auto-min-total-sites", type=int, default=500)
    parser.add_argument("--auto-min-plotted-sites", type=int, default=None)
    parser.add_argument("--auto-min-center-protection", type=float, default=0.0)
    parser.add_argument("--auto-per-lineage", type=int, default=1)
    parser.add_argument("--control-tfs", default="CTCF")
    parser.add_argument("--site-selection", choices=["all", "lineage-top"], default="all", help="Deprecated compatibility option; motif-centered paper figures use all selected motif hits.")
    parser.add_argument("--max-selected-sites", type=int, default=25, help="Deprecated compatibility option retained for older tests.")
    args = parser.parse_args(argv)

    manifest = pd.read_csv(args.manifest, sep="\t")
    groups = [group.strip() for group in args.groups.split(",") if group.strip()]
    site_dir = Path(args.tf_site_dir)
    summary = read_site_summary(Path(args.site_summary) if args.site_summary else site_dir / "motif_centered_site_summary.tsv")
    requested_tfs = discover_tfs(site_dir) if args.tfs == "auto" else [tf.strip() for tf in args.tfs.split(",") if tf.strip()]

    group_bigwigs = {}
    pass_mask = manifest["passes_filters"].map(truthy) if "passes_filters" in manifest.columns else pd.Series(True, index=manifest.index)
    for group in groups:
        rows = manifest[(manifest["group"] == group) & pass_mask]
        if rows.empty:
            continue
        if args.signal_column not in rows.columns:
            raise SystemExit(f"Manifest does not contain signal column: {args.signal_column}")
        bigwig = Path(str(rows.iloc[0][args.signal_column]))
        if bigwig.exists():
            group_bigwigs[group] = bigwig

    profiles: dict[tuple[str, str], list[float]] = {}
    counts: dict[str, int] = {}
    records = []
    screen_rows = []
    for tf in requested_tfs:
        sites = load_sites(site_path(site_dir, tf))
        counts[tf] = len(sites)
        for group, bigwig in group_bigwigs.items():
            profile = mean_profile(bigwig, sites, args.flank)
            profiles[(tf, group)] = profile
            for offset, value in zip(range(-args.flank, args.flank), profile):
                records.append({"tf": tf, "group": group, "offset_bp": offset, "cutsite_cpm": value, "n_sites": len(sites), "total_motif_hits": _summary_int(summary, tf, "total_motif_hits", len(sites)), "signal_column": args.signal_column, "site_coordinate": "motif_center", "lineage": lineage_for_tf(tf, summary)})
        row = score_tf(tf, profiles, groups, lineage_for_tf(tf, summary), args.protection_center_half_width, args.protection_flank_inner, args.protection_flank_outer)
        row["plotted_sites"] = len(sites)
        row["n_sites"] = len(sites)
        row["total_motif_hits"] = _summary_int(summary, tf, "total_motif_hits", len(sites))
        row["selection_method"] = summary.get(tf, {}).get("selection_method", "input_bed")
        row["score_min_selected"] = summary.get(tf, {}).get("score_min_selected", "")
        row["motif_prefixes"] = summary.get(tf, {}).get("motif_prefixes", "")
        screen_rows.append(row)

    screen = pd.DataFrame(screen_rows)
    if args.tfs == "auto":
        min_plotted = args.auto_min_plotted_sites if args.auto_min_plotted_sites is not None else args.auto_min_sites
        selected = select_tfs(screen, args.auto_min_total_sites, min_plotted, args.auto_per_lineage, [tf.strip() for tf in args.control_tfs.split(",") if tf.strip()], args.auto_min_center_protection)
    else:
        selected = requested_tfs
    screen["selected_for_figure"] = screen["tf"].isin(selected)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    profile_table = pd.DataFrame(records)
    profile_table = profile_table[profile_table["tf"].isin(selected)]
    profile_table.to_csv(out_prefix.with_suffix(".tsv"), sep="\t", index=False)
    screen_path = Path(args.screen_output) if args.screen_output else out_prefix.with_name(out_prefix.name + "_screen.tsv")
    screen.to_csv(screen_path, sep="\t", index=False)

    if not selected:
        raise SystemExit("No TFs passed the auto-selection criteria; inspect the screen TSV or lower --auto-min-sites.")
    plot_profiles(out_prefix, selected, groups, profiles, counts, summary, args.flank, "Cut-site signal (CPM)", False, args.protection_center_half_width, args.protection_flank_inner, args.protection_flank_outer)

    if args.footprint_like_output:
        protection_prefix = Path(args.footprint_like_output)
        protection_records = []
        for tf in selected:
            for group in groups:
                profile = profiles.get((tf, group))
                if profile is None:
                    continue
                protected = protection_profile(profile, args.protection_center_half_width, args.protection_flank_inner, args.protection_flank_outer)
                for offset, value in zip(range(-args.flank, args.flank), protected):
                    protection_records.append({"tf": tf, "group": group, "offset_bp": offset, "protection_score": value, "n_sites": counts.get(tf, 0), "total_motif_hits": _summary_int(summary, tf, "total_motif_hits", counts.get(tf, 0)), "signal_column": args.signal_column, "site_coordinate": "motif_center", "lineage": lineage_for_tf(tf, summary)})
        pd.DataFrame(protection_records).to_csv(protection_prefix.with_suffix(".tsv"), sep="\t", index=False)
        plot_profiles(protection_prefix, selected, groups, profiles, counts, summary, args.flank, "Protection score (flank - center)", True, args.protection_center_half_width, args.protection_flank_inner, args.protection_flank_outer)

    print(f"Wrote {out_prefix.with_suffix('.png')}, {out_prefix.with_suffix('.tsv')}, and {screen_path}")
    if args.footprint_like_output:
        print(f"Wrote {Path(args.footprint_like_output).with_suffix('.png')} and {Path(args.footprint_like_output).with_suffix('.tsv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
