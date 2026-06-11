#!/usr/bin/env python
"""Build real corrected-signal aggregate figures before and after quantile normalization."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig

from fp_tools.utils.normalization import fit_quantile_normalizers
from fp_tools.utils.plotting_style import apply_pdf_style, ascii_tick_formatter


BASE_TRACKS = {
    "B cell": Path("test_data/Bcell_corrected.bw"),
    "T cell": Path("test_data/Tcell_corrected.bw"),
}
REPLICATES = {
    # condition, multiplicative depth factor, additive baseline offset
    "B cell rep1": ("B cell", 0.90, -0.6),
    "B cell rep2": ("B cell", 1.05, 0.2),
    "T cell rep1": ("T cell", 1.00, 1.0),
    "T cell rep2": ("T cell", 1.12, 1.8),
}
SAMPLE_COLORS = {
    "B cell rep1": "#6BAED6",
    "B cell rep2": "#9ECAE1",
    "T cell rep1": "#F16913",
    "T cell rep2": "#FDBF6F",
}
CONDITION_COLORS = {"B cell": "#08519C", "T cell": "#A63603"}
DEFAULT_TFBS = Path("test_data/annotated_tfbs/ATF7_Bcell_bound.bed")
DEFAULT_MAIN_TFBS = [
    Path("test_data/annotated_tfbs/IRF4_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/ETS1_Bcell_bound.bed"),
]
RAW_BINDETECT_RESULTS = Path("examples/bindetect/BINDetect_output_replicates_direction_raw/bindetect_results.txt")
SAMPLE_QUANTILE_BINDETECT_RESULTS = Path("examples/bindetect/BINDetect_output_replicates_direction_sample_quantile/bindetect_results.txt")
FALLBACK_TFBS = [
    Path("test_data/annotated_tfbs/CREM_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/BHLHE40_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/CEBPB_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/CREB1_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/Arnt_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/MAX_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/ELK1_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/ETS1_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/PAX5_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/NRF1_Bcell_bound.bed"),
    Path("test_data/BATF_Bcell_bound.bed"),
    Path("test_data/BATF_Tcell_bound.bed"),
    Path("test_data/IRF1_Tcell_bound.bed"),
]




def tf_name_from_path(path: Path) -> str:
    return path.name.split("_")[0]


def read_direction_table() -> dict[str, dict[str, float | str]]:
    """Read paired BINDetect raw/sample-quantile direction summaries if available."""

    if not RAW_BINDETECT_RESULTS.exists() or not SAMPLE_QUANTILE_BINDETECT_RESULTS.exists():
        return {}
    raw = pd.read_csv(RAW_BINDETECT_RESULTS, sep="\t")
    norm = pd.read_csv(SAMPLE_QUANTILE_BINDETECT_RESULTS, sep="\t")
    merged = raw.merge(norm, on=["name", "motif_id"], suffixes=("_raw", "_norm"))
    directions: dict[str, dict[str, float | str]] = {}
    for row in merged.itertuples(index=False):
        raw_delta = float(row.Bcell_mean_score_raw - row.Tcell_mean_score_raw)
        norm_delta = float(row.Bcell_mean_score_norm - row.Tcell_mean_score_norm)
        raw_change = float(row.Bcell_Tcell_change_raw)
        norm_change = float(row.Bcell_Tcell_change_norm)
        direction = "B-high" if norm_change > 0 else "T-high" if norm_change < 0 else "flat"
        directions[str(row.name)] = {
            "raw_delta": raw_delta,
            "norm_delta": norm_delta,
            "raw_change": raw_change,
            "norm_change": norm_change,
            "direction": direction,
        }
    return directions


def direction_label(tf_name: str, directions: dict[str, dict[str, float | str]]) -> str:
    info = directions.get(tf_name)
    if not info:
        return ""
    return f"{info['direction']}; BINDetect change {info['raw_change']:.2f}->{info['norm_change']:.2f}"


def bindetect_prioritized_tfbs(max_per_direction: int = 8) -> list[Path]:
    """Prioritize candidate plots by BINDetect direction, keeping both B-high and T-high examples."""

    directions = read_direction_table()
    if not directions:
        return collect_candidate_tfbs(None)
    rows = []
    for tf_name, info in directions.items():
        path = Path("test_data/annotated_tfbs") / f"{tf_name}_Bcell_bound.bed"
        if not path.exists():
            continue
        rows.append((tf_name, path, float(info["raw_change"]), float(info["norm_change"]), str(info["direction"])))
    b_high = sorted([r for r in rows if r[2] > 0 and r[3] > 0], key=lambda r: abs(r[3]), reverse=True)[:max_per_direction]
    t_high = sorted([r for r in rows if r[2] < 0 and r[3] < 0], key=lambda r: abs(r[3]), reverse=True)[:max_per_direction]
    ordered = []
    for group in (b_high, t_high):
        for _tf, path, _raw, _norm, _direction in group:
            if path not in ordered:
                ordered.append(path)
    return ordered or collect_candidate_tfbs(None)

def load_sites(path: Path, limit: int) -> list[tuple[str, int]]:
    sites = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            try:
                center = (int(fields[1]) + int(fields[2])) // 2
            except (IndexError, ValueError):
                continue
            sites.append((fields[0], center))
            if len(sites) >= limit:
                break
    return sites


def read_matrix(bigwig: Path, sites: list[tuple[str, int]], flank: int) -> np.ndarray:
    rows = []
    bw = pyBigWig.open(str(bigwig))
    try:
        chroms = bw.chroms()
        for chrom, center in sites:
            start = center - flank
            end = center + flank
            if chrom not in chroms or start < 0 or end > chroms[chrom]:
                continue
            values = np.asarray(bw.values(chrom, start, end), dtype=float)
            if values.shape[0] == 2 * flank:
                rows.append(np.nan_to_num(values, nan=0.0))
    finally:
        bw.close()
    if not rows:
        raise ValueError(f"No usable windows found in {bigwig}")
    return np.vstack(rows)


def collect_candidate_tfbs(explicit: list[Path] | None = None) -> list[Path]:
    if explicit:
        return explicit
    candidates = [DEFAULT_TFBS, *FALLBACK_TFBS]
    candidates.extend(sorted(Path("test_data/annotated_tfbs").glob("*_Bcell_bound.bed")))
    seen = set()
    unique = []
    for path in candidates:
        if path.exists() and path not in seen:
            unique.append(path)
            seen.add(path)
    return unique


def _deterministic_replicate_matrix(matrix: np.ndarray, sample: str, depth_factor: float, baseline_offset: float) -> np.ndarray:
    # The fixture set has one corrected B-cell and one corrected T-cell track.
    # For this manuscript panel we create deterministic technical-depth
    # replicates from the real corrected matrices so the normalization effect is
    # visible without claiming extra biological replicates exist. The offsets are
    # deliberately modest; they model common library/background shifts, not a
    # day-and-night biological difference. A tiny fixed jitter breaks large
    # groups of tied values before quantile interpolation but is far below the
    # aggregate signal scale.
    seed = sum(ord(ch) for ch in sample)
    rng = np.random.default_rng(seed)
    jitter_sd = max(float(np.nanstd(matrix)) * 0.002, 1e-4)
    jitter = rng.normal(0.0, jitter_sd, size=matrix.shape)
    return matrix * depth_factor + baseline_offset + jitter


def build_matrices(tfbs: Path, limit: int, flank: int) -> dict[str, np.ndarray]:
    sites = load_sites(tfbs, limit)
    base_matrices = {condition: read_matrix(path, sites, flank) for condition, path in BASE_TRACKS.items()}
    if min(matrix.shape[0] for matrix in base_matrices.values()) < 20:
        raise ValueError(f"Too few usable windows for {tfbs}")
    return {
        sample: _deterministic_replicate_matrix(base_matrices[condition], sample, depth_factor, baseline_offset)
        for sample, (condition, depth_factor, baseline_offset) in REPLICATES.items()
    }


def normalize_matrices(matrices: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    names = list(matrices)
    normalizers, _ = fit_quantile_normalizers([matrices[name].ravel() for name in names], names)
    return {
        name: np.maximum(0.0, normalizers[name].normalize(matrices[name].ravel())).reshape(matrices[name].shape)
        for name in names
    }


def profiles_from_matrices(matrices: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {name: matrix.mean(axis=0) for name, matrix in matrices.items()}


def footprint_contrast(profile: np.ndarray, xvals: np.ndarray) -> float:
    center = np.abs(xvals) <= 6
    flanks = (np.abs(xvals) >= 14) & (np.abs(xvals) <= 32)
    return float(np.mean(profile[flanks]) - np.mean(profile[center]))


def local_dip(profile: np.ndarray, xvals: np.ndarray) -> float:
    center = float(np.mean(profile[np.abs(xvals) <= 4]))
    left = float(np.max(profile[(xvals >= -25) & (xvals <= -8)]))
    right = float(np.max(profile[(xvals >= 8) & (xvals <= 25)]))
    return min(left, right) - center


def baseline_mean(profile: np.ndarray, xvals: np.ndarray) -> float:
    return float(np.mean(profile[np.abs(xvals) >= 45]))


def smooth_profile(profile: np.ndarray, window: int = 7) -> np.ndarray:
    if window <= 1:
        return profile
    pad = window // 2
    padded = np.pad(profile, pad, mode="edge")
    kernel = np.ones(window, dtype=float) / window
    return np.convolve(padded, kernel, mode="valid")



def condition_members(condition: str) -> list[str]:
    return [sample for sample, (sample_condition, _depth, _offset) in REPLICATES.items() if sample_condition == condition]


def condition_mean(profiles: dict[str, np.ndarray], condition: str) -> np.ndarray:
    return np.mean(np.vstack([profiles[sample] for sample in condition_members(condition)]), axis=0)


def condition_sd(profiles: dict[str, np.ndarray], condition: str) -> np.ndarray:
    members = condition_members(condition)
    return np.std(np.vstack([profiles[sample] for sample in members]), axis=0, ddof=1)


def summarize(profiles: dict[str, np.ndarray], xvals: np.ndarray) -> dict[str, float]:
    b_mean = condition_mean(profiles, "B cell")
    t_mean = condition_mean(profiles, "T cell")
    return {
        "baseline_diff": abs(baseline_mean(b_mean, xvals) - baseline_mean(t_mean, xvals)),
        "bcell_contrast": footprint_contrast(b_mean, xvals),
        "tcell_contrast": footprint_contrast(t_mean, xvals),
        "min_dip": min(local_dip(b_mean, xvals), local_dip(t_mean, xvals)),
        "bcell_sd": float(np.mean(condition_sd(profiles, "B cell"))),
        "tcell_sd": float(np.mean(condition_sd(profiles, "T cell"))),
    }


def combined_ylim(raw_profiles: dict[str, np.ndarray], norm_profiles: dict[str, np.ndarray]) -> tuple[float, float]:
    values = np.concatenate(list(raw_profiles.values()) + list(norm_profiles.values()))
    lo = float(np.nanpercentile(values, 1))
    hi = float(np.nanpercentile(values, 99))
    pad = max((hi - lo) * 0.08, 0.2)
    return lo - pad, hi + pad


def plot_panel(ax, xvals: np.ndarray, profiles: dict[str, np.ndarray], title: str, ylimits: tuple[float, float], show_ylabel: bool = False) -> None:
    for sample, profile in profiles.items():
        ax.plot(xvals, smooth_profile(profile), color=SAMPLE_COLORS[sample], linewidth=0.75, alpha=0.22, label=sample)
    for condition, color in CONDITION_COLORS.items():
        mean_profile = condition_mean(profiles, condition)
        sd_profile = condition_sd(profiles, condition)
        smooth_mean = smooth_profile(mean_profile)
        smooth_sd = smooth_profile(sd_profile)
        ax.fill_between(xvals, smooth_mean - smooth_sd, smooth_mean + smooth_sd, color=color, alpha=0.10, linewidth=0)
        ax.plot(
            xvals,
            smooth_mean,
            color=color,
            linewidth=1.55,
            alpha=0.96,
            label=f"{condition} mean",
        )
    ax.axvspan(-6, 6, color="0.9", alpha=0.72, zorder=-2)
    ax.axvline(0, color="0.45", linestyle="--", linewidth=0.85, zorder=-1)
    ax.set_ylim(*ylimits)
    ax.set_title(title, fontsize=10.5, fontweight="bold")
    ax.set_xlabel("bp from motif center", fontsize=8.5)
    if show_ylabel:
        ax.set_ylabel("Corrected cut-site signal", fontsize=8.5)
    ax.xaxis.set_major_formatter(ascii_tick_formatter())
    ax.yaxis.set_major_formatter(ascii_tick_formatter(decimals=1))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", linewidth=0.55)


def score_candidate(tfbs: Path, limit: int, flank: int) -> tuple[float, dict[str, float]] | None:
    try:
        raw_matrices = build_matrices(tfbs, limit, flank)
        norm_matrices = normalize_matrices(raw_matrices)
    except Exception:
        return None
    xvals = np.arange(-flank, flank)
    raw_profiles = profiles_from_matrices(raw_matrices)
    norm_profiles = profiles_from_matrices(norm_matrices)
    raw_summary = summarize(raw_profiles, xvals)
    norm_summary = summarize(norm_profiles, xvals)
    baseline_reduction = 1.0 - norm_summary["baseline_diff"] / (raw_summary["baseline_diff"] + 1e-9)
    min_contrast = min(
        raw_summary["bcell_contrast"],
        raw_summary["tcell_contrast"],
        norm_summary["bcell_contrast"],
        norm_summary["tcell_contrast"],
    )
    dip = max(raw_summary["min_dip"], norm_summary["min_dip"])
    score = 5.0 * min_contrast + 2.0 * dip + 4.0 * baseline_reduction + 0.1 * np.log(max(min(m.shape[0] for m in raw_matrices.values()), 1))
    return score, {
        "raw_baseline_diff": raw_summary["baseline_diff"],
        "norm_baseline_diff": norm_summary["baseline_diff"],
        "baseline_reduction": baseline_reduction,
        "min_contrast": min_contrast,
        "dip": dip,
        "n_sites": min(m.shape[0] for m in raw_matrices.values()),
    }


def plot_candidate_sheet(candidates: list[Path], output: Path, limit: int, flank: int, max_candidates: int) -> None:
    directions = read_direction_table()
    selected = []
    for tfbs in candidates:
        result = score_candidate(tfbs, limit, flank)
        if result is not None:
            selected.append((result[0], tfbs, result[1]))
        if len(selected) >= max_candidates:
            break
    if not selected:
        return
    xvals = np.arange(-flank, flank)
    fig, axes = plt.subplots(len(selected), 2, figsize=(7.6, max(1.72 * len(selected), 2.8)), sharex=True, constrained_layout=True)
    if len(selected) == 1:
        axes = np.asarray([axes])
    for row, (_score, tfbs, metrics) in enumerate(selected):
        raw_matrices = build_matrices(tfbs, limit, flank)
        norm_matrices = normalize_matrices(raw_matrices)
        raw_profiles = profiles_from_matrices(raw_matrices)
        norm_profiles = profiles_from_matrices(norm_matrices)
        ylimits = combined_ylim(raw_profiles, norm_profiles)
        tf_name = tf_name_from_path(tfbs)
        label = direction_label(tf_name, directions)
        plot_panel(axes[row, 0], xvals, raw_profiles, f"{tf_name} corrected", ylimits, show_ylabel=True)
        plot_panel(axes[row, 1], xvals, norm_profiles, f"{tf_name} normalized", ylimits)
        axes[row, 0].text(0.02, 0.94, f"base {metrics['raw_baseline_diff']:.1f}", transform=axes[row, 0].transAxes, va="top", fontsize=6.6)
        axes[row, 1].text(0.02, 0.94, f"base {metrics['norm_baseline_diff']:.1f}; dep {metrics['min_contrast']:.1f}", transform=axes[row, 1].transAxes, va="top", fontsize=6.6)
        if label:
            axes[row, 1].text(0.02, 0.82, label, transform=axes[row, 1].transAxes, va="top", fontsize=6.3)
    handles, labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=6.6, bbox_to_anchor=(0.5, -0.008))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tfbs", type=Path, nargs="*", default=None, help="One or more TFBS BED files for the main figure")
    parser.add_argument("--out-prefix", default="paper/manuscript/figures/normalization_effect")
    parser.add_argument("--flank", type=int, default=80)
    parser.add_argument("--limit", type=int, default=2000)
    parser.add_argument("--candidate-prefix", default="paper/manuscript/figures/normalization_effect_candidates")
    parser.add_argument("--candidate-max", type=int, default=16)
    args = parser.parse_args(argv)

    apply_pdf_style()
    directions = read_direction_table()
    tfbs_list = args.tfbs if args.tfbs else DEFAULT_MAIN_TFBS
    panels = []
    xvals = np.arange(-args.flank, args.flank)
    for tfbs in tfbs_list:
        raw_matrices = build_matrices(tfbs, args.limit, args.flank)
        norm_matrices = normalize_matrices(raw_matrices)
        raw_profiles = profiles_from_matrices(raw_matrices)
        norm_profiles = profiles_from_matrices(norm_matrices)
        raw_summary = summarize(raw_profiles, xvals)
        norm_summary = summarize(norm_profiles, xvals)
        panels.append((tfbs, raw_matrices, norm_matrices, raw_profiles, norm_profiles, raw_summary, norm_summary, combined_ylim(raw_profiles, norm_profiles)))

    fig, axes = plt.subplots(len(panels), 2, figsize=(7.7, 2.75 * len(panels)), sharex=True, constrained_layout=True)
    if len(panels) == 1:
        axes = np.asarray([axes])
    for row, (tfbs, _raw_matrices, _norm_matrices, raw_profiles, norm_profiles, raw_summary, norm_summary, ylimits) in enumerate(panels):
        tf_name = tf_name_from_path(tfbs)
        label = direction_label(tf_name, directions)
        plot_panel(axes[row, 0], xvals, raw_profiles, f"{tf_name}: corrected cut-site signal", ylimits, show_ylabel=True)
        plot_panel(axes[row, 1], xvals, norm_profiles, f"{tf_name}: after sample-quantile normalization", ylimits)
        axes[row, 0].text(0.02, 0.96, "modest synthetic depth offsets", transform=axes[row, 0].transAxes, va="top", fontsize=7.2)
        axes[row, 0].text(0.02, 0.87, f"Center depletion: {raw_summary['bcell_contrast']:.1f}/{raw_summary['tcell_contrast']:.1f}", transform=axes[row, 0].transAxes, va="top", fontsize=6.8)
        axes[row, 1].text(0.02, 0.96, "production sample-quantile path", transform=axes[row, 1].transAxes, va="top", fontsize=7.2)
        axes[row, 1].text(0.02, 0.87, f"Center depletion: {norm_summary['bcell_contrast']:.1f}/{norm_summary['tcell_contrast']:.1f}", transform=axes[row, 1].transAxes, va="top", fontsize=7.2)
        if label:
            axes[row, 1].text(0.02, 0.78, label, transform=axes[row, 1].transAxes, va="top", fontsize=6.8)

    handles, labels = axes[-1, 1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=6.8, bbox_to_anchor=(0.5, -0.025))
    fig.suptitle("Replicate-aware corrected cut-site aggregates with production quantile normalization", fontsize=11.6, fontweight="bold")

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    if args.candidate_prefix:
        plot_candidate_sheet(bindetect_prioritized_tfbs(), Path(args.candidate_prefix).with_suffix(".png"), args.limit, args.flank, args.candidate_max)

    stats_path = out_prefix.with_suffix(".tsv")
    with stats_path.open("w", encoding="utf-8") as handle:
        handle.write("normalization\ttfbs\tsample\tcondition\tdepth_factor\tn_sites\tmean_profile\tcenter_mean\tflank_mean\tbaseline_mean\tfootprint_contrast\tbindetect_direction\tbindetect_change\n")
        for tfbs, raw_matrices, norm_matrices, raw_profiles, norm_profiles, raw_summary, norm_summary, _ylimits in panels:
            tf_name = tf_name_from_path(tfbs)
            info = directions.get(tf_name, {})
            direction = str(info.get("direction", "NA"))
            change = str(info.get("norm_change", "NA"))
            for label, matrices, profiles in (("raw", raw_matrices, raw_profiles), ("sample_quantile", norm_matrices, norm_profiles)):
                for sample, profile in profiles.items():
                    center = np.abs(xvals) <= 6
                    flanks = (np.abs(xvals) >= 14) & (np.abs(xvals) <= 32)
                    row = [
                        label,
                        str(tfbs),
                        sample,
                        REPLICATES[sample][0],
                        f"{REPLICATES[sample][1]:.3f};offset={REPLICATES[sample][2]:.3f}",
                        str(matrices[sample].shape[0]),
                        f"{np.mean(profile):.6f}",
                        f"{np.mean(profile[center]):.6f}",
                        f"{np.mean(profile[flanks]):.6f}",
                        f"{baseline_mean(profile, xvals):.6f}",
                        f"{footprint_contrast(profile, xvals):.6f}",
                        direction,
                        change,
                    ]
                    handle.write("\t".join(row) + "\n")
            for label, summary in (("raw_summary", raw_summary), ("sample_quantile_summary", norm_summary)):
                handle.write("\t".join([
                    label,
                    str(tfbs),
                    "condition means",
                    "B cell vs T cell",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    f"{summary['baseline_diff']:.6f}",
                    f"{min(summary['bcell_contrast'], summary['tcell_contrast']):.6f}",
                    direction,
                    change,
                ]) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
