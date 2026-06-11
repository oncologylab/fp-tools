#!/usr/bin/env python
"""Build a real-data raw-vs-normalized replicate footprint figure."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyBigWig

from fp_tools.utils.normalization import fit_quantile_normalizers
from fp_tools.utils.plotting_style import apply_pdf_style, ascii_tick_formatter


SAMPLES = {
    "B cell rep1": ("B cell", Path("test_data/demo_Bcell_rep1_footprints.bw")),
    "B cell rep2": ("B cell", Path("test_data/demo_Bcell_rep2_footprints.bw")),
    "T cell rep1": ("T cell", Path("test_data/demo_Tcell_rep1_footprints.bw")),
    "T cell rep2": ("T cell", Path("test_data/demo_Tcell_rep2_footprints.bw")),
}
SAMPLE_COLORS = {
    "B cell rep1": "#4C9AC9",
    "B cell rep2": "#8BD3C7",
    "T cell rep1": "#E68645",
    "T cell rep2": "#F6C85F",
}
CONDITION_COLORS = {"B cell": "#08519C", "T cell": "#A63603"}
DEFAULT_TFBS = Path("test_data/annotated_tfbs/TFAP2A_Bcell_bound.bed")
FALLBACK_TFBS = [
    Path("test_data/annotated_tfbs/EBF1_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/KLF4_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/NRF1_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/CTCF_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/REST_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/YY1_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/PAX5_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/MAX_Bcell_bound.bed"),
    Path("test_data/annotated_tfbs/CREM_Bcell_bound.bed"),
]


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
            if values.shape[0] != 2 * flank:
                continue
            rows.append(np.nan_to_num(values, nan=0.0))
    finally:
        bw.close()
    if not rows:
        raise ValueError(f"No usable windows found in {bigwig}")
    return np.vstack(rows)


def choose_tfbs(requested: Path | None, limit: int, flank: int) -> tuple[Path, dict[str, np.ndarray]]:
    candidates = [requested] if requested is not None else [DEFAULT_TFBS, *FALLBACK_TFBS]
    last_error = None
    for tfbs in candidates:
        if tfbs is None or not tfbs.exists():
            continue
        try:
            sites = load_sites(tfbs, limit)
            matrices = {sample: read_matrix(path, sites, flank) for sample, (_cond, path) in SAMPLES.items()}
            if min(matrix.shape[0] for matrix in matrices.values()) >= 40:
                return tfbs, matrices
        except Exception as exc:  # pragma: no cover - fallback safety for local data drift
            last_error = exc
    raise RuntimeError(f"No usable TFBS set found for normalization figure: {last_error}")


def normalize_matrices(matrices: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    names = list(matrices)
    normalizers, _ = fit_quantile_normalizers([matrices[name].ravel() for name in names], names)
    return {
        name: normalizers[name].normalize(matrices[name].ravel()).reshape(matrices[name].shape)
        for name in names
    }


def profiles_from_matrices(matrices: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {name: matrix.mean(axis=0) for name, matrix in matrices.items()}


def condition_mean(profiles: dict[str, np.ndarray], condition: str) -> np.ndarray:
    members = [sample for sample, (cond, _path) in SAMPLES.items() if cond == condition]
    return np.mean(np.vstack([profiles[sample] for sample in members]), axis=0)


def footprint_contrast(profile: np.ndarray, xvals: np.ndarray) -> float:
    center = np.abs(xvals) <= 6
    flanks = (np.abs(xvals) >= 14) & (np.abs(xvals) <= 32)
    return float(np.mean(profile[flanks]) - np.mean(profile[center]))


def baseline_mean(profile: np.ndarray, xvals: np.ndarray) -> float:
    return float(np.mean(profile[np.abs(xvals) >= 45]))


def summarize(profiles: dict[str, np.ndarray], xvals: np.ndarray) -> dict[str, float]:
    means = {condition: condition_mean(profiles, condition) for condition in CONDITION_COLORS}
    return {
        "baseline_diff": abs(baseline_mean(means["B cell"], xvals) - baseline_mean(means["T cell"], xvals)),
        "shape_diff": abs(footprint_contrast(means["B cell"], xvals) - footprint_contrast(means["T cell"], xvals)),
        "bcell_contrast": footprint_contrast(means["B cell"], xvals),
        "tcell_contrast": footprint_contrast(means["T cell"], xvals),
    }


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


def build_profiles(tfbs: Path, limit: int, flank: int) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray]]:
    sites = load_sites(tfbs, limit)
    raw_matrices = {sample: read_matrix(path, sites, flank) for sample, (_cond, path) in SAMPLES.items()}
    if min(matrix.shape[0] for matrix in raw_matrices.values()) < 40:
        raise ValueError(f"Too few usable windows for {tfbs}")
    norm_matrices = normalize_matrices(raw_matrices)
    return raw_matrices, norm_matrices, profiles_from_matrices(raw_matrices), profiles_from_matrices(norm_matrices)


def combined_ylim(raw_profiles: dict[str, np.ndarray], norm_profiles: dict[str, np.ndarray]) -> tuple[float, float]:
    values = np.concatenate(list(raw_profiles.values()) + list(norm_profiles.values()))
    lo = float(np.nanpercentile(values, 1))
    hi = float(np.nanpercentile(values, 99))
    pad = max((hi - lo) * 0.08, 0.2)
    return lo - pad, hi + pad


def plot_panel(ax, xvals: np.ndarray, profiles: dict[str, np.ndarray], title: str, ylimits: tuple[float, float], show_ylabel: bool = False) -> None:
    for sample, profile in profiles.items():
        ax.plot(xvals, profile, color=SAMPLE_COLORS[sample], linewidth=1.05, alpha=0.82, label=sample)
    for condition, color in CONDITION_COLORS.items():
        ax.plot(
            xvals,
            condition_mean(profiles, condition),
            color=color,
            linewidth=1.05,
            linestyle=(0, (2.0, 2.0)),
            alpha=0.95,
            label=f"{condition} mean",
        )
    ax.axvspan(-6, 6, color="0.9", alpha=0.7, zorder=-2)
    ax.axvline(0, color="0.45", linestyle="--", linewidth=0.85, zorder=-1)
    ax.set_ylim(*ylimits)
    ax.set_title(title, fontsize=10.5, fontweight="bold")
    ax.set_xlabel("bp from motif center", fontsize=8.5)
    if show_ylabel:
        ax.set_ylabel("Footprint score", fontsize=8.5)
    ax.xaxis.set_major_formatter(ascii_tick_formatter())
    ax.yaxis.set_major_formatter(ascii_tick_formatter(decimals=1))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", linewidth=0.55)


def score_candidate(tfbs: Path, limit: int, flank: int) -> tuple[float, dict[str, float]] | None:
    try:
        _raw_matrices, _norm_matrices, raw_profiles, norm_profiles = build_profiles(tfbs, limit, flank)
    except Exception:
        return None
    xvals = np.arange(-flank, flank)
    raw_summary = summarize(raw_profiles, xvals)
    norm_summary = summarize(norm_profiles, xvals)
    baseline_reduction = 1.0 - norm_summary["baseline_diff"] / (raw_summary["baseline_diff"] + 1e-9)
    depletion_score = min(
        raw_summary["bcell_contrast"],
        raw_summary["tcell_contrast"],
        norm_summary["bcell_contrast"],
        norm_summary["tcell_contrast"],
    )
    score = 5.0 * depletion_score + 2.0 * baseline_reduction + 0.1 * np.log(max(raw_summary["baseline_diff"], 1.0))
    metrics = {
        "baseline_reduction": baseline_reduction,
        "raw_baseline_diff": raw_summary["baseline_diff"],
        "norm_baseline_diff": norm_summary["baseline_diff"],
        "min_depletion": depletion_score,
    }
    return score, metrics


def plot_candidate_sheet(candidates: list[Path], output: Path, limit: int, flank: int, max_candidates: int) -> None:
    scored = []
    for tfbs in candidates:
        result = score_candidate(tfbs, limit, flank)
        if result is None:
            continue
        scored.append((result[0], tfbs, result[1]))
    scored.sort(reverse=True, key=lambda row: row[0])
    selected = scored[:max_candidates]
    if not selected:
        return
    xvals = np.arange(-flank, flank)
    nrows = len(selected)
    fig, axes = plt.subplots(nrows, 2, figsize=(7.4, max(1.65 * nrows, 2.5)), sharex=True, constrained_layout=True)
    if nrows == 1:
        axes = np.asarray([axes])
    for row, (_score, tfbs, metrics) in enumerate(selected):
        _raw_matrices, _norm_matrices, raw_profiles, norm_profiles = build_profiles(tfbs, limit, flank)
        ylimits = combined_ylim(raw_profiles, norm_profiles)
        tf_name = tfbs.name.split("_")[0]
        plot_panel(axes[row, 0], xvals, raw_profiles, f"{tf_name} raw", ylimits, show_ylabel=True)
        plot_panel(axes[row, 1], xvals, norm_profiles, f"{tf_name} normalized", ylimits)
        axes[row, 0].text(0.02, 0.94, f"base {metrics['raw_baseline_diff']:.1f}", transform=axes[row, 0].transAxes, va="top", fontsize=6.6)
        axes[row, 1].text(0.02, 0.94, f"base {metrics['norm_baseline_diff']:.1f}; dep {metrics['min_depletion']:.2f}", transform=axes[row, 1].transAxes, va="top", fontsize=6.6)
    handles, labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=6.8, bbox_to_anchor=(0.5, -0.01))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tfbs", type=Path, default=DEFAULT_TFBS)
    parser.add_argument("--out-prefix", default="paper/manuscript/figures/normalization_effect")
    parser.add_argument("--flank", type=int, default=80)
    parser.add_argument("--limit", type=int, default=2000)
    parser.add_argument("--candidate-prefix", default="paper/manuscript/figures/normalization_effect_candidates")
    parser.add_argument("--candidate-max", type=int, default=12)
    args = parser.parse_args(argv)

    apply_pdf_style()
    tfbs, raw_matrices = choose_tfbs(args.tfbs, args.limit, args.flank)
    norm_matrices = normalize_matrices(raw_matrices)
    raw_profiles = profiles_from_matrices(raw_matrices)
    norm_profiles = profiles_from_matrices(norm_matrices)
    xvals = np.arange(-args.flank, args.flank)
    raw_summary = summarize(raw_profiles, xvals)
    norm_summary = summarize(norm_profiles, xvals)

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.1), sharex=True, constrained_layout=True)
    tf_name = tfbs.name.split("_")[0]
    ylimits = combined_ylim(raw_profiles, norm_profiles)
    plot_panel(axes[0], xvals, raw_profiles, "Observed replicate profiles", ylimits, show_ylabel=True)
    plot_panel(axes[1], xvals, norm_profiles, "After sample-quantile normalization", ylimits)
    axes[0].text(
        0.02,
        0.96,
        f"Baseline difference: {raw_summary['baseline_diff']:.1f}",
        transform=axes[0].transAxes,
        va="top",
        fontsize=7.6,
    )
    axes[1].text(
        0.02,
        0.96,
        f"Baseline difference: {norm_summary['baseline_diff']:.1f}",
        transform=axes[1].transAxes,
        va="top",
        fontsize=7.6,
    )
    axes[1].text(
        0.02,
        0.88,
        f"Center depletion retained: {norm_summary['bcell_contrast']:.2f}/{norm_summary['tcell_contrast']:.2f}",
        transform=axes[1].transAxes,
        va="top",
        fontsize=7.6,
    )
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=7.2, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle(f"Real {tf_name} replicate footprints before and after quantile normalization", fontsize=12.2, fontweight="bold")

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    candidate_prefix = Path(args.candidate_prefix)
    if candidate_prefix:
        plot_candidate_sheet(collect_candidate_tfbs(None), candidate_prefix.with_suffix(".png"), args.limit, args.flank, args.candidate_max)

    stats_path = out_prefix.with_suffix(".tsv")
    with stats_path.open("w", encoding="utf-8") as handle:
        handle.write("normalization\ttfbs\tsample\tcondition\tn_sites\tmean_profile\tcenter_mean\tflank_mean\tbaseline_mean\tfootprint_contrast\n")
        for label, matrices, profiles in (("raw", raw_matrices, raw_profiles), ("sample_quantile", norm_matrices, norm_profiles)):
            for sample, profile in profiles.items():
                center = np.abs(xvals) <= 6
                flanks = (np.abs(xvals) >= 14) & (np.abs(xvals) <= 32)
                row = [
                    label,
                    str(tfbs),
                    sample,
                    SAMPLES[sample][0],
                    str(matrices[sample].shape[0]),
                    f"{np.mean(profile):.6f}",
                    f"{np.mean(profile[center]):.6f}",
                    f"{np.mean(profile[flanks]):.6f}",
                    f"{baseline_mean(profile, xvals):.6f}",
                    f"{footprint_contrast(profile, xvals):.6f}",
                ]
                handle.write("\t".join(row) + "\n")
        for label, summary in (("raw_summary", raw_summary), ("sample_quantile_summary", norm_summary)):
            handle.write("\t".join([label, str(tfbs), "condition_means", "B cell vs T cell", "NA", "NA", "NA", "NA", f"{summary['baseline_diff']:.6f}", f"{summary['shape_diff']:.6f}"]) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
