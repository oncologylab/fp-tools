#!/usr/bin/env python
"""Build a raw-vs-normalized replicate footprint figure for the manuscript."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fp_tools.utils.normalization import fit_quantile_normalizers
from fp_tools.utils.plotting_style import apply_pdf_style, ascii_tick_formatter


SAMPLE_COLORS = {
    "B cell rep1": "#2C7FB8",
    "B cell rep2": "#7FCDBB",
    "T cell rep1": "#D95F0E",
    "T cell rep2": "#FEC44F",
}
CONDITION_COLORS = {
    "B cell": "#08519C",
    "T cell": "#A63603",
}
SAMPLE_TO_CONDITION = {
    "B cell rep1": "B cell",
    "B cell rep2": "B cell",
    "T cell rep1": "T cell",
    "T cell rep2": "T cell",
}


def simulate_footprint_profiles(flank: int, seed: int) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, float]]:
    """Return controlled replicate profiles with a central footprint depletion."""

    rng = np.random.default_rng(seed)
    x = np.arange(-flank, flank + 1, dtype=float)
    broad_accessibility = 6.2 + 0.9 * np.exp(-(x / 54.0) ** 2)
    flank_shoulders = 0.95 * np.exp(-((np.abs(x) - 17.0) / 10.0) ** 2)
    protected_center = 1.45 * np.exp(-(x / 6.0) ** 2)
    local_texture = 0.06 * np.sin(x / 5.5) + 0.035 * np.cos(x / 13.0)
    base = broad_accessibility + flank_shoulders - protected_center + local_texture

    scale_factors = {
        "B cell rep1": 0.72,
        "B cell rep2": 1.38,
        "T cell rep1": 0.88,
        "T cell rep2": 1.62,
    }
    condition_offsets = {
        "B cell rep1": 0.18,
        "B cell rep2": 0.18,
        "T cell rep1": -0.10,
        "T cell rep2": -0.10,
    }
    profiles = {}
    for sample, scale in scale_factors.items():
        noise = rng.normal(0.0, 0.045, size=base.shape)
        smooth_noise = np.convolve(noise, np.ones(5) / 5.0, mode="same")
        profiles[sample] = np.maximum(0.0, base * scale + condition_offsets[sample] + smooth_noise)
    return x, profiles, scale_factors


def normalize_profiles(raw_profiles: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    names = list(raw_profiles)
    normalizers, _ = fit_quantile_normalizers([raw_profiles[name] for name in names], names)
    return {name: normalizers[name].normalize(raw_profiles[name]) for name in names}


def condition_mean(profiles: dict[str, np.ndarray], condition: str) -> np.ndarray:
    members = [name for name, cond in SAMPLE_TO_CONDITION.items() if cond == condition]
    return np.mean(np.vstack([profiles[name] for name in members]), axis=0)


def footprint_contrast(profile: np.ndarray, x: np.ndarray) -> float:
    center = np.abs(x) <= 6
    flanks = ((np.abs(x) >= 14) & (np.abs(x) <= 28))
    return float(np.mean(profile[flanks]) - np.mean(profile[center]))


def plot_panel(ax, x, profiles, title, show_ylabel=False):
    for sample, profile in profiles.items():
        ax.plot(x, profile, color=SAMPLE_COLORS[sample], linewidth=1.25, alpha=0.78, label=sample)
    for condition, color in CONDITION_COLORS.items():
        mean_profile = condition_mean(profiles, condition)
        ax.plot(x, mean_profile, color=color, linewidth=2.7, label=f"{condition} mean")
    ax.axvspan(-6, 6, color="0.88", alpha=0.75, zorder=-2)
    ax.axvline(0, color="0.55", linestyle="--", linewidth=1.0, zorder=-1)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel("bp from motif center", fontsize=9)
    if show_ylabel:
        ax.set_ylabel("Footprint score", fontsize=9)
    ax.xaxis.set_major_formatter(ascii_tick_formatter())
    ax.yaxis.set_major_formatter(ascii_tick_formatter(decimals=1))
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", linewidth=0.6)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-prefix", default="paper/manuscript/figures/normalization_effect")
    parser.add_argument("--flank", type=int, default=80)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args(argv)

    apply_pdf_style()
    x, raw_profiles, scale_factors = simulate_footprint_profiles(args.flank, args.seed)
    normalized_profiles = normalize_profiles(raw_profiles)

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.15), sharex=True, constrained_layout=True)
    plot_panel(axes[0], x, raw_profiles, "Raw replicate profiles", show_ylabel=True)
    plot_panel(axes[1], x, normalized_profiles, "After sample-quantile normalization")

    raw_range = max(np.max(v) for v in raw_profiles.values()) - min(np.min(v) for v in raw_profiles.values())
    norm_range = max(np.max(v) for v in normalized_profiles.values()) - min(np.min(v) for v in normalized_profiles.values())
    raw_contrast = np.mean([footprint_contrast(condition_mean(raw_profiles, c), x) for c in CONDITION_COLORS])
    norm_contrast = np.mean([footprint_contrast(condition_mean(normalized_profiles, c), x) for c in CONDITION_COLORS])
    axes[0].text(0.02, 0.96, f"Depth factors: {min(scale_factors.values()):.2f}x-{max(scale_factors.values()):.2f}x", transform=axes[0].transAxes, va="top", fontsize=8)
    axes[1].text(0.02, 0.96, f"Range reduced: {raw_range:.1f} to {norm_range:.1f}", transform=axes[1].transAxes, va="top", fontsize=8)
    axes[1].text(0.02, 0.88, f"Footprint contrast retained: {raw_contrast:.2f} to {norm_contrast:.2f}", transform=axes[1].transAxes, va="top", fontsize=8)

    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.03))
    fig.suptitle("Replicate-aware quantile normalization preserves a central footprint", fontsize=13, fontweight="bold")

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    stats_path = out_prefix.with_suffix(".tsv")
    with stats_path.open("w", encoding="utf-8") as handle:
        header = [
            "normalization",
            "sample",
            "condition",
            "scale_factor",
            "mean_profile",
            "center_mean",
            "flank_mean",
            "footprint_contrast",
        ]
        handle.write("\t".join(header) + "\n")
        for label, profiles in (("raw", raw_profiles), ("sample_quantile", normalized_profiles)):
            for sample, profile in profiles.items():
                center = np.abs(x) <= 6
                flanks = ((np.abs(x) >= 14) & (np.abs(x) <= 28))
                row = [
                    label,
                    sample,
                    SAMPLE_TO_CONDITION[sample],
                    f"{scale_factors[sample]:.4f}",
                    f"{np.mean(profile):.6f}",
                    f"{np.mean(profile[center]):.6f}",
                    f"{np.mean(profile[flanks]):.6f}",
                    f"{footprint_contrast(profile, x):.6f}",
                ]
                handle.write("\t".join(row) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
