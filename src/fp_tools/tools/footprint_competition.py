#!/usr/bin/env python
"""Competition-aware overlap decomposition of multiscale footprint signal.

Reads a multiscale NPZ sidecar produced by ``score-footprints
--output-multiscale-npz`` and decomposes each region's depletion signal into a
short TF-scale footprint band and a wide nucleosome-scale band. Where both bands
are simultaneously strong the components compete; this tool partitions the
positionwise signal into TF-only, nucleosome-only, and shared (competing)
fractions so overlapping footprint components can be attributed explicitly.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fp_tools.utils.multiscale import load_multiscale_npz


def select_band(scales: np.ndarray, low: float, high: float) -> np.ndarray:
    """Return a boolean mask of scales within the inclusive [low, high] band."""

    return (scales >= low) & (scales <= high)


def _band_profile(region_tensor: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Collapse the selected scale rows to a positionwise profile by max."""

    if not mask.any():
        return np.zeros(region_tensor.shape[1], dtype=float)
    band = region_tensor[mask, :]
    profile = np.nanmax(band, axis=0)
    return np.nan_to_num(profile, nan=0.0)


def decompose_region(
    region_tensor: np.ndarray,
    scales: np.ndarray,
    tf_band: tuple[float, float],
    nuc_band: tuple[float, float],
) -> dict[str, float]:
    """Decompose one region tensor into TF-only, nucleosome-only, and shared signal."""

    tf_mask = select_band(scales, tf_band[0], tf_band[1])
    nuc_mask = select_band(scales, nuc_band[0], nuc_band[1])
    tf_profile = np.clip(_band_profile(region_tensor, tf_mask), 0.0, None)
    nuc_profile = np.clip(_band_profile(region_tensor, nuc_mask), 0.0, None)

    shared = np.minimum(tf_profile, nuc_profile)
    tf_only = tf_profile - shared
    nuc_only = nuc_profile - shared

    tf_only_auc = float(tf_only.sum())
    nuc_only_auc = float(nuc_only.sum())
    shared_auc = float(shared.sum())
    total = tf_only_auc + nuc_only_auc + shared_auc

    competition_index = shared_auc / total if total > 0 else 0.0
    if total <= 0:
        dominant = "none"
    elif shared_auc >= max(tf_only_auc, nuc_only_auc):
        dominant = "competing"
    elif tf_only_auc >= nuc_only_auc:
        dominant = "tf"
    else:
        dominant = "nucleosome"

    length = region_tensor.shape[1]
    peak_pos = int(np.argmax(shared)) if length and shared.any() else -1

    return {
        "length": int(length),
        "tf_scale_score": float(tf_profile.mean()) if length else 0.0,
        "nucleosome_scale_score": float(nuc_profile.mean()) if length else 0.0,
        "tf_only_auc": tf_only_auc,
        "nucleosome_only_auc": nuc_only_auc,
        "shared_auc": shared_auc,
        "total_auc": total,
        "competition_index": competition_index,
        "dominant_component": dominant,
        "competition_peak_offset": peak_pos - int(length // 2) if peak_pos >= 0 else 0,
    }


def build_competition_report(
    npz: str | Path,
    output: str | Path,
    summary_output: str | Path | None = None,
    figure_output: str | Path | None = None,
    tf_band: tuple[float, float] = (3.0, 30.0),
    nuc_band: tuple[float, float] = (120.0, 200.0),
) -> "tuple":
    """Build per-region competition decomposition and an optional summary/figure."""

    import pandas as pd

    if tf_band[1] >= nuc_band[0]:
        raise ValueError(
            f"TF-scale band {tf_band} must end below the nucleosome-scale band {nuc_band}."
        )

    data = load_multiscale_npz(str(npz))
    tensor = data["tensor"]
    scales = data["scales"].astype(float)
    offsets = data["offsets"].astype(int)
    chroms = data["chroms"]
    starts = data["starts"].astype(int)
    ends = data["ends"].astype(int)

    if not select_band(scales, *tf_band).any():
        raise ValueError(f"No scales fall in the TF-scale band {tf_band}; scales={scales.tolist()}.")
    if not select_band(scales, *nuc_band).any():
        raise ValueError(
            f"No scales fall in the nucleosome-scale band {nuc_band}; scales={scales.tolist()}."
        )

    rows = []
    for idx in range(len(offsets) - 1):
        region_tensor = tensor[:, offsets[idx] : offsets[idx + 1]]
        stats = decompose_region(region_tensor, scales, tf_band, nuc_band)
        rows.append(
            {
                "chrom": str(chroms[idx]),
                "start": int(starts[idx]),
                "end": int(ends[idx]),
                **stats,
            }
        )

    report = pd.DataFrame(rows)
    summary = pd.DataFrame(
        [
            {
                "n_regions": len(report),
                "tf_dominant": int((report["dominant_component"] == "tf").sum()) if len(report) else 0,
                "nucleosome_dominant": int((report["dominant_component"] == "nucleosome").sum())
                if len(report)
                else 0,
                "competing": int((report["dominant_component"] == "competing").sum()) if len(report) else 0,
                "mean_competition_index": float(report["competition_index"].mean()) if len(report) else 0.0,
                "median_competition_index": float(report["competition_index"].median()) if len(report) else 0.0,
                "tf_band_low": tf_band[0],
                "tf_band_high": tf_band[1],
                "nucleosome_band_low": nuc_band[0],
                "nucleosome_band_high": nuc_band[1],
            }
        ]
    )

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out, sep="\t", index=False)
    if summary_output is not None:
        summary_out = Path(summary_output)
        summary_out.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(summary_out, sep="\t", index=False)
    if figure_output is not None:
        plot_competition_report(report, figure_output)
    return report, summary


def plot_competition_report(report, output: str | Path) -> Path:
    """Plot competition index distribution and TF vs nucleosome component scores."""

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), constrained_layout=True)

    ax = axes[0]
    if len(report):
        ax.hist(report["competition_index"], bins=20, color="tab:purple", alpha=0.85)
    ax.set_xlabel("competition index (shared / total)")
    ax.set_ylabel("regions")
    ax.set_title("Footprint competition distribution")
    ax.grid(alpha=0.25)

    ax = axes[1]
    colors = {
        "tf": "tab:blue",
        "nucleosome": "tab:orange",
        "competing": "tab:purple",
        "none": "0.6",
    }
    if len(report):
        point_colors = [colors.get(value, "0.6") for value in report["dominant_component"]]
        ax.scatter(report["tf_scale_score"], report["nucleosome_scale_score"], s=18, alpha=0.8, c=point_colors)
        lim = max(report["tf_scale_score"].max(), report["nucleosome_scale_score"].max(), 1e-9)
        ax.plot([0, lim], [0, lim], color="0.5", linewidth=0.8, linestyle="--")
    ax.set_xlabel("TF-scale score")
    ax.set_ylabel("nucleosome-scale score")
    ax.set_title("TF vs nucleosome component")
    ax.grid(alpha=0.25)

    fig.savefig(out, dpi=300 if out.suffix.lower() == ".png" else None)
    plt.close(fig)
    return out


def _parse_band(value: str) -> tuple[float, float]:
    parts = value.replace(":", ",").split(",")
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("Band must be given as low,high (e.g. 3,30).")
    low, high = float(parts[0]), float(parts[1])
    if low > high:
        raise argparse.ArgumentTypeError("Band low must be <= high.")
    return low, high


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npz", required=True, help="Multiscale NPZ sidecar from score-footprints.")
    parser.add_argument("--out", required=True, help="Per-region competition decomposition TSV.")
    parser.add_argument("--summary-out", help="Optional dataset-level summary TSV.")
    parser.add_argument("--figure-out", help="Optional PDF/SVG/PNG figure.")
    parser.add_argument(
        "--tf-band",
        type=_parse_band,
        default=(3.0, 30.0),
        help="TF-scale band as low,high in bp (default 3,30).",
    )
    parser.add_argument(
        "--nucleosome-band",
        type=_parse_band,
        default=(120.0, 200.0),
        help="Nucleosome-scale band as low,high in bp (default 120,200).",
    )
    args = parser.parse_args(argv)

    report, summary = build_competition_report(
        args.npz,
        args.out,
        summary_output=args.summary_out,
        figure_output=args.figure_out,
        tf_band=args.tf_band,
        nuc_band=args.nucleosome_band,
    )
    print(f"wrote {len(report)} competition decomposition rows to {args.out}")
    if args.summary_out:
        print(f"wrote competition summary to {args.summary_out}")
    if args.figure_out:
        print(f"wrote competition figure to {args.figure_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
