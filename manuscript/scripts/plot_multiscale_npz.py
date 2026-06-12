#!/usr/bin/env python3
"""Create a multi-panel figure from fp-tools multiscale NPZ output."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from fp_tools.utils.multiscale import aggregate_multiscale_tensor, load_multiscale_npz


def _position_axis(width: int) -> np.ndarray:
    return np.arange(width) - int(width // 2)


def plot_multiscale_npz(npz_path: str | Path, out_prefix: str | Path, title: str = "fp-tools multiscale footprint summary") -> list[Path]:
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    data = load_multiscale_npz(str(npz_path))
    aggregate = aggregate_multiscale_tensor(data, align="center")
    scales = data["scales"].astype(int)
    x = _position_axis(aggregate.shape[1])

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2), constrained_layout=True)
    fig.suptitle(title)

    heat_ax = axes[0, 0]
    im = heat_ax.imshow(
        aggregate,
        aspect="auto",
        cmap="viridis",
        interpolation="nearest",
        extent=[x[0] if len(x) else 0, x[-1] if len(x) else 0, len(scales) - 0.5, -0.5],
    )
    heat_ax.set_title("Scale-by-position aggregate")
    heat_ax.set_xlabel("bp from region center")
    heat_ax.set_ylabel("scale (bp)")
    heat_ax.set_yticks(range(len(scales)))
    heat_ax.set_yticklabels([str(scale) for scale in scales])
    fig.colorbar(im, ax=heat_ax, fraction=0.046, pad=0.04, label="depletion score")

    profile_ax = axes[0, 1]
    if aggregate.size:
        profile_ax.plot(x, np.nanmax(aggregate, axis=0), color="black", linewidth=1.4, label="max across scales")
        profile_ax.plot(x, np.nanmean(aggregate, axis=0), color="tab:blue", linewidth=1.0, label="mean across scales")
    profile_ax.set_title("Collapsed profile")
    profile_ax.set_xlabel("bp from region center")
    profile_ax.set_ylabel("depletion score")
    profile_ax.grid(alpha=0.25)
    profile_ax.legend(frameon=False, fontsize=8)

    scale_ax = axes[1, 0]
    scale_means = np.nanmean(aggregate, axis=1) if aggregate.size else np.zeros(len(scales))
    scale_ax.bar(range(len(scales)), scale_means, color="tab:green")
    scale_ax.set_title("Mean score by scale")
    scale_ax.set_xlabel("scale (bp)")
    scale_ax.set_ylabel("mean depletion score")
    scale_ax.set_xticks(range(len(scales)))
    scale_ax.set_xticklabels([str(scale) for scale in scales], rotation=45, ha="right")
    scale_ax.grid(axis="y", alpha=0.25)

    meta_ax = axes[1, 1]
    meta_ax.axis("off")
    lengths = np.diff(data["offsets"])
    lines = [
        f"NPZ: {Path(npz_path).name}",
        f"Regions: {len(lengths)}",
        f"Scales: {', '.join(str(scale) for scale in scales)}",
        f"Positions: {int(np.sum(lengths))}",
        f"Median region width: {float(np.median(lengths)):.1f} bp" if len(lengths) else "Median region width: NA",
        f"Summary bigWig method: {data['summary_method'].item()}",
    ]
    meta_ax.text(0, 1, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=9)
    meta_ax.set_title("Tensor metadata")

    outputs = []
    for suffix in ("pdf", "svg", "png"):
        path = out_prefix.with_suffix(f".{suffix}")
        fig.savefig(path, dpi=300 if suffix == "png" else None)
        outputs.append(path)
    plt.close(fig)
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--multiscale-npz", required=True, help="NPZ sidecar from call-footprints --output-multiscale-npz.")
    parser.add_argument("--out-prefix", default="manuscript/figures/figure_multiscale_summary")
    parser.add_argument("--title", default="fp-tools multiscale footprint summary")
    args = parser.parse_args()

    outputs = plot_multiscale_npz(args.multiscale_npz, args.out_prefix, title=args.title)
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
