#!/usr/bin/env python
"""Plot real 2-vs-2 differential reports before and after normalization."""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from figure_style import apply_style, bold_all_text  # noqa: E402


DEFAULT_NONE_REPORT = Path(
    "data/public/processed/buenrostro_atac_replicates/fp_tools/"
    "detect_tf_binding_jaspar2026_vertebrates_norm_none/diff_footprints_Bcell_Tcell.html"
)
DEFAULT_SAMPLE_QUANTILE_REPORT = Path(
    "data/public/processed/buenrostro_atac_replicates/fp_tools/"
    "detect_tf_binding_jaspar2026_vertebrates_norm_sample_quantile/diff_footprints_Bcell_Tcell.html"
)
DEFAULT_MOTIFS = [
    "BACH2_MA1101.3",
    "IRF4_MA1419.2",
]
CONDITION_COLORS = {
    "Bcell": "#1f77b4",
    "Tcell": "#d62728",
}
GROUP_COLORS = {
    "Bcell_up": "#2563eb",
    "Tcell_up": "#dc2626",
    "n.s.": "#8a94a6",
}
SAMPLE_STYLES = {
    "Bcell_rep1": (0, (1.0, 1.2)),
    "Bcell_rep2": (0, (3.0, 1.4)),
    "Tcell_rep1": (0, (1.0, 1.2)),
    "Tcell_rep2": (0, (3.0, 1.4)),
}


def load_report_payload(report_html: Path) -> dict:
    text = report_html.read_text(encoding="utf-8")
    match = re.search(r'reportPayloadB64="([^"]+)"', text)
    if match is None:
        raise ValueError(f"Could not find reportPayloadB64 in {report_html}")
    return json.loads(gzip.decompress(base64.b64decode(match.group(1))).decode("utf-8"))


def motif_label(motif: dict) -> str:
    motif_id = motif.get("motif_id", "")
    name = motif.get("name", motif.get("prefix", "motif"))
    return f"{name} ({motif_id})" if motif_id else name


def center_minus_flank(profile: np.ndarray, xvals: np.ndarray, center_bp: int = 10, edge_bp: int = 20) -> float:
    center = profile[np.abs(xvals) <= center_bp]
    flanks = np.concatenate([profile[xvals <= xvals.min() + edge_bp], profile[xvals >= xvals.max() - edge_bp]])
    if center.size == 0 or flanks.size == 0:
        return float("nan")
    return float(np.nanmean(center) - np.nanmean(flanks))


def get_motif(payload: dict, prefix: str) -> dict:
    motifs = payload.get("aggregate", {}).get("motifs", [])
    motif = next((item for item in motifs if item.get("prefix") == prefix), None)
    if motif is None:
        raise ValueError(f"Motif {prefix} missing from aggregate payload")
    return motif


def html_like_ylim(motif: dict) -> tuple[float, float]:
    values = []
    for condition in motif.get("conditions", []):
        values.extend(condition.get("profile", []))
        for sample in condition.get("samples", []):
            values.extend(sample.get("profile", []))
    finite = np.asarray([value for value in values if np.isfinite(value)], dtype=float)
    if finite.size == 0:
        return (-1.0, 1.0)
    ymin = min(float(np.nanmin(finite)), 0.0)
    ymax = max(float(np.nanmax(finite)), 1e-9)
    min_pad = max(1e-4, abs(ymax) * 0.05, abs(ymin) * 0.05)
    pad = max((ymax - ymin if ymax != ymin else 1.0) * 0.22, min_pad)
    return ymin - pad, ymax + pad


def volcano_limits(*payloads: dict) -> tuple[tuple[float, float], tuple[float, float]]:
    changes = []
    neglog = []
    for payload in payloads:
        for point in payload.get("points", []):
            changes.append(float(point.get("change", 0.0)))
            neglog.append(float(point.get("neglog10p", 0.0)))
    if not changes:
        return (-1.0, 1.0), (0.0, 1.0)
    x_abs = max(1.0, abs(min(changes)), abs(max(changes))) * 1.1
    y_max = max(1.0, max(neglog)) * 1.08
    return (-x_abs, x_abs), (0.0, y_max)


def plot_volcano(ax, payload: dict, title: str, xlim: tuple[float, float], ylim: tuple[float, float]) -> dict[str, object]:
    points = pd.DataFrame(payload.get("points", []))
    if points.empty:
        raise ValueError(f"No volcano points found for {title}")
    points["change"] = pd.to_numeric(points["change"], errors="coerce")
    points["neglog10p"] = pd.to_numeric(points["neglog10p"], errors="coerce")
    points = points.dropna(subset=["change", "neglog10p", "group"])

    for group in ["n.s.", "Bcell_up", "Tcell_up"]:
        subset = points[points["group"] == group]
        if subset.empty:
            continue
        ax.scatter(
            subset["change"],
            subset["neglog10p"],
            s=8 if group == "n.s." else 13,
            c=GROUP_COLORS.get(group, "#777777"),
            alpha=0.32 if group == "n.s." else 0.76,
            linewidths=0,
            label=group,
            rasterized=True,
        )

    ax.axvline(0, color="0.35", linewidth=0.75)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.grid(color="0.90", linewidth=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(title, fontsize=9.0, pad=4)
    ax.set_xlabel("Differential footprint score")
    ax.set_ylabel("-log10(P)")

    counts = points["group"].value_counts()
    b_up = int(counts.get("Bcell_up", 0))
    t_up = int(counts.get("Tcell_up", 0))
    median_change = float(points["change"].median())
    ax.tick_params(axis="both", labelsize=7.0, pad=1.5)
    bold_all_text(ax)
    return {
        "normalization": title,
        "n_motifs": len(points),
        "bcell_up": b_up,
        "tcell_up": t_up,
        "median_change": median_change,
        "mean_change": float(points["change"].mean()),
    }


def plot_motif_panel(ax, xvals: np.ndarray, motif: dict, mode: str, ylim: tuple[float, float], ylabel: str | None) -> list[dict[str, object]]:
    rows = []
    for condition in motif.get("conditions", []):
        condition_name = str(condition["name"])
        color = CONDITION_COLORS.get(condition_name, "#555555")
        mean_profile = np.asarray(condition["profile"], dtype=float)
        ax.plot(
            xvals,
            mean_profile,
            color=color,
            linewidth=1.85,
            alpha=0.98,
            label=f"{condition_name} mean",
            zorder=2,
        )
        rows.append(
            {
                "normalization": mode,
                "motif": motif.get("name", ""),
                "motif_id": motif.get("motif_id", ""),
                "output_prefix": motif.get("prefix", ""),
                "profile_type": "condition_mean",
                "condition": condition_name,
                "sample": condition_name,
                "n_sites": motif.get("n_sites", ""),
                "center_minus_flank": center_minus_flank(mean_profile, xvals),
            }
        )
    ax.axvline(0, color="0.35", linewidth=0.75, alpha=0.85)
    ax.axhline(0, color="0.82", linewidth=0.55, zorder=0)
    ax.set_ylim(*ylim)
    ax.set_title(mode, fontsize=8.4, pad=3)
    if ylabel:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=6.8, pad=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", linewidth=0.5)
    bold_all_text(ax)
    return rows


def plot_normalization_comparison(none_payload: dict, sample_quantile_payload: dict, motifs: list[str], out_prefix: Path) -> None:
    xvals = np.asarray(none_payload.get("aggregate", {}).get("x", []), dtype=float)
    xvals_norm = np.asarray(sample_quantile_payload.get("aggregate", {}).get("x", []), dtype=float)
    if not np.array_equal(xvals, xvals_norm):
        raise ValueError("Normalization reports use different aggregate x axes")

    apply_style(base_size=8)
    fig = plt.figure(figsize=(7.4, 5.6))
    grid = fig.add_gridspec(len(motifs) + 1, 2, height_ratios=[1.05] + [1.0] * len(motifs), hspace=0.56, wspace=0.24)
    volcano_axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])]
    profile_axes = np.asarray([[fig.add_subplot(grid[row + 1, col]) for col in range(2)] for row in range(len(motifs))])
    summary_rows: list[dict[str, object]] = []
    volcano_rows: list[dict[str, object]] = []

    volcano_xlim, volcano_ylim = volcano_limits(none_payload, sample_quantile_payload)
    volcano_rows.append(plot_volcano(volcano_axes[0], none_payload, "No normalization", volcano_xlim, volcano_ylim))
    volcano_rows.append(plot_volcano(volcano_axes[1], sample_quantile_payload, "Sample-quantile", volcano_xlim, volcano_ylim))
    handles, labels = volcano_axes[1].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    legend_labels = [label for label in ["Bcell_up", "Tcell_up", "n.s."] if label in by_label]
    fig.legend(
        [by_label[label] for label in legend_labels],
        legend_labels,
        loc="upper center",
        frameon=False,
        fontsize=6.8,
        markerscale=1.25,
        bbox_to_anchor=(0.52, 0.932),
        ncol=3,
        columnspacing=1.1,
        handletextpad=0.35,
    )
    for row, prefix in enumerate(motifs):
        none_motif = get_motif(none_payload, prefix)
        norm_motif = get_motif(sample_quantile_payload, prefix)
        title = f"{motif_label(none_motif)}, n={int(none_motif.get('n_sites', 0)):,}"
        profile_axes[row, 0].text(
            -0.02,
            1.12,
            title,
            transform=profile_axes[row, 0].transAxes,
            fontsize=8.0,
            fontweight="bold",
            ha="left",
            va="bottom",
            clip_on=False,
        )
        summary_rows.extend(plot_motif_panel(profile_axes[row, 0], xvals, none_motif, "No normalization", html_like_ylim(none_motif), ylabel="Signal (a.u.)"))
        summary_rows.extend(plot_motif_panel(profile_axes[row, 1], xvals, norm_motif, "Sample-quantile", html_like_ylim(norm_motif), ylabel=None))

    for ax in profile_axes[-1, :]:
        ax.set_xlabel("Distance from motif center (bp)")

    handles, labels = profile_axes[0, 1].get_legend_handles_labels()
    unique = {}
    for handle, label in zip(handles, labels):
        unique.setdefault(label, handle)
    fig.legend(
        unique.values(),
        unique.keys(),
        loc="lower center",
        ncol=6,
        frameon=False,
        fontsize=6.5,
        bbox_to_anchor=(0.5, -0.006),
        columnspacing=0.85,
        handlelength=1.9,
        handletextpad=0.35,
    )
    fig.text(0.012, 0.988, "A", fontsize=10.5, fontweight="bold", va="top")
    fig.text(0.012, 0.718, "B", fontsize=10.5, fontweight="bold", va="top")
    fig.suptitle("Effect of sample-level quantile normalization on differential footprint reports", y=0.992, fontsize=9.5)
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.075, top=0.895, hspace=0.62, wspace=0.24)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight", dpi=600)
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight", dpi=600)
    pd.DataFrame(summary_rows).to_csv(out_prefix.with_suffix(".tsv"), sep="\t", index=False)
    pd.DataFrame(volcano_rows).to_csv(out_prefix.with_name(out_prefix.name + "_volcano.tsv"), sep="\t", index=False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--none-report-html", type=Path, default=DEFAULT_NONE_REPORT)
    parser.add_argument("--sample-quantile-report-html", type=Path, default=DEFAULT_SAMPLE_QUANTILE_REPORT)
    parser.add_argument("--out-prefix", type=Path, default=Path("manuscript/figures/normalization_effect"))
    parser.add_argument("--motifs", nargs="+", default=DEFAULT_MOTIFS)
    args = parser.parse_args(argv)

    plot_normalization_comparison(
        load_report_payload(args.none_report_html),
        load_report_payload(args.sample_quantile_report_html),
        args.motifs,
        args.out_prefix,
    )
    print(f"Wrote {args.out_prefix.with_suffix('.png')}")
    print(f"Wrote {args.out_prefix.with_suffix('.pdf')}")
    print(f"Wrote {args.out_prefix.with_suffix('.svg')}")
    print(f"Wrote {args.out_prefix.with_suffix('.tsv')}")
    print(f"Wrote {args.out_prefix.with_name(args.out_prefix.name + '_volcano.tsv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
