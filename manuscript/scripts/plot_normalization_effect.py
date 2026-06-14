#!/usr/bin/env python
"""Plot real 2-vs-2 aggregate profiles before and after normalization."""

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
    "JUNB_MA1140.3",
    "ATF7_MA0834.2",
    "IRF4_MA1419.2",
]
CONDITION_COLORS = {
    "Bcell": "#1f77b4",
    "Tcell": "#d62728",
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


def smooth(profile: np.ndarray, window: int = 5) -> np.ndarray:
    if window <= 1:
        return profile
    pad = window // 2
    kernel = np.ones(window, dtype=float) / window
    return np.convolve(np.pad(profile, pad, mode="edge"), kernel, mode="valid")


def get_motif(payload: dict, prefix: str) -> dict:
    motifs = payload.get("aggregate", {}).get("motifs", [])
    motif = next((item for item in motifs if item.get("prefix") == prefix), None)
    if motif is None:
        raise ValueError(f"Motif {prefix} missing from aggregate payload")
    return motif


def collect_values(*motifs: dict) -> np.ndarray:
    values = []
    for motif in motifs:
        for condition in motif.get("conditions", []):
            values.extend(condition.get("profile", []))
            for sample in condition.get("samples", []):
                values.extend(sample.get("profile", []))
    return np.asarray(values, dtype=float)


def plot_motif_panel(ax, xvals: np.ndarray, motif: dict, mode: str, ylim: tuple[float, float], show_ylabel: bool) -> list[dict[str, object]]:
    rows = []
    for condition in motif.get("conditions", []):
        condition_name = str(condition["name"])
        color = CONDITION_COLORS.get(condition_name, "#555555")
        for sample in condition.get("samples", []):
            sample_name = str(sample.get("name", "sample"))
            profile = np.asarray(sample["profile"], dtype=float)
            ax.plot(
                xvals,
                smooth(profile),
                color=color,
                linestyle=SAMPLE_STYLES.get(sample_name, "solid"),
                linewidth=0.85,
                alpha=0.55,
                label=sample_name,
                zorder=1,
            )
            rows.append(
                {
                    "normalization": mode,
                    "motif": motif.get("name", ""),
                    "motif_id": motif.get("motif_id", ""),
                    "output_prefix": motif.get("prefix", ""),
                    "profile_type": "replicate",
                    "condition": condition_name,
                    "sample": sample_name,
                    "n_sites": motif.get("n_sites", ""),
                    "center_minus_flank": center_minus_flank(profile, xvals),
                }
            )
        mean_profile = np.asarray(condition["profile"], dtype=float)
        ax.plot(
            xvals,
            smooth(mean_profile),
            color=color,
            linewidth=1.95,
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
    ax.axvline(0, color="0.35", linewidth=0.8, alpha=0.85)
    ax.axhline(0, color="0.82", linewidth=0.55, zorder=0)
    ax.set_ylim(*ylim)
    ax.set_title(mode, fontsize=8.3)
    if show_ylabel:
        ax.set_ylabel("Cut-site signal")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", linewidth=0.55)
    bold_all_text(ax)
    return rows


def plot_normalization_comparison(none_payload: dict, sample_quantile_payload: dict, motifs: list[str], out_prefix: Path) -> None:
    xvals = np.asarray(none_payload.get("aggregate", {}).get("x", []), dtype=float)
    xvals_norm = np.asarray(sample_quantile_payload.get("aggregate", {}).get("x", []), dtype=float)
    if not np.array_equal(xvals, xvals_norm):
        raise ValueError("Normalization reports use different aggregate x axes")

    apply_style(base_size=8)
    fig, axes = plt.subplots(len(motifs), 2, figsize=(7.6, 1.82 * len(motifs)), sharex=True)
    if len(motifs) == 1:
        axes = np.asarray([axes])
    summary_rows: list[dict[str, object]] = []

    for row, prefix in enumerate(motifs):
        none_motif = get_motif(none_payload, prefix)
        norm_motif = get_motif(sample_quantile_payload, prefix)
        values = collect_values(none_motif, norm_motif)
        lo = float(np.nanpercentile(values, 1))
        hi = float(np.nanpercentile(values, 99))
        pad = max((hi - lo) * 0.08, 0.01)
        ylim = (lo - pad, hi + pad)

        title = f"{motif_label(none_motif)}, n={int(none_motif.get('n_sites', 0)):,}"
        axes[row, 0].text(-0.12, 1.06, title, transform=axes[row, 0].transAxes, fontsize=8.6, fontweight="bold")
        summary_rows.extend(plot_motif_panel(axes[row, 0], xvals, none_motif, "No normalization", ylim, show_ylabel=True))
        summary_rows.extend(plot_motif_panel(axes[row, 1], xvals, norm_motif, "Sample-quantile", ylim, show_ylabel=False))

    for ax in axes[-1, :]:
        ax.set_xlabel("Distance from motif center (bp)")

    handles, labels = axes[0, 1].get_legend_handles_labels()
    unique = {}
    for handle, label in zip(handles, labels):
        unique.setdefault(label, handle)
    fig.legend(
        unique.values(),
        unique.keys(),
        loc="lower center",
        ncol=6,
        frameon=False,
        fontsize=6.7,
        bbox_to_anchor=(0.5, -0.012),
    )
    fig.suptitle("Real 2-vs-2 corrected cut-site aggregates before and after normalization", y=0.995, fontsize=9.5)
    fig.tight_layout(rect=(0, 0.035, 1, 0.975))

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=350, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight")
    pd.DataFrame(summary_rows).to_csv(out_prefix.with_suffix(".tsv"), sep="\t", index=False)


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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
