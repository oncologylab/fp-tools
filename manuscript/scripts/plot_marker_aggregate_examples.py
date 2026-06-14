#!/usr/bin/env python
"""Plot marker aggregate examples from a diff-footprints HTML report."""

from __future__ import annotations

import argparse
import base64
import gzip
import json
import math
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


DEFAULT_MOTIFS = [
    "BACH2_MA1101.3",
    "JUNB_MA1140.3",
    "ATF7_MA0834.2",
    "IRF4_MA1419.2",
]

DEFAULT_ROLES = {
    "BACH2_MA1101.3": "B/GM12878-biased shared regulator",
    "JUNB_MA1140.3": "T-cell activation/AP-1 regulator",
    "ATF7_MA0834.2": "T-cell-deeper aggregate example",
    "IRF4_MA1419.2": "TCR-response regulator; modest aggregate",
}

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


def center_flank_score(profile: np.ndarray, xvals: np.ndarray, center_bp: int = 10, edge_bp: int = 20) -> float:
    center = profile[np.abs(xvals) <= center_bp]
    left = profile[xvals <= xvals.min() + edge_bp]
    right = profile[xvals >= xvals.max() - edge_bp]
    flanks = np.concatenate([left, right])
    if center.size == 0 or flanks.size == 0:
        return float("nan")
    return float(np.nanmean(center) - np.nanmean(flanks))


def motif_label(motif: dict) -> str:
    motif_id = motif.get("motif_id", "")
    return f"{motif.get('name', motif.get('prefix', 'motif'))} ({motif_id})" if motif_id else motif.get("name", motif.get("prefix", "motif"))


def collect_summary_rows(motif: dict, xvals: np.ndarray, role: str) -> list[dict[str, object]]:
    rows = []
    condition_scores = {}
    for condition in motif.get("conditions", []):
        condition_name = str(condition["name"])
        profile = np.asarray(condition["profile"], dtype=float)
        condition_score = center_flank_score(profile, xvals)
        condition_scores[condition_name] = condition_score
        rows.append(
            {
                "motif": motif.get("name", ""),
                "motif_id": motif.get("motif_id", ""),
                "output_prefix": motif.get("prefix", ""),
                "role": role,
                "profile_type": "condition_mean",
                "condition": condition_name,
                "sample": condition_name,
                "n_sites": motif.get("n_sites", ""),
                "center_minus_flank": condition_score,
                "interpretation": "more negative means stronger footprint-like center depletion",
            }
        )
        for sample in condition.get("samples", []):
            sample_profile = np.asarray(sample["profile"], dtype=float)
            rows.append(
                {
                    "motif": motif.get("name", ""),
                    "motif_id": motif.get("motif_id", ""),
                    "output_prefix": motif.get("prefix", ""),
                    "role": role,
                    "profile_type": "replicate",
                    "condition": condition_name,
                    "sample": sample.get("name", ""),
                    "n_sites": motif.get("n_sites", ""),
                    "center_minus_flank": center_flank_score(sample_profile, xvals),
                    "interpretation": "more negative means stronger footprint-like center depletion",
                }
            )
    if {"Bcell", "Tcell"} <= set(condition_scores):
        rows.append(
            {
                "motif": motif.get("name", ""),
                "motif_id": motif.get("motif_id", ""),
                "output_prefix": motif.get("prefix", ""),
                "role": role,
                "profile_type": "condition_delta",
                "condition": "Bcell_minus_Tcell",
                "sample": "Bcell_minus_Tcell",
                "n_sites": motif.get("n_sites", ""),
                "center_minus_flank": condition_scores["Bcell"] - condition_scores["Tcell"],
                "interpretation": "negative delta means stronger B-cell depletion; positive delta means stronger T-cell depletion",
            }
        )
    return rows


def plot_marker_aggregates(payload: dict, motifs: list[str], roles: dict[str, str], out_prefix: Path) -> None:
    aggregate = payload.get("aggregate", {})
    xvals = np.asarray(aggregate.get("x", []), dtype=float)
    motifs_by_prefix = {motif["prefix"]: motif for motif in aggregate.get("motifs", [])}
    missing = [prefix for prefix in motifs if prefix not in motifs_by_prefix]
    if missing:
        raise ValueError(f"Motifs missing from aggregate payload: {', '.join(missing)}")

    apply_style(base_size=9)
    ncols = 2
    nrows = math.ceil(len(motifs) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(7.4, 2.8 * nrows), sharex=True)
    axes = np.asarray(axes).reshape(-1)
    summary_rows: list[dict[str, object]] = []

    for ax, prefix in zip(axes, motifs):
        motif = motifs_by_prefix[prefix]
        role = roles.get(prefix, "")
        summary_rows.extend(collect_summary_rows(motif, xvals, role))
        for condition in motif.get("conditions", []):
            condition_name = str(condition["name"])
            color = CONDITION_COLORS.get(condition_name, "#555555")
            for sample in condition.get("samples", []):
                sample_name = str(sample.get("name", "sample"))
                ax.plot(
                    xvals,
                    np.asarray(sample["profile"], dtype=float),
                    color=color,
                    linestyle=SAMPLE_STYLES.get(sample_name, "solid"),
                    linewidth=0.85,
                    alpha=0.52,
                    label=sample_name,
                    zorder=1,
                )
            ax.plot(
                xvals,
                np.asarray(condition["profile"], dtype=float),
                color=color,
                linewidth=1.7,
                alpha=0.98,
                label=f"{condition_name} mean",
                zorder=2,
            )
        ax.axvline(0, color="0.35", linewidth=0.8, alpha=0.8)
        ax.set_title(f"{motif_label(motif)}\n{role}, n={int(motif.get('n_sites', 0)):,} sites", fontsize=8.5)
        ax.set_ylabel("Normalized cut-site signal")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="0.9", linewidth=0.6)
        bold_all_text(ax)

    for ax in axes[len(motifs) :]:
        ax.axis("off")
    for ax in axes[-ncols:]:
        ax.set_xlabel("Distance from motif center (bp)")
    handles, labels = axes[0].get_legend_handles_labels()
    unique = {}
    for handle, label in zip(handles, labels):
        unique.setdefault(label, handle)
    fig.legend(unique.values(), unique.keys(), loc="lower center", ncol=6, frameon=False, fontsize=6.8, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("B-cell-biased and T-cell activation-associated aggregate examples", y=1.0)
    fig.tight_layout(rect=(0, 0.035, 1, 0.98))

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight")
    pd.DataFrame(summary_rows).to_csv(out_prefix.with_suffix(".tsv"), sep="\t", index=False)


def parse_roles(role_args: list[str]) -> dict[str, str]:
    roles = dict(DEFAULT_ROLES)
    for item in role_args:
        if "=" not in item:
            raise ValueError(f"Role must use PREFIX=text format: {item}")
        prefix, text = item.split("=", 1)
        roles[prefix] = text
    return roles


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-html", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument("--motifs", nargs="+", default=DEFAULT_MOTIFS)
    parser.add_argument("--role", action="append", default=[], help="Override a panel role label as PREFIX=text.")
    args = parser.parse_args(argv)

    payload = load_report_payload(args.report_html)
    plot_marker_aggregates(payload, args.motifs, parse_roles(args.role), args.out_prefix)
    print(f"Wrote {args.out_prefix.with_suffix('.png')}")
    print(f"Wrote {args.out_prefix.with_suffix('.pdf')}")
    print(f"Wrote {args.out_prefix.with_suffix('.svg')}")
    print(f"Wrote {args.out_prefix.with_suffix('.tsv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
