#!/usr/bin/env python
"""Plot de novo motif validation from the Buenrostro replicate experiment."""

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

REPO_ROOT = SCRIPT_DIR.parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from figure_style import apply_style, bold_all_text  # noqa: E402
from fp_tools.utils.motifs import MotifList  # noqa: E402


DEFAULT_VALIDATION_DIR = (
    REPO_ROOT
    / "data/public/processed/buenrostro_atac_replicates/fp_tools/denovo_motif_validation"
)
DEFAULT_JASPAR = (
    REPO_ROOT
    / "data/public/raw/jaspar/2026/JASPAR2026_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
)
CONDITION_COLORS = {"Bcell": "#1f77b4", "Tcell": "#d62728"}
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


def jaspar_name_map(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    motifs = MotifList().from_file(str(path))
    return {motif.id: motif.name for motif in motifs}


def load_top_tomtom(validation_dir: Path, id_to_name: dict[str, str]) -> pd.DataFrame:
    rows = []
    for direction, label in [
        ("Bcell_vs_Tcell_streme", "B-cell candidate enrichment"),
        ("Tcell_vs_Bcell_streme", "T-cell candidate enrichment"),
    ]:
        path = validation_dir / "motifs" / direction / "motif_summary.tsv"
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        tomtom = df[df["source"].eq("Tomtom")].copy()
        if tomtom.empty:
            continue
        tomtom["q_value"] = pd.to_numeric(tomtom["q_value"], errors="coerce")
        tomtom = tomtom.sort_values(["motif_id", "q_value"]).drop_duplicates("motif_id")
        for _, row in tomtom.head(4).iterrows():
            target_id = str(row["target_id"])
            rows.append(
                {
                    "direction": label,
                    "de_novo_motif": str(row["motif_id"]),
                    "consensus": str(row["consensus"]),
                    "target_id": target_id,
                    "target_name": id_to_name.get(target_id, target_id),
                    "q_value": float(row["q_value"]),
                }
            )
    return pd.DataFrame(rows)


def load_result_counts(validation_dir: Path) -> pd.DataFrame:
    rows = []
    labels = {
        "denovo_only": "de novo only",
        "jaspar2026_plus_denovo": "JASPAR2026 + de novo",
        "restricted_jaspar": "restricted JASPAR",
        "restricted_jaspar_plus_denovo": "restricted + de novo",
    }
    for key, label in labels.items():
        path = validation_dir / "diff_footprints" / key / "diff_footprints_results.txt"
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        pvals = pd.to_numeric(df["Bcell_Tcell_pvalue"], errors="coerce")
        highlighted = df["Bcell_Tcell_highlighted"].astype(str).eq("True")
        denovo = df["name"].astype(str).str.contains("_denovo_", regex=False)
        rows.append(
            {
                "result_set": label,
                "n_tested": len(df),
                "n_significant": int((pvals < 0.05).sum()),
                "n_highlighted": int(highlighted.sum()),
                "n_significant_de_novo": int(((pvals < 0.05) & denovo).sum()),
                "n_highlighted_de_novo": int((highlighted & denovo).sum()),
            }
        )
    return pd.DataFrame(rows)


def center_flank_score(profile: np.ndarray, xvals: np.ndarray, center_bp: int = 10, edge_bp: int = 20) -> float:
    center = profile[np.abs(xvals) <= center_bp]
    left = profile[xvals <= xvals.min() + edge_bp]
    right = profile[xvals >= xvals.max() - edge_bp]
    flanks = np.concatenate([left, right])
    if center.size == 0 or flanks.size == 0:
        return float("nan")
    return float(np.nanmean(center) - np.nanmean(flanks))


def select_aggregate_motifs(payload: dict, n: int = 4) -> list[dict]:
    motifs = payload.get("aggregate", {}).get("motifs", [])
    scored = []
    xvals = np.asarray(payload.get("aggregate", {}).get("x", []), dtype=float)
    for motif in motifs:
        condition_scores = {}
        for condition in motif.get("conditions", []):
            condition_scores[condition["name"]] = center_flank_score(
                np.asarray(condition["profile"], dtype=float), xvals
            )
        if {"Bcell", "Tcell"} <= set(condition_scores):
            depletion_delta = abs(condition_scores["Bcell"] - condition_scores["Tcell"])
        else:
            depletion_delta = 0.0
        scored.append(
            (
                depletion_delta,
                abs(float(motif.get("change", 0.0))),
                -math.log10(max(float(motif.get("pvalue", 1.0)), 1e-300)),
                motif,
            )
        )
    scored.sort(reverse=True, key=lambda item: item[:3])
    return [item[-1] for item in scored[:n]]


def plot_validation(validation_dir: Path, jaspar: Path, out_prefix: Path) -> None:
    motif_sets = pd.read_csv(validation_dir / "motifs" / "motif_set_summary.tsv", sep="\t")
    id_to_name = jaspar_name_map(jaspar)
    tomtom = load_top_tomtom(validation_dir, id_to_name)
    counts = load_result_counts(validation_dir)
    payload = load_report_payload(
        validation_dir
        / "diff_footprints"
        / "denovo_only"
        / "diff_footprints_Bcell_Tcell.html"
    )
    xvals = np.asarray(payload.get("aggregate", {}).get("x", []), dtype=float)
    aggregate_motifs = select_aggregate_motifs(payload)

    apply_style(base_size=8.5)
    fig = plt.figure(figsize=(7.4, 8.8))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.05, 2.0], hspace=0.55, wspace=0.34)

    ax_sets = fig.add_subplot(gs[0, 0])
    set_order = [
        "de_novo_only",
        "jaspar2026_full",
        "jaspar2026_plus_denovo",
        "jaspar2026_restricted",
        "jaspar2026_restricted_plus_denovo",
    ]
    set_labels = {
        "de_novo_only": "de novo",
        "jaspar2026_full": "JASPAR2026",
        "jaspar2026_plus_denovo": "JASPAR + de novo",
        "jaspar2026_restricted": "restricted",
        "jaspar2026_restricted_plus_denovo": "restricted + de novo",
    }
    motif_sets = motif_sets.set_index("motif_set").loc[set_order].reset_index()
    colors = ["#6a3d9a", "#4c78a8", "#72b7b2", "#9aa0a6", "#f58518"]
    ax_sets.bar(
        [set_labels[m] for m in motif_sets["motif_set"]],
        motif_sets["n_motifs"].astype(int),
        color=colors,
        edgecolor="0.25",
        linewidth=0.5,
    )
    ax_sets.set_ylabel("Motifs tested")
    ax_sets.set_title("A. Motif-set construction")
    ax_sets.tick_params(axis="x", rotation=35)
    ax_sets.spines[["top", "right"]].set_visible(False)
    bold_all_text(ax_sets)

    ax_counts = fig.add_subplot(gs[0, 1])
    count_plot = counts[counts["result_set"].isin(["restricted JASPAR", "restricted + de novo"])].copy()
    x = np.arange(len(count_plot))
    ax_counts.bar(x - 0.18, count_plot["n_highlighted"], width=0.34, label="Highlighted", color="#4c78a8")
    ax_counts.bar(x + 0.18, count_plot["n_highlighted_de_novo"], width=0.34, label="de novo highlighted", color="#f58518")
    ax_counts.set_xticks(x, count_plot["result_set"], rotation=20)
    ax_counts.set_ylabel("Motif families")
    ax_counts.set_title("B. Rescue with database supplement")
    ax_counts.legend(frameon=False, fontsize=6.8)
    ax_counts.spines[["top", "right"]].set_visible(False)
    bold_all_text(ax_counts)

    ax_tomtom = fig.add_subplot(gs[1, :])
    if not tomtom.empty:
        tomtom = tomtom.sort_values("q_value").head(8).copy()
        tomtom["label"] = tomtom["de_novo_motif"] + " -> " + tomtom["target_name"] + " (" + tomtom["target_id"] + ")"
        y = np.arange(len(tomtom))[::-1]
        ax_tomtom.barh(y, -np.log10(tomtom["q_value"].clip(lower=1e-300)), color="#6a3d9a")
        ax_tomtom.set_yticks(y, tomtom["label"])
        ax_tomtom.set_xlabel("-log10 Tomtom q-value")
        ax_tomtom.set_title("C. Best known-motif matches for discovered motifs")
    else:
        ax_tomtom.text(0.5, 0.5, "No Tomtom matches found", ha="center", va="center")
    ax_tomtom.spines[["top", "right"]].set_visible(False)
    bold_all_text(ax_tomtom)

    sub_gs = gs[2, :].subgridspec(2, 2, hspace=0.42, wspace=0.32)
    summary_rows: list[dict[str, object]] = []
    for idx, motif in enumerate(aggregate_motifs):
        ax = fig.add_subplot(sub_gs[idx // 2, idx % 2])
        for condition in motif.get("conditions", []):
            condition_name = str(condition["name"])
            color = CONDITION_COLORS.get(condition_name, "#555555")
            for sample in condition.get("samples", []):
                sample_name = str(sample.get("name", "sample"))
                profile = np.asarray(sample["profile"], dtype=float)
                ax.plot(
                    xvals,
                    profile,
                    color=color,
                    linestyle=SAMPLE_STYLES.get(sample_name, "solid"),
                    linewidth=0.7,
                    alpha=0.5,
                    label=sample_name,
                )
                summary_rows.append(
                    {
                        "panel": "aggregate",
                        "motif": motif.get("name", ""),
                        "motif_id": motif.get("motif_id", ""),
                        "sample": sample_name,
                        "condition": condition_name,
                        "center_minus_flank": center_flank_score(profile, xvals),
                    }
                )
            ax.plot(
                xvals,
                np.asarray(condition["profile"], dtype=float),
                color=color,
                linewidth=1.5,
                alpha=0.98,
                label=f"{condition_name} mean",
            )
        ax.axvline(0, color="0.35", linewidth=0.75)
        ax.set_title(
            f"{chr(ord('D') + idx)}. {motif.get('name')} ({motif.get('motif_id')})\n"
            f"n={int(motif.get('n_sites', 0)):,}; p={float(motif.get('pvalue', 1.0)):.1e}",
            fontsize=7.8,
        )
        ax.set_xlabel("Distance from motif center (bp)")
        ax.set_ylabel("Normalized cut-site signal")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="0.9", linewidth=0.5)
        bold_all_text(ax)

    handles, labels = fig.axes[-1].get_legend_handles_labels()
    unique = {}
    for handle, label in zip(handles, labels):
        unique.setdefault(label, handle)
    fig.legend(unique.values(), unique.keys(), loc="lower center", ncol=6, frameon=False, fontsize=6.7)
    fig.suptitle("De novo motif validation in B-cell versus T-cell replicate ATAC-seq", y=0.995)
    fig.tight_layout(rect=(0, 0.035, 1, 0.985))

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".svg"), bbox_inches="tight")

    out_table = out_prefix.with_suffix(".tsv")
    with out_table.open("w", encoding="utf-8") as handle:
        handle.write("# motif_set_summary\n")
        motif_sets.to_csv(handle, sep="\t", index=False)
        handle.write("\n# differential_result_counts\n")
        counts.to_csv(handle, sep="\t", index=False)
        handle.write("\n# top_tomtom_matches\n")
        tomtom.to_csv(handle, sep="\t", index=False)
        handle.write("\n# aggregate_center_minus_flank\n")
        pd.DataFrame(summary_rows).to_csv(handle, sep="\t", index=False)
    print(f"Wrote {out_prefix.with_suffix('.png')}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-dir", type=Path, default=DEFAULT_VALIDATION_DIR)
    parser.add_argument("--jaspar", type=Path, default=DEFAULT_JASPAR)
    parser.add_argument("--out-prefix", type=Path, default=REPO_ROOT / "manuscript/figures/denovo_motif_validation")
    args = parser.parse_args(argv)
    plot_validation(args.validation_dir, args.jaspar, args.out_prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
