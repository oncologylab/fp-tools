#!/usr/bin/env python
"""Plot de novo motif validation from the Buenrostro replicate experiment."""

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
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

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
SELECTED_AGGREGATES = [
    "Bcell_denovo_5_Bcell_denovo_5_5-GATGAGTCA",
    "Tcell_denovo_4_Tcell_denovo_4_4-GAHGYGGAA",
    "Tcell_denovo_6_Tcell_denovo_6_6-AGGAAGTSACTGA",
    "Tcell_denovo_1_Tcell_denovo_1_1-ACAGTTTCCT",
]
AGGREGATE_ROLES = {
    "Bcell_denovo_5_Bcell_denovo_5_5-GATGAGTCA": "B-cell-deeper de novo footprint",
    "Tcell_denovo_4_Tcell_denovo_4_4-GAHGYGGAA": "T-cell-deeper de novo footprint",
    "Tcell_denovo_6_Tcell_denovo_6_6-AGGAAGTSACTGA": "T-cell-deeper de novo footprint",
    "Tcell_denovo_1_Tcell_denovo_1_1-ACAGTTTCCT": "weaker T-cell-deeper footprint",
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


def load_streme_motifs(validation_dir: Path, id_to_name: dict[str, str]) -> pd.DataFrame:
    rows = []
    for direction, label in [
        ("Bcell_vs_Tcell_streme", "B-cell candidates"),
        ("Tcell_vs_Bcell_streme", "T-cell candidates"),
    ]:
        path = validation_dir / "motifs" / direction / "motif_summary.tsv"
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        discovered = df[df["source"].eq("MEME")].copy()
        tomtom = df[df["source"].eq("Tomtom")].copy()
        tomtom["q_value"] = pd.to_numeric(tomtom.get("q_value"), errors="coerce")
        best_by_motif = tomtom.sort_values("q_value").drop_duplicates("motif_id").set_index("motif_id")
        for _, row in discovered.iterrows():
            motif_id = str(row["motif_id"])
            best = best_by_motif.loc[motif_id] if motif_id in best_by_motif.index else None
            if best is not None and pd.notna(best["q_value"]):
                target_id = str(best["target_id"])
                q_value = float(best["q_value"])
                target_name = id_to_name.get(target_id, target_id)
                match = f"{target_name} ({target_id})" if q_value <= 0.05 else "no confident match"
            else:
                target_id = ""
                target_name = ""
                q_value = float("nan")
                match = "no confident match"
            rows.append(
                {
                    "direction": label,
                    "de_novo_motif": motif_id,
                    "consensus": str(row["consensus"]),
                    "sites": int(float(row["sites"])),
                    "e_value": float(row["e_value"]),
                    "target_id": target_id,
                    "target_name": target_name,
                    "tomtom_q_value": q_value,
                    "tomtom_label": match,
                    "confident_tomtom": bool(pd.notna(q_value) and q_value <= 0.05),
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


def load_rescued_denovo(validation_dir: Path) -> pd.DataFrame:
    path = validation_dir / "diff_footprints" / "restricted_jaspar_plus_denovo" / "diff_footprints_results.txt"
    df = pd.read_csv(path, sep="\t")
    denovo = df[df["name"].astype(str).str.contains("_denovo_", regex=False)].copy()
    denovo["Bcell_Tcell_pvalue"] = pd.to_numeric(denovo["Bcell_Tcell_pvalue"], errors="coerce")
    denovo["Bcell_Tcell_change"] = pd.to_numeric(denovo["Bcell_Tcell_change"], errors="coerce")
    denovo["highlighted"] = denovo["Bcell_Tcell_highlighted"].astype(str).eq("True")
    denovo["direction"] = np.where(denovo["Bcell_Tcell_change"] >= 0, "B-cell higher", "T-cell higher")
    return denovo.sort_values(["highlighted", "Bcell_Tcell_pvalue"], ascending=[False, True])


def center_flank_score(profile: np.ndarray, xvals: np.ndarray, center_bp: int = 10, edge_bp: int = 20) -> float:
    center = profile[np.abs(xvals) <= center_bp]
    left = profile[xvals <= xvals.min() + edge_bp]
    right = profile[xvals >= xvals.max() - edge_bp]
    flanks = np.concatenate([left, right])
    if center.size == 0 or flanks.size == 0:
        return float("nan")
    return float(np.nanmean(center) - np.nanmean(flanks))


def motif_lookup(payload: dict) -> dict[str, dict]:
    return {motif["prefix"]: motif for motif in payload.get("aggregate", {}).get("motifs", [])}


def draw_box(ax, xy, text, width=0.78, height=0.105, fc="#f8fbff", ec="#4c78a8"):
    x, y = xy
    box = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.016,rounding_size=0.018",
        facecolor=fc,
        edgecolor=ec,
        linewidth=0.9,
        transform=ax.transAxes,
    )
    ax.add_patch(box)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", transform=ax.transAxes, fontsize=7.4, fontweight="bold")


def draw_arrow(ax, start, end):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=10,
            linewidth=0.9,
            color="0.3",
            transform=ax.transAxes,
        )
    )


def plot_workflow(ax):
    ax.axis("off")
    ax.set_title("A. Two ways de novo motifs enter fp-tools", loc="left", pad=8)
    draw_box(ax, (0.08, 0.78), "Corrected ATAC\ncut-site tracks", fc="#eef5ff")
    draw_box(ax, (0.08, 0.60), "call-footprints\nrank local candidates", fc="#eefaf0", ec="#2f7d32")
    draw_box(ax, (0.08, 0.42), "candidate FASTA\n+ STREME", fc="#f5f0ff", ec="#6a3d9a")
    draw_box(ax, (0.08, 0.22), "Mode 1: de novo only\ndiff-footprints", fc="#fff7ec", ec="#f58518")
    draw_box(ax, (0.08, 0.04), "Mode 2: database + de novo\nsupplement missing motifs", fc="#fff7ec", ec="#f58518")
    for y1, y2 in [(0.78, 0.705), (0.60, 0.525), (0.42, 0.325), (0.22, 0.145)]:
        draw_arrow(ax, (0.47, y1), (0.47, y2))


def plot_discovery_table(ax, streme: pd.DataFrame, rescued: pd.DataFrame):
    ax.axis("off")
    ax.set_title("B. De novo-only discovery produces testable motif families", loc="left", pad=8)
    selected_names = set(rescued[rescued["highlighted"]]["name"].astype(str))
    rows = []
    for _, row in streme.iterrows():
        source_prefix = "Bcell" if row["direction"].startswith("B-cell") else "Tcell"
        short_id = str(row["de_novo_motif"]).split("-", 1)[0]
        prefix = f"{source_prefix}_denovo_{short_id}"
        priority = 3 if row["confident_tomtom"] else 0
        if prefix in selected_names:
            priority += 2
        if prefix in {"Bcell_denovo_5", "Tcell_denovo_4", "Tcell_denovo_6", "Tcell_denovo_1"}:
            priority += 1
        if priority > 0:
            item = row.copy()
            item["priority"] = priority
            rows.append(item)
    if not rows:
        rows = [row for _, row in streme.head(6).iterrows()]
    display = pd.DataFrame(rows).sort_values(["priority", "confident_tomtom", "direction"], ascending=[False, False, True]).head(7)
    y0 = 0.88
    col_x = [0.00, 0.27, 0.50, 0.68]
    headers = ["source", "consensus", "sites", "Tomtom annotation"]
    for x, h in zip(col_x, headers):
        ax.text(x, y0, h, transform=ax.transAxes, fontsize=7.2, fontweight="bold", va="top")
    ax.plot([0, 1], [y0 - 0.035, y0 - 0.035], color="0.72", linewidth=0.7, transform=ax.transAxes)
    y = y0 - 0.085
    summary_rows = []
    for _, row in display.iterrows():
        source = "B candidates" if row["direction"].startswith("B-cell") else "T candidates"
        q = row["tomtom_q_value"]
        tomtom = row["tomtom_label"]
        if pd.notna(q) and row["confident_tomtom"]:
            tomtom = f"{tomtom}; q={q:.1e}"
        ax.text(col_x[0], y, source, transform=ax.transAxes, fontsize=6.8, va="top")
        ax.text(col_x[1], y, row["consensus"], transform=ax.transAxes, fontsize=6.8, va="top", family="monospace")
        ax.text(col_x[2], y, f"{int(row['sites']):,}", transform=ax.transAxes, fontsize=6.8, va="top")
        ax.text(col_x[3], y, tomtom, transform=ax.transAxes, fontsize=6.8, va="top")
        summary_rows.append(row.to_dict())
        y -= 0.105
    ax.text(
        0.0,
        0.02,
        "Tomtom labels are motif-similarity annotations, not definitive TF identity calls.",
        transform=ax.transAxes,
        fontsize=6.6,
        color="0.35",
        va="bottom",
    )
    return summary_rows


def plot_rescue(ax, counts: pd.DataFrame, rescued: pd.DataFrame):
    ax.set_title("C. Supplement mode rescues extra differential motif families", loc="left", pad=8)
    count_plot = counts[counts["result_set"].isin(["restricted JASPAR", "restricted + de novo"])].copy()
    x = np.arange(len(count_plot))
    ax.bar(x, count_plot["n_highlighted"], color=["#9aa0a6", "#f58518"], edgecolor="0.25", linewidth=0.5)
    ax.set_xticks(x, ["restricted\nJASPAR", "restricted\n+ de novo"])
    ax.set_ylabel("Highlighted families")
    ax.set_ylim(0, max(count_plot["n_highlighted"]) * 1.32)
    for xi, value in zip(x, count_plot["n_highlighted"]):
        ax.text(xi, value + 2.0, f"{int(value)}", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    plus = count_plot[count_plot["result_set"].eq("restricted + de novo")].iloc[0]
    base = count_plot[count_plot["result_set"].eq("restricted JASPAR")].iloc[0]
    ax.text(1.0, plus["n_highlighted"] + 9, f"+{int(plus['n_highlighted'] - base['n_highlighted'])} total\n{int(plus['n_highlighted_de_novo'])} de novo", ha="center", va="bottom", fontsize=7.0, color="#a45100", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", linewidth=0.5)
    highlighted = rescued[rescued["highlighted"]].head(6)
    text = "Highlighted de novo motifs:\n" + "\n".join(
        f"{r['name'].replace('_denovo_', ' dn')}: change {float(r['Bcell_Tcell_change']):+.2f}"
        for _, r in highlighted.iterrows()
    )
    ax.text(
        0.47,
        0.94,
        text,
        transform=ax.transAxes,
        fontsize=6.25,
        va="top",
        bbox={"boxstyle": "round,pad=0.24", "facecolor": "white", "edgecolor": "0.82", "linewidth": 0.5},
    )
    bold_all_text(ax)


def plot_aggregate_panel(ax, motif: dict, xvals: np.ndarray, panel_label: str, summary_rows: list[dict[str, object]]):
    condition_scores = {}
    for condition in motif.get("conditions", []):
        condition_name = str(condition["name"])
        color = CONDITION_COLORS.get(condition_name, "#555555")
        for sample in condition.get("samples", []):
            sample_name = str(sample.get("name", "sample"))
            profile = np.asarray(sample["profile"], dtype=float)
            sample_score = center_flank_score(profile, xvals)
            ax.plot(
                xvals,
                profile,
                color=color,
                linestyle=SAMPLE_STYLES.get(sample_name, "solid"),
                linewidth=0.72,
                alpha=0.48,
                label=sample_name,
                zorder=1,
            )
            summary_rows.append(
                {
                    "panel": "aggregate",
                    "motif": motif.get("name", ""),
                    "motif_id": motif.get("motif_id", ""),
                    "sample": sample_name,
                    "condition": condition_name,
                    "center_minus_flank": sample_score,
                    "selection_role": AGGREGATE_ROLES.get(motif.get("prefix", ""), ""),
                }
            )
        profile = np.asarray(condition["profile"], dtype=float)
        condition_scores[condition_name] = center_flank_score(profile, xvals)
        ax.plot(xvals, profile, color=color, linewidth=1.65, alpha=0.98, label=f"{condition_name} mean", zorder=2)
    stronger = min(condition_scores, key=condition_scores.get) if condition_scores else ""
    score_text = "; ".join(f"{k} {v:.3f}" for k, v in condition_scores.items())
    ax.axvline(0, color="0.35", linewidth=0.75)
    ax.axhline(0, color="0.78", linewidth=0.55, zorder=0)
    ax.set_title(
        f"{panel_label}. {motif.get('name')} ({motif.get('motif_id')})\n"
        f"{AGGREGATE_ROLES.get(motif.get('prefix', ''), '')}; n={int(motif.get('n_sites', 0)):,}; {stronger} deeper",
        fontsize=7.3,
    )
    ax.text(0.02, 0.04, f"center-flank: {score_text}", transform=ax.transAxes, fontsize=5.8, color="0.35", va="bottom")
    ax.set_xlabel("Distance from motif center (bp)")
    ax.set_ylabel("Normalized cut-site signal")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="0.9", linewidth=0.5)
    bold_all_text(ax)


def plot_validation(validation_dir: Path, jaspar: Path, out_prefix: Path) -> None:
    motif_sets = pd.read_csv(validation_dir / "motifs" / "motif_set_summary.tsv", sep="\t")
    id_to_name = jaspar_name_map(jaspar)
    streme = load_streme_motifs(validation_dir, id_to_name)
    counts = load_result_counts(validation_dir)
    rescued = load_rescued_denovo(validation_dir)
    payload = load_report_payload(
        validation_dir
        / "diff_footprints"
        / "denovo_only"
        / "diff_footprints_Bcell_Tcell.html"
    )
    xvals = np.asarray(payload.get("aggregate", {}).get("x", []), dtype=float)
    motifs_by_prefix = motif_lookup(payload)
    missing = [prefix for prefix in SELECTED_AGGREGATES if prefix not in motifs_by_prefix]
    if missing:
        raise ValueError(f"Selected aggregate motifs missing from payload: {', '.join(missing)}")

    apply_style(base_size=8.2)
    fig = plt.figure(figsize=(7.5, 8.9))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.15, 1.12, 1.35, 1.35], hspace=0.58, wspace=0.36)

    ax_workflow = fig.add_subplot(gs[0:2, 0])
    plot_workflow(ax_workflow)

    ax_table = fig.add_subplot(gs[0, 1])
    discovery_rows = plot_discovery_table(ax_table, streme, rescued)

    ax_rescue = fig.add_subplot(gs[1, 1])
    plot_rescue(ax_rescue, counts, rescued)

    summary_rows: list[dict[str, object]] = []
    panel_labels = ["D", "E", "F", "G"]
    for idx, prefix in enumerate(SELECTED_AGGREGATES):
        ax = fig.add_subplot(gs[2 + idx // 2, idx % 2])
        plot_aggregate_panel(ax, motifs_by_prefix[prefix], xvals, panel_labels[idx], summary_rows)

    handles, labels = fig.axes[-1].get_legend_handles_labels()
    unique = {}
    for handle, label in zip(handles, labels):
        unique.setdefault(label, handle)
    fig.legend(unique.values(), unique.keys(), loc="lower center", ncol=6, frameon=False, fontsize=6.6)
    fig.suptitle("De novo motif discovery validates standalone and supplement modes", y=0.995)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.955, bottom=0.075, hspace=0.74, wspace=0.38)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    svg_path = out_prefix.with_suffix(".svg")
    fig.savefig(svg_path, bbox_inches="tight")
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n")

    out_table = out_prefix.with_suffix(".tsv")
    with out_table.open("w", encoding="utf-8") as handle:
        handle.write("# motif_set_summary\n")
        motif_sets.to_csv(handle, sep="\t", index=False)
        handle.write("\n# differential_result_counts\n")
        counts.to_csv(handle, sep="\t", index=False)
        handle.write("\n# discovered_motifs\n")
        streme.to_csv(handle, sep="\t", index=False)
        handle.write("\n# displayed_discovery_rows\n")
        pd.DataFrame(discovery_rows).to_csv(handle, sep="\t", index=False)
        handle.write("\n# highlighted_de_novo_in_restricted_plus_denovo\n")
        rescued[rescued["highlighted"]].to_csv(handle, sep="\t", index=False)
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
