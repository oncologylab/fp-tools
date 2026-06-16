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


def plot_discovery_table(ax, streme: pd.DataFrame, rescued: pd.DataFrame):
    ax.axis("off")
    ax.set_title("B. Representative discovered motifs", loc="left", pad=7, fontsize=8.4)
    selected_names = set(rescued[rescued["highlighted"]]["name"].astype(str))
    rows = []
    representative = {"Bcell_denovo_5", "Bcell_denovo_1", "Tcell_denovo_6", "Tcell_denovo_5"}
    for _, row in streme.iterrows():
        source_prefix = "Bcell" if row["direction"].startswith("B-cell") else "Tcell"
        short_id = str(row["de_novo_motif"]).split("-", 1)[0]
        prefix = f"{source_prefix}_denovo_{short_id}"
        priority = 3 if row["confident_tomtom"] else 0
        if prefix in selected_names:
            priority += 2
        if prefix in representative:
            priority += 10
        if priority > 0:
            item = row.copy()
            item["prefix"] = prefix
            item["priority"] = priority
            rows.append(item)
    if not rows:
        rows = [row for _, row in streme.head(6).iterrows()]
    display = pd.DataFrame(rows).sort_values(["priority", "confident_tomtom", "direction"], ascending=[False, False, True]).head(4)
    y0 = 0.84
    col_x = [0.00, 0.24, 0.49, 0.63]
    headers = ["candidate set", "consensus", "sites", "Tomtom match"]
    for x, h in zip(col_x, headers):
        ax.text(x, y0, h, transform=ax.transAxes, fontsize=7.0, fontweight="bold", va="top")
    ax.plot([0, 1], [y0 - 0.035, y0 - 0.035], color="0.72", linewidth=0.7, transform=ax.transAxes)
    y = y0 - 0.115
    summary_rows = []
    for _, row in display.iterrows():
        source = "B-cell enriched" if row["direction"].startswith("B-cell") else "T-cell enriched"
        q = row["tomtom_q_value"]
        tomtom = row["tomtom_label"]
        if pd.notna(q) and row["confident_tomtom"]:
            tomtom = f"{tomtom.split('(', 1)[0].strip()}-like; q={q:.1e}"
        ax.text(col_x[0], y, source, transform=ax.transAxes, fontsize=6.7, va="top")
        ax.text(col_x[1], y, row["consensus"], transform=ax.transAxes, fontsize=6.7, va="top", family="monospace")
        ax.text(col_x[2], y, f"{int(row['sites']):,}", transform=ax.transAxes, fontsize=6.7, va="top")
        ax.text(col_x[3], y, tomtom, transform=ax.transAxes, fontsize=6.7, va="top")
        summary_rows.append(row.to_dict())
        y -= 0.145
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

def add_rescued_annotations(rescued: pd.DataFrame, streme: pd.DataFrame) -> pd.DataFrame:
    annotations = []
    for _, row in streme.iterrows():
        source_prefix = "Bcell" if row["direction"].startswith("B-cell") else "Tcell"
        short_id = str(row["de_novo_motif"]).split("-", 1)[0]
        annotations.append(
            {
                "name": f"{source_prefix}_denovo_{short_id}",
                "consensus": row["consensus"],
                "streme_sites": row["sites"],
                "tomtom_label": row["tomtom_label"],
                "tomtom_q_value": row["tomtom_q_value"],
                "confident_tomtom": row["confident_tomtom"],
                "candidate_direction": row["direction"],
            }
        )
    annotated = rescued.merge(pd.DataFrame(annotations), on="name", how="left")
    annotated["display_label"] = annotated.apply(rescued_display_label, axis=1)
    return annotated


def rescued_display_label(row: pd.Series) -> str:
    name = str(row.get("name", ""))
    motif_number = name.replace("Bcell_denovo_", "B dn").replace("Tcell_denovo_", "T dn")
    tomtom = str(row.get("tomtom_label", ""))
    if tomtom and tomtom != "no confident match" and tomtom != "nan":
        tomtom = tomtom.split("(", 1)[0].strip()
        return f"{motif_number}: {tomtom}-like"
    consensus = str(row.get("consensus", ""))
    return f"{motif_number}: {consensus}"


def plot_rescued_motifs(ax, counts: pd.DataFrame, rescued: pd.DataFrame, streme: pd.DataFrame) -> pd.DataFrame:
    count_plot = counts[counts["result_set"].isin(["restricted JASPAR", "restricted + de novo"])].copy()
    base = count_plot[count_plot["result_set"].eq("restricted JASPAR")].iloc[0]
    plus = count_plot[count_plot["result_set"].eq("restricted + de novo")].iloc[0]
    ax.set_title("A. Database supplement highlights additional motif families", loc="left", pad=7, fontsize=8.4)
    annotated = add_rescued_annotations(rescued[rescued["highlighted"]].copy(), streme)
    annotated = annotated.sort_values("Bcell_Tcell_change", ascending=True).reset_index(drop=True)
    ax.axis("off")

    tile_data = [
        ("restricted JASPAR", int(base["n_highlighted"]), "#f2f5f9"),
        ("restricted + de novo", int(plus["n_highlighted"]), "#edf7ef"),
        ("de novo-derived", int(plus["n_highlighted_de_novo"]), "#eef4fb"),
    ]
    for idx, (label, value, color) in enumerate(tile_data):
        x0 = 0.01 + idx * 0.325
        ax.add_patch(
            plt.Rectangle(
                (x0, 0.66),
                0.295,
                0.22,
                transform=ax.transAxes,
                facecolor=color,
                edgecolor="0.70",
                linewidth=0.65,
            )
        )
        ax.text(x0 + 0.018, 0.825, label, transform=ax.transAxes, fontsize=6.8, fontweight="bold", va="top")
        ax.text(x0 + 0.018, 0.705, f"{value}", transform=ax.transAxes, fontsize=14.0, fontweight="bold", va="bottom")
    ax.text(0.02, 0.58, "Highlighted de novo-derived motifs in the restricted-database sensitivity test", transform=ax.transAxes, fontsize=6.8, fontweight="bold")

    x = annotated["Bcell_Tcell_change"].to_numpy(dtype=float)
    max_abs = max(float(np.nanmax(np.abs(x))), 1e-6)
    center_x = 0.50
    scale = 0.34 / max_abs
    y_positions = np.linspace(0.48, 0.08, len(annotated))
    ax.plot([center_x, center_x], [0.045, 0.515], color="0.45", linewidth=0.7, transform=ax.transAxes, zorder=0)
    ax.text(0.13, 0.525, "T-cell higher", color=CONDITION_COLORS["Tcell"], transform=ax.transAxes, fontsize=6.7, fontweight="bold")
    ax.text(0.87, 0.525, "B-cell higher", color=CONDITION_COLORS["Bcell"], transform=ax.transAxes, fontsize=6.7, fontweight="bold", ha="right")
    for ypos, (_, row) in zip(y_positions, annotated.iterrows()):
        change = float(row["Bcell_Tcell_change"])
        xpos = center_x + change * scale
        color = CONDITION_COLORS["Bcell"] if change >= 0 else CONDITION_COLORS["Tcell"]
        ax.plot([center_x, xpos], [ypos, ypos], color=color, linewidth=1.35, alpha=0.74, transform=ax.transAxes)
        ax.scatter([xpos], [ypos], s=35, color=color, edgecolor="white", linewidth=0.5, transform=ax.transAxes, zorder=3)
        label_x = 0.02 if xpos >= center_x else 0.60
        ax.text(label_x, ypos, str(row["display_label"]), transform=ax.transAxes, fontsize=6.25, va="center")
        ax.text(0.94, ypos, f"{change:+.2f}", transform=ax.transAxes, fontsize=6.1, va="center", ha="right", color=color)
    bold_all_text(ax)
    return annotated


def aggregate_title(motif: dict) -> tuple[str, str]:
    prefix = str(motif.get("prefix", ""))
    if prefix.startswith("Bcell_denovo_5"):
        return "BATF-like", "B-cell candidates"
    if prefix.startswith("Tcell_denovo_6"):
        return "IKZF2-like", "T-cell candidates"
    if prefix.startswith("Tcell_denovo_4"):
        return "T-cell de novo 4", "T-cell candidates"
    if prefix.startswith("Tcell_denovo_1"):
        return "T-cell de novo 1", "T-cell candidates"
    source = "T-cell candidates" if prefix.startswith("Tcell") else "B-cell candidates"
    return str(motif.get("name", "de novo motif")).replace("_", " "), source


def plot_aggregate_panel(
    ax,
    motif: dict,
    xvals: np.ndarray,
    panel_label: str,
    summary_rows: list[dict[str, object]],
    *,
    show_xlabel: bool,
    show_ylabel: bool,
):
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
                linewidth=0.86,
                alpha=0.55,
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
        ax.plot(xvals, profile, color=color, linewidth=1.45, alpha=0.98, label=f"{condition_name} mean", zorder=2)
    stronger = min(condition_scores, key=condition_scores.get) if condition_scores else ""
    ax.axvline(0, color="0.35", linewidth=0.75)
    ax.axhline(0, color="0.78", linewidth=0.55, zorder=0)
    family, source = aggregate_title(motif)
    ax.set_title(
        f"{panel_label}. {family} ({source})\n"
        f"n={int(motif.get('n_sites', 0)):,}; deeper: {stronger}",
        fontsize=7.0,
    )
    ax.set_xlabel("Distance from motif center (bp)" if show_xlabel else "")
    ax.set_ylabel("Normalized cut-site signal" if show_ylabel else "")
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

    apply_style(base_size=8.0)
    fig = plt.figure(figsize=(7.6, 9.15))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.08, 0.70, 1.28, 1.28], hspace=0.58, wspace=0.34)

    ax_rescue = fig.add_subplot(gs[0, :])
    rescued_rows = plot_rescued_motifs(ax_rescue, counts, rescued, streme)

    ax_table = fig.add_subplot(gs[1, :])
    discovery_rows = plot_discovery_table(ax_table, streme, rescued)

    summary_rows: list[dict[str, object]] = []
    panel_labels = ["C", "D", "E", "F"]
    for idx, prefix in enumerate(SELECTED_AGGREGATES):
        ax = fig.add_subplot(gs[2 + idx // 2, idx % 2])
        plot_aggregate_panel(
            ax,
            motifs_by_prefix[prefix],
            xvals,
            panel_labels[idx],
            summary_rows,
            show_xlabel=idx >= 2,
            show_ylabel=idx % 2 == 0,
        )

    handles, labels = fig.axes[-1].get_legend_handles_labels()
    unique = {}
    for handle, label in zip(handles, labels):
        unique.setdefault(label, handle)
    fig.legend(unique.values(), unique.keys(), loc="lower center", ncol=6, frameon=False, fontsize=6.5)
    fig.suptitle("De novo motifs add testable differential-footprint signal", y=0.992, fontsize=10.5)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.94, bottom=0.064, hspace=0.62, wspace=0.34)

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
        rescued_rows.to_csv(handle, sep="\t", index=False)
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
