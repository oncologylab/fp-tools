#!/usr/bin/env python3
"""Create story figures and a written summary for nutrient-stress TF programs."""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import textwrap
import zipfile
from datetime import date
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PROGRAMS = {
    "AP-1": ["FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND"],
    "ISR/ATF": ["ATF3", "ATF4"],
    "HNF4": ["HNF4A", "HNF4G"],
    "TEAD": ["TEAD3", "TEAD4"],
    "NRF2/redox": ["NFE2L2"],
    "FOXC/ZEB": ["FOXC1", "ZEB1"],
    "KLF": ["KLF2", "KLF15"],
}

TF_ORDER = [
    "FOSL1",
    "FOSL2",
    "FOS",
    "FOSB",
    "JUN",
    "JUNB",
    "JUND",
    "ATF3",
    "ATF4",
    "HNF4A",
    "HNF4G",
    "TEAD3",
    "TEAD4",
    "NFE2L2",
    "FOXC1",
    "ZEB1",
    "KLF15",
    "KLF2",
]

EXPANDED_TF_ORDER = [
    *TF_ORDER,
    "JDP2",
    "BATF",
    "BACH2",
    "ATF7",
    "NFE2L1",
    "MAFG",
    "MAFK",
    "MAF",
    "NRF1",
    "KLF7",
    "KLF10",
    "KLF12",
    "KLF14",
    "KLF16",
    "CDX2",
    "OVOL2",
    "EGR1",
    "EGR3",
    "SREBF2",
    "MLX",
    "MLXIPL",
    "RXRA",
    "RXRB",
    "NR2F1",
    "NR2F6",
]

EXPANDED_PROGRAMS = {
    "AP-1/ATF": ["FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "JDP2", "BATF", "ATF3", "ATF4", "ATF7"],
    "NRF/MAF/redox": ["NFE2L2", "NFE2L1", "MAFG", "MAFK", "MAF", "NRF1"],
    "KLF/SP": ["KLF2", "KLF7", "KLF10", "KLF12", "KLF14", "KLF15", "KLF16", "SP1", "SP2", "SP3", "SP4"],
    "HNF/CDX lineage": ["HNF4A", "HNF4G", "CDX2", "OVOL2"],
    "TEAD/plasticity": ["TEAD3", "TEAD4", "FOXC1", "ZEB1"],
    "Metabolic regulators": ["SREBF2", "MLX", "MLXIPL", "RXRA", "RXRB", "NR2F1", "NR2F6"],
    "Immediate early": ["EGR1", "EGR3", "BACH2", "CREB5", "CREB3L4"],
}

ADDITIONAL_TF_STORIES = {
    "JDP2": "AP-1/ATF-related stress motif partner",
    "BATF": "AP-1-family immune/stress-response motif program",
    "BACH2": "AP-1/MAF-related oxidative and stress-response candidate",
    "ATF7": "ATF-family stress-response candidate",
    "NFE2L1": "NRF-family proteostasis/redox candidate",
    "MAFG": "small-MAF partner in NRF/redox motifs",
    "MAFK": "small-MAF partner in NRF/redox motifs",
    "MAF": "MAF-family stress and differentiation motif program",
    "NRF1": "mitochondrial/proteostasis regulatory candidate",
    "KLF7": "KLF metabolic/plasticity candidate",
    "KLF10": "TGF-beta and stress-responsive KLF candidate",
    "KLF12": "KLF-family plasticity/metabolic candidate",
    "KLF14": "KLF-family metabolic candidate",
    "KLF16": "KLF-family regulatory candidate",
    "CDX2": "intestinal/epithelial lineage-state candidate",
    "OVOL2": "epithelial-state and EMT-restraint candidate",
    "EGR1": "immediate-early stress-response candidate",
    "EGR3": "immediate-early stress-response candidate",
    "SREBF2": "lipid/cholesterol metabolic adaptation candidate",
    "MLX": "carbohydrate and nutrient-sensing network candidate",
    "MLXIPL": "glucose/lipid metabolic regulator candidate",
    "RXRA": "nuclear-receptor metabolic state candidate",
    "RXRB": "nuclear-receptor metabolic state candidate",
    "NR2F1": "nuclear-receptor dormancy/plasticity candidate",
    "NR2F6": "nuclear-receptor stress/plasticity candidate",
}

CORE_TF_SET = set(TF_ORDER)

STRESS_ORDER = ["FBS", "Glc", "Met.Cys", "Gln.Arg", "Gln", "Arg", "BCAA", "Trp", "Lys"]
CELL_ORDER = ["HPAFII", "AsPC1", "Panc1"]

CELL_CONTEXT = {
    "HPAFII": "epithelial-like",
    "AsPC1": "intermediate",
    "Panc1": "mesenchymal-like",
}

PLOT_FONT = "Nimbus Sans"

TARGET_ROWS = [
    {
        "axis": "FOSL1 / AP-1",
        "primary_context": "Panc1 and AsPC1",
        "evidence": "Most recurrent pan-stress footprint program; largest effects in Panc1.",
        "suggested_test": "FOSL1/AP-1 perturbation under glucose, BCAA, Trp, Lys, and Met/Cys stress.",
    },
    {
        "axis": "ATF3 / ATF4",
        "primary_context": "all three cell lines",
        "evidence": "ISR-like nutrient-stress module with recurrent footprint support.",
        "suggested_test": "ATF3/ATF4 knockdown or ISR modulation under amino-acid stress.",
    },
    {
        "axis": "HNF4A / HNF4G",
        "primary_context": "HPAFII and AsPC1",
        "evidence": "Epithelial/classical lineage signal, often reduced under stress.",
        "suggested_test": "HNF4A/HNF4G rescue or CRISPRa to test lineage stabilization.",
    },
    {
        "axis": "TEAD3 / TEAD4",
        "primary_context": "Panc1 and selected HPAFII contexts",
        "evidence": "YAP/TEAD-like plasticity response, strongest in mesenchymal-like context.",
        "suggested_test": "TEAD4 knockdown or TEAD inhibition under glucose, Trp, and FBS stress.",
    },
    {
        "axis": "NFE2L2 / NRF2",
        "primary_context": "Panc1-enriched redox axis",
        "evidence": "Redox/ferroptosis-buffering motif program under amino-acid stress.",
        "suggested_test": "NFE2L2 perturbation with ROS, GSH, lipid-peroxidation, and survival readouts.",
    },
    {
        "axis": "FOXC1 / ZEB1",
        "primary_context": "plasticity-state validation",
        "evidence": "EMT/plasticity-associated footprint layer across several stresses.",
        "suggested_test": "Use as state reporters and dependency tests after AP-1 or TEAD perturbation.",
    },
    {
        "axis": "KLF15 / KLF2",
        "primary_context": "HPAFII and AsPC1",
        "evidence": "Context-dependent metabolic/stress-response marker, especially serum and amino-acid stress.",
        "suggested_test": "Treat KLF15 first as a mechanistic readout, then test perturbation in epithelial-like cells.",
    },
]

FIRST_WAVE_TARGETS = [
    {
        "Target": "TEAD4 / TEAD3",
        "Best setting": "Panc1 0_Lys / 0_Trp",
        "Footprint strength": 0.36,
        "RNA change": -0.07,
        "Priority": "High",
        "Main meaning": "Plasticity-linked YAP/TEAD branch; useful translational anchor for EMT-blocking tests.",
    },
    {
        "Target": "FOSL1 / AP-1",
        "Best setting": "Panc1 0_Trp",
        "Footprint strength": 0.84,
        "RNA change": -1.12,
        "Priority": "High anchor",
        "Main meaning": "Strongest positive-control stress-remodeling axis; broad AP-1 footprint gain in Panc1.",
    },
    {
        "Target": "FOXC1",
        "Best setting": "Panc1 0_Lys",
        "Footprint strength": 0.31,
        "RNA change": -0.32,
        "Priority": "High",
        "Main meaning": "Plasticity-linked candidate that is less overused than ZEB1 and fits EMT-state biology.",
    },
    {
        "Target": "ATF4 / ISR",
        "Best setting": "Panc1 stress panel",
        "Footprint strength": 0.79,
        "RNA change": -0.42,
        "Priority": "High anchor",
        "Main meaning": "Integrated-stress-response branch; likely nutrient-stress survival regulator.",
    },
    {
        "Target": "MAFG",
        "Best setting": "Panc1 0_BCAA / 0_Trp",
        "Footprint strength": 0.70,
        "RNA change": -0.10,
        "Priority": "High novelty",
        "Main meaning": "Small-MAF redox cofactor candidate; strong footprint signal despite modest RNA change.",
    },
    {
        "Target": "NFE2L1",
        "Best setting": "Panc1 0_BCAA / 0_Trp",
        "Footprint strength": 0.70,
        "RNA change": -0.12,
        "Priority": "High novelty",
        "Main meaning": "Proteostasis/redox candidate that tracks the stressed mesenchymal-like state.",
    },
    {
        "Target": "BACH2",
        "Best setting": "Panc1 0_Trp",
        "Footprint strength": 0.78,
        "RNA change": 1.05,
        "Priority": "High novelty",
        "Main meaning": "BACH/MAF-related stress candidate with less established PDAC EMT literature.",
    },
]

TARGET_HANDLING = [
    {
        "Factor": "TEAD4 / TEAD3",
        "Best use": "KO or CRISPRi",
        "Reason": "Candidate plasticity drivers; useful for testing whether nutrient stress reinforces a TEAD-linked mesenchymal state.",
    },
    {
        "Factor": "FOSL1 / AP-1",
        "Best use": "KO or CRISPRi anchor",
        "Reason": "Strong positive-control axis for stress-remodeled chromatin and migration/plasticity programs.",
    },
    {
        "Factor": "FOXC1",
        "Best use": "KO or CRISPRi",
        "Reason": "Plasticity-linked candidate that may reduce invasive or EMT-like outputs.",
    },
    {
        "Factor": "MAFG / NFE2L1 / BACH2",
        "Best use": "CRISPRi first",
        "Reason": "Novel redox/proteostasis branch; partial repression is safer than full KO for broad homeostatic factors.",
    },
    {
        "Factor": "HNF4A / HNF4G / CDX2",
        "Best use": "Do not KO for EMT blocking",
        "Reason": "These look like epithelial/classical identity factors whose loss may promote rather than block transition.",
    },
    {
        "Factor": "OVOL2",
        "Best use": "Rescue or CRISPRa",
        "Reason": "May act as an epithelial brake in HPAFII under Arg/Gln.Arg stress.",
    },
    {
        "Factor": "ZEB1",
        "Best use": "Lower-priority marker",
        "Reason": "Canonical EMT factor, but not the dominant stress-induced driver in these tables.",
    },
]

NOVELTY_TARGETS = [
    {
        "Candidate": "MAFG",
        "Where strongest": "Panc1 BCAA/Trp/Lys/Glc",
        "Why interesting": "Small-MAF redox cofactor with strong footprint-dominant signal.",
    },
    {
        "Candidate": "NFE2L1",
        "Where strongest": "Panc1 BCAA/Trp/Lys/Glc",
        "Why interesting": "Proteostasis/redox branch that mirrors MAFG in the mesenchymal-like stress state.",
    },
    {
        "Candidate": "BACH2",
        "Where strongest": "Panc1 Trp; AsPC1 glucose/Gln.Arg/Met.Cys",
        "Why interesting": "BACH/MAF-related oxidative-stress candidate with less direct PDAC EMT precedent.",
    },
    {
        "Candidate": "FOXC1",
        "Where strongest": "Panc1 and HPAFII/AsPC1 selected stress states",
        "Why interesting": "Plasticity-linked factor that is less saturated than ZEB1 as a follow-up target.",
    },
    {
        "Candidate": "KLF12",
        "Where strongest": "AsPC1 0_BCAA",
        "Why interesting": "RNA rises while motif footprint falls, suggesting a compensation loop.",
    },
    {
        "Candidate": "KLF14",
        "Where strongest": "AsPC1 0_BCAA",
        "Why interesting": "Similar BCAA-sensitive compensation pattern to KLF12.",
    },
]

REFERENCES = [
    (
        "Aggrey-Fynn et al., 2025",
        "Therapeutic targeting of FOSL1 and RELA-dependent transcriptional mechanisms to suppress pancreatic cancer metastasis.",
        "Cell Death & Disease.",
        "https://www.nature.com/articles/s41419-025-07810-x",
    ),
    (
        "FOSL1 pancreatic cancer study, 2025",
        "FOSL1 drives malignant progression of pancreatic cancer by regulating stemness, invasion, metastasis, and drug resistance.",
        "Journal of Translational Medicine.",
        "https://link.springer.com/article/10.1186/s12967-025-06304-w",
    ),
    (
        "Palam et al., 2015",
        "The integrated stress response is critical for gemcitabine resistance in pancreatic cancer cells.",
        "PMC article.",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC4632294/",
    ),
    (
        "ISR review, 2021",
        "Targeting the Integrated Stress Response in Cancer Therapy.",
        "Frontiers in Pharmacology.",
        "https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2021.747837/full",
    ),
    (
        "Kloesch et al., 2022",
        "A GATA6-centred gene regulatory network involving HNFs and DeltaNp63 controls plasticity and immune escape in pancreatic cancer.",
        "Gut.",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC9733634/",
    ),
    (
        "Brunton et al., 2020",
        "HNF4A and GATA6 loss reveals therapeutically actionable subtypes in pancreatic cancer.",
        "Cell Reports.",
        "https://pubmed.ncbi.nlm.nih.gov/32402285/",
    ),
    (
        "Holden et al., 2023",
        "An allosteric pan-TEAD inhibitor blocks oncogenic YAP/TAZ signaling.",
        "Nature Chemical Biology.",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC10293011/",
    ),
    (
        "NRF2 ferroptosis review, 2024",
        "NFE2L2 and ferroptosis resistance in cancer therapy.",
        "Cancer Drug Resistance.",
        "https://www.oaepublish.com/articles/cdr.2024.123",
    ),
    (
        "Krebs et al., 2017",
        "The EMT-activator Zeb1 is a key factor for cell plasticity and promotes metastasis in pancreatic cancer.",
        "Nature Cell Biology.",
        "https://pubmed.ncbi.nlm.nih.gov/28414315/",
    ),
    (
        "ZEB1 review",
        "ZEB1 and the miR-200 family are central players in EMT, invasion, and metastasis.",
        "PMC review.",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC3837326/",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--handoff-dir",
        default="data/public/processed/nutrient_tf_handoff_20260629",
        help="Directory created by build_nutrient_tf_handoff_package.py.",
    )
    parser.add_argument("--outdir", default=None, help="Output directory for story figures.")
    return parser.parse_args()


def save_fig(fig: mpl.figure.Figure, outdir: Path, stem: str) -> None:
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(outdir / f"{stem}.{suffix}", bbox_inches="tight", dpi=300)
    plt.close(fig)


def wrap_label(text: str, width: int = 24) -> str:
    return "\n".join(textwrap.wrap(str(text), width=width, break_long_words=False))


def contrast_text_color(value: float, cmap_name: str, norm: mpl.colors.Normalize) -> str:
    if pd.isna(value):
        return "black"
    rgba = plt.get_cmap(cmap_name)(norm(value))
    luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
    return "white" if luminance < 0.48 else "black"


def style_axis_text(ax: mpl.axes.Axes, size: int = 9) -> None:
    for item in [ax.title, ax.xaxis.label, ax.yaxis.label, *ax.get_xticklabels(), *ax.get_yticklabels()]:
        item.set_fontsize(max(size, item.get_fontsize()))
        item.set_fontweight("bold")
        item.set_fontfamily(PLOT_FONT)
        item.set_color("black")


def style_colorbar(cbar: mpl.colorbar.Colorbar, size: int = 9) -> None:
    cbar.ax.yaxis.label.set_fontsize(size)
    cbar.ax.yaxis.label.set_fontweight("bold")
    cbar.ax.yaxis.label.set_fontfamily(PLOT_FONT)
    for tick in cbar.ax.get_yticklabels():
        tick.set_fontsize(size)
        tick.set_fontweight("bold")
        tick.set_fontfamily(PLOT_FONT)


def load_data(handoff: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary = pd.read_csv(handoff / "tables" / "tf_condition_summary.tsv", sep="\t")
    recurrence = pd.read_csv(handoff / "tables" / "cross_cellline_recurrence.tsv", sep="\t")
    top = pd.read_csv(handoff / "tables" / "top_candidates_by_cellline_condition.tsv", sep="\t")
    for df in (summary, recurrence, top):
        if "stress_type" in df.columns:
            df["stress_type"] = df["stress_type"].replace({"Gln": "Gln.Arg", "Arg": "Gln.Arg"})
    return summary, recurrence, top


def condition_order(df: pd.DataFrame) -> list[str]:
    def key(cond: str) -> tuple[int, float, str]:
        row = df[df["condition"].eq(cond)].iloc[0]
        stress = row["stress_type"]
        stress_idx = STRESS_ORDER.index(stress) if stress in STRESS_ORDER else 99
        dose = row["dose"]
        dose_key = -float(dose) if pd.notna(dose) else 0
        return (stress_idx, dose_key, cond)

    return sorted(df["condition"].unique(), key=key)


def selected_tf_table(summary: pd.DataFrame, tf_order: list[str] | None = None) -> pd.DataFrame:
    if tf_order is None:
        tf_order = TF_ORDER
    present = [tf for tf in tf_order if tf in set(summary["tf_gene"])]
    sub = summary[summary["tf_gene"].isin(present)].copy()
    sub["tf_gene"] = pd.Categorical(sub["tf_gene"], present, ordered=True)
    sub["cell_line"] = pd.Categorical(sub["cell_line"], CELL_ORDER, ordered=True)
    sub["rna_abs"] = sub["rna_log2fc_max_abs"].abs()
    return sub.sort_values(["cell_line", "tf_gene", "condition"])


def plot_tf_evidence_map(
    summary: pd.DataFrame,
    outdir: Path,
    tf_order: list[str],
    stem: str,
    title: str,
    footprint_cap: float = 0.35,
    rna_cap: float = 2.5,
    figsize: tuple[float, float] | None = None,
    marker_scale: tuple[float, float] = (24, 230),
) -> pd.DataFrame:
    source = selected_tf_table(summary, tf_order)
    present_tfs = [str(tf) for tf in source["tf_gene"].cat.categories]
    if figsize is None:
        figsize = (18.0, max(12.0, 2.5 + 0.23 * len(present_tfs) * len(CELL_ORDER)))
    fig, axes = plt.subplots(3, 1, figsize=figsize, sharey=True)
    cmap = plt.get_cmap("coolwarm")
    rna_cmap = cmap
    norm = mpl.colors.TwoSlopeNorm(vmin=-footprint_cap, vcenter=0, vmax=footprint_cap)
    rna_norm = mpl.colors.TwoSlopeNorm(vmin=-rna_cap, vcenter=0, vmax=rna_cap)
    size_min, size_max = marker_scale
    fdr_scaled = source["neg_log10_best_fdr"].clip(0, 220)
    if fdr_scaled.max() > fdr_scaled.min():
        source["_plot_size"] = size_min + (fdr_scaled - fdr_scaled.min()) / (fdr_scaled.max() - fdr_scaled.min()) * (size_max - size_min)
    else:
        source["_plot_size"] = 80
    for row_idx, (ax, cell) in enumerate(zip(axes, CELL_ORDER)):
        cell_df = source[source["cell_line"].eq(cell)]
        conds = condition_order(cell_df)
        x_map = {c: i for i, c in enumerate(conds)}
        y_map = {tf: i for i, tf in enumerate(present_tfs)}
        xs = cell_df["condition"].map(x_map).astype(float)
        ys = cell_df["tf_gene"].map(y_map).astype(float)
        ax.scatter(
            xs - 0.13,
            ys,
            c=cell_df["delta_fp_at_max_abs"],
            cmap=cmap,
            norm=norm,
            s=cell_df["_plot_size"],
            edgecolor="black",
            linewidth=0.25,
            alpha=0.95,
        )
        ax.scatter(
            xs + 0.20,
            ys,
            c=cell_df["rna_log2fc_max_abs"],
            cmap=rna_cmap,
            norm=rna_norm,
            s=82,
            marker="s",
            edgecolor="black",
            linewidth=0.2,
            alpha=0.95,
        )
        ax.set_xlim(-0.6, len(conds) - 0.4)
        ax.set_ylim(len(present_tfs) - 0.4, -0.6)
        ax.set_xticks(range(len(conds)))
        ax.set_xticklabels(conds, rotation=45, ha="right", fontsize=9, fontweight="bold")
        ax.set_yticks(range(len(present_tfs)))
        ax.set_yticklabels(present_tfs, fontsize=9, fontweight="bold")
        ax.grid(axis="both", color="#e5e5e5", linewidth=0.5)
        ax.set_axisbelow(True)
        for spine in ax.spines.values():
            spine.set_color("#cccccc")
        ax.set_title(f"{cell} ({CELL_CONTEXT[cell]})", loc="left", fontweight="bold")
        if row_idx == len(CELL_ORDER) - 1:
            ax.set_xlabel("Nutrient-stress condition vs 10_FBS_Ctrl")
        style_axis_text(ax, 9)
    fig.supylabel("TF or TF-family member", x=0.03)
    cax1 = fig.add_axes([0.89, 0.58, 0.016, 0.28])
    cax2 = fig.add_axes([0.94, 0.58, 0.016, 0.28])
    style_colorbar(fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax1, label=f"Footprint delta_fp\n(capped at +/-{footprint_cap:g})", extend="both"), 9)
    style_colorbar(fig.colorbar(mpl.cm.ScalarMappable(norm=rna_norm, cmap=rna_cmap), cax=cax2, label=f"RNA log2FC\n(capped at +/-{rna_cap:g})", extend="both"), 9)
    marker_handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor="#d9d9d9", markeredgecolor="black", markersize=8, label="Footprint"),
        mpl.lines.Line2D([], [], marker="s", linestyle="", markerfacecolor="#d9d9d9", markeredgecolor="black", markersize=7, label="RNA"),
    ]
    fig.legend(
        handles=marker_handles,
        title="Marker type",
        loc="center left",
        bbox_to_anchor=(0.89, 0.45),
        frameon=False,
        borderaxespad=0,
    )
    size_handles = [
        mpl.lines.Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markerfacecolor="#d9d9d9",
            markeredgecolor="black",
            markeredgewidth=0.35,
            markersize=ms,
            label=label,
        )
        for ms, label in [(4.5, "lower"), (7.5, "medium"), (11.0, "higher")]
    ]
    fig.legend(
        handles=size_handles,
        title="Footprint FDR evidence\ncircle size = -log10(FDR)",
        loc="center left",
        bbox_to_anchor=(0.89, 0.24),
        frameon=False,
        labelspacing=1.2,
        borderaxespad=0,
    )
    fig.suptitle(title, fontsize=18, fontweight="bold", y=0.985)
    fig.subplots_adjust(left=0.08, right=0.86, top=0.95, bottom=0.06, hspace=0.36)
    save_fig(fig, outdir, stem)
    return source.drop(columns=["_plot_size"], errors="ignore")


def plot_integrated_dotplot(summary: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    return plot_tf_evidence_map(
        summary,
        outdir,
        TF_ORDER,
        "fig1_integrated_tf_evidence_map",
        "Integrated TF evidence across nutrient stress",
        footprint_cap=0.35,
        rna_cap=2.5,
        figsize=(17.5, 15),
    )


def plot_expanded_tf_evidence_map(summary: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    return plot_tf_evidence_map(
        summary,
        outdir,
        EXPANDED_TF_ORDER,
        "fig6_expanded_tf_evidence_map",
        "Expanded TF evidence map across nutrient stress",
        footprint_cap=0.35,
        rna_cap=2.5,
        figsize=(19.5, 28),
        marker_scale=(14, 145),
    )


def build_program_table(summary: pd.DataFrame, programs: dict[str, list[str]] | None = None) -> pd.DataFrame:
    if programs is None:
        programs = PROGRAMS
    rows = []
    for (cell, cond, stress, dose), g0 in summary.groupby(["cell_line", "condition", "stress_type", "dose"], dropna=False):
        for program, genes in programs.items():
            g = g0[g0["tf_gene"].isin(genes)]
            if g.empty:
                continue
            idx = g["max_abs_delta_fp"].idxmax()
            best = g.loc[idx]
            rows.append(
                {
                    "cell_line": cell,
                    "condition": cond,
                    "stress_type": stress,
                    "dose": dose,
                    "program": program,
                    "dominant_tf": best["tf_gene"],
                    "signed_delta_fp": best["delta_fp_at_max_abs"],
                    "max_abs_delta_fp": best["max_abs_delta_fp"],
                    "best_fdr": g["best_fdr"].min(),
                    "rna_log2fc": best["rna_log2fc_max_abs"],
                    "evidence_class": best["dominant_evidence_class"],
                }
            )
    return pd.DataFrame(rows)


def plot_program_heatmap(program_df: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    fig, axes = plt.subplots(1, 3, figsize=(32, 8), sharey=True)
    vmax = np.nanpercentile(program_df["signed_delta_fp"].abs(), 98)
    cmap = "coolwarm"
    norm = mpl.colors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    for ax, cell in zip(axes, CELL_ORDER):
        g = program_df[program_df["cell_line"].eq(cell)]
        conds = condition_order(g)
        mat = g.pivot_table(index="program", columns="condition", values="signed_delta_fp", aggfunc="first").reindex(list(PROGRAMS), columns=conds)
        labels = g.pivot_table(index="program", columns="condition", values="dominant_tf", aggfunc="first").reindex(list(PROGRAMS), columns=conds).fillna("")
        sns.heatmap(
            mat,
            ax=ax,
            cmap=cmap,
            center=0,
            vmin=-vmax,
            vmax=vmax,
            cbar=cell == CELL_ORDER[-1],
            cbar_kws={"label": "dominant program delta_fp"},
            linewidths=0.6,
            linecolor="white",
        )
        for y, prog in enumerate(mat.index):
            for x, cond in enumerate(mat.columns):
                label = labels.loc[prog, cond]
                if label:
                    val = mat.loc[prog, cond]
                    ax.text(
                        x + 0.5,
                        y + 0.5,
                        label,
                        ha="center",
                        va="center",
                        fontsize=9,
                        color=contrast_text_color(val, cmap, norm),
                        fontweight="bold",
                    )
        ax.set_title(f"{cell} ({CELL_CONTEXT[cell]})", fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9, fontweight="bold")
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=9, fontweight="bold")
        style_axis_text(ax, 9)
    fig.suptitle("Dominant nutrient-stress TF programs", fontsize=18, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save_fig(fig, outdir, "fig2_program_level_stress_heatmap")
    return program_df


def plot_expanded_program_heatmap(program_df: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    fig, axes = plt.subplots(1, 3, figsize=(34, 10), sharey=True)
    vmax = min(0.35, max(0.2, np.nanpercentile(program_df["signed_delta_fp"].abs(), 95)))
    cmap = "coolwarm"
    norm = mpl.colors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    for ax, cell in zip(axes, CELL_ORDER):
        g = program_df[program_df["cell_line"].eq(cell)]
        conds = condition_order(g)
        mat = (
            g.pivot_table(index="program", columns="condition", values="signed_delta_fp", aggfunc="first")
            .reindex(list(EXPANDED_PROGRAMS), columns=conds)
        )
        labels = (
            g.pivot_table(index="program", columns="condition", values="dominant_tf", aggfunc="first")
            .reindex(list(EXPANDED_PROGRAMS), columns=conds)
            .fillna("")
        )
        sns.heatmap(
            mat,
            ax=ax,
            cmap=cmap,
            center=0,
            vmin=-vmax,
            vmax=vmax,
            cbar=cell == CELL_ORDER[-1],
            cbar_kws={"label": f"dominant TF delta_fp\n(capped at +/-{vmax:g})", "extend": "both"},
            linewidths=0.6,
            linecolor="white",
        )
        for y, prog in enumerate(mat.index):
            for x, cond in enumerate(mat.columns):
                label = labels.loc[prog, cond]
                if label:
                    val = mat.loc[prog, cond]
                    ax.text(
                        x + 0.5,
                        y + 0.5,
                        label,
                        ha="center",
                        va="center",
                        fontsize=9,
                        color=contrast_text_color(val, cmap, norm),
                        fontweight="bold",
                    )
        ax.set_title(f"{cell} ({CELL_CONTEXT[cell]})", fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9, fontweight="bold")
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=9, fontweight="bold")
        style_axis_text(ax, 9)
    fig.suptitle("Expanded nutrient-stress TF program map", fontsize=18, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    save_fig(fig, outdir, "fig7_expanded_program_heatmap")
    return program_df


def plot_recurrence(program_df: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    sig = program_df[program_df["best_fdr"] <= 0.05].copy()
    rows = []
    for (program, stress), g in sig.groupby(["program", "stress_type"]):
        rows.append(
            {
                "program": program,
                "stress_type": stress,
                "n_cell_lines": g["cell_line"].nunique(),
                "max_abs_delta_fp": g["max_abs_delta_fp"].max(),
                "cell_lines": ";".join(sorted(g["cell_line"].unique())),
                "dominant_tfs": ";".join(sorted(g["dominant_tf"].unique())),
            }
        )
    recur = pd.DataFrame(rows)
    stresses = [s for s in ["FBS", "Glc", "Met.Cys", "Gln.Arg", "BCAA", "Trp", "Lys"] if s in recur["stress_type"].unique()]
    mat = recur.pivot_table(index="program", columns="stress_type", values="n_cell_lines", aggfunc="max").reindex(list(PROGRAMS), columns=stresses)
    lab = recur.pivot_table(index="program", columns="stress_type", values="max_abs_delta_fp", aggfunc="max").reindex(list(PROGRAMS), columns=stresses)
    fig, ax = plt.subplots(figsize=(11, 7.5))
    cmap = "YlGnBu"
    norm = mpl.colors.Normalize(vmin=0, vmax=3)
    sns.heatmap(
        mat,
        ax=ax,
        cmap=cmap,
        vmin=0,
        vmax=3,
        linewidths=0.8,
        linecolor="white",
        cbar_kws={"label": "cell lines with FDR <= 0.05"},
    )
    for y, prog in enumerate(mat.index):
        for x, stress in enumerate(mat.columns):
            if pd.notna(mat.loc[prog, stress]):
                val = mat.loc[prog, stress]
                ax.text(
                    x + 0.5,
                    y + 0.5,
                    f"{int(val)}\n|dFP| {lab.loc[prog, stress]:.2f}",
                    ha="center",
                    va="center",
                    fontsize=9,
                    color=contrast_text_color(val, cmap, norm),
                    fontweight="bold",
                )
    ax.set_xlabel("Stress type")
    ax.set_ylabel("TF program")
    ax.set_title("Cross-cell-line recurrence of nutrient-stress TF programs", fontweight="bold", fontsize=16)
    style_axis_text(ax, 9)
    fig.tight_layout()
    save_fig(fig, outdir, "fig3_cross_cellline_recurrence_matrix")
    return recur


def plot_lineage_model(program_df: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    lineage = (
        program_df.groupby(["cell_line", "program"], as_index=False)
        .agg(max_abs_delta_fp=("max_abs_delta_fp", "max"), max_abs_rna=("rna_log2fc", lambda x: np.nanmax(np.abs(x))), dominant_tf=("dominant_tf", lambda x: x.value_counts().index[0]))
    )
    mat = lineage.pivot(index="cell_line", columns="program", values="max_abs_delta_fp").reindex(CELL_ORDER, columns=list(PROGRAMS))
    labels = lineage.pivot(index="cell_line", columns="program", values="dominant_tf").reindex(CELL_ORDER, columns=list(PROGRAMS)).fillna("")
    fig, ax = plt.subplots(figsize=(12.5, 4.6))
    cmap = "magma_r"
    norm = mpl.colors.Normalize(vmin=np.nanmin(mat.values), vmax=np.nanmax(mat.values))
    sns.heatmap(
        mat,
        ax=ax,
        cmap=cmap,
        linewidths=0.8,
        linecolor="white",
        cbar_kws={"label": "max |delta_fp|"},
    )
    for y, cell in enumerate(mat.index):
        for x, prog in enumerate(mat.columns):
            val = mat.loc[cell, prog]
            ax.text(
                x + 0.5,
                y + 0.5,
                labels.loc[cell, prog],
                ha="center",
                va="center",
                fontsize=9,
                color=contrast_text_color(val, cmap, norm),
                fontweight="bold",
            )
    ax.set_title("Lineage-state interpretation of nutrient-stress TF programs", fontweight="bold", fontsize=16)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_yticklabels([f"{c}\n{CELL_CONTEXT[c]}" for c in mat.index], rotation=0)
    style_axis_text(ax, 9)
    fig.tight_layout(rect=[0, 0, 1, 1])
    save_fig(fig, outdir, "fig4_lineage_state_model")
    return lineage


def extract_mode_source(summary: pd.DataFrame) -> pd.DataFrame:
    mode_specs = [
        {
            "mode": "HPAFII epithelial buffering",
            "cell_line": "HPAFII",
            "genes": ["HNF4A", "HNF4G", "CDX2", "OVOL2", "TEAD3", "FOXC1"],
            "conditions": ["10_Arg", "0_Gln.Arg", "5_Gln.Arg", "10_Gln.Arg", "25_Gln.Arg"],
        },
        {
            "mode": "AsPC1 BCAA compensation",
            "cell_line": "AsPC1",
            "genes": ["KLF7", "KLF10", "KLF12", "KLF14", "NRF1", "ATF7", "JDP2", "CDX2"],
            "conditions": ["0_BCAA"],
        },
        {
            "mode": "Panc1 mesenchymal lock-in",
            "cell_line": "Panc1",
            "genes": ["FOSL1", "FOS", "JUNB", "ATF3", "JDP2", "MAFG", "NFE2L1", "NFE2L2", "TEAD4"],
            "conditions": ["0_BCAA", "0_Trp", "0_Lys", "0_Glc", "6_Met.Cys"],
        },
    ]
    rows = []
    for spec in mode_specs:
        sub = summary[
            summary["cell_line"].eq(spec["cell_line"])
            & summary["tf_gene"].isin(spec["genes"])
            & summary["condition"].isin(spec["conditions"])
        ].copy()
        sub["mode"] = spec["mode"]
        sub["gene_order"] = sub["tf_gene"].map({g: i for i, g in enumerate(spec["genes"])})
        sub["condition_order"] = sub["condition"].map({c: i for i, c in enumerate(spec["conditions"])})
        rows.append(sub)
    return pd.concat(rows, ignore_index=True).sort_values(["mode", "gene_order", "condition_order"])


def plot_three_adaptation_modes(summary: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    source = extract_mode_source(summary)
    specs = [
        ("HPAFII epithelial buffering", "HNF4/CDX loss with OVOL2 gain"),
        ("AsPC1 BCAA compensation", "KLF RNA-up / footprint-loss state"),
        ("Panc1 mesenchymal lock-in", "AP-1 + ISR + NRF/MAF + TEAD"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(24, 8), gridspec_kw={"width_ratios": [1.35, 0.72, 1.4]})
    fp_norm = mpl.colors.TwoSlopeNorm(vmin=-0.55, vcenter=0, vmax=0.85)
    rna_norm = mpl.colors.TwoSlopeNorm(vmin=-3.0, vcenter=0, vmax=3.0)
    cmap = "coolwarm"
    for ax, (mode, subtitle) in zip(axes, specs):
        g = source[source["mode"].eq(mode)].copy()
        genes = list(g.sort_values("gene_order")["tf_gene"].drop_duplicates())
        conds = list(g.sort_values("condition_order")["condition"].drop_duplicates())
        x_map = {c: i for i, c in enumerate(conds)}
        y_map = {tf: i for i, tf in enumerate(genes)}
        xs = g["condition"].map(x_map).astype(float)
        ys = g["tf_gene"].map(y_map).astype(float)
        ax.scatter(
            xs - 0.14,
            ys,
            c=g["delta_fp_at_max_abs"],
            cmap=cmap,
            norm=fp_norm,
            s=170,
            marker="o",
            edgecolor="black",
            linewidth=0.35,
        )
        ax.scatter(
            xs + 0.18,
            ys,
            c=g["rna_log2fc_max_abs"],
            cmap=cmap,
            norm=rna_norm,
            s=150,
            marker="s",
            edgecolor="black",
            linewidth=0.35,
        )
        for _, row in g.iterrows():
            if row["best_fdr"] <= 0.05 and row["max_abs_delta_fp"] >= 0.45:
                ax.text(
                    x_map[row["condition"]],
                    y_map[row["tf_gene"]] - 0.34,
                    f"{row['max_abs_delta_fp']:.2f}",
                    ha="center",
                    va="center",
                    fontsize=9,
                    fontweight="bold",
                    color="black",
                )
        ax.set_xticks(range(len(conds)))
        ax.set_xticklabels(conds, rotation=45, ha="right", fontsize=9, fontweight="bold")
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes, fontsize=9, fontweight="bold")
        ax.set_xlim(-0.65, len(conds) - 0.35)
        ax.set_ylim(len(genes) - 0.6, -0.6)
        ax.grid(axis="both", color="#e5e5e5", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_title(f"{mode}\n{subtitle}", fontsize=14, fontweight="bold")
        style_axis_text(ax, 9)
    cax1 = fig.add_axes([0.89, 0.55, 0.015, 0.28])
    cax2 = fig.add_axes([0.94, 0.55, 0.015, 0.28])
    style_colorbar(fig.colorbar(mpl.cm.ScalarMappable(norm=fp_norm, cmap=cmap), cax=cax1, label="Footprint delta_fp", extend="both"), 9)
    style_colorbar(fig.colorbar(mpl.cm.ScalarMappable(norm=rna_norm, cmap=cmap), cax=cax2, label="RNA log2FC", extend="both"), 9)
    handles = [
        mpl.lines.Line2D([], [], marker="o", linestyle="", markerfacecolor="#d9d9d9", markeredgecolor="black", markersize=9, label="Footprint"),
        mpl.lines.Line2D([], [], marker="s", linestyle="", markerfacecolor="#d9d9d9", markeredgecolor="black", markersize=8, label="RNA"),
    ]
    fig.legend(handles=handles, title="Marker type", loc="center left", bbox_to_anchor=(0.89, 0.36), frameon=False)
    fig.suptitle("Three nutrient-stress adaptation modes", fontsize=18, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.06, right=0.86, top=0.84, bottom=0.16, wspace=0.32)
    save_fig(fig, outdir, "fig5_three_adaptation_modes")
    return source.drop(columns=["gene_order", "condition_order"], errors="ignore")


def plot_target_table(summary: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    rows = []
    for row in TARGET_ROWS:
        genes = []
        for token in row["axis"].replace("/", " ").split():
            token = token.strip()
            if token in summary["tf_gene"].unique():
                genes.append(token)
        g = summary[summary["tf_gene"].isin(genes)]
        if g.empty:
            best = {"cell_line": "", "condition": "", "max_abs_delta_fp": np.nan, "rna_log2fc_max_abs": np.nan, "best_fdr": np.nan}
        else:
            best = g.loc[g["max_abs_delta_fp"].idxmax()].to_dict()
        rows.append(
            {
                "TF axis": row["axis"],
                "Strongest context": f"{best.get('cell_line', '')} {best.get('condition', '')}".strip(),
                "Max |dFP|": best.get("max_abs_delta_fp", np.nan),
                "RNA log2FC": best.get("rna_log2fc_max_abs", np.nan),
                "Best FDR": best.get("best_fdr", np.nan),
                "Interpretation": row["evidence"],
            }
        )
    table = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(16, 7.5))
    ax.axis("off")
    display = table.copy()
    display["Max |dFP|"] = display["Max |dFP|"].map(lambda x: "" if pd.isna(x) else f"{x:.2f}")
    display["RNA log2FC"] = display["RNA log2FC"].map(lambda x: "" if pd.isna(x) else f"{x:.2f}")
    display["Best FDR"] = display["Best FDR"].map(lambda x: "" if pd.isna(x) else f"{x:.1e}")
    for col in ["Interpretation"]:
        display[col] = display[col].map(lambda x: wrap_label(x, 42))
    mpl_table = ax.table(cellText=display.values, colLabels=display.columns, cellLoc="left", colLoc="left", loc="center")
    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(9)
    mpl_table.scale(1, 2.05)
    widths = [0.15, 0.17, 0.10, 0.10, 0.12, 0.36]
    for (r, c), cell in mpl_table.get_celld().items():
        cell.set_edgecolor("#dddddd")
        cell.set_text_props(weight="bold", color="black", fontfamily=PLOT_FONT, fontsize=9)
        if c < len(widths):
            cell.set_width(widths[c])
        if r == 0:
            cell.set_facecolor("#efefef")
        else:
            cell.set_facecolor("#ffffff" if r % 2 else "#f8f8f8")
    ax.set_title("Prioritized nutrient-stress TF target axes", fontweight="bold", fontsize=16, pad=18)
    save_fig(fig, outdir, "fig5_prioritized_tf_target_table")
    return table


def plot_additional_candidate_table(summary: pd.DataFrame, outdir: Path) -> pd.DataFrame:
    rows = []
    for tf, story in ADDITIONAL_TF_STORIES.items():
        g = summary[summary["tf_gene"].eq(tf)]
        if g.empty:
            continue
        best = g.loc[g["max_abs_delta_fp"].idxmax()]
        rows.append(
            {
                "TF": tf,
                "Strongest context": f"{best['cell_line']} {best['condition']}",
                "Max |dFP|": best["max_abs_delta_fp"],
                "Signed dFP": best["delta_fp_at_max_abs"],
                "RNA log2FC": best["rna_log2fc_max_abs"],
                "Best FDR": best["best_fdr"],
                "Possible story": story,
            }
        )
    table = pd.DataFrame(rows).sort_values(["Max |dFP|", "Best FDR"], ascending=[False, True]).head(22)
    fig, ax = plt.subplots(figsize=(18, 11))
    ax.axis("off")
    display = table.copy()
    display["Max |dFP|"] = display["Max |dFP|"].map(lambda x: f"{x:.2f}")
    display["Signed dFP"] = display["Signed dFP"].map(lambda x: f"{x:.2f}")
    display["RNA log2FC"] = display["RNA log2FC"].map(lambda x: f"{x:.2f}")
    display["Best FDR"] = display["Best FDR"].map(lambda x: f"{x:.1e}")
    display["Possible story"] = display["Possible story"].map(lambda x: wrap_label(x, 36))
    mpl_table = ax.table(cellText=display.values, colLabels=display.columns, cellLoc="left", colLoc="left", loc="center")
    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(9)
    mpl_table.scale(1, 1.78)
    widths = [0.07, 0.15, 0.075, 0.075, 0.08, 0.10, 0.45]
    for (r, c), cell in mpl_table.get_celld().items():
        cell.set_edgecolor("#dddddd")
        cell.set_text_props(weight="bold", color="black", fontfamily=PLOT_FONT, fontsize=9)
        if c < len(widths):
            cell.set_width(widths[c])
        if r == 0:
            cell.set_facecolor("#efefef")
        else:
            cell.set_facecolor("#ffffff" if r % 2 else "#f8f8f8")
    ax.set_title("Additional TF candidates from expanded nutrient-stress scan", fontweight="bold", fontsize=16, pad=18)
    save_fig(fig, outdir, "fig8_additional_candidate_tf_table")
    return table


def format_story_table(table: pd.DataFrame, columns: list[str], max_rows: int | None = None) -> str:
    """Return a Markdown table that Pandoc converts to a native DOCX table."""
    if max_rows is not None:
        table = table.head(max_rows)
    display = table.loc[:, columns].copy()
    for col in display.columns:
        if col in {"Max |dFP|", "Signed dFP", "RNA log2FC", "Footprint strength", "Footprint direction", "RNA change"}:
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else f"{float(x):.2f}")
        elif col in {"Best FDR", "FDR"}:
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else f"{float(x):.1e}")
        else:
            display[col] = display[col].map(lambda x: "" if pd.isna(x) else str(x))
    display.columns = [col.replace("|", "\\|") for col in display.columns]

    def clean_cell(value: object) -> str:
        text = str(value).replace("\n", " ").replace("|", "\\|")
        return " ".join(text.split())

    header = [clean_cell(col) for col in display.columns]
    rows = [[clean_cell(value) for value in row] for row in display.to_numpy()]
    lines = [
        "| " + " | ".join(header) + " |",
        "| " + " | ".join(["---"] * len(header)) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def write_markdown(outdir: Path, source_tables: dict[str, pd.DataFrame]) -> Path:
    stamp = date.today().strftime("%Y%m%d")
    md = outdir / f"nutrient_tfs_{stamp}.md"
    first_wave = source_tables["first_wave_table"].copy()
    handling = source_tables["target_handling_table"].copy()
    novelty = source_tables["novelty_table"].copy()
    first_wave_md = format_story_table(
        first_wave,
        ["Target", "Best setting", "Footprint strength", "RNA change", "Priority", "Main meaning"],
    )
    handling_md = format_story_table(
        handling,
        ["Factor", "Best use", "Reason"],
    )
    novelty_md = format_story_table(
        novelty,
        ["Candidate", "Where strongest", "Why interesting"],
    )
    text = f"""# Nutrient-Stress TF Programs in Pancreatic Cancer Cell Lines

## Plain-language summary

These figures tell a more specific story than “AP-1 goes up under stress.” Nutrient stress appears to push pancreatic cancer cells into different regulatory states depending on where the cells start. HPAFII behaves like the most epithelial-like line, Panc1 behaves like the most mesenchymal-like and stress-adapted line, and AsPC1 sits between them.

The clearest model is three adaptation modes. First, Panc1 shows a **mesenchymal lock-in** pattern: AP-1, ATF/ISR, NRF-small-MAF/redox, and TEAD programs move together at the footprint level across BCAA, tryptophan, lysine, glucose, and methionine/cysteine stress. Second, AsPC1 shows a **BCAA compensation** pattern: KLF-family RNA rises while KLF-family motif footprints fall, suggesting compensation rather than straightforward activation. Third, HPAFII shows **epithelial buffering**: HNF4A, HNF4G, and CDX2 weaken under arginine or glutamine/arginine stress, while OVOL2 rises, consistent with an attempted epithelial restraint.

The footprint values below are motif-associated differential footprint signals from fp-tools. They should be read as evidence that a motif-centered regulatory program changes under nutrient stress, not as proof that one specific TF protein binds every site. The RNA values are log2 fold changes from the RUVr k=20 corrected DESeq2-normalized RNA matrix.

![Integrated TF evidence map](fig1_integrated_tf_evidence_map.png)

**Integrated TF evidence map.** Circles show footprint changes and squares show RNA changes. Red means higher in nutrient stress than in `10_FBS_Ctrl`; blue means lower. Larger circles have stronger footprint evidence. AP-1 is the broadest repeated signal, but the map also shows lineage-state and redox/plasticity branches.

![Program heatmap](fig2_program_level_stress_heatmap.png)

**Program-level stress response heatmap.** Each tile summarizes the strongest TF in a broader program: AP-1 marks the main stress response, HNF4 marks epithelial/classical identity, ATF marks nutrient-stress signaling, NRF2 marks redox protection, and TEAD/FOXC/ZEB mark plasticity.

![Recurrence matrix](fig3_cross_cellline_recurrence_matrix.png)

**Cross-cell-line recurrence.** AP-1 is the clearest repeated program, but ATF/ISR, TEAD, NRF2/redox, and plasticity programs also recur in selected stress types.

![Lineage model](fig4_lineage_state_model.png)

**Lineage-state model.** The three cell lines appear to respond from different starting states. HPAFII keeps the strongest epithelial/classical context, AsPC1 looks intermediate, and Panc1 shows the strongest AP-1, TEAD, NRF2, and plasticity-associated changes.

![Three adaptation modes](fig5_three_adaptation_modes.png)

**Three nutrient-stress adaptation modes.** HPAFII shows HNF4A/HNF4G/CDX2 loss together with OVOL2 gain under Arg/Gln.Arg stress. AsPC1 shows BCAA-sensitive KLF compensation, with RNA and footprint moving in opposite directions for several KLF-family factors. Panc1 shows a coordinated AP-1, ATF/ISR, NRF-small-MAF/redox, and TEAD footprint program across multiple severe nutrient stresses.

## First-wave target priorities

The strongest EMT-blocking target logic is not simply to knock out famous EMT genes. The data support two groups. The first group contains anchor targets with strong literature support or clear plasticity biology, such as TEAD4/TEAD3, FOSL1/AP-1, FOXC1, and ATF4/ISR. The second group contains higher-novelty redox or proteostasis candidates, especially MAFG, NFE2L1, and BACH2.

**Table 1. First-wave EMT-blocking target priorities.**

{first_wave_md}

**Table 2. Target handling logic.** Some factors are better treated as knockout or CRISPRi candidates, while others are better treated as rescue or lineage-state markers.

{handling_md}

**Table 3. Higher-novelty candidates.** These factors are less obvious than AP-1 or TEAD but may explain why nutrient stress reinforces survival and plasticity in stressed pancreatic cancer cells.

{novelty_md}

![Expanded TF evidence map](fig6_expanded_tf_evidence_map.png)

**Expanded TF evidence map.** This figure adds more TFs and motif-family members so that the AP-1/ATF, NRF/MAF, KLF, epithelial-lineage, metabolic nuclear-receptor, and immediate-early branches can be inspected together.

![Expanded program heatmap](fig7_expanded_program_heatmap.png)

**Expanded program heatmap.** This view groups the additional TFs into broader modules, including NRF/MAF redox partners, KLF-family metabolic regulators, epithelial-lineage regulators, and nuclear-receptor-like metabolic programs.

## How the story fits existing biology

The AP-1 result is the strongest repeated anchor because it appears across cell lines and has the largest footprint changes, especially in Panc1. This fits pancreatic cancer literature linking AP-1 and FOSL1 to proliferation, invasion, metastasis, drug resistance, and FOSL1/RELA-dependent transcriptional mechanisms [1,2]. RNA and footprint changes do not always move in the same direction, which is expected for AP-1 because activity depends on dimer partners, upstream signaling, and protein-level regulation.

ATF3 and ATF4 fit the same nutrient-stress story from a survival angle. They point to the integrated stress response, a pathway that helps cells survive amino-acid limitation, ER stress, and other stress conditions. In pancreatic cancer models, ISR inhibition or ATF4 depletion can make cells more sensitive to apoptosis [3,4].

HNF4A and HNF4G are the clearest epithelial/classical lineage-state signals. They are most visible in HPAFII and AsPC1, which fits studies linking HNF factors and GATA6-centered networks to classical PDAC identity, lineage plasticity, and metabolic state [5,6]. If HNF4A, HNF4G, and CDX2 decrease under stress, that likely means epithelial/classical identity is being weakened or remodeled. OVOL2 is different: in HPAFII it may represent an epithelial brake that rises while other epithelial factors decline.

TEAD3 and TEAD4 add a plasticity layer. TEAD inhibitors such as GNE-7883 show that YAP/TAZ-TEAD transcriptional output can be targeted pharmacologically and can reduce accessibility at TEAD motifs in YAP/TAZ-dependent models [7]. In these data, TEAD is most relevant for Panc1 and selected stress contexts.

The redox branch should be more central than in the earlier draft. NFE2L2/NRF2 is a known oxidative-stress regulator linked to therapy resistance and ferroptosis protection in cancer [8], but the deeper tables also point to MAFG, NFE2L1, MAFK, MAF, and BACH2. This makes the NRF-small-MAF/proteostasis branch one of the strongest high-novelty patterns, especially in Panc1.

FOXC1 and ZEB1 are best used as plasticity markers and secondary perturbation candidates. ZEB1 is a known EMT and plasticity regulator in pancreatic cancer metastasis [9,10]. In this analysis, these factors help connect nutrient stress to lineage plasticity, rather than serving as the dominant pan-stress signal.

The KLF-family result is most interesting in AsPC1 BCAA withdrawal. KLF10, KLF12, and KLF14 show RNA increases while motif-associated footprint signals fall. That pattern may mean the cells are trying to compensate transcriptionally after losing a prior KLF-linked chromatin state. This is why AsPC1 should be described as an intermediate compensation state, not only as an intermediate lineage state.

## References

"""
    for i, ref in enumerate(REFERENCES, start=1):
        text += f"{i}. {ref[0]}. {ref[1]} {ref[2]} {ref[3]}\n"
    md.write_text(text, encoding="utf-8")
    return md


def convert_docx(md: Path, docx: Path) -> None:
    subprocess.run(["pandoc", md.name, "-o", docx.name, "--standalone"], check=True, cwd=md.parent)


def zip_outputs(outdir: Path) -> Path:
    stamp = date.today().strftime("%Y%m%d")
    zip_path = outdir / f"nutrient_tfs_{stamp}.zip"
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for item in sorted(outdir.rglob("*")):
            if item == zip_path or item.is_dir() or item.suffix == ".zip":
                continue
            if item.name.startswith("nutrient_tf_story_writeup"):
                continue
            if item.name.startswith("fp_tools_nutrient_tf_story_figures"):
                continue
            zf.write(item, item.relative_to(outdir))
    return zip_path


def main() -> None:
    args = parse_args()
    handoff = Path(args.handoff_dir)
    outdir = Path(args.outdir) if args.outdir else handoff / "story_figures"
    source_dir = outdir / "story_figure_source_tables"
    outdir.mkdir(parents=True, exist_ok=True)
    source_dir.mkdir(parents=True, exist_ok=True)

    mpl.rcParams.update(
        {
            "font.family": PLOT_FONT,
            "font.sans-serif": ["Nimbus Sans", "Arial", "Helvetica", "DejaVu Sans"],
            "font.weight": "bold",
            "axes.labelweight": "bold",
            "axes.titleweight": "bold",
            "font.size": 10,
            "axes.labelcolor": "black",
            "axes.edgecolor": "black",
            "text.color": "black",
            "xtick.color": "black",
            "ytick.color": "black",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    sns.set_theme(style="whitegrid", font=PLOT_FONT)

    summary, recurrence, _ = load_data(handoff)
    dot_source = plot_integrated_dotplot(summary, outdir)
    program_df = build_program_table(summary)
    program_source = plot_program_heatmap(program_df, outdir)
    recurrence_source = plot_recurrence(program_df, outdir)
    lineage_source = plot_lineage_model(program_df, outdir)
    mode_source = plot_three_adaptation_modes(summary, outdir)
    target_source = plot_target_table(summary, outdir)
    expanded_tf_source = plot_expanded_tf_evidence_map(summary, outdir)
    expanded_program_df = build_program_table(summary, EXPANDED_PROGRAMS)
    expanded_program_source = plot_expanded_program_heatmap(expanded_program_df, outdir)
    additional_target_source = plot_additional_candidate_table(summary, outdir)
    first_wave_source = pd.DataFrame(FIRST_WAVE_TARGETS)
    target_handling_source = pd.DataFrame(TARGET_HANDLING)
    novelty_source = pd.DataFrame(NOVELTY_TARGETS)

    dot_source.to_csv(source_dir / "fig1_integrated_tf_evidence_map_source.tsv", sep="\t", index=False)
    program_source.to_csv(source_dir / "fig2_program_level_stress_heatmap_source.tsv", sep="\t", index=False)
    recurrence_source.to_csv(source_dir / "fig3_cross_cellline_recurrence_source.tsv", sep="\t", index=False)
    lineage_source.to_csv(source_dir / "fig4_lineage_state_model_source.tsv", sep="\t", index=False)
    mode_source.to_csv(source_dir / "fig5_three_adaptation_modes_source.tsv", sep="\t", index=False)
    target_source.to_csv(source_dir / "fig5_prioritized_tf_target_table_source.tsv", sep="\t", index=False)
    expanded_tf_source.to_csv(source_dir / "fig6_expanded_tf_evidence_map_source.tsv", sep="\t", index=False)
    expanded_program_source.to_csv(source_dir / "fig7_expanded_program_heatmap_source.tsv", sep="\t", index=False)
    additional_target_source.to_csv(source_dir / "fig8_additional_candidate_tf_table_source.tsv", sep="\t", index=False)
    first_wave_source.to_csv(source_dir / "table1_first_wave_targets_source.tsv", sep="\t", index=False)
    target_handling_source.to_csv(source_dir / "table2_target_handling_source.tsv", sep="\t", index=False)
    novelty_source.to_csv(source_dir / "table3_novelty_targets_source.tsv", sep="\t", index=False)
    recurrence.to_csv(source_dir / "full_cross_cellline_recurrence_input.tsv", sep="\t", index=False)

    md = write_markdown(
        outdir,
        {
            "first_wave_table": first_wave_source,
            "target_handling_table": target_handling_source,
            "novelty_table": novelty_source,
        },
    )
    stamp = date.today().strftime("%Y%m%d")
    docx = outdir / f"nutrient_tfs_{stamp}.docx"
    convert_docx(md, docx)
    zip_path = zip_outputs(outdir)
    print(f"Story figure directory: {outdir}")
    print(f"Write-up Markdown: {md}")
    print(f"Write-up DOCX: {docx}")
    print(f"Zip: {zip_path}")


if __name__ == "__main__":
    main()
