"""Streamlit GUI for fp-tools.

This module is an isolated wrapper around the packaged commands. Direct CLI
usage remains primary. The GUI supports direct form-driven runs, YAML load/save,
and batch editing while using the same normalized config model as the optional
``run-yaml-workflow --config ...`` path.
"""

from __future__ import annotations

import json
import os
import re
import shlex
from html import escape
from importlib import resources
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

from fp_tools import __version__
from fp_tools.gui_config import (
    GUI_ENUM_CHOICES,
    canonical_tool_name,
    config_to_yaml_text,
    load_yaml_config,
    make_single_config,
    normalize_config,
    parse_yaml_text,
    validate_gui_config,
)
from fp_tools.gui_jobs import default_gui_run_dir, launch_config_async, materialize_run_config, refresh_run_status

GUI_EXAMPLE_PACKAGE = "fp_tools.resources.gui_configs"

PAGE_OPTIONS = [
    "Home",
    "Run History",
    "atac-correct",
    "call-footprints",
    "match-motifs",
    "diff-footprints",
    "normalize-bigwig",
    "plot-aggregate",
    "bulk-footprinting",
    "review-multi-comparisons",
    "discover-motifs",
    "summarize-motifs",
    "pseudobulk-fragments",
    "find-signature-fp",
    "sc-footprinting",
    "Config",
]

NAV_GROUPS = [
    ("Overview", ["Home", "Run History"]),
    ("Core workflow", ["atac-correct", "call-footprints", "match-motifs", "diff-footprints"]),
    ("Workflow and interface", ["bulk-footprinting", "sc-footprinting", "review-multi-comparisons"]),
    ("Signals and reports", ["normalize-bigwig", "plot-aggregate"]),
    ("De Novo Motif Discovery", ["discover-motifs", "summarize-motifs"]),
    ("Single-cell ATAC-seq", ["pseudobulk-fragments", "find-signature-fp"]),
    ("Configuration", ["Config"]),
]

GENERIC_TOOL_DEFAULTS: dict[str, dict[str, Any]] = {
    "bulk-footprinting": {
        "sample_id": "bulk_footprinting_run",
        "sample_table": "",
        "comparison_table": "",
        "genome": "",
        "blacklist": "",
        "outdir": "",
        "motif_db": "jaspar2026_vertebrates",
        "plot_aggregate": "all",
        "review_format": "auto",
        "cores": 4,
    },
    "review-multi-comparisons": {
        "sample_id": "comparison_browser_run",
        "inputs": [],
        "labels": [],
        "output_dir": "",
        "output_html": "",
        "layout": "custom",
        "title": "Review multiple differential footprint comparisons",
    },
    "match-motifs": {
        "sample_id": "match_motifs_run",
        "signals": "",
        "genome": "",
        "peaks": "",
        "peak_header": "",
        "outdir": "",
        "cond_names": "Bcell",
        "motif_db": "jaspar2026_vertebrates",
        "skip_excel": True,
    },
    "normalize-bigwig": {
        "sample_id": "normalize_bigwig_run",
        "bigwigs": "",
        "background": "",
        "outdir": "",
        "method": "background-scale",
        "stat": "q95",
        "target": "median",
        "chrom_sizes": "",
        "workers": 2,
    },
    "discover-motifs": {
        "sample_id": "motif_discovery_run",
        "candidates": "",
        "fasta": "",
        "genome": "",
        "outdir": "",
        "method": "streme",
        "known_motif_db": "jaspar2026_vertebrates",
        "known_motifs": "",
        "script": "",
        "execute": False,
        "runtime": "auto",
    },
    "summarize-motifs": {
        "sample_id": "motif_summary_run",
        "meme_txt": "",
        "tomtom_tsv": "",
        "out_tsv": "",
        "out_html": "",
        "title": "fp-tools motif summary",
    },
    "pseudobulk-fragments": {
        "sample_id": "pseudobulk_fragments_run",
        "fragments": "",
        "annotations": "",
        "group_by": "cell_type",
        "outdir": "",
        "min_cells": 1,
        "min_fragments": 1,
        "index_output": True,
    },
    "find-signature-fp": {
        "sample_id": "signature_fp_run",
        "fragments": "",
        "annotations": "",
        "h5ad": "",
        "tf_site_dir": "",
        "outdir": "",
        "markers": "STAT6,FOSB,CEBPA,IRF8,RELA,ZNF683,NR4A1,SMAD3",
        "summary_output_prefix": "single_cell_footprinting",
        "max_motifs": 25,
    },
    "sc-footprinting": {
        "sample_id": "pseudobulk_footprints_run",
        "fragments": "",
        "annotations": "",
        "h5ad": "",
        "group_by": "cell_type",
        "outdir": "",
        "genome_sizes": "",
        "genome": "",
        "peaks": "",
        "motif_db": "jaspar2026_vertebrates",
        "dry_run": False,
    },
}

LIST_TEXT_FIELDS = {
    "signals",
    "bigwigs",
    "motifs",
    "inputs",
    "labels",
    "input_html",
    "cond_names",
    "sample_names",
    "sample_dirs",
    "match_dir",
    "aggregate_signals",
    "region_labels",
    "known_motifs",
    "TFBS",
    "read_shift",
    "markers",
}

GUI_FIELD_LABELS = {
    "TFBS": "Region BED files",
    "annotations": "Cell annotations",
    "background": "Background regions BED",
    "bigwigs": "Signal bigWig files",
    "candidates": "Candidate regions BED",
    "chrom_sizes": "Chromosome sizes",
    "comparison_table": "Comparisons TSV (optional)",
    "cond_names": "Condition names",
    "fasta": "Candidate sequences FASTA (optional)",
    "fragments": "Fragments file",
    "genome": "Genome FASTA",
    "genome_sizes": "Chromosome sizes",
    "group_by": "Annotation column",
    "h5ad": "Annotated h5ad file",
    "input_html": "Report HTML files",
    "inputs": "Comparison result folders",
    "known_motif_db": "Known motif database",
    "known_motifs": "Known motif files (optional)",
    "markers": "Marker motifs (one per line)",
    "meme_txt": "MEME motif file",
    "motif_db": "Motif database",
    "dry_run": "Validate configuration only",
    "out_html": "Output HTML",
    "out_tsv": "Output TSV",
    "outdir": "Output directory",
    "output_dir": "Output directory",
    "output_html": "Output HTML (optional)",
    "peak_header": "Peak annotation header (optional)",
    "peaks": "Accessible regions BED",
    "sample_table": "Samples TSV",
    "signals": "Footprint bigWig files",
    "tf_site_dir": "Motif-site directory",
    "tomtom_tsv": "TOMTOM matches TSV",
}

GUI_TOOL_DESCRIPTIONS = {
    "bulk-footprinting": "Run the complete bulk ATAC-seq footprinting workflow from BAM and peak files.",
    "match-motifs": "Scan accessible regions and summarize motif-associated footprint scores.",
    "review-multi-comparisons": "Combine completed differential-footprint results into one interactive report.",
    "normalize-bigwig": "Normalize multiple bigWig tracks over shared background regions.",
    "discover-motifs": "Discover motifs from candidate footprint regions.",
    "summarize-motifs": "Summarize discovered motifs and known-motif matches.",
    "pseudobulk-fragments": "Group single-cell ATAC fragments using cell annotations.",
    "find-signature-fp": "Identify and visualize cell-type footprint signatures.",
    "sc-footprinting": "Run the complete single-cell ATAC-seq footprinting workflow.",
}

GUI_ADVANCED_FIELDS = {
    "bulk-footprinting": {"blacklist", "plot_aggregate", "review_format", "cores"},
    "review-multi-comparisons": {"layout", "title", "output_html"},
    "match-motifs": {"peak_header", "skip_excel"},
    "normalize-bigwig": {"chrom_sizes", "workers", "stat", "target"},
    "discover-motifs": {"known_motifs", "script", "execute", "runtime"},
    "pseudobulk-fragments": {"min_cells", "min_fragments", "index_output"},
    "find-signature-fp": {"max_motifs", "summary_output_prefix"},
    "sc-footprinting": {"dry_run"},
}

GUI_FIELD_HELP: dict[str, dict[str, str]] = {
    "bulk-footprinting": {
        "normalization": "Normalization applied during differential footprinting.",
        "plot_aggregate": "Motifs included in aggregate profiles, or off to omit profiles.",
        "review_format": "Bundle, standalone HTML, automatic selection, or no combined report.",
    },
    "normalize-bigwig": {
        "method": "Background scaling, robust z-score transformation, or no transformation.",
        "stat": "Use median, iqr, or a quantile such as q90 or q95.",
        "target": "Across-sample statistic used as the scaling target.",
    },
    "discover-motifs": {
        "method": "MEME Suite discovery program.",
        "runtime": "Managed, system, or container source for the optional MEME tools.",
    },
}


def main() -> None:
    st.set_page_config(page_title="fp-tools GUI", layout="wide")
    _apply_page_style()
    _ensure_session_config()

    run_dir = Path(st.session_state.gui_run_dir).expanduser()

    page = _current_page_from_query()
    _sync_config_for_page(page)
    _render_sidebar_header()
    _render_sidebar_run_dir_controls()
    if not _tutorial_visible() and st.sidebar.button("Show guided tutorial", key="show_tutorial_button", width="stretch"):
        st.session_state.show_tutorial = True
        st.session_state.hide_tutorial = False
    _render_sidebar_nav(page)
    _render_config_update_notice()

    if _tutorial_visible():
        _render_tutorial_panel()

    if page == "Home":
        _render_home(run_dir)
    elif page == "Run History":
        _render_run_history(run_dir)
    elif page == "atac-correct":
        _render_atacorrect_page(run_dir)
    elif page == "call-footprints":
        _render_footprintscores_page(run_dir)
    elif page == "diff-footprints":
        _render_diff_footprints_page(run_dir)
    elif page == "plot-aggregate":
        _render_plotaggregate_page(run_dir)
    elif page == "Config":
        _render_config_page(run_dir)
    else:
        _render_generic_tool_page(run_dir, page)

def _apply_page_style() -> None:
    st.markdown(
        """
        <style>
        :root {
            color-scheme: light;
            --fp-bg: #f3f6fa;
            --fp-surface: #ffffff;
            --fp-surface-soft: #f9fafb;
            --fp-border: #cfd8e3;
            --fp-border-soft: #dbe3ec;
            --fp-text: #111827;
            --fp-text-muted: #4b5563;
            --fp-accent: #173b73;
            --fp-accent-ink: #102a5c;
            --fp-accent-soft: #eaf2ff;
            --fp-accent-hover: #102f61;
            --fp-hover: #1e293b;
            --fp-radius-control: 8px;
            --fp-radius-card: 8px;
            --fp-shadow: 0 8px 24px rgba(15, 23, 42, 0.07);
            --fp-font-body: 1rem;
            --fp-font-label: 0.92rem;
            --fp-control-height: 2.7rem;
        }
        html, body, [class*="css"], [data-testid="stAppViewContainer"] {
            font-family: Arial, -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
        }
        [data-testid="stAppViewContainer"] {
            background: var(--fp-bg);
            color: var(--fp-text);
        }
        #MainMenu {visibility: hidden;}
        footer {visibility: hidden;}
        header[data-testid="stHeader"] {
            height: 0;
            background: transparent;
            border: 0;
        }
        [data-testid="stToolbar"],
        [data-testid="stDecoration"] {
            display: none !important;
        }
        [data-testid="stSidebar"] {
            background: #0f172a;
            border-right: 1px solid #1e293b;
            min-width: 314px !important;
            max-width: 314px !important;
        }
        [data-testid="stSidebarHeader"],
        [data-testid="stSidebarCollapseButton"] {
            display: none !important;
        }
        [data-testid="stSidebar"] [data-testid="stSidebarContent"] {
            padding-top: 0.85rem;
            padding-left: 0.92rem;
            padding-right: 0.92rem;
        }
        [data-testid="stSidebar"] [data-testid="stVerticalBlock"] {
            gap: 0.3rem;
        }
        [data-testid="stSidebar"] .stMarkdown,
        [data-testid="stSidebar"] label,
        [data-testid="stSidebar"] p,
        [data-testid="stSidebar"] span,
        [data-testid="stSidebar"] .stCaption,
        [data-testid="stSidebar"] .stRadio label p {
            color: #e5edf6 !important;
        }
        [data-testid="stSidebar"] .stCaption {
            color: #aebbd0 !important;
        }
        [data-testid="stSidebar"] input,
        [data-testid="stSidebar"] textarea,
        [data-testid="stSidebar"] [data-baseweb="input"] input {
            color: var(--fp-text) !important;
            -webkit-text-fill-color: var(--fp-text) !important;
            background: var(--fp-surface) !important;
        }
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button {
            background: transparent !important;
            color: #e5edf6 !important;
            border: 1px solid transparent !important;
            font-weight: 800 !important;
            border-radius: var(--fp-radius-control) !important;
            min-height: 2.08rem !important;
            box-shadow: none !important;
            justify-content: flex-start !important;
            text-align: left !important;
            padding: 0.34rem 0.66rem !important;
        }
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button:hover {
            background: var(--fp-hover) !important;
            color: #ffffff !important;
            border-color: #334155 !important;
        }
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button:disabled,
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button:disabled:hover {
            background: #1d4ed8 !important;
            border-color: #2563eb !important;
            color: #ffffff !important;
            opacity: 1 !important;
            cursor: default !important;
            box-shadow: inset 4px 0 0 #72e0b2 !important;
        }
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button p,
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button span,
        [data-testid="stSidebar"] [data-baseweb="input"] > div,
        [data-testid="stSidebar"] [data-testid="stTextInputRootElement"] > div {
            color: #e5edf6 !important;
            background: transparent !important;
        }
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button *,
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button [data-testid="stMarkdownContainer"],
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button [data-testid="stMarkdownContainer"] * {
            color: #e5edf6 !important;
            opacity: 1 !important;
        }
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button:disabled p,
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button:disabled span,
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button:disabled * {
            color: #ffffff !important;
        }
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button p {
            width: 100% !important;
            text-align: left !important;
            font-size: 0.98rem !important;
            line-height: 1.18 !important;
        }
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button > div,
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button > div > span,
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button [data-testid="stMarkdownContainer"],
        [data-testid="stSidebar"] [class*="st-key-nav_"] .stButton > button [data-testid="stMarkdownContainer"] p {
            width: 100% !important;
            justify-content: flex-start !important;
            text-align: left !important;
        }
        [data-testid="stSidebar"] [data-baseweb="input"] {
            background: var(--fp-surface) !important;
        }
        [data-testid="stSidebar"] [data-testid="stExpander"] {
            background: #101d33 !important;
            border: 1px solid #2b3a55 !important;
            box-shadow: none !important;
        }
        [data-testid="stSidebar"] [data-testid="stExpander"] summary,
        [data-testid="stSidebar"] [data-testid="stExpander"] summary *,
        [data-testid="stSidebar"] [data-testid="stExpander"] [data-testid="stMarkdownContainer"] *,
        [data-testid="stSidebar"] [data-testid="stExpander"] label {
            color: #e5edf6 !important;
        }
        [data-testid="stTextInputRootElement"] > div,
        [data-baseweb="base-input"],
        textarea,
        [data-baseweb="select"] > div {
            background: var(--fp-surface) !important;
            border-color: #aebccc !important;
            border-radius: var(--fp-radius-control) !important;
            box-shadow: none !important;
        }
        [data-testid="stTextInputRootElement"],
        [data-testid="stTextAreaRootElement"],
        [data-testid="stNumberInputContainer"],
        [data-baseweb="select"] > div {
            background: var(--fp-surface) !important;
            border: 1px solid #aebccc !important;
            border-radius: var(--fp-radius-control) !important;
            box-shadow: none !important;
        }
        [data-testid="stTextInputRootElement"]:focus-within,
        [data-testid="stTextAreaRootElement"]:focus-within,
        [data-testid="stNumberInputContainer"]:focus-within,
        [data-baseweb="base-input"]:focus-within,
        [data-baseweb="select"] > div:focus-within,
        textarea:focus {
            border-color: #2563eb !important;
            box-shadow: 0 0 0 2px rgba(37, 99, 235, 0.14) !important;
        }
        [data-testid="stTextInputRootElement"] input,
        textarea {
            font-size: 0.95rem !important;
        }
        [data-testid="stForm"],
        [data-testid="stExpander"],
        [data-testid="stAlert"],
        [data-testid="stCodeBlock"],
        [data-testid="stDataFrame"],
        [data-testid="stTable"] {
            border-radius: var(--fp-radius-card) !important;
        }
        [data-testid="stForm"],
        [data-testid="stExpander"],
        [data-testid="stAlert"] {
            background: var(--fp-surface) !important;
            border: 1px solid var(--fp-border-soft) !important;
        }
        [data-testid="stForm"] {
            padding: 1rem 1.05rem !important;
            box-shadow: var(--fp-shadow) !important;
        }
        [data-testid="stExpander"],
        [data-testid="stAlert"] {
            padding: 0.32rem 0.5rem !important;
            box-shadow: none !important;
        }
        [data-testid="stMain"] [data-testid="stForm"] [data-testid="stVerticalBlock"] {
            gap: 0.7rem !important;
        }
        [data-testid="stMain"] [data-testid="stExpander"] summary {
            min-height: var(--fp-control-height) !important;
            padding: 0 0.45rem !important;
        }
        [data-testid="stMain"] [data-testid="stExpander"] summary p {
            color: var(--fp-text) !important;
            font-size: var(--fp-font-label) !important;
            font-weight: 700 !important;
        }
        [data-testid="stFileUploaderDropzone"] {
            background: var(--fp-surface-soft) !important;
            border: 1px dashed #aebccc !important;
            border-radius: var(--fp-radius-control) !important;
        }
        [data-testid="stMarkdownContainer"] h1,
        [data-testid="stMarkdownContainer"] h2,
        [data-testid="stMarkdownContainer"] h3,
        .stTitle,
        .stHeader,
        .stSubheader {
            letter-spacing: -0.02em;
            color: var(--fp-text);
        }
        [data-testid="stMarkdownContainer"] h1 {
            font-size: 1.8rem;
            letter-spacing: 0;
            margin-bottom: 0.2rem;
        }
        [data-testid="stMarkdownContainer"] h2 {
            font-size: 1.25rem;
            letter-spacing: 0;
        }
        [data-testid="stMarkdownContainer"] h3 {
            font-size: 1rem;
            letter-spacing: 0;
        }
        .stTitle {
            font-weight: 760 !important;
        }
        .stButton > button,
        .stDownloadButton > button,
        .stFormSubmitButton > button {
            background: var(--fp-accent) !important;
            color: #ffffff !important;
            border: 1px solid var(--fp-accent) !important;
            border-radius: var(--fp-radius-control) !important;
            min-height: 2.65rem !important;
            padding: 0.35rem 1rem !important;
            font-weight: 700 !important;
            box-shadow: var(--fp-shadow) !important;
        }
        .stButton > button:hover,
        .stDownloadButton > button:hover,
        .stFormSubmitButton > button:hover {
            background: var(--fp-accent-hover) !important;
            border-color: var(--fp-accent-hover) !important;
        }
        .stButton > button:disabled,
        .stButton > button:disabled:hover,
        .stDownloadButton > button:disabled,
        .stDownloadButton > button:disabled:hover,
        .stFormSubmitButton > button:disabled,
        .stFormSubmitButton > button:disabled:hover {
            background: #e5e7eb !important;
            border-color: #cbd5e1 !important;
            color: #64748b !important;
            box-shadow: none !important;
            cursor: not-allowed !important;
            opacity: 1 !important;
        }
        [data-testid="stMain"] [data-baseweb="radio"] [aria-checked="true"] > div,
        [data-testid="stMain"] [data-baseweb="checkbox"] [aria-checked="true"] > div {
            background-color: #2563eb !important;
            border-color: #2563eb !important;
        }
        [data-testid="stCodeBlock"] pre,
        code {
            border-radius: var(--fp-radius-control) !important;
        }
        [data-testid="stDataFrame"] {
            background: var(--fp-surface) !important;
            border: 1px solid var(--fp-border-soft) !important;
            box-shadow: var(--fp-shadow) !important;
        }
        .block-container {
            width: 100% !important;
            max-width: 1840px !important;
            padding-top: 0.75rem;
            padding-bottom: 1.15rem;
            padding-left: clamp(1rem, 1.7vw, 2.2rem);
            padding-right: clamp(1rem, 1.7vw, 2.2rem);
        }
        [data-testid="stMain"],
        [data-testid="stMain"] [data-testid="stMainBlockContainer"],
        [data-testid="stMain"] [data-testid="stHorizontalBlock"],
        [data-testid="stMain"] [data-testid="stColumn"],
        [data-testid="stMain"] [data-testid="stColumn"] > [data-testid="stVerticalBlock"],
        [data-testid="stMain"] [data-testid="stElementContainer"] {
            min-width: 0 !important;
            max-width: 100% !important;
        }
        [data-testid="stMain"] [data-testid="stCode"],
        [data-testid="stMain"] [data-testid="stCodeBlock"] {
            min-width: 0 !important;
            max-width: 100% !important;
            overflow: hidden !important;
        }
        [data-testid="stMain"] [data-testid="stCode"] pre,
        [data-testid="stMain"] [data-testid="stCodeBlock"] pre {
            min-width: 0 !important;
            max-width: 100% !important;
            overflow-x: auto !important;
        }
        [data-testid="stMain"] label,
        [data-testid="stMain"] label *,
        [data-testid="stMain"] [data-testid="stWidgetLabel"],
        [data-testid="stMain"] [data-testid="stWidgetLabel"] *,
        [data-testid="stMain"] [role="radiogroup"] label,
        [data-testid="stMain"] [role="radiogroup"] label *,
        [data-testid="stMain"] [data-baseweb="checkbox"] *,
        [data-testid="stMain"] [data-baseweb="radio"] *,
        [data-testid="stMain"] [data-baseweb="select"] *,
        [data-testid="stMain"] [data-testid="stFileUploader"] * {
            color: var(--fp-text) !important;
            opacity: 1 !important;
            -webkit-text-fill-color: var(--fp-text) !important;
            visibility: visible !important;
        }
        [data-testid="stMain"] [data-testid="stWidgetLabel"] p,
        [data-testid="stMain"] [role="radiogroup"] label p,
        [data-testid="stMain"] [data-baseweb="checkbox"] label p,
        [data-testid="stMain"] [data-baseweb="radio"] label p {
            font-size: var(--fp-font-label) !important;
            font-weight: 700 !important;
            line-height: 1.35 !important;
        }
        [data-testid="stMain"] [data-testid="stTextInputRootElement"],
        [data-testid="stMain"] [data-testid="stNumberInputContainer"],
        [data-testid="stMain"] [data-baseweb="select"] > div {
            min-height: var(--fp-control-height) !important;
        }
        [data-testid="stMain"] input,
        [data-testid="stMain"] textarea,
        [data-testid="stMain"] [data-baseweb="select"] * {
            font-size: 0.95rem !important;
        }
        [data-testid="stMain"] p,
        [data-testid="stMain"] li {
            line-height: 1.45;
        }
        [data-testid="stMain"] input,
        [data-testid="stMain"] textarea,
        [data-testid="stMain"] [data-baseweb="input"] input {
            color: var(--fp-text) !important;
            opacity: 1 !important;
            -webkit-text-fill-color: var(--fp-text) !important;
        }
        [data-testid="stMain"] input::placeholder,
        [data-testid="stMain"] textarea::placeholder {
            color: #64748b !important;
            opacity: 1 !important;
            -webkit-text-fill-color: #64748b !important;
        }
        .fp-sidebar-brand {
            color: #ffffff;
            padding: 0.18rem 0.08rem 0.28rem;
            margin: 0;
        }
        .fp-sidebar-brand-title {
            font-size: 1.2rem;
            font-weight: 800;
            letter-spacing: 0;
            line-height: 1.1;
        }
        .fp-nav-group {
            margin: 0;
            display: block;
            box-sizing: border-box;
            padding: 0.72rem 0 0.46rem;
            color: #93a4b8;
            font-size: 0.72rem;
            font-weight: 800;
            letter-spacing: 0.06em;
            text-transform: uppercase;
            line-height: 1.1;
        }
        [data-testid="stMarkdown"]:has(.fp-nav-group)
        [data-testid="stMarkdownContainer"] {
            margin-bottom: 0 !important;
        }
        .fp-workspace-meta {
            color: #aebbd0;
            font-size: 0.76rem;
            line-height: 1.35;
            margin-bottom: 0.28rem;
        }
        .fp-workspace-path {
            background: #101d33;
            border: 1px solid #334155;
            border-radius: 6px;
            color: #e5edf6;
            font-family: Arial, -apple-system, BlinkMacSystemFont, "SF Pro Text", "Helvetica Neue", sans-serif;
            font-size: 0.76rem;
            line-height: 1.35;
            margin: 0.18rem 0 0.48rem;
            max-width: 100%;
            overflow-wrap: anywhere;
            padding: 0.42rem 0.5rem;
            user-select: text;
            white-space: normal;
            word-break: break-word;
        }
        [data-testid="stSidebar"] .st-key-show_tutorial_button .stButton > button {
            font-size: 0.78rem !important;
            font-weight: 700 !important;
            min-height: 1.95rem !important;
            padding: 0.2rem 0.55rem !important;
        }
        .fp-hero {
            background: #ffffff;
            border: 1px solid var(--fp-border-soft);
            border-radius: 8px;
            padding: 0.95rem 1.05rem;
            box-shadow: var(--fp-shadow);
            margin-bottom: 0.68rem;
        }
        .fp-kicker {
            color: #173b73;
            font-size: 0.72rem;
            font-weight: 800;
            letter-spacing: 0.08em;
            text-transform: uppercase;
            margin-bottom: 0.35rem;
        }
        .fp-hero h1 {
            font-size: clamp(1.58rem, 1.45vw, 2.05rem);
            margin: 0 0 0.25rem;
            letter-spacing: 0;
            max-width: 1080px;
        }
        .fp-hero p {
            color: var(--fp-text-muted);
            font-size: var(--fp-font-body);
            line-height: 1.38;
            max-width: 980px;
            margin: 0;
        }
        .fp-pill-row {
            display: flex;
            flex-wrap: wrap;
            gap: 0.38rem;
            margin-top: 0.72rem;
        }
        .fp-pill {
            border: 1px solid #dbe3ef;
            background: #ffffff;
            border-radius: 999px;
            padding: 0.22rem 0.5rem;
            font-size: 0.76rem;
            font-weight: 800;
            color: #334155;
        }
        .fp-card-grid {
            display: grid;
            grid-template-columns: repeat(4, minmax(0, 1fr));
            gap: 0.55rem;
            margin: 0.55rem 0 0.68rem;
        }
        .fp-gui-card {
            background: var(--fp-surface);
            border: 1px solid var(--fp-border-soft);
            border-radius: 8px;
            padding: 0.66rem 0.72rem;
            box-shadow: var(--fp-shadow);
        }
        .fp-gui-card h3 {
            margin: 0 0 0.28rem;
            font-size: 0.98rem;
            letter-spacing: 0;
        }
        .fp-gui-card p {
            margin: 0;
            color: var(--fp-text-muted);
            font-size: 0.86rem;
            line-height: 1.3;
        }
        .fp-section-title {
            display: flex;
            align-items: center;
            justify-content: space-between;
            margin: 0.64rem 0 0.3rem;
            font-size: 0.98rem;
            font-weight: 800;
            color: var(--fp-text);
        }
        .fp-example-grid {
            display: grid;
            grid-template-columns: repeat(4, minmax(0, 1fr));
            gap: 0.42rem 0.55rem;
            margin-bottom: 0.55rem;
        }
        .fp-example-item {
            background: #ffffff;
            border: 1px solid var(--fp-border-soft);
            border-radius: 8px;
            padding: 0.42rem 0.55rem;
            color: #334155;
            font-family: Arial, sans-serif;
            font-size: 0.78rem;
            font-weight: 700;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }
        .fp-summary-strip {
            display: grid;
            grid-template-columns: repeat(3, minmax(0, 1fr));
            gap: 0.55rem;
            margin: 0.55rem 0 0.35rem;
        }
        .fp-summary-metric {
            background: #ffffff;
            border: 1px solid var(--fp-border-soft);
            border-radius: 8px;
            padding: 0.5rem 0.65rem;
        }
        .fp-summary-label {
            color: #64748b;
            font-size: 0.72rem;
            font-weight: 800;
            text-transform: uppercase;
            letter-spacing: 0.04em;
        }
        .fp-summary-value {
            color: #111827;
            font-size: 1.18rem;
            font-weight: 800;
            margin-top: 0.15rem;
        }
        .fp-tool-shell {
            display: grid;
            grid-template-columns: minmax(0, 1.35fr) minmax(360px, 0.78fr);
            gap: clamp(0.8rem, 1.2vw, 1.4rem);
            align-items: start;
        }
        .fp-page-heading {
            margin: 0 0 0.68rem;
        }
        .fp-page-heading h1 {
            font-size: clamp(1.58rem, 1.45vw, 2.05rem);
            margin: 0 0 0.25rem;
        }
        .fp-page-heading p {
            color: var(--fp-text-muted);
            font-size: var(--fp-font-body);
            line-height: 1.38;
            margin: 0;
        }
        .fp-run-card {
            padding: 0.08rem 0.08rem 0;
        }
        .fp-run-card h3 {
            font-size: 1rem;
            margin: 0 0 0.45rem;
        }
        .fp-run-summary {
            display: grid;
            grid-template-columns: repeat(3, minmax(0, 1fr));
            gap: 0.4rem;
            margin: 0.45rem 0;
        }
        .fp-run-summary div {
            border: 1px solid var(--fp-border-soft);
            border-radius: 8px;
            padding: 0.42rem 0.5rem;
            background: #f8fafc;
        }
        .fp-run-summary span {
            display: block;
            color: #64748b;
            font-size: 0.68rem;
            font-weight: 800;
            text-transform: uppercase;
        }
        .fp-run-summary strong {
            display: block;
            color: #111827;
            font-size: 0.92rem;
            margin-top: 0.08rem;
            overflow-wrap: normal;
            word-break: normal;
            hyphens: none;
        }
        [data-testid="stMain"] [class*="st-key-fp_run_panel_"] {
            background: var(--fp-surface) !important;
            border: 1px solid var(--fp-border-soft) !important;
            border-radius: var(--fp-radius-card) !important;
            box-shadow: var(--fp-shadow) !important;
        }
        .fp-validation-errors {
            box-sizing: border-box;
            width: 100%;
            min-width: 0;
            max-width: 100%;
            margin: 0.35rem 0 0.55rem;
            padding-left: 1.2rem;
            white-space: normal;
            overflow-x: hidden;
        }
        .fp-validation-errors li,
        .fp-validation-errors li * {
            box-sizing: border-box;
            min-width: 0;
            max-width: 100%;
            margin: 0.24rem 0;
            overflow-wrap: anywhere !important;
            word-break: normal !important;
            white-space: normal !important;
        }
        .fp-tutorial-panel {
            background: #eaf2ff;
            border: 1px solid #bdd7ff;
            border-radius: 8px;
            color: #173b73;
            padding: 0.68rem 0.78rem;
            margin: 0 0 0.68rem;
        }
        .fp-tutorial-panel h2 {
            margin: 0 0 0.22rem;
            font-size: 1rem;
        }
        .fp-tutorial-panel p,
        .fp-tutorial-panel li {
            font-size: 0.86rem;
            line-height: 1.28;
        }
        .fp-tutorial-panel ol {
            display: grid;
            grid-template-columns: repeat(4, minmax(0, 1fr));
            gap: 0.55rem;
            margin: 0.42rem 0 0;
            padding-left: 1.2rem;
        }
        @media (max-width: 1500px) {
            .fp-run-summary {
                grid-template-columns: repeat(2, minmax(0, 1fr));
            }
            .fp-run-summary div:first-child {
                grid-column: 1 / -1;
            }
            [data-testid="stMain"] [class*="st-key-fp_tool_shell_"] > div > [data-testid="stHorizontalBlock"] {
                align-items: stretch !important;
                flex-direction: column !important;
                gap: 0.72rem !important;
            }
            [data-testid="stMain"] [class*="st-key-fp_tool_shell_"] > div > [data-testid="stHorizontalBlock"] > [data-testid="stColumn"] {
                flex: 1 1 auto !important;
                min-width: 0 !important;
                width: 100% !important;
            }
        }
        @media (max-width: 900px) {
            header[data-testid="stHeader"] {
                height: 3.4rem !important;
                pointer-events: none !important;
            }
            [data-testid="stToolbar"] {
                align-items: center !important;
                display: flex !important;
                height: 2.75rem !important;
                left: 0.55rem !important;
                pointer-events: auto !important;
                position: fixed !important;
                right: auto !important;
                top: 0.45rem !important;
                visibility: visible !important;
                width: 2.75rem !important;
                z-index: 1000000 !important;
            }
            [data-testid="stToolbar"] > div,
            [data-testid="stToolbar"] > div > div,
            [data-testid="stToolbar"] > div > div > div:has([data-testid="stExpandSidebarButton"]) {
                align-items: center !important;
                display: flex !important;
                height: 2.75rem !important;
                margin: 0 !important;
                padding: 0 !important;
                width: 2.75rem !important;
            }
            [data-testid="stToolbarActions"],
            [data-testid="stAppDeployButton"],
            [data-testid="stMainMenu"] {
                display: none !important;
            }
            [data-testid="stExpandSidebarButton"] {
                align-items: center !important;
                background: #0f172a !important;
                border: 1px solid #334155 !important;
                border-radius: 8px !important;
                box-shadow: 0 4px 12px rgba(15, 23, 42, 0.18) !important;
                display: flex !important;
                height: 2.75rem !important;
                justify-content: center !important;
                min-height: 2.75rem !important;
                min-width: 2.75rem !important;
                padding: 0 !important;
                visibility: visible !important;
                width: 2.75rem !important;
            }
            [data-testid="stExpandSidebarButton"] span {
                color: #ffffff !important;
            }
            [data-testid="stExpandSidebarButton"]::after {
                clip-path: inset(50%);
                content: "Open navigation";
                height: 1px;
                overflow: hidden;
                position: absolute;
                white-space: nowrap;
                width: 1px;
            }
            [data-testid="stSidebar"] {
                min-width: 300px !important;
                max-width: 300px !important;
            }
            [data-testid="stSidebar"][aria-expanded="true"] {
                transition: none !important;
                transform: none !important;
                width: 300px !important;
            }
            [data-testid="stSidebar"][aria-expanded="false"] {
                transition: none !important;
                transform: translateX(-300px) !important;
            }
            [data-testid="stMain"] .block-container {
                padding-top: 4.25rem !important;
            }
            [data-testid="stSidebarHeader"] {
                align-items: center !important;
                background: #0f172a !important;
                display: flex !important;
                height: 3rem !important;
                justify-content: flex-end !important;
                min-height: 3rem !important;
                padding: 0.18rem 0.28rem !important;
                position: sticky !important;
                top: 0 !important;
                z-index: 10 !important;
            }
            [data-testid="stSidebarCollapseButton"] {
                display: flex !important;
            }
            [data-testid="stSidebarCollapseButton"] button {
                align-items: center !important;
                background: #1e293b !important;
                border: 1px solid #475569 !important;
                border-radius: 8px !important;
                display: flex !important;
                height: 2.5rem !important;
                justify-content: center !important;
                min-height: 2.5rem !important;
                min-width: 2.5rem !important;
                width: 2.5rem !important;
            }
            [data-testid="stSidebarCollapseButton"] button span {
                color: #ffffff !important;
            }
            [data-testid="stSidebarCollapseButton"] button::after {
                clip-path: inset(50%);
                content: "Close navigation";
                height: 1px;
                overflow: hidden;
                position: absolute;
                white-space: nowrap;
                width: 1px;
            }
            .fp-card-grid, .fp-example-grid, .fp-summary-strip, .fp-tool-shell, .fp-tutorial-panel ol {
                grid-template-columns: 1fr;
            }
            [data-testid="stMain"] [data-testid="stHorizontalBlock"] {
                align-items: stretch !important;
                flex-direction: column !important;
                gap: 0.72rem !important;
                width: 100% !important;
            }
            [data-testid="stMain"] [data-testid="stHorizontalBlock"] > [data-testid="stColumn"] {
                flex: 1 1 auto !important;
                min-width: 0 !important;
                width: 100% !important;
            }
            .fp-run-card {
                position: static;
            }
            .fp-run-summary {
                grid-template-columns: 1fr;
            }
            .fp-run-summary div:first-child {
                grid-column: auto;
            }
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def _ensure_session_config() -> None:
    if "config_revision" not in st.session_state:
        st.session_state.config_revision = 0
    if "current_config" not in st.session_state:
        st.session_state.current_config = make_single_config(
            "atac-correct",
            {"bams": [], "genome": "", "peaks": "", "blacklist": "", "outdir": "", "cores": 1},
            job_id="run",
        )
    if "gui_run_dir" not in st.session_state:
        st.session_state.gui_run_dir = os.environ.get("FP_TOOLS_GUI_RUN_DIR", str(default_gui_run_dir()))
    if "touched_tools" not in st.session_state:
        st.session_state.touched_tools = set()


def _tool_layout(label: str):
    """Keep form and run controls together, then stack them on narrower screens."""

    with st.container(key=f"fp_tool_shell_{label.replace('-', '_')}"):
        return st.columns([0.67, 0.33], gap="large")


def _config_widget_key(name: str) -> str:
    revision = int(st.session_state.get("config_revision", 0))
    return f"cfg_{revision}_{name}"


def _render_config_update_notice() -> None:
    message = st.session_state.pop("config_update_notice", "")
    if message:
        st.toast(str(message))


def _tool_pages() -> set[str]:
    return set(PAGE_OPTIONS).difference({"Home", "Run History", "Config"})


def _sync_config_for_page(page: str) -> None:
    if page not in _tool_pages():
        return
    current_tool = _current_config_tool()
    canonical_page = canonical_tool_name(page)
    if current_tool != canonical_page:
        _set_config(
            _default_config_for_tool(canonical_page),
            rerun=False,
            notify=False,
        )


def _current_config_tool() -> str:
    try:
        normalized = normalize_config(st.session_state.current_config)
    except Exception:
        return ""
    for section in ("samples", "comparisons"):
        for item in normalized.get(section, []):
            tool = str(item.get("tool", "")).strip()
            if tool:
                return canonical_tool_name(tool)
    return ""


def _default_config_for_tool(tool: str) -> dict[str, Any]:
    tool = canonical_tool_name(tool)
    defaults: dict[str, Any]
    if tool == "atac-correct":
        defaults = {"bams": [], "genome": "", "peaks": "", "blacklist": "", "outdir": "", "cores": 1}
    elif tool == "call-footprints":
        defaults = {"signal": "", "regions": "", "output": "", "score": "footprint", "cores": 1}
    elif tool == "diff-footprints":
        defaults = {
            "comparison_axis": "conditions",
            "motifs": "",
            "motif_db": "jaspar2026_vertebrates",
            "signals": [],
            "genome": "",
            "peaks": "",
            "peak_header": "",
            "outdir": "",
            "cond_names": ["Bcell"],
            "cores": 1,
            "skip_excel": False,
        }
    elif tool == "plot-aggregate":
        defaults = {"TFBS": [], "signals": [], "output": "", "grid": "", "output_aggregated_scores": "", "output_aggregated_signals": ""}
    else:
        raw_defaults = GENERIC_TOOL_DEFAULTS.get(tool, {"sample_id": f"{tool.replace('-', '_')}_run"})
        defaults = _drop_single_meta(raw_defaults)
    return make_single_config(tool, defaults, job_id=f"{tool.replace('-', '_')}_run")


def _current_page_from_query() -> str:
    raw_page = st.query_params.get("page", "Home")
    if isinstance(raw_page, list):
        raw_page = raw_page[0] if raw_page else "Home"
    page = str(raw_page)
    if page not in PAGE_OPTIONS:
        session_page = st.session_state.get("gui_page")
        page = str(session_page) if session_page in PAGE_OPTIONS else "Home"
    st.session_state.gui_page = page
    return page


def _render_sidebar_header() -> None:
    st.sidebar.markdown(
        """
        <div class="fp-sidebar-brand">
          <div class="fp-sidebar-brand-title">fp-tools</div>
        </div>
        """,
        unsafe_allow_html=True,
    )


def _render_sidebar_nav(active_page: str) -> None:
    for group, pages in NAV_GROUPS:
        st.sidebar.markdown(f'<div class="fp-nav-group">{escape(group)}</div>', unsafe_allow_html=True)
        for page in pages:
            key = "nav_" + page.replace(" ", "_").replace("-", "_")
            if st.sidebar.button(page, key=key, width="stretch", disabled=page == active_page):
                st.session_state.gui_page = page
                st.query_params["page"] = page
                st.rerun()


def _render_home(run_dir: Path) -> None:
    st.markdown(
        """
        <section class="fp-hero">
          <h1>fp-tools</h1>
          <p>Choose a workflow, enter your files, and run ATAC-seq footprinting.</p>
        </section>
        """,
        unsafe_allow_html=True,
    )
    columns = st.columns(3)
    actions = [
        ("Bulk workflow", "BAM and peak inputs", "bulk-footprinting"),
        ("Single-cell workflow", "Fragments and cell annotations", "sc-footprinting"),
        ("Load YAML", "Open or edit a saved workflow", "Config"),
    ]
    for column, (title, description, page) in zip(columns, actions):
        with column:
            st.markdown(f"### {title}\n{description}")
            if st.button(f"Open {title}", key=f"home_{page}", width="stretch"):
                st.session_state.gui_page = page
                st.query_params["page"] = page
                st.rerun()


def _render_home_config_snapshot() -> None:
    normalized = normalize_config(st.session_state.current_config)
    st.markdown(
        f"""
        <div class="fp-section-title">Current config snapshot</div>
        <div class="fp-summary-strip">
          <div class="fp-summary-metric"><div class="fp-summary-label">Run mode</div><div class="fp-summary-value">{escape(str(normalized["run_mode"]))}</div></div>
          <div class="fp-summary-metric"><div class="fp-summary-label">Sample jobs</div><div class="fp-summary-value">{len(normalized["samples"])}</div></div>
          <div class="fp-summary-metric"><div class="fp-summary-label">Comparisons</div><div class="fp-summary-value">{len(normalized["comparisons"])}</div></div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    with st.expander("Preview runnable YAML", expanded=False):
        st.code(config_to_yaml_text(normalized), language="yaml")


def _render_tutorial_panel() -> None:
    st.markdown(
        """
        <section class="fp-tutorial-panel">
          <h2>Guided workflow</h2>
          <p>Use the sidebar, load an example, check the YAML, then launch.</p>
          <ol>
            <li><b>Choose command.</b> Pick a tool from the sidebar.</li>
            <li><b>Load example.</b> Use bundled YAML or paste paths.</li>
            <li><b>Review YAML.</b> Confirm the config in the run panel.</li>
            <li><b>Inspect outputs.</b> Use Run History for logs and files.</li>
          </ol>
        </section>
        """,
        unsafe_allow_html=True,
    )
    if st.button("Hide guided tutorial", key="hide_tutorial_button"):
        st.session_state.show_tutorial = False
        st.session_state.hide_tutorial = True
        st.rerun()


def _render_page_heading(title: str, description: str = "") -> None:
    desc_html = f"<p>{escape(description)}</p>" if description else ""
    st.markdown(
        f"""
        <section class="fp-hero fp-page-heading">
          <h1>{escape(title)}</h1>
          {desc_html}
        </section>
        """,
        unsafe_allow_html=True,
    )


def _render_run_history(run_dir: Path) -> None:
    _render_page_heading("Run History", "Inspect background jobs, commands, logs, child runs, and detected outputs.")
    rows = []
    if run_dir.exists():
        for child in sorted(run_dir.iterdir(), reverse=True):
            status_path = child / "status.json"
            if status_path.exists():
                try:
                    status = refresh_run_status(child) or {}
                except json.JSONDecodeError:
                    continue
                rows.append(
                    {
                        "run_dir": str(child),
                        "tool": status.get("tool", ""),
                        "job_id": status.get("job_id", child.name),
                        "status": status.get("status", ""),
                        "exit_code": str(status.get("exit_code", "")),
                        "started_at": status.get("started_at", ""),
                        "finished_at": status.get("finished_at", ""),
                    }
                )
    if not rows:
        st.info("No runs yet. Start a workflow to see its status, logs, and outputs here.")
        left, right = st.columns(2)
        with left:
            if st.button("Start bulk workflow", width="stretch"):
                st.session_state.gui_page = "bulk-footprinting"
                st.query_params["page"] = "bulk-footprinting"
                st.rerun()
        with right:
            if st.button("Start single-cell workflow", width="stretch"):
                st.session_state.gui_page = "sc-footprinting"
                st.query_params["page"] = "sc-footprinting"
                st.rerun()
        return

    history = pd.DataFrame(rows)
    st.dataframe(history, width="stretch", hide_index=True)
    selected = st.selectbox("Inspect run", options=[""] + history["run_dir"].tolist())
    if not selected:
        return
    selected_path = Path(selected)
    st.code(
        (selected_path / "command.txt").read_text(encoding="utf-8")
        if (selected_path / "command.txt").exists()
        else "",
        language="bash",
    )
    batch_index = selected_path / "batch_index.tsv"
    if batch_index.exists():
        try:
            batch_df = pd.read_csv(batch_index, sep="\t")
            if "exit_code" in batch_df.columns:
                batch_df["exit_code"] = batch_df["exit_code"].fillna("").astype(str)
            st.dataframe(batch_df, width="stretch", hide_index=True)
        except Exception:
            st.text(batch_index.read_text(encoding="utf-8"))
    cols = st.columns(2)
    with cols[0]:
        if (selected_path / "launcher_stdout.log").exists():
            st.text_area("launcher stdout", value=(selected_path / "launcher_stdout.log").read_text(encoding="utf-8"), height=180)
    with cols[1]:
        if (selected_path / "launcher_stderr.log").exists():
            st.text_area("launcher stderr", value=(selected_path / "launcher_stderr.log").read_text(encoding="utf-8"), height=180)

    child_dirs = sorted([path for path in selected_path.iterdir() if path.is_dir() and (path / "status.json").exists()])
    if child_dirs:
        child_choice = st.selectbox("Inspect child job", options=[""] + [path.name for path in child_dirs], key=f"child_{selected}")
        if child_choice:
            child_path = selected_path / child_choice
            child_status = json.loads((child_path / "status.json").read_text(encoding="utf-8"))
            st.json(child_status)
            outputs = _discover_outputs(child_path)
            if outputs:
                st.subheader("Detected outputs")
                output_rows = [{"path": str(path), "exists": path.exists(), "size_bytes": path.stat().st_size if path.exists() else ""} for path in outputs]
                st.dataframe(pd.DataFrame(output_rows), width="stretch", hide_index=True)
            cols = st.columns(2)
            with cols[0]:
                if (child_path / "stdout.log").exists():
                    st.text_area("child stdout", value=(child_path / "stdout.log").read_text(encoding="utf-8"), height=260, key=f"child_stdout_{child_choice}")
            with cols[1]:
                if (child_path / "stderr.log").exists():
                    st.text_area("child stderr", value=(child_path / "stderr.log").read_text(encoding="utf-8"), height=260, key=f"child_stderr_{child_choice}")


def _render_atacorrect_page(run_dir: Path) -> None:
    _render_page_heading("atac-correct", "Bias-correct ATAC-seq cut-site signal from BAM, genome, peaks, and blacklist inputs.")
    form_col, run_col = _tool_layout("atac-correct")
    with form_col:
        _render_page_loader("atac-correct")
        mode_options = ["Single run", "Batch sample list"]
        mode_default = _config_form_mode("atac-correct")
        mode = st.radio(
            "Mode",
            mode_options,
            index=mode_options.index(mode_default),
            horizontal=True,
            key=_config_widget_key("at_mode"),
        )
        single = _current_single_params("atac-correct")
        batch_rows = _current_sample_rows(
            "atac-correct",
            default_rows=[{"sample_id": "sample1", "bams": "", "genome": "", "peaks": "", "blacklist": "", "outdir": "", "cores": 1}],
        )
        if mode == "Single run":
            with st.form(_config_widget_key("ataccorrect_single_form")):
                left, right = st.columns(2)
                with left:
                    bams_value = single.get("bams", "")
                    if isinstance(bams_value, list):
                        bams_value = "\n".join(str(value) for value in bams_value)
                    bams = st.text_area(
                        "BAM files (one per line)",
                        value=str(bams_value),
                        height=82,
                        help="One BAM path per line",
                        key=_config_widget_key("atacorrect_bams"),
                    )
                    genome = st.text_input(
                        "Genome FASTA",
                        value=str(single.get("genome", "")),
                        key=_config_widget_key("atacorrect_genome"),
                    )
                    peaks = st.text_input(
                        "Peaks BED",
                        value=str(single.get("peaks", "")),
                        key=_config_widget_key("atacorrect_peaks"),
                    )
                with right:
                    blacklist = st.text_input(
                        "Blacklist BED (optional)",
                        value=str(single.get("blacklist", "")),
                        key=_config_widget_key("atacorrect_blacklist"),
                    )
                    outdir = st.text_input(
                        "Output directory",
                        value=str(single.get("outdir", "")),
                        key=_config_widget_key("atacorrect_outdir"),
                    )
                    cores = st.number_input(
                        "Cores",
                        min_value=1,
                        value=int(single.get("cores", 1)),
                        step=1,
                        key=_config_widget_key("atacorrect_cores"),
                    )
                submitted = st.form_submit_button("Update page config")
            if submitted:
                _set_config(
                    _updated_single_config(
                        "atac-correct",
                        {
                            "bams": [line.strip() for line in bams.splitlines() if line.strip()],
                            "genome": genome,
                            "peaks": peaks,
                            "blacklist": blacklist,
                            "outdir": outdir,
                            "cores": int(cores),
                        },
                        job_id="atacorrect_run",
                    )
                )
        else:
            rows = _data_editor("atac-correct sample list", batch_rows, key="atacorrect_batch_editor")
            if st.button("Update page config from sample list", key="atacorrect_batch_set"):
                _set_config(
                    {
                        "version": 1,
                        "run_mode": "batch",
                        "defaults": {},
                        "samples": [{"tool": "atac-correct", **row} for row in rows],
                        "comparisons": [],
                    }
                )
    with run_col:
        _render_run_controls(run_dir, label="atacorrect")


def _render_footprintscores_page(run_dir: Path) -> None:
    _render_page_heading("call-footprints", "Score footprint signal over genomic regions from corrected signal tracks.")
    form_col, run_col = _tool_layout("call-footprints")
    with form_col:
        _render_page_loader("call-footprints")
        mode_options = ["Single run", "Batch sample list"]
        mode_default = _config_form_mode("call-footprints")
        mode = st.radio(
            "Mode",
            mode_options,
            index=mode_options.index(mode_default),
            horizontal=True,
            key=_config_widget_key("fs_mode"),
        )
        single = _current_single_params("call-footprints")
        batch_rows = _current_sample_rows(
            "call-footprints",
            default_rows=[{"sample_id": "sample1", "signal": "", "regions": "", "output": "", "score": "footprint", "cores": 1}],
        )
        if mode == "Single run":
            with st.form(_config_widget_key("footprintscores_single_form")):
                left, right = st.columns([0.68, 0.32])
                with left:
                    signal = st.text_input(
                        "Signal bigWig",
                        value=str(single.get("signal", "")),
                        key=_config_widget_key("footprintscores_signal"),
                    )
                    regions = st.text_input(
                        "Regions BED",
                        value=str(single.get("regions", "")),
                        key=_config_widget_key("footprintscores_regions"),
                    )
                    output = st.text_input(
                        "Output bigWig",
                        value=str(single.get("output", "")),
                        key=_config_widget_key("footprintscores_output"),
                    )
                with right:
                    score_values = ["footprint", "sum", "mean", "none"]
                    score_default = str(single.get("score", "footprint"))
                    score = st.selectbox(
                        "Score",
                        score_values,
                        index=score_values.index(score_default) if score_default in score_values else 0,
                        key=_config_widget_key("footprintscores_score"),
                    )
                    cores = st.number_input(
                        "Cores",
                        min_value=1,
                        value=int(single.get("cores", 1)),
                        step=1,
                        key=_config_widget_key("fs_single_cores"),
                    )
                submitted = st.form_submit_button("Update page config")
            if submitted:
                _set_config(
                    _updated_single_config(
                        "call-footprints",
                        {
                            "signal": signal,
                            "regions": regions,
                            "output": output,
                            "score": score,
                            "cores": int(cores),
                        },
                        job_id="footprintscores_run",
                    )
                )
        else:
            rows = _data_editor("call-footprints sample list", batch_rows, key="footprintscores_batch_editor")
            if st.button("Update page config from sample list", key="footprintscores_batch_set"):
                _set_config(
                    {
                        "version": 1,
                        "run_mode": "batch",
                        "defaults": {},
                        "samples": [{"tool": "call-footprints", **row} for row in rows],
                        "comparisons": [],
                    }
                )
    with run_col:
        _render_run_controls(run_dir, label="footprintscores")


def _render_diff_footprints_page(run_dir: Path) -> None:
    _render_page_heading("diff-footprints", "Run motif-aware differential footprint detection across conditions or comparisons.")
    form_col, run_col = _tool_layout("diff-footprints")
    with form_col:
        _render_page_loader("diff-footprints")
        mode_options = [
            "Single condition",
            "Batch single-condition list",
            "Batch comparison list",
        ]
        mode_default = _config_form_mode("diff-footprints")
        mode = st.radio(
            "Comparison setup",
            mode_options,
            index=mode_options.index(mode_default),
            horizontal=True,
            format_func={
                "Single condition": "One comparison",
                "Batch single-condition list": "Batch samples",
                "Batch comparison list": "Batch comparisons",
            }.get,
            key=_config_widget_key("diff_footprints_mode"),
        )
        single = _current_single_params("diff-footprints")
        sample_rows = _current_sample_rows(
            "diff-footprints",
            default_rows=[
                {
                    "sample_id": "sample1",
                    "motifs": "",
                    "motif_db": "jaspar2026_vertebrates",
                    "signals": "",
                    "genome": "",
                    "peaks": "",
                    "peak_header": "",
                    "outdir": "",
                    "cond_names": "Sample1",
                    "cores": 1,
                    "skip_excel": False,
                }
            ],
        )
        comparison_rows = _current_comparison_rows(
            "diff-footprints",
            default_rows=[
                {
                    "comparison_id": "bcell_vs_tcell",
                    "motifs": "",
                    "motif_db": "jaspar2026_vertebrates",
                    "signals": "",
                    "cond_names": "Bcell,Bcell,Tcell,Tcell",
                    "genome": "",
                    "peaks": "",
                    "peak_header": "",
                    "outdir": "",
                    "cores": 1,
                    "skip_excel": False,
                }
            ],
        )
        if mode == "Single condition":
            comparison_axis = st.selectbox(
                "Comparison axis",
                ["conditions", "regions"],
                index=0 if str(single.get("comparison_axis", "conditions")) == "conditions" else 1,
                help="Compare samples or compare region sets measured in the same sample(s).",
                key=_config_widget_key("diff_footprints_comparison_axis"),
            )
            with st.form(_config_widget_key("diff_footprints_single_form")):
                left, right = st.columns(2)
                with left:
                    motifs_value = single.get("motifs", "")
                    if isinstance(motifs_value, list):
                        motifs_value = _join_multi(motifs_value)
                    motifs = st.text_area(
                        "Motif files (one per line)",
                        value=str(motifs_value),
                        height=76,
                        key=_config_widget_key("diff_footprints_motifs"),
                    )
                    motif_db = st.text_input(
                        "Motif database",
                        value=str(single.get("motif_db", "jaspar2026_vertebrates")),
                        key=_config_widget_key("diff_footprints_motif_db"),
                    )
                    genome = st.text_input(
                        "Genome FASTA",
                        value=str(single.get("genome", "")),
                        key=_config_widget_key("diff_footprints_genome"),
                    )
                    peaks = st.text_input(
                        "Peaks BED",
                        value=str(single.get("peaks", "")),
                        key=_config_widget_key("diff_footprints_peaks"),
                    )
                with right:
                    peak_header = st.text_input(
                        "Peak annotation header (optional)",
                        value=str(single.get("peak_header", "")),
                        key=_config_widget_key("diff_footprints_peak_header"),
                    )
                    outdir = st.text_input(
                        "Output directory",
                        value=str(single.get("outdir", "")),
                        key=_config_widget_key("diff_footprints_outdir"),
                    )
                    cores = st.number_input(
                        "Cores",
                        min_value=1,
                        value=int(single.get("cores", 1)),
                        step=1,
                        key=_config_widget_key("diff_footprints_single_cores"),
                    )
                    skip_excel = st.checkbox(
                        "Skip Excel",
                        value=bool(single.get("skip_excel", False)),
                        key=_config_widget_key("diff_footprints_skip_excel"),
                    )
                signals = st.text_area(
                    "Footprint bigWig files",
                    value=_join_multi(single.get("signals", [])),
                    height=94,
                    key=_config_widget_key("diff_footprints_signals"),
                )
                cond_names = _join_multi(single.get("cond_names", ["Bcell"]))
                regions = _join_multi(single.get("regions", []))
                region_labels = ",".join(single.get("region_labels", []))
                region_strata_column = int(single.get("region_strata_column", 0) or 0)
                if comparison_axis == "conditions":
                    cond_names = st.text_area(
                        "Condition names",
                        value=cond_names,
                        height=76,
                        key=_config_widget_key("diff_footprints_cond_names"),
                    )
                else:
                    regions = st.text_area(
                        "Region BED files",
                        value=regions,
                        height=76,
                        help="One region set per line.",
                        key=_config_widget_key("diff_footprints_regions"),
                    )
                    region_labels = st.text_input(
                        "Region labels (optional)",
                        value=region_labels,
                        help="Comma-separated labels in the same order as the BED files.",
                        key=_config_widget_key("diff_footprints_region_labels"),
                    )
                    region_strata_column = st.number_input(
                        "Matching-stratum BED column (0 = none)",
                        min_value=0,
                        value=region_strata_column,
                        step=1,
                        key=_config_widget_key("diff_footprints_region_strata_column"),
                    )
                submitted = st.form_submit_button("Update page config")
            if submitted:
                _set_config(
                    _updated_single_config(
                        "diff-footprints",
                        {
                            "comparison_axis": comparison_axis,
                            "motifs": _split_multi(motifs),
                            "motif_db": motif_db,
                            "signals": _split_multi(signals),
                            "genome": genome,
                            "peaks": peaks,
                            "peak_header": peak_header,
                            "outdir": outdir,
                            "cond_names": _split_multi(cond_names),
                            "regions": _split_multi(regions),
                            "region_labels": _split_multi(region_labels),
                            "region_strata_column": int(region_strata_column) or None,
                            "cores": int(cores),
                            "skip_excel": bool(skip_excel),
                        },
                        job_id="diff_footprints_single",
                    )
                )
        elif mode == "Batch single-condition list":
            rows = _data_editor("diff-footprints single-condition sample list", sample_rows, key="diff_footprints_sample_editor")
            if st.button("Update page config from single-condition list", key="diff_footprints_sample_set"):
                _set_config(
                    {
                        "version": 1,
                        "run_mode": "batch",
                        "defaults": {},
                        "samples": [
                            {
                                "tool": "diff-footprints",
                                **row,
                                "signals": _split_multi(str(row.get("signals", ""))),
                                "cond_names": _split_multi(str(row.get("cond_names", ""))),
                                "skip_excel": _as_bool(row.get("skip_excel", False)),
                            }
                            for row in rows
                        ],
                        "comparisons": [],
                    }
                )
        else:
            rows = _data_editor("diff-footprints comparison list", comparison_rows, key="diff_footprints_comparison_editor")
            if st.button("Update page config from comparison list", key="diff_footprints_comparison_set"):
                _set_config(
                    {
                        "version": 1,
                        "run_mode": "batch",
                        "defaults": {},
                        "samples": [],
                        "comparisons": [
                            {
                                "tool": "diff-footprints",
                                **row,
                                "signals": _split_multi(str(row.get("signals", ""))),
                                "cond_names": _split_multi(str(row.get("cond_names", ""))),
                                "skip_excel": _as_bool(row.get("skip_excel", False)),
                            }
                            for row in rows
                        ],
                    }
                )
    with run_col:
        _render_run_controls(run_dir, label="diff-footprints")


def _render_plotaggregate_page(run_dir: Path) -> None:
    _render_page_heading("plot-aggregate", "Create aggregate footprint plots and optional aggregate score tables.")
    form_col, run_col = _tool_layout("plot-aggregate")
    with form_col:
        _render_page_loader("plot-aggregate")
        mode_options = ["Single run", "Batch sample list"]
        mode_default = _config_form_mode("plot-aggregate")
        mode = st.radio(
            "Mode",
            mode_options,
            index=mode_options.index(mode_default),
            horizontal=True,
            key=_config_widget_key("pa_mode"),
        )
        single = _current_single_params("plot-aggregate")
        batch_rows = _current_sample_rows(
            "plot-aggregate",
            default_rows=[
                {
                    "sample_id": "panel1",
                    "TFBS": "",
                    "signals": "",
                    "output": "",
                    "grid": "",
                    "output_aggregated_scores": "",
                    "output_aggregated_signals": "",
                }
            ],
        )
        if mode == "Single run":
            with st.form(_config_widget_key("plotaggregate_single_form")):
                left, right = st.columns(2)
                with left:
                    tfbs = st.text_area(
                        "Region BED files",
                        value=_join_multi(single.get("TFBS", [])),
                        height=92,
                        key=_config_widget_key("plotaggregate_tfbs"),
                    )
                    signals = st.text_area(
                        "Signal bigWig files",
                        value=_join_multi(single.get("signals", [])),
                        height=92,
                        key=_config_widget_key("plotaggregate_signals"),
                    )
                with right:
                    output = st.text_input(
                        "Output PDF",
                        value=str(single.get("output", "")),
                        key=_config_widget_key("plotaggregate_output"),
                    )
                    grid = st.text_input(
                        "Panel grid (optional; for example, 2x5)",
                        value=str(single.get("grid", "")),
                        key=_config_widget_key("plotaggregate_grid"),
                    )
                    score_csv = st.text_input(
                        "Aggregated score CSV (optional)",
                        value=str(single.get("output_aggregated_scores", "")),
                        key=_config_widget_key("plotaggregate_score_csv"),
                    )
                    signal_csv = st.text_input(
                        "Aggregated signal CSV (optional)",
                        value=str(single.get("output_aggregated_signals", "")),
                        key=_config_widget_key("plotaggregate_signal_csv"),
                    )
                submitted = st.form_submit_button("Update page config")
            if submitted:
                config: dict[str, Any] = {
                    "TFBS": _split_multi(tfbs),
                    "signals": _split_multi(signals),
                    "output": output,
                    "grid": grid.strip(),
                    "output_aggregated_scores": score_csv.strip(),
                    "output_aggregated_signals": signal_csv.strip(),
                }
                _set_config(_updated_single_config("plot-aggregate", config, job_id="plotaggregate_run"))
        else:
            rows = _data_editor("plot-aggregate panel list", batch_rows, key="plotaggregate_batch_editor")
            if st.button("Update page config from panel list", key="plotaggregate_batch_set"):
                _set_config(
                    {
                        "version": 1,
                        "run_mode": "batch",
                        "defaults": {},
                        "samples": [
                            {
                                "tool": "plot-aggregate",
                                **row,
                                "TFBS": _split_multi(str(row.get("TFBS", ""))),
                                "signals": _split_multi(str(row.get("signals", ""))),
                            }
                            for row in rows
                        ],
                        "comparisons": [],
                    }
                )
    with run_col:
        _render_run_controls(run_dir, label="plotaggregate")


def _render_generic_tool_page(run_dir: Path, tool: str) -> None:
    tool = canonical_tool_name(tool)
    _render_page_heading(tool, GUI_TOOL_DESCRIPTIONS.get(tool, "Configure and run this fp-tools command."))
    form_col, run_col = _tool_layout(tool)
    with form_col:
        _render_page_loader(tool)
        defaults = GENERIC_TOOL_DEFAULTS.get(tool, {"sample_id": f"{tool.replace('-', '_')}_run"})
        current = {**_drop_single_meta(defaults), **_current_single_params(tool)}
        with st.form(_config_widget_key(f"{tool}_generic_form")):
            edited: dict[str, Any] = {}
            core_simple: list[tuple[str, Any]] = []
            core_wide: list[tuple[str, Any]] = []
            advanced: list[tuple[str, Any]] = []
            advanced_keys = GUI_ADVANCED_FIELDS.get(tool, set())
            for key, default_value in current.items():
                if key in {"tool", "sample_id", "job_id", "comparison_id", "extra_args"}:
                    continue
                if key in advanced_keys:
                    advanced.append((key, default_value))
                elif key in LIST_TEXT_FIELDS or isinstance(default_value, list):
                    core_wide.append((key, default_value))
                else:
                    core_simple.append((key, default_value))
            columns = st.columns(2)
            for idx, (key, default_value) in enumerate(core_simple):
                with columns[idx % 2]:
                    edited[key] = _render_generic_field(tool, key, default_value)
            for key, default_value in core_wide:
                edited[key] = _render_generic_field(tool, key, default_value)
            with st.expander("Advanced options", expanded=False):
                advanced_columns = st.columns(2)
                for idx, (key, default_value) in enumerate(advanced):
                    with advanced_columns[idx % 2]:
                        edited[key] = _render_generic_field(tool, key, default_value)
                extra_args = st.text_input(
                    "Additional CLI arguments",
                    value=_format_extra_args(current.get("extra_args", [])),
                    help="Optional arguments not represented by the fields above.",
                    key=_config_widget_key(f"{tool}_extra_args"),
                )
            submitted = st.form_submit_button("Update page config")
        if submitted:
            config = _prepare_generic_params(tool, edited)
            if extra_args.strip():
                config["extra_args"] = _parse_extra_args(extra_args)
            elif "extra_args" in current:
                config["extra_args"] = []
            rendered_fields = set(edited).union({"extra_args"})
            merged = dict(current)
            for key in rendered_fields:
                if key in config:
                    merged[key] = config[key]
                else:
                    merged.pop(key, None)
            _set_config(
                _updated_single_config(
                    tool,
                    merged,
                    job_id=f"{tool.replace('-', '_')}_run",
                )
            )
    with run_col:
        _render_run_controls(run_dir, label=tool.replace("-", "_"))


def _render_config_page(run_dir: Path) -> None:
    _render_page_heading("Config", "Download, load, edit, save, and run YAML workflow configs.")
    form_col, run_col = _tool_layout("config")
    with form_col:
        st.download_button(
            "Download current YAML",
            data=config_to_yaml_text(st.session_state.current_config),
            file_name="fp_tools_config.yml",
            mime="text/yaml",
            width="stretch",
        )
        load_left, load_right = st.columns(2)
        with load_left:
            uploader = st.file_uploader("Load YAML file", type=["yml", "yaml"], key="config_uploader")
            if uploader is not None and st.button("Apply uploaded YAML", key="config_apply_upload"):
                _set_config(parse_yaml_text(uploader.getvalue().decode("utf-8")))
        with load_right:
            load_path = st.text_input("Load config from path", key="config_load_path")
            if st.button("Load YAML from path", key="config_load_path_btn") and load_path.strip():
                _set_config(normalize_config(load_yaml_config(load_path.strip())))

        yaml_text = st.text_area(
            "Current YAML",
            value=config_to_yaml_text(st.session_state.current_config),
            height=360,
            key=_config_widget_key("config_yaml_text"),
        )
        edit_left, edit_right = st.columns(2)
        with edit_left:
            if st.button("Apply YAML text", key="config_apply_text"):
                _set_config(parse_yaml_text(yaml_text))
        with edit_right:
            save_path = st.text_input("Save current YAML to path", key="config_save_path")
            if st.button("Save YAML", key="config_save_btn") and save_path.strip():
                path = Path(save_path.strip()).expanduser()
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(yaml_text, encoding="utf-8")
                st.success(f"Saved config to {path}")
    with run_col:
        _render_run_controls(run_dir, label="config")


def _render_page_loader(tool: str) -> None:
    tool = canonical_tool_name(tool)
    with st.expander(f"Load {tool} config", expanded=False):
        example_files = _example_files_for_tool(tool)
        if example_files:
            example_choice = st.selectbox(
                "Example YAML",
                options=[""] + [path.name for path in example_files],
                key=f"{tool}_example_select",
            )
            if st.button("Load example", key=f"{tool}_load_example") and example_choice:
                example = next(path for path in example_files if path.name == example_choice)
                _set_config(parse_yaml_text(example.read_text(encoding="utf-8")))
        upload = st.file_uploader("Upload YAML", type=["yml", "yaml"], key=f"{tool}_uploader")
        if upload is not None and st.button("Apply uploaded YAML", key=f"{tool}_apply_upload"):
            _set_config(parse_yaml_text(upload.getvalue().decode("utf-8")))
        path_text = st.text_input("Config path", key=f"{tool}_config_path")
        if st.button("Load YAML from path", key=f"{tool}_load_path") and path_text.strip():
            _set_config(normalize_config(load_yaml_config(path_text.strip())))


def _validation_errors_markup(messages: list[str]) -> str:
    items = "".join(f"<li>{escape(message)}</li>" for message in messages)
    return (
        '<ul class="fp-validation-errors" '
        f'aria-label="Configuration validation errors">{items}</ul>'
    )


def _friendly_validation_message(message: str) -> str:
    """Translate validator diagnostics into concise GUI guidance."""

    text = re.sub(r"^[A-Za-z0-9_.-]+:\s*", "", str(message), count=1)
    text = text.replace("file does not exist:", "File not found:")
    text = text.replace("directory does not exist:", "Folder not found:")
    text = text.replace("is required", "is required.")
    text = text.replace("'sample_table'", "Samples TSV")
    text = text.replace("'comparison_table'", "Comparisons TSV")
    text = text.replace("'extra_args'", "Additional CLI arguments")
    for field, label in GUI_FIELD_LABELS.items():
        text = text.replace(f"'{field}'", label)
    return text[:1].upper() + text[1:] if text else "Check this configuration."


def _render_run_controls(run_dir: Path, label: str) -> None:
    normalized = normalize_config(st.session_state.current_config)
    validation_errors = validate_gui_config(normalized)
    has_validation_errors = bool(validation_errors)
    tool = _current_config_tool() or "none"
    touched = tool in set(st.session_state.get("touched_tools", set()))
    friendly_errors = [_friendly_validation_message(message) for message in validation_errors]
    with st.container(border=True, key=f"fp_run_panel_{label}"):
        st.markdown(
            f"""
            <div class="fp-run-card">
              <h3>Run</h3>
              <div class="fp-run-summary">
                <div><span>Tool</span><strong>{escape(tool)}</strong></div>
                <div><span>Samples</span><strong>{len(normalized["samples"])}</strong></div>
                <div><span>Comparisons</span><strong>{len(normalized["comparisons"])}</strong></div>
              </div>
            </div>
            """,
            unsafe_allow_html=True,
        )
        if validation_errors and touched:
            st.error("Check the highlighted setup details before starting.")
            st.markdown(_validation_errors_markup(friendly_errors), unsafe_allow_html=True)
        elif validation_errors:
            st.info("Add the required inputs, then update the page config.")
        else:
            st.success("Config is ready to run.")
        with st.expander("Preview runnable YAML", expanded=False):
            st.code(config_to_yaml_text(normalized), language="yaml")
        if st.button(
            "Start run",
            key=f"run_{label}",
            width="stretch",
            disabled=has_validation_errors,
        ):
            if validation_errors:
                st.error("Run not started. Fix the config errors above first.")
                return
            run_dir_path, config_path = materialize_run_config(
                normalized,
                run_root=run_dir,
                label=label,
            )
            _status_path, pid = launch_config_async(config_path, run_dir_path, label)
            st.write(f"Run folder: `{run_dir_path}`")
            st.success(f"Run started in background (pid {pid}). Open Run History to monitor logs and status.")


def _show_current_summary() -> None:
    normalized = normalize_config(st.session_state.current_config)
    st.subheader("Current config summary")
    st.write(
        {
            "run_mode": normalized["run_mode"],
            "sample_jobs": len(normalized["samples"]),
            "comparison_jobs": len(normalized["comparisons"]),
        }
    )
    st.code(config_to_yaml_text(normalized), language="yaml")


def _render_sidebar_run_dir_controls() -> None:
    with st.sidebar.expander("Workspace", expanded=False):
        current_run_dir = str(Path(st.session_state.gui_run_dir).expanduser())
        st.markdown(
            f"""
            <div class="fp-workspace-meta">fp-tools v{escape(__version__)}</div>
            <div class="fp-workspace-meta">Current run directory</div>
            <div class="fp-workspace-path">{escape(current_run_dir)}</div>
            """,
            unsafe_allow_html=True,
        )
        run_dir_input = st.text_input("GUI run dir", value=str(st.session_state.gui_run_dir), key="sidebar_run_dir")
        if st.button("Apply run dir", key="sidebar_apply_run_dir", width="stretch"):
            path = Path(run_dir_input).expanduser()
            path.mkdir(parents=True, exist_ok=True)
            st.session_state.gui_run_dir = str(path)
            st.success(f"Run dir set to {path}")


def _discover_outputs(child_run_dir: Path) -> list[Path]:
    config_path = child_run_dir / "config.yml"
    if not config_path.exists():
        return []
    try:
        config = normalize_config(load_yaml_config(config_path))
    except Exception:
        return []

    outputs: list[Path] = []
    items = config["samples"] or config["comparisons"]
    if not items:
        return []
    item = items[0]
    tool = canonical_tool_name(str(item.get("tool", "")))
    base = Path.cwd()

    def add_path(value: Any) -> None:
        text = str(value).strip()
        if not text:
            return
        path = Path(text).expanduser()
        if not path.is_absolute():
            path = (base / path).resolve()
        outputs.append(path)

    if tool == "atac-correct":
        outdir = str(item.get("outdir", "")).strip()
        if outdir:
            outdir_path = Path(outdir).expanduser()
            if not outdir_path.is_absolute():
                outdir_path = (base / outdir_path).resolve()
            outputs.append(outdir_path)
            if outdir_path.exists():
                outputs.extend(sorted(path for path in outdir_path.iterdir() if path.is_file()))
    elif tool == "call-footprints":
        add_path(item.get("output", ""))
    elif tool == "plot-aggregate":
        for key in ("output", "output_aggregated_scores", "output_aggregated_signals", "output_csv"):
            add_path(item.get(key, ""))
    elif tool == "diff-footprints":
        outdir = str(item.get("outdir", "")).strip()
        if outdir:
            outdir_path = Path(outdir).expanduser()
            if not outdir_path.is_absolute():
                outdir_path = (base / outdir_path).resolve()
            outputs.append(outdir_path)
            if outdir_path.exists():
                preferred = [
                    "diff_footprints_results.txt",
                    "diff_footprints_figures.pdf",
                    "diff_footprints_clusters.pdf",
                    "diff_footprints_results_skewness_report.pdf",
                ]
                outputs.extend(outdir_path / name for name in preferred if (outdir_path / name).exists())
                outputs.extend(sorted(path for path in outdir_path.glob("diff_footprints_*.html")))
    elif tool == "normalize-bigwig":
        add_path(item.get("outdir", ""))
    elif tool == "discover-motifs":
        add_path(item.get("outdir", ""))
        add_path(item.get("script", ""))
    elif tool == "summarize-motifs":
        for key in ("out_tsv", "out_html"):
            add_path(item.get(key, ""))
    elif tool in {"pseudobulk-fragments", "find-signature-fp", "sc-footprinting", "bulk-footprinting"}:
        add_path(item.get("outdir", ""))
    elif tool == "review-multi-comparisons":
        add_path(item.get("output_dir", ""))
        add_path(item.get("output_html", ""))

    seen: set[str] = set()
    deduped: list[Path] = []
    for path in outputs:
        key = str(path)
        if key not in seen:
            seen.add(key)
            deduped.append(path)
    return deduped


def _prepare_generic_params(tool: str, params: dict[str, Any]) -> dict[str, Any]:
    prepared: dict[str, Any] = {}
    for key, value in params.items():
        if key in LIST_TEXT_FIELDS:
            prepared[key] = _split_multi(str(value))
        elif isinstance(value, str) and value.strip() == "":
            continue
        else:
            prepared[key] = value
    if tool == "plot-aggregate" and prepared.get("input_html"):
        prepared.pop("manifest", None)
    return prepared


def _render_generic_field(tool: str, key: str, default_value: Any) -> Any:
    value = default_value
    help_text = GUI_FIELD_HELP.get(tool, {}).get(key)
    if isinstance(value, list):
        value = _join_multi(value)
    if isinstance(value, bool):
        return st.checkbox(
            _human_label(key),
            value=value,
            key=_config_widget_key(f"{tool}_{key}"),
            help=help_text,
        )
    choices = GUI_ENUM_CHOICES.get(tool, {}).get(key)
    if choices:
        selected = str(value or choices[0])
        index = choices.index(selected) if selected in choices else 0
        return st.selectbox(
            _human_label(key),
            choices,
            index=index,
            key=_config_widget_key(f"{tool}_{key}"),
            help=help_text,
        )
    if key in LIST_TEXT_FIELDS:
        return st.text_area(
            _human_label(key),
            value=str(value or ""),
            height=88,
            key=_config_widget_key(f"{tool}_{key}"),
            help=help_text,
        )
    if isinstance(value, int):
        return int(
            st.number_input(
                _human_label(key),
                min_value=0,
                value=int(value),
                step=1,
                key=_config_widget_key(f"{tool}_{key}"),
                help=help_text,
            )
        )
    return st.text_input(
        _human_label(key),
        value=str(value or ""),
        key=_config_widget_key(f"{tool}_{key}"),
        help=help_text,
    )


def _human_label(key: str) -> str:
    if key in GUI_FIELD_LABELS:
        return GUI_FIELD_LABELS[key]
    replacements = {"h5ad": "h5ad", "id": "ID", "tsv": "TSV", "html": "HTML", "pdf": "PDF"}
    words = [replacements.get(word.lower(), word.capitalize()) for word in key.split("_")]
    return " ".join(words)


def _set_config(
    config: dict[str, Any],
    *,
    rerun: bool = True,
    notify: bool = True,
) -> None:
    st.session_state.current_config = normalize_config(config)
    if notify:
        touched = set(st.session_state.get("touched_tools", set()))
        current_tool = _current_config_tool()
        if current_tool:
            touched.add(current_tool)
        st.session_state.touched_tools = touched
    st.session_state.config_revision = int(st.session_state.get("config_revision", 0)) + 1
    if notify:
        st.session_state.config_update_notice = "Current config updated."
    if rerun:
        st.rerun()


def _data_editor(title: str, default_rows: list[dict[str, Any]], key: str) -> list[dict[str, Any]]:
    st.subheader(title)
    df = pd.DataFrame(default_rows)
    edited = st.data_editor(
        df,
        num_rows="dynamic",
        width="stretch",
        key=_config_widget_key(key),
    )
    cleaned = edited.fillna("").to_dict(orient="records")
    return [row for row in cleaned if any(str(value).strip() for value in row.values())]


def _current_config() -> dict[str, Any]:
    return normalize_config(st.session_state.current_config)


def _config_form_mode(tool: str) -> str:
    tool = canonical_tool_name(tool)
    config = _current_config()
    if tool == "diff-footprints":
        if config["comparisons"]:
            return "Batch comparison list"
        if config["run_mode"] == "batch" or len(config["samples"]) > 1:
            return "Batch single-condition list"
        return "Single condition"
    if config["run_mode"] == "batch" or len(config["samples"]) > 1:
        return "Batch sample list"
    return "Single run"


def _updated_single_config(
    tool: str,
    params: dict[str, Any],
    *,
    job_id: str,
) -> dict[str, Any]:
    tool = canonical_tool_name(tool)
    config = _current_config()
    if not config["comparisons"] and len(config["samples"]) == 1:
        item = dict(config["samples"][0])
        if canonical_tool_name(str(item.get("tool", ""))) == tool:
            item.update(params)
            item["tool"] = tool
            return {
                **config,
                "run_mode": "single",
                "samples": [item],
                "comparisons": [],
            }
    return make_single_config(tool, params, job_id=job_id)


def _current_single_params(tool: str) -> dict[str, Any]:
    tool = canonical_tool_name(tool)
    config = _current_config()
    if config["comparisons"]:
        return {}
    matching = [item for item in config["samples"] if canonical_tool_name(str(item.get("tool", ""))) == tool]
    if len(matching) != 1:
        return {}
    return _drop_single_meta({**config["defaults"], **matching[0]})


def _current_sample_rows(tool: str, default_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    tool = canonical_tool_name(tool)
    config = _current_config()
    matching = [
        _prepare_row_for_editor({**config["defaults"], **item})
        for item in config["samples"]
        if canonical_tool_name(str(item.get("tool", ""))) == tool
    ]
    return matching or default_rows


def _current_comparison_rows(tool: str, default_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    tool = canonical_tool_name(tool)
    config = _current_config()
    matching = [
        _prepare_row_for_editor({**config["defaults"], **item})
        for item in config["comparisons"]
        if canonical_tool_name(str(item.get("tool", ""))) == tool
    ]
    return matching or default_rows


def _prepare_row_for_editor(item: dict[str, Any]) -> dict[str, Any]:
    row = _drop_editor_meta(item)
    for key, value in list(row.items()):
        if isinstance(value, list):
            row[key] = ",".join(str(v) for v in value)
    return row


def _drop_single_meta(item: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in item.items()
        if key not in {"tool", "sample_id", "comparison_id", "job_id", "label", "name", "description"}
    }


def _drop_editor_meta(item: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in item.items()
        if key not in {"tool", "job_id", "label", "name", "description"}
    }


def _all_example_files() -> list[Any]:
    root = resources.files(GUI_EXAMPLE_PACKAGE)
    return sorted(
        (path for path in root.iterdir() if path.is_file() and path.name.endswith(".yml")),
        key=lambda path: path.name,
    )


def _example_files_for_tool(tool: str) -> list[Any]:
    tool = canonical_tool_name(tool)
    prefixes = {
        "atac-correct": ["atacorrect", "atac_correct"],
        "call-footprints": ["call_footprints", "footprintscores"],
        "diff-footprints": ["diff_footprints"],
        "plot-aggregate": ["plotaggregate", "plot_aggregate"],
    }.get(tool, [tool.replace("-", "_")])
    files: list[Any] = []
    candidates = _all_example_files()
    for prefix in prefixes:
        files.extend(path for path in candidates if path.name.startswith(f"{prefix}_"))
    return sorted({path.name: path for path in files}.values(), key=lambda path: path.name)


def _parse_extra_args(text: str) -> list[str]:
    lexer = shlex.shlex(str(text), posix=True)
    lexer.whitespace_split = True
    lexer.commenters = ""
    lexer.escape = ""
    return list(lexer)


def _format_extra_args(value: Any) -> str:
    if isinstance(value, str):
        return value
    return shlex.join([str(item) for item in (value or [])])


def _tutorial_visible() -> bool:
    if st.session_state.get("hide_tutorial"):
        return False
    return bool(st.session_state.get("show_tutorial")) or os.environ.get("FP_TOOLS_GUI_TUTORIAL") == "1"


def _render_tutorial_overlay() -> None:
    _render_tutorial_panel()


def _split_multi(text: str) -> list[str]:
    raw = str(text).replace(";", "\n").replace(",", "\n").replace("|", "\n")
    return [line.strip() for line in raw.splitlines() if line.strip()]


def _join_multi(values: list[Any]) -> str:
    return "\n".join(str(value) for value in values if str(value).strip())


def _as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    return text in {"1", "true", "yes", "y", "on"}


if __name__ == "__main__":
    main()
