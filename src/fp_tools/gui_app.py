"""Streamlit GUI for fp-tools.

This module is an isolated wrapper around the packaged commands. Direct CLI
usage remains primary. The GUI supports direct form-driven runs, YAML load/save,
and batch editing while using the same normalized config model as the optional
``run-yaml-workflow --config ...`` path.
"""

from __future__ import annotations

import json
import os
from html import escape
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

from fp_tools import __version__
from fp_tools.gui_config import (
    canonical_tool_name,
    config_to_yaml_text,
    load_yaml_config,
    make_single_config,
    normalize_config,
    parse_yaml_text,
    validate_config,
)
from fp_tools.gui_jobs import default_gui_run_dir, launch_config_async, materialize_run_config, refresh_run_status

GUI_EXAMPLE_DIR = Path("examples/gui_configs")

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
        "reads_table": "reads.tsv",
        "sample_table": "",
        "comparison_table": "comparisons.tsv",
        "genome": "hg38",
        "outdir": "results/bulk_footprinting",
        "motif_db": "jaspar2026_vertebrates",
        "plot_aggregate": "all",
        "review_format": "auto",
        "cores": 4,
        "runtime": "auto",
    },
    "review-multi-comparisons": {
        "sample_id": "comparison_browser_run",
        "inputs": ["results/bulk_footprinting/comparisons"],
        "output_dir": "results/bulk_footprinting/reports/review_multi_comparisons",
    },
    "match-motifs": {
        "sample_id": "match_motifs_run",
        "signals": "test_data/Bcell_footprints.bw",
        "genome": "test_data/genome.fa.gz",
        "peaks": "test_data/merged_peaks_annotated.bed",
        "peak_header": "test_data/merged_peaks_annotated_header.txt",
        "outdir": "examples/gui_demo_outputs/match_motifs_single",
        "cond_names": "Bcell",
        "motif_db": "jaspar2026_vertebrates",
        "skip_excel": True,
    },
    "normalize-bigwig": {
        "sample_id": "normalize_bigwig_run",
        "bigwigs": "test_data/Bcell_corrected.bw\ntest_data/Tcell_corrected.bw",
        "background": "test_data/merged_peaks.bed",
        "outdir": "examples/gui_demo_outputs/normalize_bigwig",
        "method": "background-scale",
        "stat": "q95",
        "target": "median",
    },
    "discover-motifs": {
        "sample_id": "motif_discovery_run",
        "candidates": "test_data/merged_peaks.bed",
        "genome": "test_data/genome.fa.gz",
        "outdir": "examples/gui_demo_outputs/motif_discovery",
        "method": "streme",
        "known_motif_db": "jaspar2026_vertebrates",
    },
    "summarize-motifs": {
        "sample_id": "motif_summary_run",
        "meme_txt": "",
        "tomtom_tsv": "",
        "out_tsv": "examples/gui_demo_outputs/motif_summary/motif_summary.tsv",
        "out_html": "examples/gui_demo_outputs/motif_summary/motif_summary.html",
        "title": "fp-tools motif summary",
    },
    "pseudobulk-fragments": {
        "sample_id": "pseudobulk_fragments_run",
        "fragments": "data/public/raw/10x_pbmc5k_scatac/atac_pbmc_5k_nextgem_fragments.tsv.gz",
        "annotations": "data/public/processed/pseudobulk_pbmc5k_scatac/pbmc5k_scprinter_broad_annotations.tsv",
        "group_by": "cell_type",
        "outdir": "examples/gui_demo_outputs/pseudobulk_fragments",
        "min_cells": 1,
        "min_fragments": 1,
        "index_output": True,
    },
    "find-signature-fp": {
        "sample_id": "signature_fp_run",
        "fragments": "data/public/raw/10x_pbmc5k_scatac/atac_pbmc_5k_nextgem_fragments.tsv.gz",
        "annotations": "data/public/processed/pseudobulk_pbmc5k_scatac/pbmc5k_scprinter_broad_annotations.tsv",
        "h5ad": "data/public/processed/pseudobulk_pbmc5k_scatac/pbmc5k_scprinter_broad.h5ad",
        "tf_site_dir": "data/public/processed/pseudobulk_pbmc5k_scatac/footprint_demo/tf_sites",
        "outdir": "examples/gui_demo_outputs/signature_fp",
        "markers": "STAT6,FOSB,CEBPA,IRF8,RELA,ZNF683,NR4A1,SMAD3",
        "summary_output_prefix": "single_cell_footprinting",
        "max_motifs": 25,
    },
    "sc-footprinting": {
        "sample_id": "pseudobulk_footprints_run",
        "fragments": "data/public/raw/10x_pbmc5k_scatac/atac_pbmc_5k_nextgem_fragments.tsv.gz",
        "annotations": "data/public/processed/pseudobulk_pbmc5k_scatac/pbmc5k_scprinter_broad_annotations.tsv",
        "h5ad": "data/public/processed/pseudobulk_pbmc5k_scatac/pbmc5k_scprinter_broad.h5ad",
        "group_by": "cell_type",
        "outdir": "examples/gui_demo_outputs/pseudobulk_footprints",
        "genome_sizes": "data/public/processed/pseudobulk_pbmc5k_scatac/hg38.chrom.sizes",
        "genome": "data/public/raw/genome/hg38.fa",
        "peaks": "data/public/raw/10x_pbmc5k_scatac/atac_pbmc_5k_snatac2_selected_bins.demo.bed",
        "motif_db": "jaspar2026_vertebrates",
        "dry_run": True,
    },
}

LIST_TEXT_FIELDS = {
    "signals",
    "bigwigs",
    "motifs",
    "input_html",
    "cond_names",
    "TFBS",
    "read_shift",
}


def main() -> None:
    st.set_page_config(page_title="fp-tools GUI", layout="wide")
    _apply_page_style()
    _ensure_session_config()

    run_dir = Path(st.session_state.gui_run_dir).expanduser()

    page = _current_page_from_query()
    _sync_config_for_page(page)
    _render_sidebar_header(run_dir)
    _render_sidebar_run_dir_controls()
    if not _tutorial_visible() and st.sidebar.button("Show guided tutorial", key="show_tutorial_button", width="stretch"):
        st.session_state.show_tutorial = True
        st.session_state.hide_tutorial = False
    _render_sidebar_nav(page)

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
        textarea {
            background: var(--fp-surface) !important;
            border-color: var(--fp-border) !important;
            border-radius: var(--fp-radius-control) !important;
            box-shadow: none !important;
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
            box-shadow: var(--fp-shadow) !important;
            padding: 0.35rem 0.55rem !important;
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
        .stDownloadButton > button {
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
        .stDownloadButton > button:hover {
            background: var(--fp-accent-hover) !important;
            border-color: var(--fp-accent-hover) !important;
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
        .fp-sidebar-brand {
            background: #101d33;
            color: #ffffff;
            border-radius: 8px;
            padding: 0.76rem 0.82rem;
            margin: 0.03rem 0 0.42rem;
            border: 1px solid #2b3a55;
        }
        .fp-sidebar-brand-title {
            font-size: 1.28rem;
            font-weight: 800;
            letter-spacing: 0;
            line-height: 1.1;
        }
        .fp-sidebar-brand-subtitle {
            margin-top: 0.18rem;
            color: #aebbd0;
            font-size: 0.78rem;
            line-height: 1.28;
        }
        .fp-nav-group {
            margin: 0.72rem 0 0.28rem;
            display: block;
            padding-bottom: 0.18rem;
            color: #93a4b8;
            font-size: 0.72rem;
            font-weight: 800;
            letter-spacing: 0.06em;
            text-transform: uppercase;
            line-height: 1.1;
        }
        .fp-run-dir-pill {
            display: block;
            color: #aebbd0;
            font-size: 0.76rem;
            line-height: 1.35;
            margin: -0.08rem 0 0.32rem;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
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
            font-size: clamp(0.95rem, 0.7vw, 1.08rem);
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
            margin: 0 0 0.55rem;
        }
        .fp-page-heading h1 {
            font-size: clamp(1.42rem, 1.14vw, 1.84rem);
            margin: 0;
        }
        .fp-page-heading p {
            color: var(--fp-text-muted);
            font-size: 0.92rem;
            line-height: 1.35;
            margin: 0.18rem 0 0;
        }
        .fp-run-card {
            background: #ffffff;
            border: 1px solid var(--fp-border-soft);
            border-radius: 8px;
            padding: 0.78rem 0.86rem;
            box-shadow: var(--fp-shadow);
            position: sticky;
            top: 0.75rem;
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
            word-break: break-word;
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
        @media (max-width: 900px) {
            [data-testid="stSidebar"] {
                min-width: 300px !important;
                max-width: 300px !important;
            }
            .fp-card-grid, .fp-example-grid, .fp-summary-strip, .fp-tool-shell, .fp-tutorial-panel ol {
                grid-template-columns: 1fr;
            }
            .fp-run-card {
                position: static;
            }
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def _ensure_session_config() -> None:
    if "current_config" not in st.session_state:
        st.session_state.current_config = make_single_config(
            "atac-correct",
            {"bams": [], "genome": "", "peaks": "", "blacklist": "", "outdir": "", "cores": 1},
            job_id="run",
        )
    if "gui_run_dir" not in st.session_state:
        st.session_state.gui_run_dir = os.environ.get("FP_TOOLS_GUI_RUN_DIR", str(default_gui_run_dir()))


def _tool_pages() -> set[str]:
    return set(PAGE_OPTIONS).difference({"Home", "Run History", "Config"})


def _sync_config_for_page(page: str) -> None:
    if page not in _tool_pages():
        return
    current_tool = _current_config_tool()
    canonical_page = canonical_tool_name(page)
    if current_tool != canonical_page:
        st.session_state.current_config = _default_config_for_tool(canonical_page)


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
        defaults = _prepare_generic_params(tool, _drop_single_meta(raw_defaults))
    return make_single_config(tool, defaults, job_id=f"{tool.replace('-', '_')}_run")


def _current_page_from_query() -> str:
    session_page = st.session_state.get("gui_page")
    if session_page in PAGE_OPTIONS:
        return str(session_page)
    raw_page = st.query_params.get("page", "Home")
    if isinstance(raw_page, list):
        raw_page = raw_page[0] if raw_page else "Home"
    page = str(raw_page)
    if page not in PAGE_OPTIONS:
        page = "Home"
    st.session_state.gui_page = page
    return page


def _render_sidebar_header(run_dir: Path) -> None:
    st.sidebar.markdown(
        f"""
        <div class="fp-sidebar-brand">
          <div class="fp-sidebar-brand-title">fp-tools</div>
          <div class="fp-sidebar-brand-subtitle">Command-first footprinting workflows · v{escape(__version__)}</div>
        </div>
        <span class="fp-run-dir-pill" title="{escape(str(run_dir))}">Run dir: {escape(str(run_dir))}</span>
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
        f"""
        <section class="fp-hero">
          <div class="fp-kicker">fp-tools GUI</div>
          <h1>Run footprint workflows from one clean control surface</h1>
          <p>Choose a command, load an example, review the YAML, and launch the same reproducible workflow used by the command line.</p>
          <div class="fp-pill-row">
            <span class="fp-pill">Public bind ready</span>
            <span class="fp-pill">YAML-first</span>
            <span class="fp-pill">Batch aware</span>
            <span class="fp-pill">Run history</span>
          </div>
        </section>
        <div class="fp-section-title">Typical workflow</div>
        <div class="fp-card-grid">
          <div class="fp-gui-card"><h3>1. Choose command</h3><p>Open bulk ATAC, motif, report, pseudobulk, or signature tools.</p></div>
          <div class="fp-gui-card"><h3>2. Load example</h3><p>Start from bundled YAML or paste your own paths.</p></div>
          <div class="fp-gui-card"><h3>3. Review YAML</h3><p>Check the exact config before running.</p></div>
          <div class="fp-gui-card"><h3>4. Inspect outputs</h3><p>Open logs, tables, bigWigs, and HTML reports.</p></div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    st.caption(
        "Remote access: launch `fp-tools-gui --host 0.0.0.0 --port 8891`, open the printed URL, "
        "and allow that TCP port in the firewall or cloud security group."
    )
    if GUI_EXAMPLE_DIR.exists():
        files = sorted(path.name for path in GUI_EXAMPLE_DIR.glob("*.yml"))
        if files:
            example_items = "\n".join(f'<div class="fp-example-item">{escape(name)}</div>' for name in files[:12])
            st.markdown(
                f"""
                <div class="fp-section-title">Example YAML configs <span>{escape(str(GUI_EXAMPLE_DIR))}</span></div>
                <div class="fp-example-grid">{example_items}</div>
                """,
                unsafe_allow_html=True,
            )
    _render_home_config_snapshot()


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
        <div class="fp-page-heading">
          <h1>{escape(title)}</h1>
          {desc_html}
        </div>
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
        st.caption("No run history yet.")
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
    form_col, run_col = st.columns([0.64, 0.36], gap="large")
    with form_col:
        _render_page_loader("atac-correct")
        mode = st.radio("Mode", ["Single run", "Batch sample list"], horizontal=True, key="at_mode")
        single = _current_single_params("atac-correct")
        batch_rows = _current_sample_rows(
            "atac-correct",
            default_rows=[{"sample_id": "sample1", "bams": "", "genome": "", "peaks": "", "blacklist": "", "outdir": "", "cores": 1}],
        )
        if mode == "Single run":
            with st.form("atacorrect_single_form"):
                left, right = st.columns(2)
                with left:
                    bams_value = single.get("bams", "")
                    if isinstance(bams_value, list):
                        bams_value = "\n".join(str(value) for value in bams_value)
                    bams = st.text_area("BAMs", value=str(bams_value), height=82, help="One BAM path per line")
                    genome = st.text_input("Genome FASTA", value=str(single.get("genome", "")))
                    peaks = st.text_input("Peaks BED", value=str(single.get("peaks", "")))
                with right:
                    blacklist = st.text_input("Blacklist BED", value=str(single.get("blacklist", "")))
                    outdir = st.text_input("Output directory", value=str(single.get("outdir", "")))
                    cores = st.number_input("Cores", min_value=1, value=int(single.get("cores", 1)), step=1)
                submitted = st.form_submit_button("Update page config")
            if submitted:
                _set_config(
                    make_single_config(
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
    form_col, run_col = st.columns([0.64, 0.36], gap="large")
    with form_col:
        _render_page_loader("call-footprints")
        mode = st.radio("Mode", ["Single run", "Batch sample list"], horizontal=True, key="fs_mode")
        single = _current_single_params("call-footprints")
        batch_rows = _current_sample_rows(
            "call-footprints",
            default_rows=[{"sample_id": "sample1", "signal": "", "regions": "", "output": "", "score": "footprint", "cores": 1}],
        )
        if mode == "Single run":
            with st.form("footprintscores_single_form"):
                left, right = st.columns([0.68, 0.32])
                with left:
                    signal = st.text_input("Signal bigWig", value=str(single.get("signal", "")))
                    regions = st.text_input("Regions BED", value=str(single.get("regions", "")))
                    output = st.text_input("Output bigWig", value=str(single.get("output", "")))
                with right:
                    score_values = ["footprint", "sum", "mean", "none"]
                    score_default = str(single.get("score", "footprint"))
                    score = st.selectbox("Score", score_values, index=score_values.index(score_default) if score_default in score_values else 0)
                    cores = st.number_input("Cores", min_value=1, value=int(single.get("cores", 1)), step=1, key="fs_single_cores")
                submitted = st.form_submit_button("Update page config")
            if submitted:
                _set_config(
                    make_single_config(
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
    form_col, run_col = st.columns([0.64, 0.36], gap="large")
    with form_col:
        _render_page_loader("diff-footprints")
        mode = st.radio(
            "Mode",
            ["Single condition", "Batch single-condition list", "Batch comparison list"],
            horizontal=True,
            key="diff_footprints_mode",
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
            with st.form("diff_footprints_single_form"):
                comparison_axis = st.selectbox(
                    "Comparison axis",
                    ["conditions", "regions"],
                    index=0 if str(single.get("comparison_axis", "conditions")) == "conditions" else 1,
                    help="Compare conditions, or compare two or more region sets in the same biological sample(s).",
                )
                left, right = st.columns(2)
                with left:
                    motifs = st.text_input("Motifs", value=str(single.get("motifs", "")))
                    motif_db = st.text_input("Motif database", value=str(single.get("motif_db", "jaspar2026_vertebrates")))
                    genome = st.text_input("Genome FASTA", value=str(single.get("genome", "")))
                    peaks = st.text_input("Peaks BED", value=str(single.get("peaks", "")))
                with right:
                    peak_header = st.text_input("Peak header", value=str(single.get("peak_header", "")))
                    outdir = st.text_input("Output directory", value=str(single.get("outdir", "")))
                    cores = st.number_input("Cores", min_value=1, value=int(single.get("cores", 1)), step=1, key="diff_footprints_single_cores")
                    skip_excel = st.checkbox("Skip Excel", value=bool(single.get("skip_excel", False)))
                signals = st.text_area("Signals", value=_join_multi(single.get("signals", [])), height=94)
                cond_names = st.text_area("Condition names", value=_join_multi(single.get("cond_names", ["Bcell"])), height=76)
                regions = st.text_area(
                    "Region-set BED files",
                    value=_join_multi(single.get("regions", [])),
                    height=76,
                    help="Required only for region comparisons; one BED path per line.",
                )
                region_labels = st.text_input(
                    "Region labels",
                    value=",".join(single.get("region_labels", [])),
                    help="Optional comma-separated labels in the same order as the BED files.",
                )
                region_strata_column = st.number_input(
                    "Matching-stratum BED column (0 = none)",
                    min_value=0,
                    value=int(single.get("region_strata_column", 0) or 0),
                    step=1,
                )
                submitted = st.form_submit_button("Update page config")
            if submitted:
                _set_config(
                    make_single_config(
                        "diff-footprints",
                        {
                            "comparison_axis": comparison_axis,
                            "motifs": motifs,
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
    form_col, run_col = st.columns([0.64, 0.36], gap="large")
    with form_col:
        _render_page_loader("plot-aggregate")
        mode = st.radio("Mode", ["Single run", "Batch sample list"], horizontal=True, key="pa_mode")
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
            with st.form("plotaggregate_single_form"):
                left, right = st.columns(2)
                with left:
                    tfbs = st.text_area("TFBS paths", value=_join_multi(single.get("TFBS", [])), height=92)
                    signals = st.text_area("Signal paths", value=_join_multi(single.get("signals", [])), height=92)
                with right:
                    output = st.text_input("Output PDF", value=str(single.get("output", "")))
                    grid = st.text_input("Grid (optional, e.g. 2x5)", value=str(single.get("grid", "")))
                    score_csv = st.text_input("Aggregated score CSV", value=str(single.get("output_aggregated_scores", "")))
                    signal_csv = st.text_input("Aggregated signal CSV", value=str(single.get("output_aggregated_signals", "")))
                submitted = st.form_submit_button("Update page config")
            if submitted:
                config: dict[str, Any] = {
                    "TFBS": _split_multi(tfbs),
                    "signals": _split_multi(signals),
                    "output": output,
                }
                if grid.strip():
                    config["grid"] = grid.strip()
                if score_csv.strip():
                    config["output_aggregated_scores"] = score_csv.strip()
                if signal_csv.strip():
                    config["output_aggregated_signals"] = signal_csv.strip()
                _set_config(make_single_config("plot-aggregate", config, job_id="plotaggregate_run"))
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
    _render_page_heading(tool, "Configure this command as YAML, then run it here or from the CLI.")
    form_col, run_col = st.columns([0.64, 0.36], gap="large")
    with form_col:
        _render_page_loader(tool)
        defaults = GENERIC_TOOL_DEFAULTS.get(tool, {"sample_id": f"{tool.replace('-', '_')}_run"})
        current = _current_single_params(tool) or defaults
        st.caption("This page writes the same YAML config used by run-yaml-workflow and the direct CLI.")
        with st.form(f"{tool}_generic_form"):
            edited: dict[str, Any] = {}
            simple_fields: list[tuple[str, Any]] = []
            wide_fields: list[tuple[str, Any]] = []
            for key, default_value in current.items():
                if key in {"tool", "sample_id", "job_id", "comparison_id"}:
                    continue
                if key in LIST_TEXT_FIELDS or isinstance(default_value, list):
                    wide_fields.append((key, default_value))
                else:
                    simple_fields.append((key, default_value))
            columns = st.columns(2)
            for idx, (key, default_value) in enumerate(simple_fields):
                with columns[idx % 2]:
                    edited[key] = _render_generic_field(tool, key, default_value)
            for key, default_value in wide_fields:
                edited[key] = _render_generic_field(tool, key, default_value)
            extra_args = st.text_input(
                "Extra CLI args",
                value=str(current.get("extra_args", "") if not isinstance(current.get("extra_args"), list) else " ".join(current.get("extra_args", []))),
            )
            submitted = st.form_submit_button("Update page config")
        if submitted:
            config = _prepare_generic_params(tool, edited)
            if extra_args.strip():
                config["extra_args"] = extra_args.split()
            _set_config(make_single_config(tool, config, job_id=str(current.get("sample_id", f"{tool.replace('-', '_')}_run"))))
    with run_col:
        _render_run_controls(run_dir, label=tool.replace("-", "_"))


def _render_config_page(run_dir: Path) -> None:
    _render_page_heading("Config", "Download, load, edit, save, and run YAML workflow configs.")
    form_col, run_col = st.columns([0.64, 0.36], gap="large")
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
            key="config_yaml_text",
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
                _set_config(normalize_config(load_yaml_config(GUI_EXAMPLE_DIR / example_choice)))
        upload = st.file_uploader("Upload YAML", type=["yml", "yaml"], key=f"{tool}_uploader")
        if upload is not None and st.button("Apply uploaded YAML", key=f"{tool}_apply_upload"):
            _set_config(parse_yaml_text(upload.getvalue().decode("utf-8")))
        path_text = st.text_input("Config path", key=f"{tool}_config_path")
        if st.button("Load YAML from path", key=f"{tool}_load_path") and path_text.strip():
            _set_config(normalize_config(load_yaml_config(path_text.strip())))


def _render_run_controls(run_dir: Path, label: str) -> None:
    normalized = normalize_config(st.session_state.current_config)
    validation_errors = validate_config(normalized)
    tool = _current_config_tool() or "none"
    with st.container(border=True):
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
        if validation_errors:
            st.error("Config needs fixes before launch.")
            for message in validation_errors:
                st.write(f"- {message}")
        else:
            st.success("Config is ready to run.")
        with st.expander("Preview runnable YAML", expanded=False):
            st.code(config_to_yaml_text(normalized), language="yaml")
        if st.button("Start run", key=f"run_{label}", width="stretch"):
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
    if isinstance(value, list):
        value = _join_multi(value)
    if isinstance(value, bool):
        return st.checkbox(_human_label(key), value=value)
    if key in LIST_TEXT_FIELDS:
        return st.text_area(_human_label(key), value=str(value or ""), height=88, key=f"{tool}_{key}")
    if isinstance(value, int):
        return int(st.number_input(_human_label(key), min_value=0, value=int(value), step=1, key=f"{tool}_{key}"))
    return st.text_input(_human_label(key), value=str(value or ""), key=f"{tool}_{key}")


def _human_label(key: str) -> str:
    return key.replace("_", " ").replace("tfbs", "TFBS").title()


def _set_config(config: dict[str, Any]) -> None:
    st.session_state.current_config = normalize_config(config)
    st.success("Current config updated.")


def _data_editor(title: str, default_rows: list[dict[str, Any]], key: str) -> list[dict[str, Any]]:
    st.subheader(title)
    df = pd.DataFrame(default_rows)
    edited = st.data_editor(df, num_rows="dynamic", width="stretch", key=key)
    cleaned = edited.fillna("").to_dict(orient="records")
    return [row for row in cleaned if any(str(value).strip() for value in row.values())]


def _current_config() -> dict[str, Any]:
    return normalize_config(st.session_state.current_config)


def _current_single_params(tool: str) -> dict[str, Any]:
    tool = canonical_tool_name(tool)
    config = _current_config()
    if config["comparisons"]:
        return {}
    matching = [item for item in config["samples"] if canonical_tool_name(str(item.get("tool", ""))) == tool]
    if len(matching) != 1:
        return {}
    return _drop_single_meta(matching[0])


def _current_sample_rows(tool: str, default_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    tool = canonical_tool_name(tool)
    config = _current_config()
    matching = [_prepare_row_for_editor(item) for item in config["samples"] if canonical_tool_name(str(item.get("tool", ""))) == tool]
    return matching or default_rows


def _current_comparison_rows(tool: str, default_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    tool = canonical_tool_name(tool)
    config = _current_config()
    matching = [_prepare_row_for_editor(item) for item in config["comparisons"] if canonical_tool_name(str(item.get("tool", ""))) == tool]
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


def _example_files_for_tool(tool: str) -> list[Path]:
    if not GUI_EXAMPLE_DIR.exists():
        return []
    tool = canonical_tool_name(tool)
    prefixes = {
        "atac-correct": ["atacorrect", "atac_correct"],
        "call-footprints": ["call_footprints", "footprintscores"],
        "diff-footprints": ["diff_footprints"],
        "plot-aggregate": ["plotaggregate", "plot_aggregate"],
    }.get(tool, [tool.replace("-", "_")])
    files: list[Path] = []
    for prefix in prefixes:
        files.extend(GUI_EXAMPLE_DIR.glob(f"{prefix}_*.yml"))
    return sorted(set(files))


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
