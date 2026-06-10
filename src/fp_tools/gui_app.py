"""Streamlit GUI for fp-tools.

This module is an isolated wrapper around the packaged commands. Direct CLI
usage remains primary. The GUI supports direct form-driven runs, YAML load/save,
and batch editing while using the same normalized config model as the optional
``fp-tools-run --config ...`` path.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pandas as pd
import streamlit as st

from fp_tools import __version__
from fp_tools.gui_config import (
    config_to_yaml_text,
    load_yaml_config,
    make_single_config,
    normalize_config,
    parse_yaml_text,
    validate_config,
)
from fp_tools.gui_jobs import default_gui_run_dir, launch_config_async, materialize_run_config, refresh_run_status

GUI_EXAMPLE_DIR = Path("examples/gui_configs")


def main() -> None:
    st.set_page_config(page_title="fp-tools GUI", layout="wide")
    _apply_page_style()
    _ensure_session_config()

    run_dir = Path(st.session_state.gui_run_dir).expanduser()

    st.sidebar.markdown("## fp-tools")
    st.sidebar.caption(f"Version {__version__}")
    st.sidebar.caption(f"Run dir: {run_dir}")
    _render_sidebar_run_dir_controls()
    page = st.sidebar.radio(
        "Navigation",
        ["Home", "Run History", "ATACorrect", "FootprintScores", "BINDetect", "PlotAggregate", "Config"],
        label_visibility="collapsed",
    )

    if page == "Home":
        _render_home(run_dir)
    elif page == "Run History":
        _render_run_history(run_dir)
    elif page == "ATACorrect":
        _render_atacorrect_page(run_dir)
    elif page == "FootprintScores":
        _render_footprintscores_page(run_dir)
    elif page == "BINDetect":
        _render_bindetect_page(run_dir)
    elif page == "PlotAggregate":
        _render_plotaggregate_page(run_dir)
    else:
        _render_config_page(run_dir)


def _apply_page_style() -> None:
    st.markdown(
        """
        <style>
        :root {
            color-scheme: light;
            --fp-bg: #f5f5f7;
            --fp-surface: #ffffff;
            --fp-surface-soft: #fbfbfd;
            --fp-border: #d2d6dc;
            --fp-border-soft: #e5e7eb;
            --fp-text: #111827;
            --fp-text-muted: #4b5563;
            --fp-accent: #111827;
            --fp-accent-hover: #1f2937;
            --fp-hover: #f3f4f6;
            --fp-radius-control: 10px;
            --fp-radius-card: 14px;
            --fp-shadow: 0 1px 3px rgba(15, 23, 42, 0.06);
        }
        html, body, [class*="css"], [data-testid="stAppViewContainer"] {
            font-family: -apple-system, BlinkMacSystemFont, "SF Pro Display", "SF Pro Text", "Helvetica Neue", Arial, sans-serif;
        }
        [data-testid="stAppViewContainer"] {
            background: var(--fp-bg);
            color: var(--fp-text);
        }
        #MainMenu {visibility: hidden;}
        footer {visibility: hidden;}
        header[data-testid="stHeader"] {
            background: rgba(255, 255, 255, 0.88);
            backdrop-filter: blur(10px);
            -webkit-backdrop-filter: blur(10px);
            border-bottom: 1px solid var(--fp-border-soft);
        }
        [data-testid="stSidebar"] {
            background: var(--fp-surface-soft);
            border-right: 1px solid var(--fp-border-soft);
        }
        [data-testid="stSidebar"] .stMarkdown,
        [data-testid="stSidebar"] label,
        [data-testid="stSidebar"] p,
        [data-testid="stSidebar"] span,
        [data-testid="stSidebar"] .stCaption,
        [data-testid="stSidebar"] .stRadio label p {
            color: var(--fp-text) !important;
        }
        [data-testid="stSidebar"] .stCaption {
            color: var(--fp-text-muted) !important;
        }
        [data-testid="stSidebar"] input,
        [data-testid="stSidebar"] textarea,
        [data-testid="stSidebar"] [data-baseweb="input"] input {
            color: var(--fp-text) !important;
            -webkit-text-fill-color: var(--fp-text) !important;
            background: var(--fp-surface) !important;
        }
        [data-testid="stSidebar"] .stRadio label p {
            font-size: 1rem;
            font-weight: 700;
            letter-spacing: 0.01em;
        }
        [data-testid="stSidebar"] .stRadio > div {
            gap: 0.42rem;
        }
        [data-testid="stSidebar"] .stRadio [role="radiogroup"] > label {
            padding: 0.5rem 0.72rem;
            border-radius: var(--fp-radius-control);
            border: 1px solid var(--fp-border-soft);
            background: var(--fp-surface);
            transition: background 120ms ease, border-color 120ms ease;
        }
        [data-testid="stSidebar"] .stRadio [role="radiogroup"] > label:hover {
            background: var(--fp-hover);
            border-color: var(--fp-border);
        }
        [data-testid="stSidebar"] .stButton > button {
            background: var(--fp-surface) !important;
            color: var(--fp-text) !important;
            border: 1px solid var(--fp-border) !important;
            font-weight: 700 !important;
            border-radius: var(--fp-radius-control) !important;
            min-height: 2.5rem !important;
            box-shadow: var(--fp-shadow) !important;
        }
        [data-testid="stSidebar"] .stButton > button:hover {
            background: var(--fp-hover) !important;
            color: var(--fp-text) !important;
            border-color: var(--fp-border) !important;
        }
        [data-testid="stSidebar"] .stButton > button p,
        [data-testid="stSidebar"] .stButton > button span,
        [data-testid="stSidebar"] [data-baseweb="input"] > div,
        [data-testid="stSidebar"] [data-testid="stTextInputRootElement"] > div {
            color: var(--fp-text) !important;
            background: var(--fp-surface) !important;
        }
        [data-testid="stSidebar"] [data-baseweb="input"] {
            background: var(--fp-surface) !important;
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
            font-size: 0.96rem !important;
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
            padding-top: 1.8rem;
            padding-bottom: 2.6rem;
            max-width: 1280px;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def _ensure_session_config() -> None:
    if "current_config" not in st.session_state:
        st.session_state.current_config = make_single_config(
            "ATACorrect",
            {"bam": "", "genome": "", "peaks": "", "blacklist": "", "outdir": "", "cores": 1},
            job_id="run",
        )
    if "gui_run_dir" not in st.session_state:
        st.session_state.gui_run_dir = os.environ.get("FP_TOOLS_GUI_RUN_DIR", str(default_gui_run_dir()))


def _render_home(run_dir: Path) -> None:
    st.title("fp-tools GUI")
    st.write("Per-user browser wrapper around the packaged fp-tools commands.")
    st.write(f"Run directory: `{run_dir}`")
    st.info(
        "Direct CLI remains primary. The GUI builds normalized configs, can save or load YAML, "
        "and can run single jobs or batch jobs without changing the core tool behavior."
    )
    if GUI_EXAMPLE_DIR.exists():
        st.subheader("Example YAML configs")
        st.write(f"Ready-to-load examples are under `{GUI_EXAMPLE_DIR}`.")
        files = sorted(path.name for path in GUI_EXAMPLE_DIR.glob("*.yml"))
        if files:
            st.code("\n".join(files), language="text")
    _show_current_summary()


def _render_run_history(run_dir: Path) -> None:
    st.title("Run History")
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
    st.title("ATACorrect")
    _render_page_loader("ATACorrect")
    mode = st.radio("Mode", ["Single run", "Batch sample list"], horizontal=True, key="at_mode")
    single = _current_single_params("ATACorrect")
    batch_rows = _current_sample_rows("ATACorrect", default_rows=[
        {"sample_id": "sample1", "bam": "", "genome": "", "peaks": "", "blacklist": "", "outdir": "", "cores": 1}
    ])
    if mode == "Single run":
        with st.form("atacorrect_single_form"):
            bam = st.text_input("BAM", value=str(single.get("bam", "")))
            genome = st.text_input("Genome FASTA", value=str(single.get("genome", "")))
            peaks = st.text_input("Peaks BED", value=str(single.get("peaks", "")))
            blacklist = st.text_input("Blacklist BED", value=str(single.get("blacklist", "")))
            outdir = st.text_input("Output directory", value=str(single.get("outdir", "")))
            cores = st.number_input("Cores", min_value=1, value=int(single.get("cores", 1)), step=1)
            submitted = st.form_submit_button("Update page config")
        if submitted:
            _set_config(
                make_single_config(
                    "ATACorrect",
                    {
                        "bam": bam,
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
        rows = _data_editor("ATACorrect sample list", batch_rows, key="atacorrect_batch_editor")
        if st.button("Update page config from sample list", key="atacorrect_batch_set"):
            _set_config(
                {
                    "version": 1,
                    "run_mode": "batch",
                    "defaults": {},
                    "samples": [{"tool": "ATACorrect", **row} for row in rows],
                    "comparisons": [],
                }
            )
    _render_run_controls(run_dir, label="atacorrect")


def _render_footprintscores_page(run_dir: Path) -> None:
    st.title("FootprintScores")
    _render_page_loader("FootprintScores")
    mode = st.radio("Mode", ["Single run", "Batch sample list"], horizontal=True, key="fs_mode")
    single = _current_single_params("FootprintScores")
    batch_rows = _current_sample_rows(
        "FootprintScores",
        default_rows=[{"sample_id": "sample1", "signal": "", "regions": "", "output": "", "score": "footprint", "cores": 1}],
    )
    if mode == "Single run":
        with st.form("footprintscores_single_form"):
            signal = st.text_input("Signal bigWig", value=str(single.get("signal", "")))
            regions = st.text_input("Regions BED", value=str(single.get("regions", "")))
            output = st.text_input("Output bigWig", value=str(single.get("output", "")))
            score_values = ["footprint", "sum", "mean", "none"]
            score_default = str(single.get("score", "footprint"))
            score = st.selectbox("Score", score_values, index=score_values.index(score_default) if score_default in score_values else 0)
            cores = st.number_input("Cores", min_value=1, value=int(single.get("cores", 1)), step=1, key="fs_single_cores")
            submitted = st.form_submit_button("Update page config")
        if submitted:
            _set_config(
                make_single_config(
                    "FootprintScores",
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
        rows = _data_editor("FootprintScores sample list", batch_rows, key="footprintscores_batch_editor")
        if st.button("Update page config from sample list", key="footprintscores_batch_set"):
            _set_config(
                {
                    "version": 1,
                    "run_mode": "batch",
                    "defaults": {},
                    "samples": [{"tool": "FootprintScores", **row} for row in rows],
                    "comparisons": [],
                }
            )
    _render_run_controls(run_dir, label="footprintscores")


def _render_bindetect_page(run_dir: Path) -> None:
    st.title("BINDetect")
    _render_page_loader("BINDetect")
    mode = st.radio(
        "Mode",
        ["Single condition", "Batch single-condition list", "Batch comparison list"],
        horizontal=True,
        key="bindetect_mode",
    )
    single = _current_single_params("BINDetect")
    sample_rows = _current_sample_rows(
        "BINDetect",
        default_rows=[
            {
                "sample_id": "sample1",
                "motifs": "",
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
        "BINDetect",
        default_rows=[
            {
                "comparison_id": "bcell_vs_tcell",
                "motifs": "",
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
        with st.form("bindetect_single_form"):
            motifs = st.text_input("Motifs", value=str(single.get("motifs", "")))
            signals = st.text_area("Signals", value=_join_multi(single.get("signals", [])))
            genome = st.text_input("Genome FASTA", value=str(single.get("genome", "")))
            peaks = st.text_input("Peaks BED", value=str(single.get("peaks", "")))
            peak_header = st.text_input("Peak header", value=str(single.get("peak_header", "")))
            outdir = st.text_input("Output directory", value=str(single.get("outdir", "")))
            cond_names = st.text_area("Condition names", value=_join_multi(single.get("cond_names", ["Bcell"])))
            cores = st.number_input("Cores", min_value=1, value=int(single.get("cores", 1)), step=1, key="bindetect_single_cores")
            skip_excel = st.checkbox("Skip Excel", value=bool(single.get("skip_excel", False)))
            submitted = st.form_submit_button("Update page config")
        if submitted:
            _set_config(
                make_single_config(
                    "BINDetect",
                    {
                        "motifs": motifs,
                        "signals": _split_multi(signals),
                        "genome": genome,
                        "peaks": peaks,
                        "peak_header": peak_header,
                        "outdir": outdir,
                        "cond_names": _split_multi(cond_names),
                        "cores": int(cores),
                        "skip_excel": bool(skip_excel),
                    },
                    job_id="bindetect_single",
                )
            )
    elif mode == "Batch single-condition list":
        rows = _data_editor("BINDetect single-condition sample list", sample_rows, key="bindetect_sample_editor")
        if st.button("Update page config from single-condition list", key="bindetect_sample_set"):
            _set_config(
                {
                    "version": 1,
                    "run_mode": "batch",
                    "defaults": {},
                    "samples": [
                        {
                            "tool": "BINDetect",
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
        rows = _data_editor("BINDetect comparison list", comparison_rows, key="bindetect_comparison_editor")
        if st.button("Update page config from comparison list", key="bindetect_comparison_set"):
            _set_config(
                {
                    "version": 1,
                    "run_mode": "batch",
                    "defaults": {},
                    "samples": [],
                    "comparisons": [
                        {
                            "tool": "BINDetect",
                            **row,
                            "signals": _split_multi(str(row.get("signals", ""))),
                            "cond_names": _split_multi(str(row.get("cond_names", ""))),
                            "skip_excel": _as_bool(row.get("skip_excel", False)),
                        }
                        for row in rows
                    ],
                }
            )
    _render_run_controls(run_dir, label="bindetect")


def _render_plotaggregate_page(run_dir: Path) -> None:
    st.title("PlotAggregate")
    _render_page_loader("PlotAggregate")
    mode = st.radio("Mode", ["Single run", "Batch sample list"], horizontal=True, key="pa_mode")
    single = _current_single_params("PlotAggregate")
    batch_rows = _current_sample_rows(
        "PlotAggregate",
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
            tfbs = st.text_area("TFBS paths", value=_join_multi(single.get("TFBS", [])))
            signals = st.text_area("Signal paths", value=_join_multi(single.get("signals", [])))
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
            _set_config(make_single_config("PlotAggregate", config, job_id="plotaggregate_run"))
    else:
        rows = _data_editor("PlotAggregate panel list", batch_rows, key="plotaggregate_batch_editor")
        if st.button("Update page config from panel list", key="plotaggregate_batch_set"):
            _set_config(
                {
                    "version": 1,
                    "run_mode": "batch",
                    "defaults": {},
                    "samples": [
                        {
                            "tool": "PlotAggregate",
                            **row,
                            "TFBS": _split_multi(str(row.get("TFBS", ""))),
                            "signals": _split_multi(str(row.get("signals", ""))),
                        }
                        for row in rows
                    ],
                    "comparisons": [],
                }
            )
    _render_run_controls(run_dir, label="plotaggregate")


def _render_config_page(run_dir: Path) -> None:
    st.title("Config")
    st.download_button(
        "Download current YAML",
        data=config_to_yaml_text(st.session_state.current_config),
        file_name="fp_tools_config.yml",
        mime="text/yaml",
        width="stretch",
    )
    uploader = st.file_uploader("Load YAML file", type=["yml", "yaml"], key="config_uploader")
    if uploader is not None and st.button("Apply uploaded YAML", key="config_apply_upload"):
        _set_config(parse_yaml_text(uploader.getvalue().decode("utf-8")))

    load_path = st.text_input("Load config from path", key="config_load_path")
    if st.button("Load YAML from path", key="config_load_path_btn") and load_path.strip():
        _set_config(normalize_config(load_yaml_config(load_path.strip())))

    yaml_text = st.text_area(
        "Current YAML",
        value=config_to_yaml_text(st.session_state.current_config),
        height=460,
        key="config_yaml_text",
    )
    col1, col2 = st.columns(2)
    with col1:
        if st.button("Apply YAML text", key="config_apply_text"):
            _set_config(parse_yaml_text(yaml_text))
    with col2:
        save_path = st.text_input("Save current YAML to path", key="config_save_path")
        if st.button("Save YAML", key="config_save_btn") and save_path.strip():
            path = Path(save_path.strip()).expanduser()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(yaml_text, encoding="utf-8")
            st.success(f"Saved config to {path}")
    _render_run_controls(run_dir, label="config")


def _render_page_loader(tool: str) -> None:
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
    st.subheader("Run")
    normalized = normalize_config(st.session_state.current_config)
    validation_errors = validate_config(normalized)
    st.code(config_to_yaml_text(normalized), language="yaml")
    if validation_errors:
        st.error("Config needs fixes before launch.")
        for message in validation_errors:
            st.write(f"- {message}")
    if st.button("Start run", key=f"run_{label}"):
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
    st.sidebar.markdown("### Run directory")
    run_dir_input = st.sidebar.text_input("GUI run dir", value=str(st.session_state.gui_run_dir), key="sidebar_run_dir")
    if st.sidebar.button("Apply run dir", key="sidebar_apply_run_dir", width="stretch"):
        path = Path(run_dir_input).expanduser()
        path.mkdir(parents=True, exist_ok=True)
        st.session_state.gui_run_dir = str(path)
        st.sidebar.success(f"Run dir set to {path}")


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
    tool = str(item.get("tool", ""))
    base = Path.cwd()

    def add_path(value: Any) -> None:
        text = str(value).strip()
        if not text:
            return
        path = Path(text).expanduser()
        if not path.is_absolute():
            path = (base / path).resolve()
        outputs.append(path)

    if tool == "ATACorrect":
        outdir = str(item.get("outdir", "")).strip()
        if outdir:
            outdir_path = Path(outdir).expanduser()
            if not outdir_path.is_absolute():
                outdir_path = (base / outdir_path).resolve()
            outputs.append(outdir_path)
            if outdir_path.exists():
                outputs.extend(sorted(path for path in outdir_path.iterdir() if path.is_file()))
    elif tool == "FootprintScores":
        add_path(item.get("output", ""))
    elif tool == "PlotAggregate":
        for key in ("output", "output_aggregated_scores", "output_aggregated_signals", "output_csv"):
            add_path(item.get(key, ""))
    elif tool == "BINDetect":
        outdir = str(item.get("outdir", "")).strip()
        if outdir:
            outdir_path = Path(outdir).expanduser()
            if not outdir_path.is_absolute():
                outdir_path = (base / outdir_path).resolve()
            outputs.append(outdir_path)
            if outdir_path.exists():
                preferred = [
                    "bindetect_results.txt",
                    "bindetect_figures.pdf",
                    "bindetect_clusters.pdf",
                    "bindetect_results_skewness_report.pdf",
                ]
                outputs.extend(outdir_path / name for name in preferred if (outdir_path / name).exists())
                outputs.extend(sorted(path for path in outdir_path.glob("bindetect_*.html")))

    seen: set[str] = set()
    deduped: list[Path] = []
    for path in outputs:
        key = str(path)
        if key not in seen:
            seen.add(key)
            deduped.append(path)
    return deduped


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
    config = _current_config()
    if config["comparisons"]:
        return {}
    matching = [item for item in config["samples"] if str(item.get("tool")) == tool]
    if len(matching) != 1:
        return {}
    return _drop_single_meta(matching[0])


def _current_sample_rows(tool: str, default_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    config = _current_config()
    matching = [_prepare_row_for_editor(item) for item in config["samples"] if str(item.get("tool")) == tool]
    return matching or default_rows


def _current_comparison_rows(tool: str, default_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    config = _current_config()
    matching = [_prepare_row_for_editor(item) for item in config["comparisons"] if str(item.get("tool")) == tool]
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
    prefix = tool.lower().replace("footprintscores", "footprintscores").replace("plotaggregate", "plotaggregate")
    return sorted(path for path in GUI_EXAMPLE_DIR.glob(f"{prefix}_*.yml"))


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
