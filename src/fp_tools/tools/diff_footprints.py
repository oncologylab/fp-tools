#!/usr/bin/env python
"""
Differential-footprint command driver for motif-associated footprint analysis.

This module handles:
- motif/scored-signal integration
- bound versus unbound site calling
- differential binding statistics between conditions
- summary tables, PDFs, and interactive HTML outputs

It also includes replicate grouping support and optional skewness report integration.
"""

import os
import sys
import argparse
import base64
import gzip
import json
import time
import glob
import random
import shutil
import copy
import tempfile
import numpy as np
import multiprocessing as mp
import itertools
import pandas as pd
import seaborn as sns
from collections import Counter
from types import SimpleNamespace
import warnings

# ML / stats
import sklearn
from sklearn import mixture
import scipy
from kneed import KneeLocator  # noqa: F401

# Plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Bio
import pysam
import pyBigWig as pybw

# Internal (fp_tools namespace)
from fp_tools.parsers import add_diff_footprints_arguments
from fp_tools.tools.diff_footprint_helpers import *
from fp_tools.tools import diff_footprint_skew_report as skewrep
from fp_tools.tools.diff_footprint_replicate_report import build_replicate_report

from fp_tools.utils.utilities import (
    check_required, check_files, make_directory, merge_dicts,
    monitor_progress, expand_dirs, check_cores, file_writer
)
from fp_tools.utils.regions import *
from fp_tools.utils.motifs import *
from fp_tools.utils.motif_databases import motif_db_table, resolve_motif_inputs
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.project_layout import (
    analysis_peaks_path,
    comparison_dir,
    corrected_bigwig_path,
    footprint_bigwig_path,
    is_project_layout,
    match_motifs_dir,
    normalized_bigwig_path,
    project_root,
    read_comparison_table,
    read_sample_table,
    sample_dir,
    samples_for_condition,
    samples_root,
)
from fp_tools.utils.plotting_style import PDF_FONT_SIZE, apply_pdf_style, apply_ascii_minus_to_figure

# tame some noisy warnings during curve fitting
from scipy.optimize import OptimizeWarning


def _benjamini_hochberg(pvalues):
    """Return BH-adjusted q-values for a 1D array-like of p-values."""

    pvals = np.asarray(pvalues, dtype=float)
    qvals = np.full(pvals.shape, np.nan, dtype=float)
    finite = np.isfinite(pvals)
    if not finite.any():
        return qvals
    clipped = np.clip(pvals[finite], 0.0, 1.0)
    order = np.argsort(clipped)
    ranked = clipped[order]
    n = float(len(ranked))
    adjusted = ranked * n / np.arange(1, len(ranked) + 1, dtype=float)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    unsorted = np.empty_like(adjusted)
    unsorted[order] = adjusted
    qvals[finite] = unsorted
    return qvals

warnings.simplefilter("ignore", OptimizeWarning)
warnings.simplefilter("ignore", RuntimeWarning)
apply_pdf_style()


def _signal_sample_stems(signals):
    """Return file-stem labels for signal paths, with duplicate stems disambiguated."""

    raw_names = [os.path.basename(os.path.splitext(str(bw))[0]) for bw in signals]
    counts = Counter()
    names = []
    for name in raw_names:
        counts[name] += 1
        names.append(name if counts[name] == 1 else f"{name}_{counts[name]}")
    return names


def _resolve_sample_names(args, default_names):
    provided = getattr(args, "sample_names", None)
    if provided is None:
        return list(default_names)
    sample_names = list(provided)
    if len(sample_names) != len(args.signals):
        raise ValueError("--sample-names must have the same length as --signals")
    duplicates = sorted(name for name, count in Counter(sample_names).items() if count > 1)
    if duplicates:
        raise ValueError("--sample-names must be unique; duplicate labels: " + ", ".join(duplicates))
    return sample_names


def _disambiguate_sample_condition_collisions(sample_names, condition_names, match_only=False):
    """Avoid duplicate sample-level and condition-level score columns."""

    if match_only:
        return list(sample_names)
    condition_set = set(condition_names or [])
    out = []
    for sample in sample_names:
        out.append(f"{sample}_sample" if sample in condition_set else sample)
    duplicates = sorted(name for name, count in Counter(out).items() if count > 1)
    if duplicates:
        raise ValueError(
            "Sample names collide with condition-derived internal labels; "
            "please provide distinct --sample-names. Duplicates: " + ", ".join(duplicates)
        )
    return out


def _find_one_file(root, patterns, label):
    matches = []
    for pattern in patterns:
        matches.extend(sorted(root.glob(pattern)))
    files = [p for p in matches if p.is_file()]
    if len(files) != 1:
        found = ", ".join(str(p) for p in files[:5]) or "none"
        raise ValueError(f"Expected exactly one {label} in {root}; found {len(files)} ({found})")
    return files[0]


def _find_optional_file(root, patterns):
    matches = []
    for pattern in patterns:
        matches.extend(sorted(root.glob(pattern)))
    files = [p for p in matches if p.is_file()]
    if len(files) > 1:
        found = ", ".join(str(p) for p in files[:5])
        raise ValueError(f"Expected at most one file in {root}; found {len(files)} ({found})")
    return files[0] if files else None


def _find_preferred_optional_file(root, patterns):
    for pattern in patterns:
        files = [p for p in sorted(root.glob(pattern)) if p.is_file()]
        if len(files) > 1:
            found = ", ".join(str(p) for p in files[:5])
            raise ValueError(f"Expected at most one file for pattern {pattern!r} in {root}; found {len(files)} ({found})")
        if files:
            return files[0]
    return None


def _sample_name_from_folder(sample_dir):
    manifest = sample_dir / "fp_tools_sample.json"
    if manifest.exists():
        try:
            import json
            data = json.loads(manifest.read_text(encoding="utf-8"))
            for key in ("sample_name", "sample", "name"):
                value = str(data.get(key, "")).strip()
                if value:
                    return value
        except Exception:
            pass
    return sample_dir.name


def _resolve_folder_inputs(args):
    """Resolve --sample-dirs/--project-dir into signal paths and match caches."""

    from pathlib import Path

    dirs = []
    if getattr(args, "sample_dirs", None):
        dirs.extend(Path(p).resolve() for p in args.sample_dirs)
    if getattr(args, "project_dir", None):
        project = Path(args.project_dir).resolve()
        discovered = []
        for child in sorted(project.iterdir()):
            if not child.is_dir():
                continue
            if (child / "match_motifs").is_dir() or (child / "motif_matches").is_dir():
                discovered.append(child)
        dirs.extend(discovered)
    if not dirs:
        return args

    if getattr(args, "signals", None):
        raise ValueError("--signals cannot be combined with --sample-dirs or --project-dir")
    seen = set()
    unique_dirs = []
    for directory in dirs:
        if directory in seen:
            continue
        seen.add(directory)
        unique_dirs.append(directory)
    signals = []
    match_dirs = []
    default_names = []
    aggregate_signals = []
    for sample_dir in unique_dirs:
        if not sample_dir.exists():
            raise FileNotFoundError(f"Sample directory does not exist: {sample_dir}")
        footprint = _find_optional_file(
            sample_dir,
            ["*_footprints.bw", "*_footprint.bw", "footprints.bw", "footprint.bw"],
        )
        match_dir = sample_dir / "match_motifs"
        if not match_dir.is_dir():
            match_dir = sample_dir / "motif_matches"
        if not match_dir.is_dir():
            raise FileNotFoundError(f"Sample directory lacks match_motifs/ output: {sample_dir}")
        sample_name = _sample_name_from_folder(sample_dir)
        signals.append(str(footprint) if footprint is not None else f"cached:{sample_name}")
        match_dirs.append(str(match_dir))
        default_names.append(sample_name)
        aggregate_signal = _find_preferred_optional_file(
            sample_dir,
            [
                "normalize/*_corrected_q95_scaled.bw",
                "atac_correct/*_corrected.bw",
                "*_corrected_q95_scaled.bw",
                "*_corrected.bw",
            ],
        )
        if aggregate_signal:
            aggregate_signals.append(str(aggregate_signal))
    args.signals = signals
    args.cached_match_dirs = match_dirs
    args.cached_without_bigwigs = any(str(signal).startswith("cached:") for signal in signals)
    args.folder_default_sample_names = default_names
    if not getattr(args, "aggregate_signals", None) and len(aggregate_signals) == len(signals):
        args.aggregate_signals = aggregate_signals
    return args


def _resolve_motif_arguments(args):
    try:
        args.motifs = resolve_motif_inputs(getattr(args, "motifs", None), getattr(args, "motif_db", None))
    except ValueError as exc:
        sys.exit(f"ERROR: {exc}")
    return args.motifs


def norm_fit(x, mean, std, scale):
    return scale * scipy.stats.norm.pdf(x, mean, std)


def _prepare_condition_metadata(args):
    """Derive condition, replicate, and comparison metadata from CLI arguments."""

    default_sample_names = list(getattr(args, "folder_default_sample_names", None) or _signal_sample_stems(args.signals))
    sample_names = _resolve_sample_names(args, default_sample_names)
    default_condition_names = default_sample_names
    if getattr(args, "match_only", False) and getattr(args, "sample_names", None) is not None:
        default_condition_names = sample_names
    args.cond_names = (
        list(default_condition_names)
        if args.cond_names is None else args.cond_names
    )
    args.outdir = os.path.abspath(args.outdir)

    orig = list(args.cond_names)
    if len(orig) != len(args.signals):
        raise ValueError("--cond-names must have the same length as --signals")
    if getattr(args, "norm_off", False):
        args.normalization = "none"

    idxs = {}
    for i, nm in enumerate(orig):
        idxs.setdefault(nm, []).append(i)
    args.cond_groups = idxs
    args.cond_names = list(idxs.keys())
    sample_names = _disambiguate_sample_condition_collisions(
        sample_names,
        args.cond_names,
        match_only=getattr(args, "match_only", False),
    )
    args.condition_replicates = {cond: len(indices) for cond, indices in idxs.items()}
    args.signal_sample_names = list(sample_names)
    args.sample_names = list(sample_names)
    args.sample_to_condition = {}
    args.condition_samples = {cond: [] for cond in args.cond_names}
    for cond, indices in idxs.items():
        for signal_idx in indices:
            sample_name = sample_names[signal_idx]
            args.sample_to_condition[sample_name] = cond
            args.condition_samples[cond].append(sample_name)

    if getattr(args, "match_only", False):
        args.comparisons = []
    elif args.time_series:
        args.comparisons = list(zip(args.cond_names[:-1], args.cond_names[1:]))
    else:
        args.comparisons = list(itertools.combinations(args.cond_names, 2))
    return args


def _existing_result_motifs(info_table, comparison, args, motif_lookup=None):
    """Build lightweight motif records from an existing diff-footprints result table."""

    c1, c2 = comparison
    base = f"{c1}_{c2}"
    required = ["output_prefix", "name", "motif_id", base + "_change", base + "_pvalue"]
    missing = [column for column in required if column not in info_table.columns]
    if missing:
        raise ValueError(
            "Existing results table is missing required column(s) for "
            f"{base}: {', '.join(missing)}"
        )

    rows = info_table.copy()
    rows[base + "_change_numeric"] = pd.to_numeric(rows[base + "_change"], errors="coerce").fillna(0.0)
    rows[base + "_pvalue_numeric"] = pd.to_numeric(rows[base + "_pvalue"], errors="coerce").fillna(1.0)
    qvalue_col = base + "_qvalue_bh"
    if qvalue_col in rows.columns:
        rows[base + "_qvalue_numeric"] = pd.to_numeric(rows[qvalue_col], errors="coerce").fillna(1.0)
    else:
        rows[base + "_qvalue_numeric"] = _benjamini_hochberg(rows[base + "_pvalue_numeric"].to_numpy())
    filtered_p = rows.loc[rows[base + "_pvalue_numeric"] > 0, base + "_pvalue_numeric"]
    pval_min = np.percentile(filtered_p, 5) if len(filtered_p) else 1.0
    change_min, change_max = np.percentile(rows[base + "_change_numeric"], [5, 95]) if len(rows) else (0.0, 0.0)

    motifs = []
    for _, row in rows.iterrows():
        prefix = str(row["output_prefix"])
        change = float(row[base + "_change_numeric"])
        pvalue = float(row[base + "_pvalue_numeric"])
        qvalue = float(row[base + "_qvalue_numeric"])
        highlighted_col = base + "_highlighted"
        if highlighted_col in row and not pd.isna(row[highlighted_col]):
            highlighted = str(row[highlighted_col]).strip().lower() in {"true", "1", "yes"}
        else:
            highlighted = (change < change_min) or (change > change_max) or (pvalue < pval_min)
        if highlighted:
            group = f"{c2}_up" if change < 0 else f"{c1}_up"
        else:
            group = "n.s."
        logo_path = os.path.join(args.outdir, prefix, prefix + ".png")
        logo = ""
        if os.path.exists(logo_path):
            with open(logo_path, "rb") as handle:
                logo = base64.b64encode(handle.read()).decode("utf-8")
        matched_motif = (motif_lookup or {}).get(prefix)
        motifs.append(
            SimpleNamespace(
                prefix=prefix,
                name=str(row.get("name", prefix)),
                id=str(row.get("motif_id", "")),
                change=change,
                pvalue=pvalue,
                qvalue=qvalue,
                logpvalue=-np.log10(max(pvalue, 1e-308)),
                highlighted=highlighted,
                group=group,
                base=logo,
                counts=getattr(matched_motif, "counts", None),
            )
        )
    return motifs


def _load_reuse_motif_lookup(args, logger):
    """Load motif matrices for vector logos during report-only regeneration."""

    motif_lookup = {}
    motif_paths = getattr(args, "motifs", None) or []
    if not motif_paths:
        return motif_lookup
    try:
        motif_list = MotifList()
        for f in expand_dirs(motif_paths):
            motif_list += MotifList().from_file(f)
        for motif in motif_list:
            motif.set_prefix(args.naming)
        motif_prefixes = [m.prefix.upper() for m in motif_list]
        name_count = Counter(motif_prefixes)
        duplicated = [k for k, v in name_count.items() if v > 1]
        motif_count = {dup: 1 for dup in duplicated}
        for motif in motif_list:
            if motif.prefix.upper() in duplicated:
                original = motif.prefix
                motif.prefix = f"{motif.prefix}_{motif_count[motif.prefix.upper()]}"
                motif_count[original.upper()] += 1
            motif_lookup[motif.prefix] = motif
    except Exception as exc:
        logger.warning(f"Could not load motif matrices for SVG logos in reuse mode: {exc}")
    return motif_lookup


def run_diff_footprints_reuse_existing_results(args):
    """Regenerate final diff-footprints reports from completed motif-level outputs."""

    _resolve_motif_arguments(args)
    _prepare_condition_metadata(args)
    logger = FpToolsLogger("diff-footprints", args.verbosity)
    logger.begin()
    parser = add_diff_footprints_arguments(argparse.ArgumentParser())
    args.cores = check_cores(args.cores, logger)
    logger.arguments_overview(parser, args)

    results_path = os.path.join(args.outdir, args.prefix + "_results.txt")
    if not os.path.exists(results_path):
        raise FileNotFoundError(
            f"--reuse-existing-results requires an existing results table: {results_path}"
        )
    logger.info(f"Reusing existing diff-footprints results table: {results_path}")
    info_table = pd.read_csv(results_path, sep="\t")

    if len(args.cond_names) < 2:
        raise ValueError("--reuse-existing-results requires at least two unique conditions")

    write_replicate_report = args.replicate_report == "on" or (
        args.replicate_report == "auto" and any(count > 1 for count in args.condition_replicates.values())
    )
    if write_replicate_report:
        report_out = args.replicate_report_out or os.path.join(args.outdir, args.prefix + "_replicate_report.tsv")
        summary_out = args.replicate_summary_out or os.path.join(args.outdir, args.prefix + "_replicate_summary.tsv")
        figure_out = args.replicate_figure_out or os.path.join(args.outdir, args.prefix + "_replicate_report.png")
        try:
            build_replicate_report(
                results_path,
                report_out,
                summary_output=summary_out,
                figure_output=figure_out,
                replicate_map=args.replicate_map,
            )
            logger.info(f"Regenerated replicate-aware report from existing results: {report_out}")
        except Exception as exc:
            logger.warning(f"Could not regenerate replicate-aware report: {exc}")

    motif_lookup = _load_reuse_motif_lookup(args, logger)

    for comparison in args.comparisons:
        c1, c2 = comparison
        base = f"{c1}_{c2}"
        logger.info(f"Regenerating interactive diff-footprints report for {c1} / {c2}")
        motifs = _existing_result_motifs(info_table, comparison, args, motif_lookup=motif_lookup)
        aggregate_data = None
        if getattr(args, "aggregate_signals", None) and getattr(args, "plot_aggregate", "off") != "off":
            try:
                aggregate_data = build_diff_footprint_aggregate_payload(motifs, info_table, comparison, args)
                if aggregate_data is None or len(aggregate_data.get("motifs", [])) == 0:
                    logger.warning(f"No reusable aggregate profiles were found for {base}; writing volcano-only HTML")
            except Exception as exc:
                logger.warning(f"Could not build aggregate payload from existing motif BEDs: {exc}")
        html_out = os.path.join(args.outdir, args.prefix + "_" + base + ".html")
        plot_interactive_diff_footprints(
            motifs,
            [c1, c2],
            html_out,
            aggregate_data=aggregate_data,
            title="Differential footprint report",
            report_label=getattr(args, "report_label", None),
        )
        logger.info(f"Wrote {html_out}")

    logger.info("Reuse mode skipped motif scanning, per-motif processing, static PDFs, and clustering.")
    logger.end()


def _read_cached_all_bed(path, peak_cols, sample_col_count=1):
    """Read a cached *_all.bed file without a header."""

    rows = []
    expected = 6 + peak_cols + sample_col_count
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < expected:
                raise ValueError(f"Cached motif BED has {len(parts)} columns but expected at least {expected}: {path}")
            rows.append(parts)
    return rows


def _cached_motif_bed_map(match_dir):
    """Map motif prefix to cached all.bed path in a match-motifs output directory."""

    paths = {}
    for bed in sorted(os.path.abspath(p) for p in glob.glob(os.path.join(match_dir, "*", "beds", "*_all.bed"))):
        prefix = os.path.basename(bed)[:-len("_all.bed")]
        paths[prefix] = bed
    if not paths:
        raise FileNotFoundError(f"No cached motif BEDs found under {match_dir}")
    return paths


def _cached_motif_logo_map(match_dir):
    """Map motif prefix to cached logo PNG path in a match-motifs output directory."""

    paths = {}
    for logo in sorted(os.path.abspath(p) for p in glob.glob(os.path.join(match_dir, "*", "*.png"))):
        prefix = os.path.splitext(os.path.basename(logo))[0]
        paths[prefix] = logo
    return paths


def _cached_result_table(match_dir):
    candidates = sorted(glob.glob(os.path.join(match_dir, "*_results.txt")))
    if not candidates:
        return None
    return candidates[0]


def _cached_distance_table(match_dir):
    candidates = sorted(glob.glob(os.path.join(match_dir, "*_distances.txt")))
    if not candidates:
        return None
    return candidates[0]


def _cached_cluster_map(match_dir):
    result_table = _cached_result_table(match_dir)
    if result_table is None:
        return {}
    try:
        table = pd.read_csv(result_table, sep="\t", usecols=["output_prefix", "cluster"])
    except Exception:
        return {}
    return dict(zip(table["output_prefix"].astype(str), table["cluster"].astype(str)))


def _copy_cached_logos(args, motif_list, logo_filenames, logger):
    cached_dirs = list(getattr(args, "cached_match_dirs", []) or [])
    if not cached_dirs:
        return 0
    logo_maps = [_cached_motif_logo_map(path) for path in cached_dirs]
    copied = 0
    for motif in motif_list:
        src = None
        for mapping in logo_maps:
            if motif.prefix in mapping:
                src = mapping[motif.prefix]
                break
        if src is None:
            continue
        dst = logo_filenames[motif.prefix]
        make_directory(os.path.dirname(dst))
        shutil.copyfile(src, dst)
        copied += 1
    if copied:
        logger.info(f"Reused {copied} cached motif logo(s) from match-motifs sample folders")
    return copied


def _write_cached_tfbs_tmp_files(args, motif_names, logger):
    """Merge per-sample match-motifs caches into the .tmp files process_tfbs expects."""

    cached_dirs = list(getattr(args, "cached_match_dirs", []) or [])
    if len(cached_dirs) != len(args.sample_names):
        raise ValueError("--sample-dirs/--project-dir cache count must match resolved samples")

    peak_cols = len(args.peak_header_list)
    build_overlap_clusters = bool(getattr(args, "static_plots", False))
    if not build_overlap_clusters:
        args.cached_cluster_map = _cached_cluster_map(cached_dirs[0])
        args.cached_distance_table = _cached_distance_table(cached_dirs[0])

    if not build_overlap_clusters and not getattr(args, "write_motif_outputs", True):
        maps = [_cached_motif_bed_map(path) for path in cached_dirs]
        common = set(maps[0])
        for mapping in maps[1:]:
            common &= set(mapping)
        missing = [name for name in motif_names if name not in common]
        if missing:
            raise ValueError("Cached match-motifs outputs are missing motif(s): " + ", ".join(missing[:10]))
        args.cached_motif_bed_maps = maps
        logger.info(
            f"Reusing cached motif-site scores from {len(cached_dirs)} sample folder(s) "
            "with parallel per-motif reads"
        )
        return {}

    compact_paths = [_match_motif_site_cache_path(path) for path in cached_dirs]
    try:
        compact_overlaps = _write_cached_tfbs_tmp_files_from_compact(
            args, motif_names, compact_paths, peak_cols, build_overlap_clusters
        )
    except Exception as exc:
        logger.warning(f"Compact motif-site cache could not be reused ({exc}); falling back to per-motif BED files")
        compact_overlaps = None
    if compact_overlaps is not None:
        logger.info(f"Reused cached motif-site scores from {len(cached_dirs)} sample folder(s) using compact motif-site cache")
        return compact_overlaps

    maps = [_cached_motif_bed_map(path) for path in cached_dirs]
    common = set(maps[0])
    for mapping in maps[1:]:
        common &= set(mapping)
    missing = [name for name in motif_names if name not in common]
    if missing:
        raise ValueError("Cached match-motifs outputs are missing motif(s): " + ", ".join(missing[:10]))

    global_tfbs = RegionList() if build_overlap_clusters else None
    for motif in motif_names:
        per_sample_rows = [_read_cached_all_bed(mapping[motif], peak_cols, sample_col_count=1) for mapping in maps]
        _write_merged_cached_rows_for_motif(args, motif, peak_cols, per_sample_rows, global_tfbs)
    logger.info(f"Reused cached motif-site scores from {len(cached_dirs)} sample folder(s) using per-motif BED files")
    if build_overlap_clusters:
        return global_tfbs.count_overlaps() if len(global_tfbs) else {}
    return {}


def _process_tfbs_from_cached_beds(TF_name, args, log2fc_params):
    """Build one cached motif temp file inside a worker, then run normal summary logic."""

    peak_cols = len(args.peak_header_list)
    maps = getattr(args, "cached_motif_bed_maps", None)
    if not maps:
        return process_tfbs(TF_name, args, log2fc_params)
    per_sample_rows = [_read_cached_all_bed(mapping[TF_name], peak_cols, sample_col_count=1) for mapping in maps]
    _write_merged_cached_rows_for_motif(args, TF_name, peak_cols, per_sample_rows, None)
    return process_tfbs(TF_name, args, log2fc_params)


def _match_cache_paths(match_dir):
    cache_dir = os.path.join(match_dir, "cache")
    return (
        os.path.join(cache_dir, "background_scores.tsv.gz"),
        os.path.join(cache_dir, "manifest.json"),
    )


def _match_motif_site_cache_path(match_dir):
    return os.path.join(match_dir, "cache", "motif_sites.tsv.gz")


def _write_match_motifs_cache(args, background, logger):
    """Persist background scores needed by folder-based diff-footprints reuse."""

    cache_tsv, manifest_json = _match_cache_paths(args.outdir)
    make_directory(os.path.dirname(cache_tsv))
    keys = background.get("keys") or []
    sample_names = list(args.sample_names)
    lengths = [len(keys)] + [len(background["sample_signal"].get(sample, [])) for sample in sample_names]
    if len(set(lengths)) != 1:
        raise ValueError(f"Cannot write match-motifs cache; background row counts differ: {lengths}")
    with gzip.open(cache_tsv, "wt", encoding="utf-8", compresslevel=1) as handle:
        handle.write("\t".join(["chrom", "start", "end", "offset"] + sample_names) + "\n")
        for row_idx, key in enumerate(keys):
            scores = [f"{float(background['sample_signal'][sample][row_idx]):.8g}" for sample in sample_names]
            handle.write("\t".join(list(key) + scores) + "\n")
    manifest = {
        "format": "fp-tools match-motifs background cache",
        "version": 1,
        "sample_names": sample_names,
        "condition_names": list(args.cond_names),
        "peak_header": list(getattr(args, "peak_header_list", [])),
        "background_rows": len(keys),
        "normalization": getattr(args, "normalization", "none"),
    }
    with open(manifest_json, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    logger.info(f"Wrote match-motifs reuse cache: {cache_tsv}")


def _write_match_motif_site_cache(args, motif_names, logger):
    """Persist per-motif site scores for fast folder-based diff-footprints reuse."""

    cache_tsv = _match_motif_site_cache_path(args.outdir)
    make_directory(os.path.dirname(cache_tsv))
    peak_cols = len(args.peak_header_list)
    expected = 6 + peak_cols + 1
    written = 0
    with gzip.open(cache_tsv, "wt", encoding="utf-8", compresslevel=1) as handle:
        header = [
            "motif",
            "TFBS_chr",
            "TFBS_start",
            "TFBS_end",
            "TFBS_name",
            "TFBS_score",
            "TFBS_strand",
        ] + list(args.peak_header_list) + ["score"]
        handle.write("\t".join(header) + "\n")
        for motif in motif_names:
            path = os.path.join(args.outdir, motif, "beds", motif + "_all.bed")
            if not os.path.exists(path):
                continue
            with open(path, "r", encoding="utf-8") as bed_handle:
                for line in bed_handle:
                    if not line.strip():
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < expected:
                        raise ValueError(f"Motif BED has {len(parts)} columns but expected at least {expected}: {path}")
                    handle.write(motif + "\t" + "\t".join(parts[:expected]) + "\n")
                    written += 1
    logger.info(f"Wrote compact motif-site reuse cache with {written} rows: {cache_tsv}")


def _load_match_motif_site_cache(match_dir):
    cache_tsv = _match_motif_site_cache_path(match_dir)
    if not os.path.exists(cache_tsv):
        return None
    rows_by_motif = {}
    with gzip.open(cache_tsv, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if len(header) < 8 or header[0] != "motif":
            raise ValueError(f"Unexpected motif-site cache header in {cache_tsv}")
        for line_no, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                raise ValueError(f"Motif-site cache line {line_no} has {len(parts)} columns but expected {len(header)}: {cache_tsv}")
            motif = parts[0]
            rows_by_motif.setdefault(motif, []).append(parts[1:])
    return rows_by_motif


def _read_next_motif_site_cache_group(handle, expected_cols):
    line = handle.readline()
    if not line:
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) < expected_cols + 1:
        raise ValueError("Motif-site cache row has fewer columns than expected")
    motif = parts[0]
    rows = [parts[1:expected_cols + 1]]
    while True:
        pos = handle.tell()
        line = handle.readline()
        if not line:
            break
        parts = line.rstrip("\n").split("\t")
        if len(parts) < expected_cols + 1:
            raise ValueError("Motif-site cache row has fewer columns than expected")
        if parts[0] != motif:
            handle.seek(pos)
            break
        rows.append(parts[1:expected_cols + 1])
    return motif, rows


def _write_cached_tfbs_tmp_files_from_compact(args, motif_names, cache_paths, peak_cols, build_overlap_clusters):
    expected_cols = 6 + peak_cols + 1
    motif_order = {motif: index for index, motif in enumerate(motif_names)}
    handles = []
    groups = []
    global_tfbs = RegionList() if build_overlap_clusters else None
    try:
        for path in cache_paths:
            if not os.path.exists(path):
                return None
            handle = gzip.open(path, "rt", encoding="utf-8")
            header = handle.readline().rstrip("\n").split("\t")
            if len(header) < expected_cols + 1 or header[0] != "motif":
                return None
            handles.append(handle)
            groups.append(_read_next_motif_site_cache_group(handle, expected_cols))

        for motif_index, motif in enumerate(motif_names):
            per_sample_rows = []
            for sample_idx, group in enumerate(groups):
                if group is None:
                    per_sample_rows.append([])
                    continue
                group_motif, rows = group
                group_index = motif_order.get(group_motif)
                if group_index is None:
                    return None
                if group_index < motif_index:
                    return None
                if group_motif == motif:
                    per_sample_rows.append(rows)
                    groups[sample_idx] = _read_next_motif_site_cache_group(handles[sample_idx], expected_cols)
                else:
                    per_sample_rows.append([])
            _write_merged_cached_rows_for_motif(args, motif, peak_cols, per_sample_rows, global_tfbs)
        if build_overlap_clusters:
            return global_tfbs.count_overlaps() if len(global_tfbs) else {}
        return {}
    finally:
        for handle in handles:
            handle.close()


def _write_merged_cached_rows_for_motif(args, motif, peak_cols, per_sample_rows, global_tfbs=None):
    lengths = {len(rows) for rows in per_sample_rows}
    if len(lengths) != 1:
        raise ValueError(f"Cached motif site count mismatch for {motif}: {sorted(lengths)}")
    out_rows = []
    for row_idx in range(len(per_sample_rows[0])):
        key = per_sample_rows[0][row_idx][:6 + peak_cols]
        scores = []
        for sample_idx, rows in enumerate(per_sample_rows):
            other_key = rows[row_idx][:6 + peak_cols]
            if other_key != key:
                raise ValueError(
                    f"Cached motif site mismatch for {motif} at row {row_idx + 1} "
                    f"between sample 1 and sample {sample_idx + 1}"
                )
            scores.append(rows[row_idx][6 + peak_cols])
        out_rows.append(key + scores)
    tmp_root = getattr(args, "tmp_tfbs_root", None) or args.outdir
    tmp_path = os.path.join(tmp_root, motif, "beds", motif + ".tmp")
    make_directory(os.path.dirname(tmp_path))
    with open(tmp_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join("\t".join(row) for row in out_rows))
        if out_rows:
            handle.write("\n")
    if global_tfbs is not None and out_rows:
        global_tfbs.extend(RegionList().from_bed(tmp_path))


def _aggregate_site_maps_from_cached_match_dirs(args):
    cached_dirs = list(getattr(args, "cached_match_dirs", []) or [])
    if not cached_dirs:
        return None
    maps = []
    for sample_name, match_dir in zip(args.sample_names, cached_dirs):
        sample_map = {}
        for bed in sorted(glob.glob(os.path.join(match_dir, "*", "beds", "*_all.bed"))):
            prefix = os.path.basename(bed)[:-len("_all.bed")]
            sample_map.setdefault(prefix, {})["all"] = os.path.abspath(bed)
        for bed in sorted(glob.glob(os.path.join(match_dir, "*", "beds", f"*_{sample_name}_bound.bed"))):
            prefix = os.path.basename(bed)[:-len(f"_{sample_name}_bound.bed")]
            sample_map.setdefault(prefix, {})["bound"] = os.path.abspath(bed)
        maps.append(sample_map)
    return maps


def _load_match_background_cache(match_dir):
    cache_tsv, manifest_json = _match_cache_paths(match_dir)
    if not os.path.exists(cache_tsv):
        raise FileNotFoundError(
            f"Missing match-motifs background cache: {cache_tsv}. "
            "Re-run match-motifs with the current fp-tools version before using --sample-dirs."
        )
    manifest = {}
    if os.path.exists(manifest_json):
        with open(manifest_json, "r", encoding="utf-8") as handle:
            manifest = json.load(handle)
    df = pd.read_csv(cache_tsv, sep="\t")
    if df.shape[1] < 5:
        raise ValueError(f"Background cache has no sample-score columns: {cache_tsv}")
    key_cols = ["chrom", "start", "end", "offset"]
    if list(df.columns[:4]) != key_cols:
        raise ValueError(f"Unexpected background cache header in {cache_tsv}")
    score_cols = list(df.columns[4:])
    return df[key_cols].astype(str), {col: df[col].to_numpy(dtype=float) for col in score_cols}, manifest


def _background_from_signal_chunk(regions, args):
    """Collect background scores from signals using the same random-position rule as scanning."""

    sample_bigwigs = {sample: pybw.open(args.signals[idx], "rb") for idx, sample in enumerate(args.sample_names)}
    background = {
        "keys": [],
        "gc": [],
        "signal": {c: [] for c in args.cond_names},
        "sample_signal": {s: [] for s in args.sample_names},
    }
    logger = FpToolsLogger("", args.verbosity, getattr(args, "log_q", None))
    rand_window = 200
    for region in regions:
        reglen = region.get_length()
        random.seed(reglen)
        rand_positions = random.sample(range(reglen), max(1, int(reglen / rand_window)))
        for pos in rand_positions:
            background["keys"].append([region.chrom, str(region.start), str(region.end), str(pos)])
        sample_footprints = {}
        for sample_name in args.sample_names:
            bw = sample_bigwigs[sample_name]
            arr = region.get_signal(bw, logger=logger, key=sample_name)
            sample_footprints[sample_name] = arr
            for pos in rand_positions:
                background["sample_signal"][sample_name].append(arr[pos])
        for condition in args.cond_names:
            rep_signals = [sample_footprints[sample_name] for sample_name in args.condition_samples[condition]]
            stacked = np.vstack(rep_signals)
            values = np.mean(stacked, axis=0)
            for pos in rand_positions:
                background["signal"][condition].append(values[pos])
    for bw in sample_bigwigs.values():
        bw.close()
    return background


def _collect_cached_background(peak_chunks, args, pool, worker_cores):
    if getattr(args, "cached_without_bigwigs", False):
        cached_dirs = list(getattr(args, "cached_match_dirs", []) or [])
        if len(cached_dirs) != len(args.sample_names):
            raise ValueError("--sample-dirs/--project-dir cache count must match resolved samples")
        first_keys = None
        sample_arrays = {}
        for match_dir, sample_name in zip(cached_dirs, args.sample_names):
            keys, score_map, _manifest = _load_match_background_cache(match_dir)
            if first_keys is None:
                first_keys = keys
            elif not keys.equals(first_keys):
                raise ValueError(f"Background cache coordinates differ between sample folders; first mismatch: {match_dir}")
            if sample_name in score_map:
                sample_arrays[sample_name] = score_map[sample_name]
            elif len(score_map) == 1:
                sample_arrays[sample_name] = next(iter(score_map.values()))
            else:
                raise ValueError(
                    f"Background cache for {match_dir} has multiple samples but no column named {sample_name}"
                )
        background = {
            "keys": first_keys.astype(str).values.tolist() if first_keys is not None else [],
            "gc": [],
            "sample_signal": sample_arrays,
            "signal": {},
        }
        for condition in args.cond_names:
            stacked = np.vstack([sample_arrays[sample_name] for sample_name in args.condition_samples[condition]])
            background["signal"][condition] = np.mean(stacked, axis=0)
        return background
    if worker_cores == 1:
        results = [_background_from_signal_chunk(chunk, args) for chunk in peak_chunks]
    else:
        tasks = [pool.apply_async(_background_from_signal_chunk, (chunk, args)) for chunk in peak_chunks]
        results = [task.get() for task in tasks]
    return merge_dicts(results)



# ----------------------------------------------------------------------------- #
def run_diff_footprints(args):
    """Run the differential-footprint pipeline from parsed CLI arguments."""
    if getattr(args, "reuse_existing_results", False):
        return run_diff_footprints_reuse_existing_results(args)

    _resolve_folder_inputs(args)
    _resolve_motif_arguments(args)
    check_required(args, ["genome", "peaks"])
    if not getattr(args, "cached_match_dirs", None):
        check_required(args, ["signals"])
    _prepare_condition_metadata(args)

    motif_output_mode = getattr(args, "motif_outputs", "auto")
    args.write_motif_outputs = (
        getattr(args, "match_only", False)
        or motif_output_mode == "full"
        or (motif_output_mode == "auto" and getattr(args, "plot_aggregate", "off") != "off")
    )

    # outputs we’ll create
    outfiles = [
        os.path.abspath(os.path.join(args.outdir, args.prefix + "_distances.txt")),
        os.path.abspath(os.path.join(args.outdir, args.prefix + "_results.txt")),
    ]
    if not getattr(args, "skip_excel", False):
        outfiles.append(os.path.abspath(os.path.join(args.outdir, args.prefix + "_results.xlsx")))
    if args.write_motif_outputs:
        states = ["bound", "unbound"]
        outfiles += [os.path.abspath(os.path.join(
            args.outdir, "*", "beds", f"*_{cond}_{state}.bed"))
            for (cond, state) in itertools.product(args.cond_names, states)]
        outfiles += [
            os.path.abspath(os.path.join(args.outdir, "*", "beds", "*_all.bed")),
            os.path.abspath(os.path.join(args.outdir, "*", "*_overview.txt")),
        ]
        if not getattr(args, "skip_excel", False):
            outfiles.append(os.path.abspath(os.path.join(args.outdir, "*", "*_overview.xlsx")))
    if getattr(args, "static_plots", False):
        outfiles += [
            os.path.abspath(os.path.join(args.outdir, args.prefix + "_figures.pdf")),
            os.path.abspath(os.path.join(args.outdir, args.prefix + "_clusters.pdf")),
        ]
    if getattr(args, "per_motif_plots", False):
        outfiles.append(os.path.abspath(os.path.join(args.outdir, "*", "plots", "*_log2fcs.pdf")))
    if getattr(args, "skew_report", False):
        outfiles.append(os.path.abspath(os.path.join(args.outdir, args.prefix + "_results_skewness_report.pdf")))

    # ------------------------------ logger/pools ------------------------------ #
    logger = FpToolsLogger("diff-footprints", args.verbosity)
    logger.begin()
    parser = add_diff_footprints_arguments(argparse.ArgumentParser())
    logger.arguments_overview(parser, args)
    logger.output_files(outfiles)

    args.cores = check_cores(args.cores, logger)
    writer_cores = max(1, int(args.cores * 0.1))
    worker_cores = max(1, args.cores - writer_cores)
    logger.debug(f"Worker cores: {worker_cores}")
    logger.debug(f"Writer cores: {writer_cores}")

    pool = mp.Pool(processes=worker_cores)
    writer_pool = mp.Pool(processes=writer_cores)

    # ------------------------------ inputs ----------------------------------- #
    logger.info("----- Processing input data -----")
    logger.info("Checking reading/writing of files")
    files_to_check = [args.motifs, args.genome, args.peaks]
    if not getattr(args, "cached_without_bigwigs", False):
        files_to_check.insert(0, args.signals)
    check_files(files_to_check, action="r")
    check_files([path for path in outfiles if "*" not in path], action="w")
    make_directory(args.outdir)

    # condition comparisons
    no_conditions = len(args.cond_names)  # NOTE: use unique conditions (not #signals)
    comparisons = args.comparisons

    # debug/fig PDFs
    if args.debug:
        debug_out = os.path.join(args.outdir, args.prefix + "_debug.pdf")
        debug_pdf = PdfPages(debug_out, keep_empty=True)

    figure_pdf = None
    cluster_pdf = None
    if getattr(args, "static_plots", False):
        fig_out = os.path.join(args.outdir, args.prefix + "_figures.pdf")
        figure_pdf = PdfPages(fig_out, keep_empty=True)
        cluster_out = os.path.join(args.outdir, args.prefix + "_clusters.pdf")
        cluster_pdf = PdfPages(cluster_out, keep_empty=True)

        plt.figure()
        plt.axis('off')
        plt.text(0.5, 0.8, "DIFF-FOOTPRINTS FIGURES", ha="center", va="center", fontsize=PDF_FONT_SIZE, fontweight="bold")
        titles = ["Raw score distributions"]
        if no_conditions > 1 and not args.norm_off:
            titles.append("Normalized score distributions")
        if args.debug:
            for (c1, c2) in comparisons:
                titles.append(f"Background log2FCs ({c1} / {c2})")
        for (c1, c2) in comparisons:
            titles.append(f"diff-footprints volcano plot ({c1} / {c2})")
        plt.text(0.1, 0.6, "\n".join([f"Page {i+2}) {t}" for i, t in enumerate(titles)]) + "\n\n", va="top", fontsize=PDF_FONT_SIZE, fontweight="bold")
        apply_ascii_minus_to_figure(plt.gcf())
        figure_pdf.savefig(bbox_inches='tight')
        plt.close()

        plt.figure()
        plt.axis('off')
        plt.text(0.5, 0.8, "DIFF-FOOTPRINTS CLUSTERS", ha="center", va="center", fontsize=PDF_FONT_SIZE, fontweight="bold")
        cluster_titles = [f"Cluster overview ({c1} / {c2})" for (c1, c2) in comparisons]
        plt.text(0.1, 0.6, "\n".join([f"Page {i+2}) {t}" for i, t in enumerate(cluster_titles)]) + "\n\n", va="top", fontsize=PDF_FONT_SIZE, fontweight="bold")
        apply_ascii_minus_to_figure(plt.gcf())
        cluster_pdf.savefig(bbox_inches='tight')
        plt.close()

    # ------------------------------ peaks ------------------------------------ #
    logger.info("Reading peaks")
    peaks = RegionList().from_bed(args.peaks)
    logger.info(f"- Found {len(peaks)} regions in input peaks")

    n_cols = len(peaks[0])
    for i, peak in enumerate(peaks):
        if len(peak) != n_cols:
            logger.error(
                f"The lines in --peaks have a varying number of columns. "
                f"Line 1 has {n_cols}, but line {i+1} has {len(peak)}."
            )
            sys.exit(1)

    peaks = peaks.merge()
    logger.info(f"- Merged to {len(peaks)} regions")
    if len(peaks) == 0:
        logger.error("Input --peaks file is empty!")
        sys.exit(1)

    peak_columns = len(peaks[0])
    logger.debug(f"--peaks have {peak_columns} columns")
    if args.peak_header is not None:
        content = open(args.peak_header, "r").read()
        args.peak_header_list = content.split()
        logger.debug(f"Peak header: {args.peak_header_list}")
        if len(args.peak_header_list) != peak_columns:
            logger.error(
                f"Length of --peak_header ({len(args.peak_header_list)}) "
                f"does not fit number of columns in --peaks ({peak_columns})."
            )
            sys.exit(1)
    else:
        args.peak_header_list = (
            ["peak_chr", "peak_start", "peak_end"] +
            [f"additional_{num+1}" for num in range(peak_columns - 3)]
        )
    logger.debug(f"Peak header list: {args.peak_header_list}")

    # boundaries vs fasta / signals
    logger.info("Checking for match between --peaks and --fasta/--signals boundaries")
    logger.info(f"- Comparing peaks to {args.genome}")
    fasta_obj = pysam.FastaFile(args.genome)
    fasta_boundaries = dict(zip(fasta_obj.references, fasta_obj.lengths))
    fasta_obj.close()
    logger.debug(f"Fasta boundaries: {fasta_boundaries}")
    peaks = peaks.apply_method(OneRegion.check_boundary, fasta_boundaries, "exit")

    if not getattr(args, "cached_without_bigwigs", False):
        for signal in args.signals:
            logger.info(f"- Comparing peaks to {signal}")
            pybw_obj = pybw.open(signal)
            pybw_header = pybw_obj.chroms()
            pybw_obj.close()
            logger.debug(f"Signal boundaries: {pybw_header}")
            peaks = peaks.apply_method(OneRegion.check_boundary, pybw_header, "exit")

    # GC content (for motif background)
    logger.info("Estimating GC content from peak sequences")
    peak_chunks = peaks.chunks(args.split)
    gc_content_pool = pool.starmap(get_gc_content, itertools.product(peak_chunks, [args.genome]))
    gc_content = np.mean(gc_content_pool)
    args.gc = gc_content
    bg = np.array([(1-args.gc)/2.0, args.gc/2.0, args.gc/2.0, (1-args.gc)/2.0])
    logger.info(f"- GC content estimated at {gc_content*100:.2f}%")

    # ------------------------------ motifs ----------------------------------- #
    logger.info("Reading motifs")
    motif_list = MotifList()
    args.motifs = expand_dirs(args.motifs)
    for f in args.motifs:
        try:
            motif_list += MotifList().from_file(f)
        except Exception as e:
            logger.error(f"Error reading motifs from '{f}'. Error: {e}")
            sys.exit(1)

    no_pfms = len(motif_list)
    logger.info(f"- Read {no_pfms} motifs")

    logger.debug("Getting motifs ready")
    motif_list.bg = bg
    for motif in motif_list:
        motif.set_prefix(args.naming)
        motif.bg = bg
        logger.spam(f"Getting pssm for motif {motif.name}")
        motif.get_pssm()

    # ensure output prefixes unique (case-insensitive)
    motif_prefixes = [m.prefix.upper() for m in motif_list]
    name_count = Counter(motif_prefixes)
    if max(name_count.values()) > 1:
        duplicated = [k for k, v in name_count.items() if v > 1]
        logger.warning("The motif output names (from --naming) are not unique.")
        logger.warning(f"These names occur >1 time: {duplicated}")
        logger.warning("They will be renamed with '_1', '_2', ...")
        motif_count = {dup: 1 for dup in duplicated}
        for i, m in enumerate(motif_list):
            if m.prefix.upper() in duplicated:
                original = m.prefix
                m.prefix = f"{m.prefix}_{motif_count[m.prefix.upper()]}"
                logger.debug(f"Renamed motif {i+1}: {original} -> {m.prefix}")
                motif_count[original.upper()] += 1

    motif_names = [m.prefix for m in motif_list]

    logger.debug("Getting match threshold per motif")
    outlist = pool.starmap(OneMotif.get_threshold, itertools.product(motif_list, [args.motif_pvalue]))
    motif_list = MotifList(outlist)
    for m in motif_list:
        logger.debug(f"Motif {m.name}: threshold {m.threshold}")

    cached_summary_mode = bool(getattr(args, "cached_match_dirs", None)) and not getattr(args, "write_motif_outputs", True)
    temp_tfbs_dir = None
    if cached_summary_mode:
        temp_tfbs_dir = tempfile.mkdtemp(prefix="fp_tools_diff_tfbs_")
        args.tmp_tfbs_root = temp_tfbs_dir
        args.aggregate_site_maps = _aggregate_site_maps_from_cached_match_dirs(args)
        logger.info("Using temporary motif-site files for cached summary mode")
    else:
        logger.info("Creating folder structure for each TF")
        for TF in motif_names:
            make_directory(os.path.join(args.outdir, TF))
            make_directory(os.path.join(args.outdir, TF, "beds"))
            make_directory(os.path.join(args.outdir, TF, "plots"))

        # logos
        logo_filenames = {m.prefix: os.path.join(args.outdir, m.prefix, m.prefix + ".png") for m in motif_list}
        copied_logos = _copy_cached_logos(args, motif_list, logo_filenames, logger) if getattr(args, "cached_match_dirs", None) else 0
        missing_logo_motifs = [m for m in motif_list if not os.path.exists(logo_filenames[m.prefix])]
        if missing_logo_motifs:
            logger.info("Plotting sequence logos for each motif" if not copied_logos else "Plotting missing sequence logos")
            task_list = [pool.apply_async(OneMotif.logo_to_file, (m, logo_filenames[m.prefix],)) for m in missing_logo_motifs]
            monitor_progress(task_list, logger)
            _ = [t.get() for t in task_list]
            logger.comment("")

        logger.debug("Getting base64 strings per motif")
        for m in motif_list:
            with open(logo_filenames[m.prefix], "rb") as png:
                m.base = base64.b64encode(png.read()).decode("utf-8")

    # --------------------- scan/cache motifs + match to signals --------------- #
    logger.comment("")
    logger.start_logger_queue()
    args.log_q = logger.queue
    cached_mode = bool(getattr(args, "cached_match_dirs", None))
    if cached_mode:
        logger.info("Reusing cached motif-site scores from match-motifs sample folders")
        writer_pool.terminate()
        writer_pool.join()
        TF_overlaps = _write_cached_tfbs_tmp_files(args, motif_names, logger)
        logger.info("Collecting cached background scores")
        background = _collect_cached_background(peak_chunks, args, pool, worker_cores)
    else:
        manager = mp.Manager()
        logger.info("Scanning for motifs and matching to signals...")

        # bed writer queues (one or more writers)
        logger.debug("Setting up writer queues")
        qs_list, writer_qs = [], {}
        TF_names_chunks = [motif_names[i::writer_cores] for i in range(writer_cores)]
        writer_tasks = []
        for TF_sub in TF_names_chunks:
            logger.debug(f"Creating writer queue for {TF_sub}")
            files = [os.path.join(args.outdir, TF, "beds", TF + ".tmp") for TF in TF_sub]
            q = manager.Queue()
            qs_list.append(q)
            writer_tasks.append(writer_pool.apply_async(file_writer, args=(q, dict(zip(TF_sub, files)), args)))
            for TF in TF_sub:
                writer_qs[TF] = q
        writer_pool.close()  # no more writer jobs

        # scan in parallel
        results = []
        if worker_cores == 1:
            logger.debug("Running with cores = 1")
            for chunk in peak_chunks:
                results.append(scan_and_score(chunk, motif_list, args, args.log_q, writer_qs))
        else:
            logger.debug("Sending jobs to worker pool")
            tlist = [pool.apply_async(scan_and_score, (chunk, motif_list, args, args.log_q, writer_qs))
                     for chunk in peak_chunks]
            monitor_progress(tlist, logger)
            results = [t.get() for t in tlist]

        logger.info("Done scanning for TFBS across regions!")
        logger.info("Waiting for bedfiles to write")

        # stop writer queues
        logger.debug("Stop all queues by inserting None")
        for q in qs_list:
            q.put((None, None))

        # wait for writers to complete
        finished = 0
        while finished == 0:
            logger.debug(f"Writer task return status: {[t.get() if t.ready() else 'NA' for t in writer_tasks]}")
            if sum([t.ready() for t in writer_tasks]) == len(writer_tasks):
                finished = 1
                return_codes = [t.get() for t in writer_tasks]
                if sum(return_codes) != 0:
                    logger.error("Bedfile writer finished with an error")
                else:
                    logger.debug("Bedfile writer(s) finished!")
            time.sleep(0.5)

        logger.debug("Joining bed_writer queues")
        for i, q in enumerate(qs_list):
            logger.debug(f"- Queue {i} (size {q.qsize()})")
        writer_pool.join()

        # ---------------------- background + normalization ----------------------- #
        logger.info("Merging results from subsets")
        background = merge_dicts([r[0] for r in results])
        TF_overlaps = merge_dicts([r[1] for r in results])
        if getattr(args, "match_only", False):
            _write_match_motifs_cache(args, background, logger)
        results = None

    # fill possible missing overlap keys
    for TF1 in motif_list:
        if TF1.prefix not in TF_overlaps:
            TF_overlaps[TF1.prefix] = 0
        for TF2 in motif_list:
            tup = (TF1.prefix, TF2.prefix)
            if tup not in TF_overlaps:
                TF_overlaps[tup] = 0

    for cond in args.cond_names:
        background["signal"][cond] = np.array(background["signal"][cond], dtype=float)
    for sample_name in args.sample_names:
        background["sample_signal"][sample_name] = np.array(background["sample_signal"][sample_name], dtype=float)

    n_bg_values = len(background["signal"][args.cond_names[0]])
    logger.debug(f"Collected {n_bg_values} background values")
    if n_bg_values < 1000:
        logger.warning(
            "Low number of background values (<1000). Bound/unbound threshold and "
            "cross-condition normalization may be unstable. Prefer the full union peak set."
        )

    # raw score distributions
    fig = plot_score_distribution([background["signal"][c] for c in args.cond_names],
                                  labels=args.cond_names, title="Raw scores per condition")
    if figure_pdf is not None:
        apply_ascii_minus_to_figure(fig)
        figure_pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # normalization
    args.norm_objects = {}
    if args.normalization == "none" or len(args.cond_names) == 1:
        for cond in args.cond_names:
            args.norm_objects[cond] = ArrayNorm("constant", popt=1.0, value_min=0, value_max=1)
        for sample_name in args.sample_names:
            args.norm_objects[sample_name] = ArrayNorm("constant", popt=1.0, value_min=0, value_max=1)
    elif args.normalization == "sample-quantile":
        logger.comment("")
        logger.info("Normalizing scores across input samples")
        lists = [background["sample_signal"][s] for s in args.sample_names]
        args.norm_objects = quantile_normalization(lists, args.sample_names, pdfpages=debug_pdf if args.debug else None, logger=logger)
        for sample_name in args.sample_names:
            original = background["sample_signal"][sample_name]
            normalized = args.norm_objects[sample_name].normalize(original)
            normalized[normalized < 0] = 0
            background["sample_signal"][sample_name] = normalized
        for cond in args.cond_names:
            stacked = np.vstack([background["sample_signal"][sample] for sample in args.condition_samples[cond]])
            background["signal"][cond] = np.mean(stacked, axis=0)
        fig = plot_score_distribution([background["signal"][c] for c in args.cond_names],
                                      labels=args.cond_names, title="Sample-quantile normalized scores per condition")
        if figure_pdf is not None:
            apply_ascii_minus_to_figure(fig)
            figure_pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
    else:
        logger.comment("")
        logger.info("Normalizing scores across conditions")
        lists = [background["signal"][c] for c in args.cond_names]
        args.norm_objects = quantile_normalization(lists, args.cond_names, pdfpages=debug_pdf if args.debug else None, logger=logger)

        for cond in args.cond_names:
            original = background["signal"][cond]
            logger.debug(f"Background nans ({cond}): {np.isnan(original).sum()}")
            normalized = args.norm_objects[cond].normalize(original)
            normalized[normalized < 0] = 0
            background["signal"][cond] = normalized
            logger.debug(f"Background nans after norm ({cond}): {np.isnan(normalized).sum()}")

        fig = plot_score_distribution([background["signal"][c] for c in args.cond_names],
                                      labels=args.cond_names, title="Condition-quantile normalized scores per condition")
        if figure_pdf is not None:
            apply_ascii_minus_to_figure(fig)
            figure_pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    # ---------------------- threshold (bound/unbound) ------------------------ #
    logger.info("Estimating bound/unbound threshold")
    bg_values = np.array([background["signal"][c] for c in args.cond_names]).flatten()
    logger.debug(f"Size of background array collected: {bg_values.size}")
    bg_values = bg_values[np.isfinite(bg_values)]
    bg_values = bg_values[~np.isclose(bg_values, 0.0)]
    logger.debug(f"Size after filtering > 0: {bg_values.size}")
    if len(bg_values) == 0:
        logger.error("All background scores are zero. Check inputs.")
        sys.exit(1)

    x_max = np.percentile(bg_values, [99])
    bg_values = bg_values[bg_values < x_max]
    logger.debug(f"Size after filtering < x_max ({x_max}): {bg_values.size}")

    log_vals = np.log(bg_values).reshape(-1, 1)
    gmm = sklearn.mixture.GaussianMixture(n_components=2, random_state=1).fit(log_vals)
    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_).flatten()
    chosen_i = np.argmax(means)
    log_params = scipy.stats.lognorm.fit(bg_values, f0=stds[chosen_i], fscale=np.exp(means[chosen_i]))

    mode = scipy.optimize.fmin(lambda x: -scipy.stats.lognorm.pdf(x, *log_params), 0, disp=False)[0]
    logger.debug(f"- Mode estimated at: {mode}")
    args.pseudo = mode / 2.0
    logger.debug(f"Pseudocount estimated at: {args.pseudo:.5f}")

    leftside_x = np.linspace(scipy.stats.lognorm(*log_params).ppf([0.01]), mode, 100)
    leftside_pdf = scipy.stats.lognorm.pdf(leftside_x, *log_params)
    leftside_x_scale = leftside_x - np.min(leftside_x)
    mirrored_x = np.concatenate([leftside_x, np.max(leftside_x) + leftside_x_scale]).flatten()
    mirrored_pdf = np.concatenate([leftside_pdf, leftside_pdf[::-1]]).flatten()
    popt, _ = scipy.optimize.curve_fit(
        lambda x, std, sc: sc * scipy.stats.norm.pdf(x, mode, std),
        mirrored_x, mirrored_pdf
    )
    norm_params = (mode, popt[0])
    threshold = round(scipy.stats.norm.ppf(1 - args.bound_pvalue, *norm_params), 5)
    args.thresholds = {c: threshold for c in args.cond_names}
    logger.stats(f"- Threshold estimated at: {threshold}")

    # ------------------ background log2fc for comparisons -------------------- #
    logger.comment("")
    log2fc_params = {}
    if len(args.cond_names) > 1:
        logger.info("Calculating background log2 fold-changes between conditions")
        for (c1, c2) in comparisons:
            logger.info(f"- {c1} / {c2}")
            s1 = np.copy(background["signal"][c1])
            s2 = np.copy(background["signal"][c2])
            included = np.logical_or(s1 > 0, s2 > 0)
            s1, s2 = s1[included], s2[included]
            log2fcs = np.log2((s1 + args.pseudo) / (s2 + args.pseudo))
            lower, upper = np.percentile(log2fcs, [1, 99])
            fit_vals = log2fcs[(log2fcs >= lower) & (log2fcs <= upper)]
            diff_dist = scipy.stats.norm
            normp = diff_dist.fit(fit_vals)
            logger.debug(f"({c1} / {c2}) Background log2fc distribution: {normp}")
            log2fc_params[(c1, c2)] = normp

            if args.debug:
                fig, ax = plt.subplots(1, 1)
                plt.hist(log2fcs, density=True, bins='auto', label=f"Background log2fc ({c1} / {c2})")
                xvals = np.linspace(plt.xlim()[0], plt.xlim()[1], 100)
                pdf = diff_dist.pdf(xvals, *normp)
                plt.plot(xvals, pdf, label="Distribution fit")
                plt.title(f"Background log2FCs ({c1} / {c2})")
                plt.xlabel("Log2 fold change"); plt.ylabel("Density")
                apply_ascii_minus_to_figure(fig)
                debug_pdf.savefig(fig, bbox_inches='tight'); plt.close()

    background = None  # free mem

    # ------------------ per-TF processing (bound/unbound, stats) ------------- #
    logger.comment("")
    logger.info("Processing scanned TFBS individually")

    info_columns = ["total_tfbs"]
    info_columns += [f"{cond}_{metric}" for cond, metric in itertools.product(args.cond_names, ["threshold", "bound", "n_replicates", "score_sd"])]
    info_columns += [f"{c1}_{c2}_{metric}" for (c1, c2), metric in itertools.product(comparisons, ["change", "pvalue", "mean_delta_fp", "mean_log2fc", "delta_fp_se", "log2fc_se"])]
    info_table = pd.DataFrame(np.zeros((len(motif_names), len(info_columns))),
                              columns=info_columns, index=motif_names)

    results = []
    process_func = _process_tfbs_from_cached_beds if getattr(args, "cached_motif_bed_maps", None) else process_tfbs
    if args.cores == 1:
        for name in motif_names:
            logger.info(f"- {name}")
            results.append(process_func(name, args, log2fc_params))
    else:
        tlist = [pool.apply_async(process_func, (name, args, log2fc_params)) for name in motif_names]
        monitor_progress(tlist, logger)
        results = [t.get() for t in tlist]

    logger.info("Concatenating results from subsets")
    info_table = pd.concat(results)

    pool.terminate()
    pool.join()
    logger.stop_logger_queue()

    matrix_out = os.path.join(args.outdir, args.prefix + "_distances.txt")
    cached_cluster_map = getattr(args, "cached_cluster_map", {}) or {}
    if cached_cluster_map and not set(motif_names).issubset(cached_cluster_map):
        cached_cluster_map = {}
    clustering = None
    if cached_cluster_map and not getattr(args, "static_plots", False):
        cached_distance = getattr(args, "cached_distance_table", None)
        if cached_distance and os.path.exists(cached_distance):
            shutil.copyfile(cached_distance, matrix_out)
        else:
            pd.DataFrame(index=motif_names, columns=motif_names).to_csv(matrix_out, sep="\t")
    else:
        clustering = RegionCluster(TF_overlaps)
        clustering.cluster(threshold=args.cluster_threshold)

        convert = {m.prefix: m.name for m in motif_list}
        for cluster in clustering.clusters:
            for name in convert:
                clustering.clusters[cluster]["cluster_name"] = clustering.clusters[cluster]["cluster_name"].replace(name, convert[name])

        clustering.write_distance_mat(matrix_out)

    logger.comment("")
    if getattr(args, "write_motif_outputs", True):
        logger.info("Writing motif summaries and motif-site files")
    else:
        logger.info("Writing motif summaries")

    names, ids = [], []
    for prefix in info_table.index:
        m = [m for m in motif_list if m.prefix == prefix]
        names.append(m[0].name); ids.append(m[0].id)
    info_table.insert(0, "output_prefix", info_table.index)
    info_table.insert(1, "name", names)
    info_table.insert(2, "motif_id", ids)

    cluster_names = []
    for name in info_table.index:
        if cached_cluster_map and name in cached_cluster_map:
            cluster_names.append(cached_cluster_map[name])
        else:
            for cluster in clustering.clusters:
                if name in clustering.clusters[cluster]["member_names"]:
                    cluster_names.append(clustering.clusters[cluster]["cluster_name"])
                    break
            else:
                cluster_names.append(str(name))
    info_table.insert(3, "cluster", cluster_names)

    info_table_clustered = info_table.groupby("cluster").mean(numeric_only=True).reset_index()

    info_table["total_tfbs"] = info_table["total_tfbs"].map(int)
    for cond in args.cond_names:
        info_table[f"{cond}_bound"] = info_table[f"{cond}_bound"].map(int)
        if f"{cond}_n_replicates" in info_table.columns:
            info_table[f"{cond}_n_replicates"] = info_table[f"{cond}_n_replicates"].map(int)

    repeated_conditions = any(count > 1 for count in args.condition_replicates.values())

    for (c1, c2) in comparisons:
        base = f"{c1}_{c2}"
        info_table[base + "_change"] = info_table[base + "_change"].astype(float).round(5)
        raw_pvals = pd.to_numeric(info_table[base + "_pvalue"], errors="coerce").fillna(1.0).astype(float)
        qvals = _benjamini_hochberg(raw_pvals.to_numpy())
        info_table[base + "_pvalue"] = raw_pvals.map("{:.5E}".format, na_action="ignore")
        info_table[base + "_qvalue_bh"] = pd.Series(qvals, index=info_table.index).map("{:.5E}".format, na_action="ignore")
        info_table[base + "_significant_fdr05"] = pd.Series(qvals <= 0.05, index=info_table.index).fillna(False).astype(bool)

        names_series = info_table["output_prefix"]
        changes = info_table[base + "_change"].astype(float)
        pvals = raw_pvals
        filtered_p = pvals[pvals > 0]
        pval_min = np.percentile(filtered_p, 5) if len(filtered_p) >= 1 else 1.0
        change_min, change_max = np.percentile(changes, [5, 95])

        for i, (chg, p) in enumerate(zip(changes, pvals)):
            # info_table.at[names_series[i], base + "_highlighted"] = (chg < change_min) or (chg > change_max) or (p < pval_min)
            name_key = names_series.iloc[i] if hasattr(names_series, "iloc") else names_series[i]
            info_table.at[name_key, f"{base}_highlighted"] = (chg < change_min) or (chg > change_max) or (p < pval_min)

    if not repeated_conditions:
        single_sample_uncertainty_cols = [f"{cond}_score_sd" for cond in args.cond_names]
        single_sample_uncertainty_cols += [
            f"{c1}_{c2}_{metric}"
            for (c1, c2), metric in itertools.product(comparisons, ["delta_fp_se", "log2fc_se"])
        ]
        info_table = info_table.drop(columns=[c for c in single_sample_uncertainty_cols if c in info_table.columns])

    diff_results_out = os.path.join(args.outdir, args.prefix + "_results.txt")
    info_table.to_csv(diff_results_out, sep="\t", index=False, header=True, na_rep="NA")
    if getattr(args, "match_only", False):
        _write_match_motif_site_cache(args, motif_names, logger)

    write_replicate_report = args.replicate_report == "on" or (
        args.replicate_report == "auto" and (repeated_conditions or args.replicate_map is not None)
    )
    if write_replicate_report and len(args.cond_names) > 1:
        report_out = args.replicate_report_out or os.path.join(args.outdir, args.prefix + "_replicate_report.tsv")
        summary_out = args.replicate_summary_out or os.path.join(args.outdir, args.prefix + "_replicate_summary.tsv")
        figure_out = args.replicate_figure_out or os.path.join(args.outdir, args.prefix + "_replicate_report.png")
        try:
            build_replicate_report(
                diff_results_out,
                report_out,
                summary_output=summary_out,
                figure_output=figure_out,
                replicate_map=args.replicate_map,
            )
            logger.info(f"Wrote replicate-aware differential-footprint report to {report_out}")
        except Exception as exc:
            logger.warning(f"Could not write replicate-aware differential-footprint report: {exc}")

    if not args.skip_excel:
        diff_results_excel = os.path.join(args.outdir, args.prefix + "_results.xlsx")
        with pd.ExcelWriter(diff_results_excel, engine='xlsxwriter') as writer:
            info_table.to_excel(writer, index=False, sheet_name="Individual motifs")
            info_table_clustered.to_excel(writer, index=False, sheet_name="Motif clusters")
            for sheet in writer.sheets:
                ws = writer.sheets[sheet]
                n_rows = ws.dim_rowmax
                n_cols = ws.dim_colmax
                ws.autofilter(0, 0, n_rows, n_cols)

    if getattr(args, "skew_report", False) and comparisons:
        skew_pdf = os.path.join(args.outdir, args.prefix + "_results_skewness_report.pdf")
        try:
            # Prefer a programmatic API if available
            if hasattr(skewrep, "generate_skew_report"):
                skewrep.generate_skew_report(
                    results_tsv=diff_results_out,
                    out_pdf=skew_pdf,
                    out_json=None,
                    skew_method="perm",
                    skew_stat="bowley",
                    n_perm=20000,
                    seed=1,
                )
                logger.info(f"Skew/shift report saved → {os.path.basename(skew_pdf)}")
            else:
                # Backward-compatible fallback to module main() style runner
                # (expects skewrep.main_from_kwargs to exist; see note below)
                if hasattr(skewrep, "main_from_kwargs"):
                    skewrep.main_from_kwargs(
                        results_tsv=diff_results_out,
                        out_pdf=skew_pdf,
                        out_json=None,
                        skew_method="perm",
                        skew_stat="bowley",
                        n_perm=20000,
                        seed=1,
                    )
                    logger.info(f"Skew/shift report saved → {os.path.basename(skew_pdf)}")
                else:
                    logger.warning(
                        "diff_footprint_skew_report has no generate_skew_report() or main_from_kwargs(); skipping PDF.")
        except Exception as e:
            logger.warning(f"Could not generate skew/shift report: {e}")

    # ------------------------------ plots ------------------------------------ #
    if no_conditions > 1:
        logger.info("Creating diff-footprints plot(s)")
        change_cols = [c for c in info_table.columns if "_change" in c]
        pvalue_cols = [c for c in info_table.columns if "_pvalue" in c]
        info_table[change_cols] = info_table[change_cols].fillna(0)
        info_table[pvalue_cols] = info_table[pvalue_cols].fillna(1)

        for (c1, c2) in comparisons:
            base = f"{c1}_{c2}"
            for m in motif_list:
                name = m.prefix
                m.change = float(info_table.at[name, base + "_change"])
                m.pvalue = float(info_table.at[name, base + "_pvalue"])
                qvalue_col = base + "_qvalue_bh"
                m.qvalue = float(info_table.at[name, qvalue_col]) if qvalue_col in info_table.columns else 1.0
                m.logpvalue = -np.log10(m.pvalue) if m.pvalue > 0 else -np.log10(1e-308)
                m.highlighted = info_table.at[name, base + "_highlighted"]
                if m.highlighted:
                    m.group = f"{c2}_up" if m.change < 0 else f"{c1}_up"
                else:
                    m.group = "n.s."
            if figure_pdf is not None and cluster_pdf is not None:
                logger.info(f"- {c1} / {c2} (static plot)")
                volcano_fig, cluster_fig = plot_diff_footprints(motif_list, clustering, [c1, c2], args)
                apply_ascii_minus_to_figure(volcano_fig)
                apply_ascii_minus_to_figure(cluster_fig)
                figure_pdf.savefig(volcano_fig, bbox_inches='tight'); plt.close(volcano_fig)
                cluster_pdf.savefig(cluster_fig, bbox_inches='tight'); plt.close(cluster_fig)

            logger.info(f"- {c1} / {c2} (interactive plot)")
            html_out = os.path.join(args.outdir, args.prefix + "_" + base + ".html")
            aggregate_data = None
            if getattr(args, "aggregate_signals", None) and getattr(args, "plot_aggregate", "off") != "off":
                try:
                    aggregate_data = build_diff_footprint_aggregate_payload(motif_list, info_table, [c1, c2], args)
                except Exception as exc:
                    logger.warning(f"Could not build aggregate payload for interactive HTML: {exc}")
            plot_interactive_diff_footprints(
                motif_list,
                [c1, c2],
                html_out,
                aggregate_data=aggregate_data,
                title="Differential footprint report",
                report_label=getattr(args, "report_label", None),
            )

    if args.debug and len(args.cond_names) > 1:
        logger.info("Plotting heatmap across conditions (debug)")
        mean_columns = [c + "_mean_score" for c in args.cond_names]
        heatmap_table = info_table[mean_columns].apply(pd.to_numeric, errors="coerce")
        heatmap_table.index = info_table["output_prefix"]
        finite_rows = np.isfinite(heatmap_table.to_numpy()).all(axis=1)
        variable_rows = heatmap_table.nunique(axis=1, dropna=False) > 1
        valid_rows = finite_rows & variable_rows.to_numpy()
        dropped = int((~valid_rows).sum())
        if dropped > 0:
            logger.warning(
                f"Skipping {dropped} motif row(s) with non-finite or zero-variance values in debug heatmap."
            )
        heatmap_table = heatmap_table.loc[valid_rows]

        if heatmap_table.empty:
            logger.warning("Debug heatmap skipped because no finite, variable motif rows remained after filtering.")
        else:
            rows, cols = heatmap_table.shape
            figsize = (7 + cols, max(10, rows / 8.0))
            cm = sns.clustermap(
                heatmap_table, figsize=figsize, z_score=0, col_cluster=False,
                yticklabels=True, xticklabels=True, cbar_pos=(0, 0, .4, .005),
                dendrogram_ratio=(0.3, 0.01), cbar_kws={"orientation": "horizontal", 'label': 'Row z-score'},
                method="single"
            )
            plt.setp(cm.ax_heatmap.get_xticklabels(), fontsize=PDF_FONT_SIZE, fontweight="bold", rotation=45, ha="right")
            plt.setp(cm.ax_heatmap.get_yticklabels(), fontsize=PDF_FONT_SIZE, fontweight="bold")
            cm.ax_col_dendrogram.set_title('Mean scores across conditions', fontsize=PDF_FONT_SIZE, fontweight="bold")
            cm.ax_heatmap.set_ylabel("Transcription factor motifs", fontsize=PDF_FONT_SIZE, fontweight="bold", rotation=270)
            plt.tight_layout()
            apply_ascii_minus_to_figure(cm.fig)
            debug_pdf.savefig(cm.fig, bbox_inches='tight'); plt.close(cm.fig)

    if args.debug:
        debug_pdf.close()
    if figure_pdf is not None:
        figure_pdf.close()
    if cluster_pdf is not None:
        cluster_pdf.close()
    if temp_tfbs_dir is not None:
        shutil.rmtree(temp_tfbs_dir, ignore_errors=True)
    logger.end()


# ----------------------------------------------------------------------------- #
def run_cli():
    parser = argparse.ArgumentParser()
    parser = add_diff_footprints_arguments(parser)
    args = parser.parse_args()
    if getattr(args, "list_motif_dbs", False):
        print(motif_db_table())
        return
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()
    run_diff_footprints(args)


def _apply_match_motifs_project_layout(args, parser):
    if not (is_project_layout(getattr(args, "layout", None)) and getattr(args, "sample_table", None)):
        return
    if not getattr(args, "outdir", None):
        parser.error("--layout project requires --outdir")
    project = project_root(getattr(args, "outdir", None))
    if not getattr(args, "sample_table", None):
        parser.error("--layout project requires --sample-table")
    samples = read_sample_table(args.sample_table)
    args.signals = [str(footprint_bigwig_path(project, row.sample)) for row in samples]
    args.sample_names = [row.sample for row in samples]
    args.cond_names = [row.condition for row in samples]
    args.sample_output_root = str(samples_root(project))
    if not getattr(args, "peaks", None):
        args.peaks = str(analysis_peaks_path(project))


def _run_project_comparison_table(args, parser):
    if not getattr(args, "outdir", None):
        parser.error("--layout project requires --outdir")
    project = project_root(getattr(args, "outdir", None))
    if not getattr(args, "sample_table", None):
        parser.error("--layout project with --comparison-table requires --sample-table")
    samples = read_sample_table(args.sample_table)
    comparisons = read_comparison_table(args.comparison_table)
    if not getattr(args, "peaks", None):
        args.peaks = str(analysis_peaks_path(project))
    for comparison in comparisons:
        cond1_samples = samples_for_condition(samples, comparison.cond1, comparison.cond1_samples)
        cond2_samples = samples_for_condition(samples, comparison.cond2, comparison.cond2_samples)
        selected = cond1_samples + cond2_samples
        comparison_args = copy.copy(args)
        comparison_args.comparison_table = None
        comparison_args.project_dir = None
        comparison_args.sample_dirs = [str(sample_dir(project, row.sample)) for row in selected]
        comparison_args.sample_names = [row.sample for row in selected]
        comparison_args.cond_names = [comparison.cond1] * len(cond1_samples) + [comparison.cond2] * len(cond2_samples)
        comparison_args.outdir = str(comparison_dir(project, comparison.comparison))
        comparison_args.aggregate_signals = [
            str(normalized_bigwig_path(project, row.sample))
            if normalized_bigwig_path(project, row.sample).exists()
            else str(corrected_bigwig_path(project, row.sample))
            for row in selected
        ]
        comparison_args.motif_outputs = "summary" if getattr(args, "motif_outputs", "auto") == "auto" else args.motif_outputs
        comparison_args.skip_excel = True if not getattr(args, "skip_excel", False) else args.skip_excel
        run_diff_footprints(comparison_args)


def match_motifs_cli():
    parser = argparse.ArgumentParser(prog="match-motifs")
    parser = add_diff_footprints_arguments(parser, command_name="match-motifs")
    args = parser.parse_args()
    if getattr(args, "list_motif_dbs", False):
        print(motif_db_table())
        return
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()
    _apply_match_motifs_project_layout(args, parser)
    if not args.signals:
        parser.error("match-motifs expects at least one --signals bigWig")
    if args.sample_names is not None and len(args.sample_names) != len(args.signals):
        parser.error("match-motifs expects one --sample-names value per --signals bigWig when provided")
    if args.cond_names is not None and len(args.cond_names) != len(args.signals):
        parser.error("match-motifs expects one --cond-names value per --signals bigWig when provided")
    if args.prefix == "diff_footprints":
        args.prefix = "motif_matches"
    args.method = "motif"
    args.match_only = True
    args.replicate_report = "off"
    args.static_plots = False
    args.per_motif_plots = False
    args.skew_report = False
    if getattr(args, "sample_output_root", None):
        if args.sample_names is None:
            parser.error("--sample-output-root requires --sample-names")
        if len(args.sample_names) != len(args.signals):
            parser.error("--sample-names must contain one value per --signals bigWig")
        cond_names = args.cond_names or args.sample_names
        for signal, sample, cond in zip(args.signals, args.sample_names, cond_names):
            sample_args = copy.copy(args)
            sample_args.signals = [signal]
            sample_args.sample_names = [sample]
            sample_args.cond_names = [cond]
            sample_args.outdir = os.path.join(os.path.abspath(args.sample_output_root), sample, "match_motifs")
            run_diff_footprints(sample_args)
        return
    run_diff_footprints(args)


def diff_footprints_cli():
    parser = argparse.ArgumentParser(prog="diff-footprints")
    parser = add_diff_footprints_arguments(parser, command_name="diff-footprints")
    args = parser.parse_args()
    if getattr(args, "list_motif_dbs", False):
        print(motif_db_table())
        return
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()
    if is_project_layout(getattr(args, "layout", None)) and getattr(args, "comparison_table", None):
        _run_project_comparison_table(args, parser)
        return
    folder_mode = bool(args.sample_dirs or args.project_dir)
    n_inputs = len(args.sample_dirs or []) if args.sample_dirs else len(args.signals or [])
    if args.project_dir and not args.sample_dirs:
        n_inputs = None
    if args.method == "motif" and not folder_mode and (not args.signals or len(args.signals) < 2):
        parser.error("diff-footprints expects at least two --signals bigWigs, or --sample-dirs/--project-dir")
    if args.signals and folder_mode:
        parser.error("diff-footprints cannot combine --signals with --sample-dirs/--project-dir")
    if n_inputs is not None and args.sample_names is not None and len(args.sample_names) != n_inputs:
        parser.error("diff-footprints expects one --sample-names value per input sample when provided")
    if n_inputs is not None and args.cond_names is not None and len(args.cond_names) != n_inputs:
        parser.error("diff-footprints expects one --cond-names value per input sample when provided")
    run_diff_footprints(args)


if __name__ == '__main__':
    run_cli()
