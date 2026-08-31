#!/usr/bin/env python

"""
plot-aggregate command for aggregate signal visualization from TFBS and bigWigs.

This implementation provides:
- aggregate plotting from explicit BED inputs or a BED directory
- optional CSV exports of aggregated signals and scores
- grid-based plot layouts with consistent subplot sizing
- fp-tools-specific layout and reporting behavior
"""

import argparse
import copy
import itertools
import os
import re
import sys
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from fp_tools.utils import bigwig as pyBigWig
from sklearn import preprocessing

from fp_tools.parsers import add_aggregate_arguments
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.multiscale import aggregate_multiscale_tensor, load_multiscale_npz
from fp_tools.utils.normalization import fit_quantile_normalizers
from fp_tools.utils.plotting_style import PDF_FONT_SIZE, apply_pdf_style, ascii_tick_formatter
from fp_tools.utils.project_layout import (
    corrected_bigwig_path,
    is_project_layout,
    match_motifs_dir,
    project_root,
    read_sample_table,
    reports_dir,
)
from fp_tools.utils.regions import OneRegion, RegionList
from fp_tools.utils.signals import fast_rolling_math
from fp_tools.utils.utilities import check_files, check_required, make_directory

PANEL_SIZE_IN = 3.05
DEFAULT_GRID_COLS = 5


def default_multiscale_output(output_path):
    """Derive the companion multiscale figure path from the main aggregate path."""

    root, ext = os.path.splitext(output_path)
    if not ext:
        ext = ".pdf"
    return f"{root}_multiscale{ext}"


def plot_multiscale_aggregate_npz(npz_path, output_path, title="Multiscale aggregate"):
    """Render a plot-aggregate companion figure from a multiscale NPZ sidecar."""

    data = load_multiscale_npz(str(npz_path))
    aggregate = aggregate_multiscale_tensor(data, align="center")
    scales = data["scales"].astype(int)
    if aggregate.size == 0:
        aggregate = np.zeros((len(scales), 1), dtype=float)

    width = aggregate.shape[1]
    center = width // 2
    xvals = np.arange(width) - center
    profile = np.nanmean(aggregate, axis=0)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(5.2, 4.6),
        gridspec_kw={"height_ratios": [2.5, 1.2]},
        constrained_layout=True,
    )
    image = axes[0].imshow(
        aggregate,
        aspect="auto",
        cmap="viridis",
        interpolation="nearest",
        extent=[xvals[0], xvals[-1], len(scales) - 0.5, -0.5],
    )
    axes[0].set_yticks(np.arange(len(scales)))
    axes[0].set_yticklabels([str(scale) for scale in scales])
    axes[0].set_ylabel("scale bp", fontsize=PDF_FONT_SIZE, fontweight="bold")
    axes[0].set_title(title, fontsize=PDF_FONT_SIZE, fontweight="bold")
    axes[0].xaxis.set_major_formatter(ascii_tick_formatter())
    fig.colorbar(image, ax=axes[0], label="depletion")

    axes[1].plot(xvals, profile, color="black", linewidth=1)
    axes[1].axvline(0, color="grey", linestyle="dashed", linewidth=0.8)
    axes[1].set_xlabel("bp from center", fontsize=PDF_FONT_SIZE, fontweight="bold")
    axes[1].set_ylabel("mean", fontsize=PDF_FONT_SIZE, fontweight="bold")
    axes[1].xaxis.set_major_formatter(ascii_tick_formatter())
    axes[1].yaxis.set_major_formatter(ascii_tick_formatter(decimals=2))

    plt.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    return output_path


def parse_grid_spec(grid_spec):
    """Parse a grid spec like '2x5' into integer rows/cols."""
    if grid_spec is None:
        return None
    match = re.fullmatch(r"\s*(\d+)\s*[xX]\s*(\d+)\s*", grid_spec)
    if match is None:
        raise ValueError("Grid must be formatted as <rows>x<cols>, e.g. 2x5.")
    rows, cols = int(match.group(1)), int(match.group(2))
    if rows < 1 or cols < 1:
        raise ValueError("Grid rows and columns must both be >= 1.")
    return rows, cols



def build_condition_groups(signal_labels, cond_names=None):
    """Return ordered condition labels and condition-to-sample mapping."""

    if cond_names is None:
        cond_names = list(signal_labels)
    if len(cond_names) != len(signal_labels):
        raise ValueError("--cond-names must have the same length as --signals")
    groups = {}
    for signal_name, cond in zip(signal_labels, cond_names):
        groups.setdefault(cond, []).append(signal_name)
    return list(groups.keys()), groups


def apply_quantile_normalization_to_signal_dict(signal_dict, region_names, regions_dict, signal_names, condition_names, condition_groups, mode, logger):
    """Apply shared diff-footprints-style quantile normalization to per-region signal arrays."""

    mode = (mode or "none").replace("_", "-")
    if mode == "none":
        return copy.deepcopy(signal_dict)
    normalized = copy.deepcopy(signal_dict)
    if mode == "sample-quantile":
        names = list(signal_names)
        arrays = [np.concatenate([signal_dict[name][reg.tup()] for rid in region_names for reg in regions_dict[rid]]) for name in names]
        norm_objects, _ = fit_quantile_normalizers(arrays, names, logger=logger)
        for name in names:
            for tup, arr in normalized[name].items():
                normalized[name][tup] = np.maximum(0.0, norm_objects[name].normalize(arr))
        return normalized
    if mode == "condition-quantile":
        arrays = []
        for cond in condition_names:
            sample_arrays = []
            for sample in condition_groups[cond]:
                sample_arrays.append(np.concatenate([signal_dict[sample][reg.tup()] for rid in region_names for reg in regions_dict[rid]]))
            arrays.append(np.mean(np.vstack(sample_arrays), axis=0))
        norm_objects, _ = fit_quantile_normalizers(arrays, condition_names, logger=logger)
        for cond in condition_names:
            for sample in condition_groups[cond]:
                for tup, arr in normalized[sample].items():
                    normalized[sample][tup] = np.maximum(0.0, norm_objects[cond].normalize(arr))
        return normalized
    raise ValueError(f"Unsupported normalization mode: {mode}")


def normalize_site_profiles(signalmat, mode="none", clip=5.0):
    """Normalize motif-centered profiles independently using their outer flanks.

    ``flank-rms`` gives each site comparable influence in an aggregate while
    retaining the sign of corrected cut-site signal. A lower-quartile scale
    floor prevents nearly signal-free sites from being amplified arbitrarily.
    """

    matrix = np.asarray(signalmat, dtype=float)
    if matrix.ndim != 2:
        raise ValueError("site profiles must be a two-dimensional matrix")
    if mode == "none" or matrix.shape[0] == 0:
        return matrix.copy()
    if mode not in {"flank-center", "flank-rms"}:
        raise ValueError(f"Unsupported site normalization mode: {mode}")

    outer_width = max(3, matrix.shape[1] // 6)
    outer = np.concatenate((matrix[:, :outer_width], matrix[:, -outer_width:]), axis=1)
    flank_mean = np.nanmean(outer, axis=1)
    centered = matrix - flank_mean[:, None]
    if mode == "flank-center":
        return centered

    scale = np.sqrt(np.nanmean(np.square(outer), axis=1))
    finite_positive = scale[np.isfinite(scale) & (scale > 0)]
    scale_floor = float(np.quantile(finite_positive, 0.25)) if len(finite_positive) else 1.0
    normalized = centered / np.maximum(np.nan_to_num(scale, nan=0.0), scale_floor)[:, None]
    if clip and clip > 0:
        normalized = np.clip(normalized, -float(clip), float(clip))
    return normalized


def profile_shape_diagnostics(signalmat, motif_width, flank):
    """Summarize motif-centered depletion and its site-level uncertainty."""

    matrix = np.asarray(signalmat, dtype=float)
    n_sites = int(matrix.shape[0])
    if n_sites == 0:
        return {
            "n_sites": 0, "mean_depletion": np.nan, "depletion_se": np.nan,
            "depletion_z": np.nan, "depletion_effect_size": np.nan,
            "shoulder_asymmetry": np.nan, "detectability": "underpowered",
        }

    half_width = max(1, int(np.ceil(float(motif_width) / 2.0)))
    center_start = max(0, int(flank) - half_width)
    center_end = min(matrix.shape[1], int(flank) + half_width + 1)
    shoulder_inner = max(half_width + 4, 16)
    shoulder_outer = min(int(flank), max(shoulder_inner + 8, 40))
    left = matrix[:, max(0, int(flank) - shoulder_outer):max(0, int(flank) - shoulder_inner)]
    right = matrix[:, min(matrix.shape[1], int(flank) + shoulder_inner + 1):min(matrix.shape[1], int(flank) + shoulder_outer + 1)]
    center = matrix[:, center_start:center_end]
    if center.shape[1] == 0 or left.shape[1] == 0 or right.shape[1] == 0:
        return {
            "n_sites": n_sites, "mean_depletion": np.nan, "depletion_se": np.nan,
            "depletion_z": np.nan, "depletion_effect_size": np.nan,
            "shoulder_asymmetry": np.nan, "detectability": "underpowered",
        }

    center_mean = np.nanmean(center, axis=1)
    left_mean = np.nanmean(left, axis=1)
    right_mean = np.nanmean(right, axis=1)
    depletion = (left_mean + right_mean) / 2.0 - center_mean
    finite = depletion[np.isfinite(depletion)]
    n_finite = int(len(finite))
    mean_depletion = float(np.mean(finite)) if n_finite else np.nan
    sd = float(np.std(finite, ddof=1)) if n_finite > 1 else np.nan
    se = sd / np.sqrt(n_finite) if n_finite > 1 and sd > 0 else np.nan
    zscore = mean_depletion / se if np.isfinite(se) and se > 0 else np.nan
    effect_size = mean_depletion / sd if np.isfinite(sd) and sd > 0 else np.nan
    asymmetry = float(np.nanmean(np.abs(left_mean - right_mean)))
    asymmetry /= max(abs(mean_depletion), np.finfo(float).eps)

    if n_finite < 30:
        status = "underpowered"
    elif mean_depletion <= 0 or not np.isfinite(zscore) or zscore < 1.96:
        status = "not detected"
    elif effect_size >= 0.5 and zscore >= 5:
        status = "strong"
    elif effect_size >= 0.2 and zscore >= 3:
        status = "detectable"
    else:
        status = "weak"
    return {
        "n_sites": n_finite,
        "mean_depletion": mean_depletion,
        "depletion_se": float(se),
        "depletion_z": float(zscore),
        "depletion_effect_size": float(effect_size),
        "shoulder_asymmetry": asymmetry,
        "detectability": status,
    }


def calculate_group_aggregates(signal_dict, regions_dict, region_names, condition_names, condition_groups, motif_widths, args):
    """Calculate aggregate mean profiles and optional replicate SD profiles."""

    sample_aggregates = {sample: {} for samples in condition_groups.values() for sample in samples}
    sample_site_se = {sample: {} for sample in sample_aggregates}
    sample_site_matrices = {sample: {} for sample in sample_aggregates}
    for sample_name in sample_aggregates:
        for region_name in region_names:
            profiles = [signal_dict[sample_name][reg.tup()] for reg in regions_dict[region_name]]
            signalmat = np.vstack(profiles) if profiles else np.zeros((0, args.width), dtype=float)
            if signalmat.shape[0] == 0:
                aggregate = np.zeros(args.width)
            else:
                max_values = np.max(signalmat, axis=1)
                upper_limit = np.percentile(max_values, [100 * args.remove_outliers])[0]
                signalmat = signalmat[max_values <= upper_limit]
                if signalmat.shape[0] == 0:
                    aggregate = np.zeros(args.width)
                else:
                    signalmat = normalize_site_profiles(
                        signalmat,
                        mode=getattr(args, "site_normalization", "none"),
                        clip=getattr(args, "site_normalization_clip", 5.0),
                    )
                    if args.log_transform:
                        signal_mat_abs = np.abs(signalmat)
                        signal_mat_log = np.log2(signal_mat_abs + 1)
                        signal_mat_log[signalmat < 0] *= -1
                        signalmat = signal_mat_log
                    if args.smooth > 1:
                        smoothed_rows = []
                        for site_profile in signalmat:
                            extended = np.pad(site_profile, args.smooth, "edge")
                            smoothed = fast_rolling_math(extended.astype("float64"), args.smooth, "mean")
                            smoothed_rows.append(smoothed[args.smooth:-args.smooth])
                        signalmat = np.vstack(smoothed_rows)
                    aggregate = np.nanmean(signalmat, axis=0)
                    if args.normalize:
                        aggregate = preprocessing.minmax_scale(aggregate)
            sample_aggregates[sample_name][region_name] = aggregate
            sample_site_matrices[sample_name][region_name] = signalmat
            if signalmat.shape[0] > 1:
                sample_site_se[sample_name][region_name] = np.nanstd(signalmat, axis=0, ddof=1) / np.sqrt(signalmat.shape[0])
            else:
                sample_site_se[sample_name][region_name] = np.full(aggregate.shape, np.nan)

    aggregate_dict = {cond: {} for cond in condition_names}
    aggregate_sd_dict = {cond: {} for cond in condition_names}
    aggregate_ci_dict = {cond: {} for cond in condition_names}
    aggregated_fp_scores = {}
    stats_rows = []
    for cond in condition_names:
        samples = condition_groups[cond]
        for region_name in region_names:
            stack = np.vstack([sample_aggregates[sample][region_name] for sample in samples])
            mean_profile = np.nanmean(stack, axis=0)
            sd_profile = np.nanstd(stack, axis=0, ddof=1) if len(samples) > 1 else np.full(mean_profile.shape, np.nan)
            aggregate_dict[cond][region_name] = mean_profile
            aggregate_sd_dict[cond][region_name] = sd_profile
            site_se_stack = np.vstack([sample_site_se[sample][region_name] for sample in samples])
            combined_site_se = np.sqrt(np.nansum(np.square(site_se_stack), axis=0)) / len(samples)
            aggregate_ci_dict[cond][region_name] = 1.96 * combined_site_se
            flank_len = args.flank
            motif_w = motif_widths[region_name]
            mid_start = max(0, int(flank_len - np.floor(motif_w / 2.0)))
            mid_end = min(len(mean_profile), int(flank_len + np.ceil(motif_w / 2.0)))
            mean_mid = np.mean(mean_profile[mid_start:mid_end]) if mid_end - mid_start > 0 else 0.0
            flank_indices = list(range(0, mid_start)) + list(range(mid_end, len(mean_profile)))
            mean_flank = np.mean(mean_profile[flank_indices]) if flank_indices else 0.0
            depletion = max(0.0, mean_flank - mean_mid)
            aggregated_fp_scores[(cond, region_name)] = mean_flank + depletion
            pooled_matrix = np.vstack([sample_site_matrices[sample][region_name] for sample in samples])
            diagnostics = profile_shape_diagnostics(pooled_matrix, motif_w, flank_len)
            stats_rows.append({
                "condition": cond,
                "regions": region_name,
                "n_replicates": len(samples),
                "mean_profile": float(np.nanmean(mean_profile)),
                "mean_profile_sd": float(np.nanmean(sd_profile)) if np.isfinite(sd_profile).any() else np.nan,
                "mean_flank": float(mean_flank),
                "mean_center": float(mean_mid),
                "aggregate_fp_score": float(aggregated_fp_scores[(cond, region_name)]),
                "site_normalization": getattr(args, "site_normalization", "none"),
                **diagnostics,
            })
    return aggregate_dict, aggregate_sd_dict, aggregate_ci_dict, aggregated_fp_scores, stats_rows


def plot_normalization_comparison(raw_aggregates, norm_aggregates, condition_names, region_names, output, title="Raw vs quantile-normalized aggregates"):
    """Write a compact raw-vs-normalized aggregate comparison figure."""

    total = len(condition_names) * len(region_names)
    fig, axes = plt.subplots(total, 2, figsize=(7.2, max(2.2, total * 1.8)), squeeze=False, constrained_layout=True)
    row = 0
    for region_name in region_names:
        for cond in condition_names:
            raw = raw_aggregates[cond][region_name]
            norm = norm_aggregates[cond][region_name]
            flank = len(raw) // 2
            xvals = np.arange(len(raw)) - flank
            axes[row, 0].plot(xvals, raw, color="0.25", linewidth=1)
            axes[row, 1].plot(xvals, norm, color="tab:blue", linewidth=1)
            axes[row, 0].set_ylabel(f"{cond}\n{region_name}", fontsize=PDF_FONT_SIZE)
            axes[row, 0].set_title("Raw" if row == 0 else "", fontsize=PDF_FONT_SIZE, fontweight="bold")
            axes[row, 1].set_title("Normalized" if row == 0 else "", fontsize=PDF_FONT_SIZE, fontweight="bold")
            for ax in axes[row]:
                ax.axvline(0, color="0.7", linewidth=0.8, linestyle="--")
                ax.set_xlim(-flank, flank)
                ax.xaxis.set_major_formatter(ascii_tick_formatter())
                ax.yaxis.set_major_formatter(ascii_tick_formatter(decimals=2))
            row += 1
    fig.suptitle(title, fontsize=PDF_FONT_SIZE, fontweight="bold")
    dpi = 300 if str(output).lower().endswith(".png") else None
    fig.savefig(output, bbox_inches="tight", dpi=dpi)
    plt.close(fig)


def _path_stem(path):
    return os.path.splitext(os.path.basename(str(path)))[0]


def _aggregate_output_format(args):
    requested = getattr(args, "format", "auto") or "auto"
    if requested != "auto":
        return requested
    ext = os.path.splitext(str(args.output or ""))[1].lower()
    return "html" if ext in {".html", ".htm"} else "pdf"


def _match_dir_rows(args):
    signals = list(args.signals or [])
    if not signals:
        raise ValueError("--signals is required with --match-dir")

    labels = list(args.signal_labels or [_path_stem(path) for path in signals])
    if len(labels) != len(signals):
        raise ValueError("--signal-labels must have the same length as --signals")

    groups_defined = bool(args.cond_names)
    conditions = list(args.cond_names or labels)
    if len(conditions) != len(signals):
        raise ValueError("--cond-names must have the same length as --signals")

    match_dirs = list(getattr(args, "sample_dirs", None) or args.match_dir or [])
    if len(match_dirs) == 1:
        match_dirs = match_dirs * len(signals)
    if len(match_dirs) != len(signals):
        raise ValueError("--match-dir must contain one directory or one directory per --signals input")

    rows = []
    for signal, label, condition, match_dir in zip(signals, labels, conditions, match_dirs):
        rows.append({
            "sample": label,
            "label": label,
            "condition": condition,
            "_groups_defined": groups_defined,
            "signal": signal,
            "match_dir": match_dir,
        })
    return rows


def _run_match_dir_html(args, logger):
    from fp_tools.tools.plot_aggregate_batch import _read_manifest, build_payload, build_payload_from_tfbs, merge_payloads, read_embedded_payload, write_html

    payloads = []
    if getattr(args, "manifest", None):
        manifest_rows = _read_manifest(args.manifest)
        for row in manifest_rows:
            if "match_dir" not in row and row.get("sample_dir"):
                row["match_dir"] = row["sample_dir"]
        payloads.append(
            build_payload(
                manifest_rows,
                flank=max(1, int(args.flank)),
                top_n=max(1, int(args.top_n)),
                normalization=args.normalization,
                motif_names=args.motifs,
                site_set=args.site_set,
            )
        )
    if getattr(args, "input_html", None):
        for path in args.input_html:
            payloads.append(read_embedded_payload(path))
    if getattr(args, "match_dir", None) or getattr(args, "sample_dirs", None):
        try:
            rows = _match_dir_rows(args)
        except ValueError as exc:
            logger.error(f"ERROR: {exc}")
            sys.exit(1)
        payloads.append(
            build_payload(
                rows,
                flank=max(1, int(args.flank)),
                top_n=max(1, int(args.top_n)),
                normalization=args.normalization,
                motif_names=args.motifs,
                site_set=args.site_set,
            )
        )
    elif getattr(args, "TFBS", None) and getattr(args, "signals", None):
        labels = list(args.signal_labels or [_path_stem(path) for path in args.signals])
        groups_defined = bool(args.cond_names)
        conditions = list(args.cond_names or labels)
        if len(labels) != len(args.signals) or len(conditions) != len(args.signals):
            logger.error("ERROR: --signal-labels and --cond-names must have the same length as --signals")
            sys.exit(1)
        tfbs_inputs = list(args.TFBS)
        if len(tfbs_inputs) == 1 and os.path.isdir(tfbs_inputs[0]):
            tfbs_inputs = sorted(str(path) for path in Path(tfbs_inputs[0]).glob("*.bed"))
            if args.top_n:
                tfbs_inputs = tfbs_inputs[: max(1, int(args.top_n))]
        if not tfbs_inputs:
            logger.error("ERROR: no BED files found for --TFBS")
            sys.exit(1)
        motif_labels = list(args.TFBS_labels or [_path_stem(path) for path in tfbs_inputs])
        payloads.append(
            build_payload_from_tfbs(
                tfbs_inputs,
                args.signals,
                labels,
                conditions,
                flank=max(1, int(args.flank)),
                normalization=args.normalization,
                motif_labels=motif_labels,
                groups_defined=groups_defined,
            )
        )
    if not payloads:
        logger.error("ERROR: --format html requires --sample-dirs/--match-dir, --manifest, --input-html, or --TFBS with --signals")
        sys.exit(1)
    payload = merge_payloads(payloads)
    make_directory(os.path.dirname(os.path.abspath(args.output)) or ".")
    write_html(payload, args.output, args.title, default_layout=args.default_layout, show_summary=not getattr(args, "hide_summary", False))
    logger.info(f"Wrote interactive aggregate HTML report to {args.output}")


def _resolve_match_dir_tfbs(args, logger):
    if not getattr(args, "match_dir", None):
        return
    from fp_tools.tools.plot_aggregate_batch import _discover_motifs, _motif_bed_path, _select_motif_prefixes

    try:
        rows = _match_dir_rows(args)
    except ValueError as exc:
        logger.error(f"ERROR: {exc}")
        sys.exit(1)

    motif_meta = {}
    motif_scores = {}
    for row in rows:
        for motif in _discover_motifs(row["match_dir"]):
            prefix = str(motif["prefix"])
            motif_meta.setdefault(prefix, motif)
            motif_scores[prefix] = max(float(motif.get("score") or 0.0), motif_scores.get(prefix, 0.0))
    ranked = sorted(motif_meta.values(), key=lambda row: (-motif_scores[str(row["prefix"])], str(row.get("name") or row["prefix"])))
    selected = _select_motif_prefixes(ranked, args.motifs, max(1, int(args.top_n)))
    first = rows[0]
    tfbs = []
    labels = []
    for prefix in selected:
        bed = _motif_bed_path(
            first["match_dir"],
            prefix,
            condition=first["condition"],
            sample=first["sample"],
            site_set=args.site_set,
        )
        if not bed.exists():
            logger.warning(f"Skipping motif {prefix}; no BED file found at {bed}")
            continue
        tfbs.append(str(bed))
        labels.append(str(motif_meta[prefix].get("name") or prefix))
    if not tfbs:
        logger.error("ERROR: no motif BED files were resolved from --match-dir")
        sys.exit(1)
    args.TFBS = tfbs
    if args.TFBS_labels is None:
        args.TFBS_labels = labels


def run_aggregate(args):
    """Create aggregate plots and optional aggregate exports."""

    from fp_tools.utils.intervals import filter_regions

    apply_pdf_style()
    logger = FpToolsLogger("plot-aggregate", args.verbosity)
    logger.begin()

    if is_project_layout(getattr(args, "layout", None)) and getattr(args, "sample_table", None):
        if not getattr(args, "outdir", None):
            logger.error("ERROR: --layout project requires --outdir")
            sys.exit(1)
        project = project_root(getattr(args, "outdir", None))
        samples = read_sample_table(args.sample_table)
        args.signals = [str(corrected_bigwig_path(project, row.sample)) for row in samples]
        args.signal_labels = [row.sample for row in samples]
        args.cond_names = [row.condition for row in samples]
        args.match_dir = [str(match_motifs_dir(project, row.sample)) for row in samples]
        args.format = "html" if getattr(args, "format", "auto") == "auto" else args.format
        if getattr(args, "output", None) == "fp-tools_aggregate.pdf":
            args.output = str(reports_dir(project) / "plot_aggregate.html")

    html_inputs = any([
        getattr(args, "match_dir", None),
        getattr(args, "sample_dirs", None),
        getattr(args, "manifest", None),
        getattr(args, "input_html", None),
        getattr(args, "TFBS", None),
    ])
    if html_inputs and _aggregate_output_format(args) == "html":
        logger.arguments_overview(add_aggregate_arguments(argparse.ArgumentParser()), args)
        _run_match_dir_html(args, logger)
        logger.end()
        return

    _resolve_match_dir_tfbs(args, logger)

    if args.TFBS is None:
        args.TFBS = []
    if len(args.TFBS) == 1 and os.path.isdir(args.TFBS[0]):
        bed_dir = args.TFBS[0]
        beds = sorted(
            os.path.join(bed_dir, name)
            for name in os.listdir(bed_dir)
            if name.endswith(".bed")
        )
        if not beds:
            logger.error(f"No .bed files found under {bed_dir}")
            sys.exit(1)
        args.TFBS = beds

    if args.output_aggregated_signals is None and args.output_csv is not None:
        args.output_aggregated_signals = args.output_csv
    if args.multiscale_npz is not None and args.output_multiscale_aggregate is None:
        args.output_multiscale_aggregate = default_multiscale_output(args.output)

    logger.arguments_overview(add_aggregate_arguments(argparse.ArgumentParser()), args)
    logger.output_files([
        args.output,
        args.output_txt,
        args.output_aggregated_signals,
        args.output_aggregated_scores,
        args.output_aggregated_stats,
        args.output_multiscale_aggregate,
    ])

    check_required(args, ["TFBS", "signals"])
    check_files([args.TFBS, args.signals, args.regions, args.whitelist, args.blacklist, args.multiscale_npz], action="r")

    out_parent_dirs = []
    for output_path in [
        args.output,
        args.output_txt,
        args.output_aggregated_signals,
        args.output_aggregated_scores,
        args.output_aggregated_stats,
        args.output_multiscale_aggregate,
    ]:
        if output_path:
            parent = os.path.dirname(os.path.abspath(output_path))
            if parent:
                out_parent_dirs.append(parent)
    for parent in sorted(set(out_parent_dirs)):
        make_directory(parent)

    check_files(
        [args.output, args.output_txt, args.output_aggregated_signals, args.output_aggregated_scores, args.output_aggregated_stats, args.output_multiscale_aggregate, args.normalization_comparison_output],
        action="w",
    )

    if args.TFBS_labels is not None and len(args.TFBS) != len(args.TFBS_labels):
        logger.error(
            f"ERROR: --TFBS and --TFBS-labels have different lengths ({len(args.TFBS)} vs. {len(args.TFBS_labels)})"
        )
        sys.exit(1)
    if args.region_labels is not None and len(args.regions) != len(args.region_labels):
        logger.error(
            f"ERROR: --regions and --region-labels have different lengths ({len(args.regions)} vs. {len(args.region_labels)})"
        )
        sys.exit(1)
    if args.signal_labels is not None and len(args.signals) != len(args.signal_labels):
        logger.error(
            f"ERROR: --signals and --signal-labels have different lengths ({len(args.signals)} vs. {len(args.signal_labels)})"
        )
        sys.exit(1)
    if args.cond_names is not None and len(args.signals) != len(args.cond_names):
        logger.error(
            f"ERROR: --signals and --cond-names have different lengths ({len(args.signals)} vs. {len(args.cond_names)})"
        )
        sys.exit(1)

    args.TFBS_labels = (
        [os.path.splitext(os.path.basename(path))[0] for path in args.TFBS]
        if args.TFBS_labels is None else args.TFBS_labels
    )
    args.region_labels = (
        [os.path.splitext(os.path.basename(path))[0] for path in args.regions]
        if args.region_labels is None else args.region_labels
    )
    args.signal_labels = (
        [os.path.splitext(os.path.basename(path))[0] for path in args.signals]
        if args.signal_labels is None else args.signal_labels
    )

    if len(set(args.TFBS_labels)) < len(args.TFBS_labels):
        logger.error("ERROR: --TFBS-labels are not allowed to contain duplicates.")
        sys.exit(1)

    control_label = args.control_label
    if control_label is not None and control_label not in args.signal_labels:
        logger.error(
            f"ERROR: --control-label '{control_label}' not found among signal-labels: {args.signal_labels}"
        )
        sys.exit(1)

    logger.info("---- Processing input ----")
    logger.info("Reading information from .bed-files")

    region_names = []
    if len(args.regions) > 0:
        logger.info("Overlapping sites to --regions")
        regions_dict = {}
        for tfbs_idx, region_idx in itertools.product(range(len(args.TFBS)), range(len(args.regions))):
            tfbs_f = args.TFBS[tfbs_idx]
            region_f = args.regions[region_idx]

            name = args.TFBS_labels[tfbs_idx] + " <OVERLAPPING> " + args.region_labels[region_idx]
            region_names.append(name)
            regions_dict[name] = filter_regions(RegionList().from_bed(tfbs_f), region_f)

            if args.negate:
                name = args.TFBS_labels[tfbs_idx] + " <NOT OVERLAPPING> " + args.region_labels[region_idx]
                region_names.append(name)
                regions_dict[name] = filter_regions(RegionList().from_bed(tfbs_f), region_f, invert=True)
    else:
        region_names = list(args.TFBS_labels)
        regions_dict = {
            args.TFBS_labels[i]: RegionList().from_bed(args.TFBS[i])
            for i in range(len(args.TFBS))
        }
        for name in region_names:
            logger.stats(f"COUNT {name}: {len(regions_dict[name])} sites")

    if len(args.whitelist) > 0 or len(args.blacklist) > 0:
        logger.info("Subsetting regions on whitelist/blacklist")
        for regions_id in regions_dict:
            logger.stats(f"Found {len(regions_dict[regions_id])} sites in {regions_id}")

            if len(args.whitelist) > 0:
                for whitelist_f in args.whitelist:
                    regions_dict[regions_id] = filter_regions(regions_dict[regions_id], whitelist_f)
                    logger.stats(f"Overlapped to whitelist -> {len(regions_dict[regions_id])}")

            if len(args.blacklist) > 0:
                for blacklist_f in args.blacklist:
                    regions_dict[regions_id] = filter_regions(regions_dict[regions_id], blacklist_f, invert=True)
                    logger.stats(f"Removed blacklist -> {len(regions_dict[regions_id])}")

    motif_widths = {}
    for regions_id, site_list in regions_dict.items():
        motif_widths[regions_id] = site_list[0].get_width() if len(site_list) > 0 else 0

    logger.info("Reading signal from bigwigs")
    args.width = args.flank * 2
    signal_dict = {}

    for signal_idx, signal_f in enumerate(args.signals):
        signal_name = args.signal_labels[signal_idx]
        signal_dict[signal_name] = {}
        pybw = pyBigWig.open(signal_f)
        boundaries = pybw.chroms()

        logger.info(f"- Reading signal from {signal_name}")
        for regions_id in regions_dict:
            original = copy.deepcopy(regions_dict[regions_id])
            regions_dict[regions_id].apply_method(OneRegion.set_width, args.width)

            invalid = [
                idx for idx, region in enumerate(regions_dict[regions_id])
                if region.check_boundary(boundaries, action="remove") is None
            ]
            for invalid_idx in reversed(invalid):
                logger.warning(
                    "Region '{reg}' ('{orig}' before flank extension) from bed regions '{rid}' is out of boundaries. Excluding.".format(
                        reg=regions_dict[regions_id][invalid_idx].pretty(),
                        orig=original[invalid_idx].pretty(),
                        rid=regions_id,
                    )
                )
                del regions_dict[regions_id][invalid_idx]

            for one_region in regions_dict[regions_id]:
                tup = one_region.tup()
                if tup not in signal_dict[signal_name]:
                    signal_dict[signal_name][tup] = one_region.get_signal(pybw, logger=logger, key=signal_name)

        pybw.close()

    sample_names = args.signal_labels
    try:
        signal_names, condition_groups = build_condition_groups(sample_names, args.cond_names)
    except ValueError as exc:
        logger.error(f"ERROR: {exc}")
        sys.exit(1)

    logger.info("Calculating aggregate signals")
    raw_aggregate_dict, raw_aggregate_sd_dict, _, _, _ = calculate_group_aggregates(
        signal_dict, regions_dict, region_names, signal_names, condition_groups, motif_widths, args
    )
    normalized_signal_dict = apply_quantile_normalization_to_signal_dict(
        signal_dict, region_names, regions_dict, sample_names, signal_names, condition_groups, args.normalization, logger
    )
    aggregate_dict, aggregate_sd_dict, aggregate_ci_dict, aggregated_fp_scores, aggregate_stats = calculate_group_aggregates(
        normalized_signal_dict, regions_dict, region_names, signal_names, condition_groups, motif_widths, args
    )
    aggregate_stats_lookup = {
        (row["condition"], row["regions"]): row for row in aggregate_stats
    }

    if args.normalization_comparison_output is not None and args.normalization != "none":
        plot_normalization_comparison(
            raw_aggregate_dict,
            aggregate_dict,
            signal_names,
            region_names,
            args.normalization_comparison_output,
            title=f"Raw vs {args.normalization} aggregates",
        )

    signal_dict = None
    normalized_signal_dict = None

    if args.output_txt is not None:
        with open(args.output_txt, "w") as handle:
            handle.write("### AGGREGATE\n")
            handle.write("# Signal\tRegions\tAggregate\n")
            for signal_name in signal_names:
                for region_name in region_names:
                    agg_txt = ",".join(f"{val:.4f}" for val in aggregate_dict[signal_name][region_name])
                    handle.write(f"{signal_name}\t{region_name}\t{agg_txt}\n")

    if args.output_aggregated_signals is not None:
        with open(args.output_aggregated_signals, "w") as handle:
            header = ["pos"]
            combos = []
            for signal_name in signal_names:
                for region_name in region_names:
                    header.append(f"{signal_name}___{region_name}")
                    combos.append((signal_name, region_name))
            handle.write(",".join(header) + "\n")

            flank = int(args.width / 2.0)
            xvals_full = np.arange(-flank, flank + 1)
            xvals_positions = np.delete(xvals_full, flank)

            for idx, pos in enumerate(xvals_positions):
                row_vals = [str(pos)]
                for signal_name, region_name in combos:
                    row_vals.append(f"{aggregate_dict[signal_name][region_name][idx]:.6f}")
                handle.write(",".join(row_vals) + "\n")

    if args.output_aggregated_scores is not None:
        with open(args.output_aggregated_scores, "w") as handle:
            handle.write(",".join(["TFBS"] + signal_names) + "\n")
            for region_name in region_names:
                row_vals = [region_name]
                for signal_name in signal_names:
                    row_vals.append(f"{aggregated_fp_scores[(signal_name, region_name)]:.6f}")
                handle.write(",".join(row_vals) + "\n")

    if args.output_aggregated_stats is not None:
        with open(args.output_aggregated_stats, "w") as handle:
            header = [
                "condition", "regions", "n_replicates", "n_sites", "mean_profile",
                "mean_profile_sd", "mean_flank", "mean_center", "aggregate_fp_score",
                "site_normalization", "mean_depletion", "depletion_se", "depletion_z",
                "depletion_effect_size", "shoulder_asymmetry", "detectability",
            ]
            handle.write(",".join(header) + "\n")
            for row in aggregate_stats:
                handle.write(",".join(str(row[col]) for col in header) + "\n")

    logger.comment("")
    logger.info("---- Plotting aggregates ----")

    if control_label is not None:
        all_plots = [
            (region_name, sig)
            for region_name in region_names
            for sig in signal_names
            if sig != control_label
        ]
        suptitle_text = f"Comparison vs {control_label}"
    else:
        all_plots = list(itertools.product(signal_names, region_names))
        suptitle_text = args.title

    total_panels = len(all_plots)

    def panel_scale_group(panel_index, combo):
        if args.share_y == "both":
            return "all"
        if control_label is not None:
            region_name, other_signal = combo
            if args.share_y == "signals":
                return ("region", region_name)
            if args.share_y == "sites":
                return ("signal", other_signal)
        else:
            signal_name, region_name = combo
            if args.share_y == "signals":
                return ("region", region_name)
            if args.share_y == "sites":
                return ("signal", signal_name)
        return ("panel", panel_index)

    scale_values = {}
    panel_scale_groups = []
    for panel_index, combo in enumerate(all_plots):
        group = panel_scale_group(panel_index, combo)
        panel_scale_groups.append(group)
        if control_label is not None:
            region_name, other_signal = combo
            values = [aggregate_dict[control_label][region_name], aggregate_dict[other_signal][region_name]]
        else:
            signal_name, region_name = combo
            values = [aggregate_dict[signal_name][region_name]]
        scale_values.setdefault(group, []).extend(values)
    scale_limits = {}
    for group, values in scale_values.items():
        finite = np.concatenate(values)
        finite = finite[np.isfinite(finite)]
        if len(finite):
            lower, upper = float(np.min(finite)), float(np.max(finite))
            value_range = upper - lower
            pad = 0.05 * value_range if value_range > 0 else 0.1
            scale_limits[group] = (lower - pad, upper + pad)
        else:
            scale_limits[group] = (-0.1, 0.1)

    grid_spec = None
    if args.grid is not None:
        try:
            grid_spec = parse_grid_spec(args.grid)
        except ValueError as exc:
            logger.error(f"ERROR: {exc}")
            sys.exit(1)

    if grid_spec is None:
        n_cols = min(total_panels, DEFAULT_GRID_COLS)
        n_rows = int(np.ceil(total_panels / DEFAULT_GRID_COLS))
    else:
        n_rows, n_cols = grid_spec
        if n_rows * n_cols < total_panels:
            logger.error(
                f"ERROR: grid {args.grid} only has room for {n_rows * n_cols} panels, "
                f"but {total_panels} are required."
            )
            sys.exit(1)
    logger.info(f"Arranging {total_panels} panels into {n_rows} rows x {n_cols} columns")

    fig, axarr = plt.subplots(
        n_rows,
        n_cols,
        figsize=(n_cols * PANEL_SIZE_IN, n_rows * PANEL_SIZE_IN),
        sharex=True,
        sharey=False,
        constrained_layout=True,
    )

    if n_rows == 1 and n_cols == 1:
        ax_matrix = np.array([[axarr]])
    elif n_rows == 1:
        ax_matrix = np.array([axarr])
    elif n_cols == 1:
        ax_matrix = np.array([[a] for a in axarr])
    else:
        ax_matrix = axarr

    flank = int(args.width / 2.0)
    xvals_full = np.arange(-flank, flank + 1)
    xvals = np.delete(xvals_full, flank)

    for idx, combo in enumerate(all_plots):
        row = idx // n_cols
        col = idx % n_cols
        ax = ax_matrix[row, col]
        ax.set_box_aspect(1)

        if control_label is not None:
            region_name, other_signal = combo
            num_sites = len(regions_dict[region_name])

            ax.plot(xvals, aggregate_dict[control_label][region_name], color="black", linewidth=1, label=control_label, zorder=1)
            if args.show_site_ci and np.isfinite(aggregate_ci_dict[control_label][region_name]).any():
                ci = aggregate_ci_dict[control_label][region_name]
                mean = aggregate_dict[control_label][region_name]
                ax.fill_between(xvals, mean - ci, mean + ci, color="black", alpha=0.10, linewidth=0, zorder=0)
            if args.show_replicate_sd and np.isfinite(aggregate_sd_dict[control_label][region_name]).any():
                sd = aggregate_sd_dict[control_label][region_name]
                mean = aggregate_dict[control_label][region_name]
                ax.fill_between(xvals, mean - sd, mean + sd, color="black", alpha=0.15, linewidth=0, zorder=0)
            ax.plot(xvals, aggregate_dict[other_signal][region_name], color="tab:red", linewidth=1, label=other_signal, zorder=2)
            if args.show_site_ci and np.isfinite(aggregate_ci_dict[other_signal][region_name]).any():
                ci = aggregate_ci_dict[other_signal][region_name]
                mean = aggregate_dict[other_signal][region_name]
                ax.fill_between(xvals, mean - ci, mean + ci, color="tab:red", alpha=0.10, linewidth=0, zorder=1)
            if args.show_replicate_sd and np.isfinite(aggregate_sd_dict[other_signal][region_name]).any():
                sd = aggregate_sd_dict[other_signal][region_name]
                mean = aggregate_dict[other_signal][region_name]
                ax.fill_between(xvals, mean - sd, mean + sd, color="tab:red", alpha=0.15, linewidth=0, zorder=1)
            ax.set_ylabel(region_name, fontsize=PDF_FONT_SIZE, fontweight="bold")
            ax.set_xlabel("bp from center", fontsize=PDF_FONT_SIZE, fontweight="bold")
            ax.text(0.98, 0.98, str(num_sites), transform=ax.transAxes, fontsize=PDF_FONT_SIZE, fontweight="bold", va="top", ha="right")
            ax.legend(loc="lower right", fontsize=PDF_FONT_SIZE, frameon=False)
        else:
            signal_name, region_name = combo
            num_sites = len(regions_dict[region_name])

            ax.plot(xvals, aggregate_dict[signal_name][region_name], color="tab:blue", linewidth=1)
            if args.show_site_ci and np.isfinite(aggregate_ci_dict[signal_name][region_name]).any():
                ci = aggregate_ci_dict[signal_name][region_name]
                mean = aggregate_dict[signal_name][region_name]
                ax.fill_between(xvals, mean - ci, mean + ci, color="tab:blue", alpha=0.16, linewidth=0)
            if args.show_replicate_sd and np.isfinite(aggregate_sd_dict[signal_name][region_name]).any():
                sd = aggregate_sd_dict[signal_name][region_name]
                mean = aggregate_dict[signal_name][region_name]
                ax.fill_between(xvals, mean - sd, mean + sd, color="tab:blue", alpha=0.18, linewidth=0)
            ax.text(0.98, 0.98, str(num_sites), transform=ax.transAxes, fontsize=PDF_FONT_SIZE, fontweight="bold", va="top", ha="right")

            if args.plot_boundaries:
                mw = motif_widths[region_name]
                mstart = -np.floor(mw / 2.0)
                mend = np.ceil(mw / 2.0) - 1
                ax.axvline(mstart, color="grey", linestyle="dashed", linewidth=1)
                ax.axvline(mend, color="grey", linestyle="dashed", linewidth=1)

            ax.set_title(signal_name, fontsize=PDF_FONT_SIZE, fontweight="bold")
            ax.set_ylabel(region_name, fontsize=PDF_FONT_SIZE, fontweight="bold")
            ax.set_xlabel("bp from center", fontsize=PDF_FONT_SIZE, fontweight="bold")

            if args.shape_diagnostics:
                diagnostics = aggregate_stats_lookup[(signal_name, region_name)]
                label = diagnostics["detectability"].replace("_", " ")
                zscore = diagnostics["depletion_z"]
                ztext = f"z={zscore:.1f}" if np.isfinite(zscore) else "z=NA"
                ax.text(
                    0.02, 0.98, f"{label}; {ztext}", transform=ax.transAxes,
                    fontsize=max(6, PDF_FONT_SIZE - 1), fontweight="bold", va="top", ha="left",
                )

        ax.set_xlim(-flank, flank)
        ax.set_ylim(*scale_limits[panel_scale_groups[idx]])
        ax.xaxis.set_major_formatter(ascii_tick_formatter())
        ax.yaxis.set_major_formatter(ascii_tick_formatter(decimals=2))
        ax.tick_params(axis="x", labelbottom=True)
        ax.tick_params(axis="y", labelleft=True)

    for extra_idx in range(total_panels, n_rows * n_cols):
        r = extra_idx // n_cols
        c = extra_idx % n_cols
        ax_matrix[r, c].axis("off")

    fig.suptitle(suptitle_text, fontsize=PDF_FONT_SIZE, fontweight="bold")
    plt.savefig(args.output, bbox_inches="tight")
    plt.close()
    if args.multiscale_npz is not None:
        logger.info("Plotting multiscale aggregate sidecar")
        plot_multiscale_aggregate_npz(
            args.multiscale_npz,
            args.output_multiscale_aggregate,
            title=f"{args.title} multiscale",
        )
    logger.end()


def main(argv=None):
    parser = add_aggregate_arguments(argparse.ArgumentParser())
    args = parser.parse_args(argv)
    if args.motif_grid:
        if len(args.input_html) != 1:
            parser.error("--motif-grid requires exactly one --input-html review report")
        from fp_tools.tools.motif_aggregate_grid import run_grid

        grid_args = argparse.Namespace(
            input_html=args.input_html[0],
            order_htmls=args.order_htmls,
            outdir=args.outdir,
            layout=args.layout,
            output=args.output,
            source_tsv=args.output_txt,
            rows_per_page=args.rows_per_page,
            flank=args.flank,
            fill_missing_profiles=args.fill_missing_profiles,
            recompute_missing_profiles=args.recompute_missing_profiles,
            cores=args.cores,
            repeat_column_labels=args.repeat_column_labels,
            title=args.title,
        )
        return run_grid(grid_args, parser)
    run_aggregate(args)
    return 0


if __name__ == "__main__":
    main()
