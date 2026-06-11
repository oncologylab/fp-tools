#!/usr/bin/env python

"""
PlotAggregate command for aggregate signal visualization from TFBS and bigWigs.

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

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyBigWig
from sklearn import preprocessing

from fp_tools.parsers import add_aggregate_arguments
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.multiscale import aggregate_multiscale_tensor, load_multiscale_npz
from fp_tools.utils.normalization import fit_quantile_normalizers
from fp_tools.utils.plotting_style import PDF_FONT_SIZE, apply_pdf_style, ascii_tick_formatter
from fp_tools.utils.regions import OneRegion, RegionList
from fp_tools.utils.signals import fast_rolling_math
from fp_tools.utils.utilities import check_files, check_required, make_directory

PANEL_SIZE_IN = 2.55
DEFAULT_GRID_COLS = 5


def default_multiscale_output(output_path):
    """Derive the companion multiscale figure path from the main aggregate path."""

    root, ext = os.path.splitext(output_path)
    if not ext:
        ext = ".pdf"
    return f"{root}_multiscale{ext}"


def plot_multiscale_aggregate_npz(npz_path, output_path, title="Multiscale aggregate"):
    """Render a PlotAggregate companion figure from a multiscale NPZ sidecar."""

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
    """Apply shared BINDetect-style quantile normalization to per-region signal arrays."""

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


def calculate_group_aggregates(signal_dict, regions_dict, region_names, condition_names, condition_groups, motif_widths, args):
    """Calculate aggregate mean profiles and optional replicate SD profiles."""

    sample_aggregates = {sample: {} for samples in condition_groups.values() for sample in samples}
    for sample_name in sample_aggregates:
        for region_name in region_names:
            signalmat = np.array([signal_dict[sample_name][reg.tup()] for reg in regions_dict[region_name]])
            if signalmat.shape[0] == 0:
                aggregate = np.zeros(args.width)
            else:
                max_values = np.max(signalmat, axis=1)
                upper_limit = np.percentile(max_values, [100 * args.remove_outliers])[0]
                signalmat = signalmat[max_values <= upper_limit]
                if signalmat.shape[0] == 0:
                    aggregate = np.zeros(args.width)
                else:
                    if args.log_transform:
                        signal_mat_abs = np.abs(signalmat)
                        signal_mat_log = np.log2(signal_mat_abs + 1)
                        signal_mat_log[signalmat < 0] *= -1
                        signalmat = signal_mat_log
                    aggregate = np.nanmean(signalmat, axis=0)
                    if args.normalize:
                        aggregate = preprocessing.minmax_scale(aggregate)
                    if args.smooth > 1:
                        agg_ext = np.pad(aggregate, args.smooth, "edge")
                        agg_smooth = fast_rolling_math(agg_ext.astype("float64"), args.smooth, "mean")
                        aggregate = agg_smooth[args.smooth:-args.smooth]
            sample_aggregates[sample_name][region_name] = aggregate

    aggregate_dict = {cond: {} for cond in condition_names}
    aggregate_sd_dict = {cond: {} for cond in condition_names}
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
            flank_len = args.flank
            motif_w = motif_widths[region_name]
            mid_start = max(0, int(flank_len - np.floor(motif_w / 2.0)))
            mid_end = min(len(mean_profile), int(flank_len + np.ceil(motif_w / 2.0)))
            mean_mid = np.mean(mean_profile[mid_start:mid_end]) if mid_end - mid_start > 0 else 0.0
            flank_indices = list(range(0, mid_start)) + list(range(mid_end, len(mean_profile)))
            mean_flank = np.mean(mean_profile[flank_indices]) if flank_indices else 0.0
            depletion = max(0.0, mean_flank - mean_mid)
            aggregated_fp_scores[(cond, region_name)] = mean_flank + depletion
            stats_rows.append({
                "condition": cond,
                "regions": region_name,
                "n_replicates": len(samples),
                "mean_profile": float(np.nanmean(mean_profile)),
                "mean_profile_sd": float(np.nanmean(sd_profile)) if np.isfinite(sd_profile).any() else np.nan,
                "mean_flank": float(mean_flank),
                "mean_center": float(mean_mid),
                "aggregate_fp_score": float(aggregated_fp_scores[(cond, region_name)]),
            })
    return aggregate_dict, aggregate_sd_dict, aggregated_fp_scores, stats_rows


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

def run_aggregate(args):
    """Create aggregate plots and optional aggregate exports."""

    # Import lazily so PlotAggregate --help does not trigger pybedtools/genomepy
    # cache initialization before argument parsing.
    import pybedtools as pb

    apply_pdf_style()
    logger = FpToolsLogger("PlotAggregate", args.verbosity)
    logger.begin()

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

            overlap = pb.BedTool(tfbs_f).intersect(pb.BedTool(region_f), u=True)
            name = args.TFBS_labels[tfbs_idx] + " <OVERLAPPING> " + args.region_labels[region_idx]
            region_names.append(name)
            regions_dict[name] = RegionList().from_bed(overlap.fn)

            if args.negate:
                overlap_neg = pb.BedTool(tfbs_f).intersect(pb.BedTool(region_f), v=True)
                name = args.TFBS_labels[tfbs_idx] + " <NOT OVERLAPPING> " + args.region_labels[region_idx]
                region_names.append(name)
                regions_dict[name] = RegionList().from_bed(overlap_neg.fn)
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
            sites = pb.BedTool(regions_dict[regions_id].as_bed(), from_string=True)
            logger.stats(f"Found {len(regions_dict[regions_id])} sites in {regions_id}")

            if len(args.whitelist) > 0:
                for whitelist_f in args.whitelist:
                    sites = sites.intersect(pb.BedTool(whitelist_f), u=True)
                    logger.stats(f"Overlapped to whitelist -> {len(sites)}")

            if len(args.blacklist) > 0:
                for blacklist_f in args.blacklist:
                    sites = sites.intersect(pb.BedTool(blacklist_f), v=True)
                    logger.stats(f"Removed blacklist -> {len(sites)}")

            regions_dict[regions_id] = RegionList().from_bed(sites.fn)

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
    raw_aggregate_dict, raw_aggregate_sd_dict, _, _ = calculate_group_aggregates(
        signal_dict, regions_dict, region_names, signal_names, condition_groups, motif_widths, args
    )
    normalized_signal_dict = apply_quantile_normalization_to_signal_dict(
        signal_dict, region_names, regions_dict, sample_names, signal_names, condition_groups, args.normalization, logger
    )
    aggregate_dict, aggregate_sd_dict, aggregated_fp_scores, aggregate_stats = calculate_group_aggregates(
        normalized_signal_dict, regions_dict, region_names, signal_names, condition_groups, motif_widths, args
    )

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

    all_values = np.concatenate([
        aggregate_dict[sig][reg]
        for sig in signal_names
        for reg in region_names
    ])
    y_min_global = np.nanmin(all_values)
    y_max_global = np.nanmax(all_values)
    y_range = y_max_global - y_min_global
    pad = 0.05 * y_range if np.isfinite(y_range) and y_range > 0 else 0.1
    y_min_global -= pad
    y_max_global += pad

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
            header = ["condition", "regions", "n_replicates", "mean_profile", "mean_profile_sd", "mean_flank", "mean_center", "aggregate_fp_score"]
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
        sharey=True,
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
            if args.show_replicate_sd and np.isfinite(aggregate_sd_dict[control_label][region_name]).any():
                sd = aggregate_sd_dict[control_label][region_name]
                mean = aggregate_dict[control_label][region_name]
                ax.fill_between(xvals, mean - sd, mean + sd, color="black", alpha=0.15, linewidth=0, zorder=0)
            ax.plot(xvals, aggregate_dict[other_signal][region_name], color="tab:red", linewidth=1, label=other_signal, zorder=2)
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

        ax.set_xlim(-flank, flank)
        if idx == 0:
            ax.set_ylim(y_min_global, y_max_global)
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


def main():
    parser = add_aggregate_arguments(argparse.ArgumentParser())
    args = parser.parse_args()
    run_aggregate(args)


if __name__ == "__main__":
    main()
