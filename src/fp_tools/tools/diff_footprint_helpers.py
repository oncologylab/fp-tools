#!/usr/bin/env python
"""
Helper functions for diff-footprints scoring, summaries, and output generation.

This module contains reusable routines for:
- score normalization
- per-motif result summaries
- static PDF plotting
- self-contained interactive HTML volcano reports
"""

import base64
import gzip
import html
import json
import random
import itertools
from datetime import datetime
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import scipy
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

# Plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from adjustText import adjust_text  # noqa: F401

# Bio
from fp_tools.utils import bigwig as pyBigWig
from fp_tools.utils.fasta import open_fasta

# Internal (fp_tools namespace)
from fp_tools.utils.regions import *
# from fp_tools.utils.utilities import fast_rolling_math, merge_dicts, file_writer
from fp_tools.utils.motifs import *
from fp_tools.utils.signals import *
from fp_tools.utils.utilities import show_worker_progress
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.normalization import ArrayNorm, fit_quantile_normalizers
from fp_tools.utils.plotting_style import PDF_FONT_SIZE, apply_pdf_style, apply_ascii_minus_to_figure

# bump open-file limit
try:
    import resource

    def bump_nofile_limit(target=4096):
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        # only raise soft up to the hard limit
        new_soft = min(int(target), int(hard))
        if soft < new_soft:
            resource.setrlimit(resource.RLIMIT_NOFILE, (new_soft, hard))

    bump_nofile_limit(4096)
except (ImportError, ValueError):
    # resource not available (e.g. non-Unix) or call failed – just skip
    pass

apply_pdf_style()


def _finite_scale(*standard_deviations):
    """Return a finite positive scale for standardized mean differences."""

    values = [abs(float(value)) for value in standard_deviations if np.isfinite(value)]
    scale = float(np.mean(values)) if values else 0.0
    return max(scale, float(np.finfo(float).eps))


def dict_to_tab(dict_list, fname, chosen_columns, header=False):
    out_str = ("\t".join(chosen_columns) + "\n") if header else ""
    out_str += "\n".join(["\t".join([str(line[c]) for c in chosen_columns]) for line in dict_list])
    out_str += "\n" if out_str else ""
    with open(fname, "w") as f:
        f.write(out_str)


def quantile_normalization(list_of_arrays, names, pdfpages=None, logger=FpToolsLogger()):
    norm_objects, diagnostics = fit_quantile_normalizers(list_of_arrays, names, logger=logger)
    array_quantiles = diagnostics["array_quantiles"]
    mean_array_quantiles = diagnostics["mean_array_quantiles"]

    if pdfpages is not None:
        fig, ax = plt.subplots()
        for i in range(len(names)):
            plt.plot(array_quantiles[i], mean_array_quantiles, label=f"Quantiles for '{names[i]}'")
        plt.title("Quantile-quantile plot", fontsize=PDF_FONT_SIZE, fontweight="bold")
        ax.set_xlabel("Value quantiles"); ax.set_ylabel("Mean quantiles")
        ax.plot([0, 1], [0, 1], transform=ax.transAxes, linestyle="dashed", color="black", label="Expected")
        plt.legend()
        apply_ascii_minus_to_figure(fig)
        pdfpages.savefig(fig, bbox_inches='tight'); plt.close()

    for i, bigwig in enumerate(names):
        xdata = array_quantiles[i]
        ydata = np.divide(mean_array_quantiles, xdata, out=np.ones_like(mean_array_quantiles), where=~np.isclose(xdata, 0.0))

        fig, ax = plt.subplots(nrows=2, ncols=1, constrained_layout=True)
        ax[0].set_xlabel("Original value"); ax[0].set_ylabel("Multiplication factor")
        ax[0].set_title(f"Multiplication needed for normalization of '{bigwig}'", fontsize=PDF_FONT_SIZE, fontweight="bold")
        ax[0].plot(xdata, ydata, color="black", linewidth=3, label="Original")
        ax[0].plot(xdata, norm_objects[bigwig].get_norm_factor(xdata), label="Norm function")
        ax[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))

        arr = np.sort(np.asarray(list_of_arrays[i], dtype=float))
        normalized = norm_objects[bigwig].normalize(arr)
        ax[1].set_title("Normalized vs. original", fontsize=PDF_FONT_SIZE, fontweight="bold")
        ax[1].plot(arr, normalized)
        ax[1].set_xlabel("Original"); ax[1].set_ylabel("Normalized values")
        max_lim = max(ax[1].get_xlim()[1], ax[1].get_ylim()[1])
        ax[1].set_xlim(0, max_lim); ax[1].set_ylim(0, max_lim)
        ax[1].plot([0, 1], [0, 1], transform=ax[1].transAxes, ls="--", color="grey")
        ax[1].grid()
        if pdfpages is not None:
            apply_ascii_minus_to_figure(fig)
            pdfpages.savefig(fig, bbox_inches="tight")
        plt.close()

    return norm_objects
def plot_score_distribution(list_of_arr, labels=None, title="Score distribution"):
    labels = labels or [f"arr_{i}" for i in range(len(list_of_arr))]
    fig, ax = plt.subplots(1, 1)
    xlim = []
    for i, arr in enumerate(list_of_arr):
        values = np.array(arr)
        x_max = np.percentile(values, [99])
        values = values[values < x_max]
        xlim.append(x_max)
        plt.hist(values, bins=100, alpha=.4, density=True, label=labels[i])
    ax.set_xlabel("Scores"); ax.set_ylabel("Density")
    ax.set_xlim(0, min(xlim))
    plt.legend(); plt.title(title, fontsize=PDF_FONT_SIZE, fontweight="bold")
    return fig


def get_gc_content(regions, fasta):
    """Mean GC fraction inside regions."""
    nuc_count = {"T": 0, "t": 0, "A": 0, "a": 0, "G": 1, "g": 1, "C": 1, "c": 1}
    gc = 0; total = 0
    fasta_obj = open_fasta(fasta)
    for region in regions:
        seq = fasta_obj.fetch(region.chrom, region.start, region.end)
        gc += sum([nuc_count.get(nuc, 0.5) for nuc in seq])
        total += region.end - region.start
    fasta_obj.close()
    return gc / float(total)


# ----------------------------------------------------------------------------- #
def scan_and_score(regions, motifs_obj, args, log_q, qs):
    """Scan motifs in regions, pull per-condition signals (averaging replicates), enqueue TFBS lines."""
    logger = FpToolsLogger("", args.verbosity, log_q)
    logger.debug("Setting up scanner/bigwigs/fasta")
    motifs_obj.setup_moods_scanner()

    # open all bigwigs as individual samples; repeated condition names define replicate groups
    sample_bigwigs = {}
    signal_to_sample = {}
    for condition, rep_idxs in args.cond_groups.items():
        files = [args.signals[i] for i in rep_idxs]
        logger.debug(f"[scan_and_score] Condition '{condition}' -> opening {files}")
        for rep_no, signal_idx in enumerate(rep_idxs, start=1):
            signal_sample_names = getattr(args, "signal_sample_names", None)
            sample_name = (
                signal_sample_names[signal_idx]
                if signal_sample_names and signal_idx < len(signal_sample_names)
                else f"{condition}_rep{rep_no}"
            )
            sample_bigwigs[sample_name] = pyBigWig.open(args.signals[signal_idx], "rb")
            signal_to_sample[signal_idx] = sample_name

    fasta_obj = open_fasta(args.genome)
    chrom_boundaries = dict(zip(fasta_obj.references, fasta_obj.lengths))

    rand_window = 200
    background_signal = {
        "keys": [],
        "gc": [],
        "signal": {c: [] for c in args.cond_names},
        "sample_signal": {s: [] for s in args.sample_names},
    }

    logger.debug("Scanning for motif occurrences")
    all_TFBS = {motif.prefix: RegionList() for motif in motifs_obj}

    # progress bar over regions (per worker)
    total_regions = len(regions)
    if tqdm is not None and show_worker_progress(args.verbosity, total_regions):
        region_iter = enumerate(
            tqdm(
                regions,
                total=total_regions,
                desc=f"scan_and_score pid={os.getpid()}",
                unit="region",
                leave=False,
            )
        )
    else:
        region_iter = enumerate(regions)

    for i, region in region_iter:
        logger.spam(f"Processing region: {region.tup()}")

        if region.end > chrom_boundaries[region.chrom]:
            logger.error(
                f"Region {region} beyond chromosome boundaries ({region.chrom}: {chrom_boundaries[region.chrom]})")
            raise Exception

        reglen = region.get_length()
        random.seed(reglen)
        rand_positions = random.sample(range(reglen), max(1, int(reglen / rand_window)))
        logger.spam(f"Random indices: {rand_positions} for len {reglen}")
        for pos in rand_positions:
            background_signal["keys"].append([region.chrom, str(region.start), str(region.end), str(pos)])

        # read signals for all samples, then summarize replicate groups per condition
        sample_footprints = {}
        for sample_name in args.sample_names:
            bw = sample_bigwigs[sample_name]
            arr = region.get_signal(bw, logger=logger, key=sample_name)
            if len(arr) == 0:
                logger.error(f"Error reading signal for '{sample_name}' in region {region}")
                raise Exception
            sample_footprints[sample_name] = arr
            for pos in rand_positions:
                background_signal["sample_signal"][sample_name].append(arr[pos])

        footprints = {}
        for condition in args.cond_names:
            rep_signals = [sample_footprints[sample_name] for sample_name in args.condition_samples[condition]]
            stacked = np.vstack(rep_signals)
            footprints[condition] = np.mean(stacked, axis=0)
            logger.spam(
                f"[scan_and_score] region {i} '{condition}': "
                f"averaged {len(rep_signals)} reps -> len {footprints[condition].shape[0]}"
            )
            for pos in rand_positions:
                background_signal["signal"][condition].append(footprints[condition][pos])

        # scan DNA sequence for motif occurrences
        seq = fasta_obj.fetch(region.chrom, region.start, region.end)
        region_TFBS = motifs_obj.scan_sequence(seq, region)

        # extend lines with peak columns and condition scores
        extra_columns = region
        for TFBS in region_TFBS:
            motif_len = TFBS.end - TFBS.start
            pos = TFBS.start - region.start + int(motif_len / 2.0)
            TFBS.extend(extra_columns)
            for sample_name in args.sample_names:
                score = sample_footprints[sample_name][pos]
                TFBS.append(f"{score:.5f}")

        for TFBS in region_TFBS:
            all_TFBS[TFBS.name].append(TFBS)

    global_TFBS = RegionList()
    for name in all_TFBS:
        all_TFBS[name] = all_TFBS[name].resolve_overlaps()
        bed_content = all_TFBS[name].as_bed()
        qs[name].put((name, bed_content))
        global_TFBS.extend(all_TFBS[name])
        all_TFBS[name] = []

    overlap = global_TFBS.count_overlaps()

    fasta_obj.close()
    for bw in sample_bigwigs.values():
        bw.close()

    logger.stop()
    logger.debug(f"Done: 'scan_and_score' finished for this chunk (time elapsed: {logger.total_time})")
    return (background_signal, overlap)


def process_tfbs(TF_name, args, log2fc_params, bed_rows=None):
    """Split into bound/unbound, write per-TF BED/overview, return TF summary row."""
    logger = FpToolsLogger("", args.verbosity, args.log_q)
    write_motif_outputs = bool(getattr(args, "write_motif_outputs", True))
    write_cache_motif_all = bool(getattr(args, "write_cache_motif_all", False))
    aggregate_site_set = (getattr(args, "aggregate_site_set", "all") or "all").replace("_", "-")
    write_aggregate_sites = bool(
        not write_motif_outputs
        and getattr(args, "aggregate_signals", None)
        and getattr(args, "plot_aggregate", "off") != "off"
    )

    tmp_root = getattr(args, "tmp_tfbs_root", None) or args.outdir
    bed_outdir = os.path.join(tmp_root, TF_name, "beds")
    filename = os.path.join(bed_outdir, TF_name + ".tmp")
    tmp_files = [] if bed_rows is not None else [filename]
    no_cond = len(args.cond_names)
    comparisons = args.comparisons
    diff_dist = scipy.stats.norm

    if args.output_peaks is not None and bed_rows is None:
        from fp_tools.utils.intervals import intersect_bed

        intersection_path = filename + ".output_peaks.bed"
        intersect_bed(filename, args.output_peaks, intersection_path)
        filename = intersection_path
        tmp_files.append(intersection_path)

    stime = datetime.now()
    header = ["TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_name", "TFBS_score", "TFBS_strand"] \
             + args.peak_header_list \
             + [f"{sample}_score" for sample in args.sample_names]
    if bed_rows is None:
        with open(filename) as f:
            bedlines = [dict(zip(header, line.rstrip().split("\t"))) for line in f.readlines()]
    else:
        bedlines = [dict(zip(header, row)) for row in bed_rows]
    n_rows = len(bedlines)
    logger.spam(f"{TF_name} - Reading took: {datetime.now() - stime}")
    if n_rows == 0:
        logger.warning(f"No TFBS found for TF {TF_name} - outputs will be empty (xlsx skipped).")

    # local: normalize, aggregate replicates, threshold, delta/log2fc
    stime = datetime.now()
    bedlines = sorted(bedlines, key=lambda line: (line["TFBS_chr"], int(line["TFBS_start"]), int(line["TFBS_end"])))
    for line in bedlines:
        for sample_name in args.sample_names:
            line[sample_name + "_score"] = float(line[sample_name + "_score"])
            if args.normalization == "sample-quantile":
                val = args.norm_objects[sample_name].normalize(line[sample_name + "_score"])
            elif args.normalization == "condition-quantile":
                cond = args.sample_to_condition[sample_name]
                val = args.norm_objects[cond].normalize(line[sample_name + "_score"])
            else:
                val = line[sample_name + "_score"]
            line[sample_name + "_score"] = round(max(0.0, float(val)), 5)

        for condition in args.cond_names:
            threshold = args.thresholds[condition]
            rep_values = np.array([line[sample + "_score"] for sample in args.condition_samples[condition]], dtype=float)
            mean_score = float(np.mean(rep_values)) if len(rep_values) else np.nan
            sd_score = float(np.std(rep_values, ddof=1)) if len(rep_values) > 1 else np.nan
            line[condition + "_score"] = round(mean_score, 5)
            line[condition + "_score_sd"] = round(sd_score, 5) if np.isfinite(sd_score) else "NA"
            line[condition + "_bound"] = 1 if line[condition + "_score"] > threshold else 0

        for (cond1, cond2) in comparisons:
            base = f"{cond1}_{cond2}"
            line[base + "_delta_fp"] = round(line[cond1 + "_score"] - line[cond2 + "_score"], 5)
            line[base + "_log2fc"] = round(np.log2((line[cond1 + "_score"] + args.pseudo) /
                                                   (line[cond2 + "_score"] + args.pseudo)), 5)

    condition_columns = [f"{cond}_score" for cond in args.cond_names]
    condition_sd_columns = [f"{cond}_score_sd" for cond in args.cond_names]
    overview_columns = header + condition_columns + condition_sd_columns + [c + "_bound" for c in args.cond_names] \
                       + [f"{c1}_{c2}_delta_fp" for (c1, c2) in comparisons] \
                       + [f"{c1}_{c2}_log2fc" for (c1, c2) in comparisons]

    bed_table = pd.DataFrame(bedlines, columns=overview_columns)
    logger.spam(f"Read table {bed_table.shape} for TF {TF_name}")

    if write_motif_outputs or write_cache_motif_all or write_aggregate_sites:
        # write *_all.bed
        outfile = os.path.join(bed_outdir, TF_name + "_all.bed")
        dict_to_tab(bedlines, outfile, header + condition_columns + condition_sd_columns)

    if write_motif_outputs or (write_aggregate_sites and aggregate_site_set == "bound"):
        # write bound/unbound per condition
        for condition in args.cond_names:
            chosen_columns = header[:-len(args.sample_names)] + [condition + "_score"]
            states = ["bound", "unbound"] if write_motif_outputs else ["bound"]
            for state in states:
                chosen_bool = 1 if state == "bound" else 0
                subset = [bl for bl in bedlines if bl[condition + "_bound"] == chosen_bool]
                outfile = os.path.join(bed_outdir, f"{TF_name}_{condition}_{state}.bed")
                dict_to_tab(subset, outfile, chosen_columns)

    if write_motif_outputs:
        # overview (txt + optional xlsx)
        overview_txt = os.path.join(args.outdir, TF_name, TF_name + "_overview.txt")
        dict_to_tab(bedlines, overview_txt, overview_columns, header=True)

        if not args.skip_excel and n_rows > 0:
            try:
                overview_excel = os.path.join(args.outdir, TF_name, TF_name + "_overview.xlsx")
                with pd.ExcelWriter(overview_excel, engine='xlsxwriter') as writer:
                    bed_table.to_excel(writer, index=False, columns=overview_columns)
                    ws = writer.sheets['Sheet1']
                    n_rows_x, n_cols_x = bed_table.shape
                    ws.autofilter(0, 0, n_rows_x, n_cols_x)
            except Exception as e:
                logger.warning(f"Could not write Excel for TF {TF_name}. Exception: {e}")

    # global summary row
    info_columns = ["total_tfbs"]
    info_columns += [f"{sample}_mean_score" for sample in args.sample_names]
    info_columns += [f"{cond}_{metric}" for cond, metric in itertools.product(args.cond_names, ["mean_score", "score_sd", "n_replicates", "bound"])]
    info_columns += [f"{c1}_{c2}_{metric}" for (c1, c2), metric in itertools.product(comparisons, ["change", "pvalue", "mean_delta_fp", "mean_log2fc", "delta_fp_se", "log2fc_se"])]
    # A single-replicate workflow may deliberately use the same sample and
    # condition label. Keep one shared mean-score field in that case.
    info_columns = list(dict.fromkeys(info_columns))
    info_table = pd.DataFrame(np.nan, columns=info_columns, index=[TF_name])

    info_table.at[TF_name, "total_tfbs"] = n_rows
    for sample_name in args.sample_names:
        # Read from the row dictionaries rather than ``bed_table``. A valid
        # single-replicate invocation can use the same label for the sample and
        # condition, which gives the DataFrame duplicate ``*_score`` columns
        # and makes column selection two-dimensional.
        sample_values = pd.to_numeric(
            pd.Series([row.get(sample_name + "_score", np.nan) for row in bedlines]),
            errors="coerce",
        )
        info_table.at[TF_name, sample_name + "_mean_score"] = (
            round(float(sample_values.mean()), 8) if len(sample_values.dropna()) else np.nan
        )
    for condition in args.cond_names:
        info_table.at[TF_name, condition + "_mean_score"] = round(float(np.mean(bed_table[condition + "_score"])), 5) if n_rows > 0 else np.nan
        sd_values = pd.to_numeric(bed_table[condition + "_score_sd"], errors="coerce") if n_rows > 0 else pd.Series(dtype=float)
        info_table.at[TF_name, condition + "_score_sd"] = round(float(np.nanmean(sd_values)), 5) if len(sd_values.dropna()) else np.nan
        info_table.at[TF_name, condition + "_n_replicates"] = args.condition_replicates.get(condition, 1)
        info_table.at[TF_name, condition + "_bound"] = int(np.sum(bed_table[condition + "_bound"].values))

    # per-comparison stats and figure
    write_per_motif_plots = bool(getattr(args, "per_motif_plots", False))
    log2fc_pdf = None
    if write_per_motif_plots:
        fig_out = os.path.join(args.outdir, TF_name, "plots", TF_name + "_log2fcs.pdf")
        log2fc_pdf = PdfPages(fig_out, keep_empty=False)

    if n_rows > 0:
        for (cond1, cond2) in comparisons:
            base = f"{cond1}_{cond2}"
            included = np.logical_or(bed_table[cond1 + "_score"].values > 0, bed_table[cond2 + "_score"].values > 0)
            subset = bed_table[included].copy()
            subset.loc[:, "peak_id"] = ["_".join([chrom, str(start), str(end)])
                                        for (chrom, start, end) in zip(subset.iloc[:, 0].values,
                                                                       subset.iloc[:, 1].values,
                                                                       subset.iloc[:, 2].values)]
            observed_log2fcs = subset.groupby('peak_id')[base + '_log2fc'].mean().reset_index()[base + "_log2fc"].values
            observed_deltas = subset.groupby('peak_id')[base + '_delta_fp'].mean().reset_index()[base + "_delta_fp"].values
            if len(observed_log2fcs) == 0:
                # Motif sites can exist while every compared score is zero.
                # There is no evidence of a differential footprint in that
                # case, and attempting a distribution fit would emit warnings
                # and propagate NaN values into comparison-wide thresholds.
                info_table.at[TF_name, base + "_mean_delta_fp"] = 0.0
                info_table.at[TF_name, base + "_mean_log2fc"] = 0.0
                info_table.at[TF_name, base + "_change"] = 0.0
                info_table.at[TF_name, base + "_pvalue"] = 1.0
                continue

            info_table.at[TF_name, base + "_mean_delta_fp"] = np.round(float(np.mean(observed_deltas)), 5)
            info_table.at[TF_name, base + "_mean_log2fc"] = np.round(float(np.mean(observed_log2fcs)), 5)
            n1 = max(1, args.condition_replicates.get(cond1, 1))
            n2 = max(1, args.condition_replicates.get(cond2, 1))
            sd1 = pd.to_numeric(subset[cond1 + "_score_sd"], errors="coerce").to_numpy(dtype=float)
            sd2 = pd.to_numeric(subset[cond2 + "_score_sd"], errors="coerce").to_numpy(dtype=float)
            mu1 = pd.to_numeric(subset[cond1 + "_score"], errors="coerce").to_numpy(dtype=float)
            mu2 = pd.to_numeric(subset[cond2 + "_score"], errors="coerce").to_numpy(dtype=float)
            if np.isfinite(sd1).any() and np.isfinite(sd2).any():
                delta_se = np.sqrt(np.nanmean((sd1 ** 2) / n1 + (sd2 ** 2) / n2))
                log2fc_se = (1.0 / np.log(2.0)) * np.sqrt(np.nanmean((sd1 ** 2) / (n1 * (mu1 + args.pseudo) ** 2) + (sd2 ** 2) / (n2 * (mu2 + args.pseudo) ** 2)))
                info_table.at[TF_name, base + "_delta_fp_se"] = np.round(float(delta_se), 5)
                info_table.at[TF_name, base + "_log2fc_se"] = np.round(float(log2fc_se), 5)

            bg_mean, bg_std = log2fc_params[(cond1, cond2)]
            obs_params = scipy.stats.norm.fit(observed_log2fcs)
            obs_mean, obs_std = obs_params
            n_obs = len(observed_log2fcs)

            if obs_mean != bg_mean:
                change = (obs_mean - bg_mean) / _finite_scale(obs_std, bg_std)
                info_table.at[TF_name, base + "_change"] = np.round(change, 5)

                np.random.seed(n_obs)
                sample_changes = []
                for _ in range(100):
                    sample = scipy.stats.norm.rvs(bg_mean, bg_std, size=n_obs)
                    sm, ss = float(np.mean(sample)), float(np.std(sample))
                    sample_changes.append((sm - bg_mean) / _finite_scale(ss, bg_std))
                ttest = scipy.stats.ttest_1samp(
                    sample_changes,
                    float(info_table.at[TF_name, base + "_change"]),
                )
                pvalue = float(ttest.pvalue)
                info_table.at[TF_name, base + "_pvalue"] = pvalue if np.isfinite(pvalue) else 1.0
            else:
                info_table.at[TF_name, base + "_change"] = 0
                info_table.at[TF_name, base + "_pvalue"] = 1

            if write_per_motif_plots:
                fig, ax = plt.subplots(1, 1)
                ax.hist(observed_log2fcs, bins='auto', label="Observed log2fcs", density=True)
                xvals = np.linspace(plt.xlim()[0], plt.xlim()[1], 100)
                ax.plot(xvals, scipy.stats.norm.pdf(xvals, *obs_params), label="Observed (fit)", color="red", ls="--")
                ax.axvline(obs_mean, color="red", label="Observed mean")
                ax.plot(xvals, scipy.stats.norm.pdf(xvals, bg_mean, bg_std), label="Background (fit)", color="black", ls="--")
                ax.axvline(bg_mean, color="black", label="Background mean")
                x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
                ax.set_aspect(((x1 - x0) / (y1 - y0)) / 1.5)
                ax.legend(); plt.xlabel("Log2 fold change"); plt.ylabel("Density")
                plt.title(f"Differential binding for \"{TF_name}\"\n({cond1} / {cond2})", fontsize=PDF_FONT_SIZE, fontweight="bold")
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                plt.tight_layout()
                apply_ascii_minus_to_figure(fig)
                log2fc_pdf.savefig(fig, bbox_inches='tight'); plt.close(fig)

    if log2fc_pdf is not None:
        log2fc_pdf.close()

    # cleanup tmp unless the caller needs the temporary motif-site files for
    # compact cache construction after per-motif summaries are computed.
    if not getattr(args, "keep_tmp_tfbs_for_cache", False):
        for fn in tmp_files:
            try:
                os.remove(fn)
            except Exception:
                logger.error(f"Could not remove temporary file {fn} (harmless).")

    return info_table


# ------------------------------ plotting utils ------------------------------ #
def plot_diff_footprints(motifs, cluster_obj, conditions, args):
    import warnings as _warnings
    _warnings.filterwarnings("ignore")

    cond1, cond2 = conditions
    n_IDS = cluster_obj.n

    diff_scores = {
        m.prefix: {
            "change": float(getattr(m, "change", 0)),
            "pvalue": float(getattr(m, "pvalue", 1)),
            "log10pvalue": -np.log10(float(getattr(m, "pvalue", 1))) if float(getattr(m, "pvalue", 1)) > 0 else -np.log10(1e-308),
            "volcano_label": m.name,
            "overview_label": f"{m.name} ({m.id})",
            "group": getattr(m, "group", "n.s.")
        }
        for m in motifs
    }

    xvalues = np.array([v["change"] for v in diff_scores.values()], dtype=float)
    yvalues = np.array([v["log10pvalue"] for v in diff_scores.values()], dtype=float)
    xvalues = np.where(np.isfinite(xvalues), xvalues, 0.0)
    yvalues = np.where(np.isfinite(yvalues), yvalues, 0.0)

    y_min = np.percentile(yvalues[yvalues < -np.log10(1e-300)], 95) if (yvalues < -np.log10(1e-300)).any() else np.percentile(yvalues, 95)
    x_min, x_max = np.percentile(xvalues, [5, 95])

    for TF, v in diff_scores.items():
        if v["change"] < x_min or v["change"] > x_max or v["log10pvalue"] > y_min:
            v["show"] = True
            v["color"] = "blue" if v["change"] < 0 else ("red" if v["change"] > 0 else "black")
        else:
            v["show"] = False
            v["color"] = "black"

    node_color = cluster_obj.node_color
    IDS = np.array(cluster_obj.names)

    # Volcano plot lives in the main diff-footprints PDF
    volcano_fig, ax1 = plt.subplots(figsize=(4.0, 4.0))
    ax1.set_title("diff-footprints volcano plot", fontsize=PDF_FONT_SIZE, fontweight="bold", pad=12)
    ax1.scatter(xvalues, yvalues, color="black", s=5)
    ylim = ax1.get_ylim(); y_extra = (ylim[1] - ylim[0]) * 0.1
    ax1.set_ylim(ylim[0], ylim[1] + y_extra)
    xlim = ax1.get_xlim(); x_extra = (xlim[1] - xlim[0]) * 0.1
    lim = np.max([abs(xlim[0]-x_extra), abs(xlim[1]+x_extra)])
    ax1.set_xlim(-lim, lim)
    x0, x1 = ax1.get_xlim(); y0, y1 = ax1.get_ylim()
    ax1.set_aspect((x1 - x0) / (y1 - y0))
    ax1.set_xlabel("Differential binding score")
    ax1.set_ylabel("-log10(pvalue)")

    # Clustering/overview plot lives in a separate PDF
    l = 10 + 7 * (n_IDS / 25)
    limit = 2**16 / 100 - 1
    l = limit if l > limit else l
    cluster_fig = plt.figure(figsize=(8, l))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 1.0], figure=cluster_fig)
    gs.update(wspace=0.30, hspace=0.05, bottom=0.04, top=0.97)
    ax2 = cluster_fig.add_subplot(gs[0, 0])
    ax3 = cluster_fig.add_subplot(gs[0, 1])

    # Dendrogram
    if len(IDS) > 1:
        dendro_dat = dendrogram(
            cluster_obj.linkage_mat, labels=list(IDS), no_labels=True,
            orientation="right", ax=ax3, above_threshold_color="black",
            link_color_func=lambda k: cluster_obj.node_color[k]
        )
        labels = dendro_dat["ivl"]
        ax3.set_xlabel("TF distance (clusters colored below threshold)")
        ax3.set_ylabel("TF clustering based on TFBS overlap", rotation=270, labelpad=20)
        x0, x1 = ax3.get_xlim(); y0, y1 = ax3.get_ylim()
        ax3.set_aspect(((x1 - x0) / (y1 - y0)) * len(IDS) / 10)
    else:
        ax3.axis('off')
        labels = IDS
    ax3.axvline(x=args.cluster_threshold, linestyle="dashed", alpha=0.5, color="grey")

    # Long scatter overview
    ax2.set_xlabel("Differential binding score\n" + f"({cond2} <-> {cond1})")
    ax2.set_ylim(0.5, len(labels) + 0.5)
    ax2.set_ylabel("Transcription factors")
    ax2.set_yticks(range(1, len(labels) + 1))
    ax2.set_yticklabels([diff_scores[TF]["overview_label"] for TF in labels])
    ax2.axvline(0, color="grey", linestyle="--")
    for y, TF in enumerate(labels):
        idx = np.where(IDS == TF)[0][0]
        score = diff_scores[TF]["change"]
        fill = "full" if diff_scores[TF]["show"] else "none"
        ax2.axhline(y + 1, color="grey", linewidth=1)
        ax2.plot(score, y + 1, marker='o', color=node_color[idx], fillstyle=fill)
        ax2.yaxis.get_ticklabels()[y].set_color(node_color[idx])

    lim2 = np.max(np.abs(ax2.get_xlim()))
    ax2.set_xlim((-lim2, lim2))
    x0, x1 = ax2.get_xlim(); y0, y1 = ax2.get_ylim()
    ax2.set_aspect(((x1 - x0) / (y1 - y0)) * n_IDS / 10)

    # label/highlight volcano
    txts = []
    for TF, v in diff_scores.items():
        ax1.scatter(v["change"], v["log10pvalue"], color=v["color"], s=4.5)
        if v["show"]:
            txts.append(
                ax1.text(
                    v["change"],
                    v["log10pvalue"],
                    v["volcano_label"],
                    fontsize=PDF_FONT_SIZE,
                    fontweight="bold",
                )
            )

    if txts:
        adjust_text(
            txts,
            ax=ax1,
            expand_points=(1.3, 1.5),
            expand_text=(1.15, 1.25),
            force_points=(0.35, 0.45),
            force_text=(0.3, 0.4),
            lim=300,
        )

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor="red", label=f"Higher in {cond1}"),
        Line2D([0], [0], marker='o', color='w', markerfacecolor="blue", label=f"Higher in {cond2}"),
    ]
    ax1.legend(handles=legend_elements, loc="lower left", framealpha=0.5)

    volcano_fig.tight_layout()
    cluster_fig.tight_layout()
    apply_ascii_minus_to_figure(volcano_fig)
    apply_ascii_minus_to_figure(cluster_fig)
    return volcano_fig, cluster_fig




def _read_bed_centers(path):
    centers = []
    if not os.path.exists(path):
        return centers
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            centers.append((fields[0], (start + end) // 2))
    return centers


def _mean_profile(bigwig_path, centers, flank, norm=None):
    profiles = []
    with pyBigWig.open(bigwig_path) as bw:
        chroms = bw.chroms()
        for chrom, center in centers:
            if chrom not in chroms:
                continue
            start = center - flank
            end = center + flank
            if start < 0 or end > chroms[chrom] or end <= start:
                continue
            values = np.asarray(bw.values(chrom, start, end, numpy=True), dtype=float)
            if values.size != flank * 2:
                continue
            values = np.nan_to_num(values, nan=0.0, posinf=0.0, neginf=0.0)
            if norm is not None:
                values = norm.normalize(values)
            profiles.append(values)
    if not profiles:
        return [0.0] * (flank * 2)
    return [round(float(v), 6) for v in np.nanmean(np.vstack(profiles), axis=0)]




class AggregateAffineNorm:
    """Sign-preserving affine scaler for aggregate cut-site profiles."""

    def __init__(self, source_center, scale, target_center):
        self.source_center = float(source_center)
        self.scale = float(scale)
        self.target_center = float(target_center)

    def normalize(self, values):
        arr = np.asarray(values, dtype=float)
        return (arr - self.source_center) * self.scale + self.target_center


class AggregateSizeFactorNorm:
    """Multiplicative size-factor scaler for aggregate cut-site profiles."""

    def __init__(self, size_factor):
        self.size_factor = float(size_factor)
        if not np.isfinite(self.size_factor) or self.size_factor <= 1e-12:
            self.size_factor = 1.0

    def normalize(self, values):
        arr = np.asarray(values, dtype=float)
        return arr / self.size_factor


def _mean_positive_signal(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    arr = arr[arr > 0]
    if arr.size == 0:
        return 1.0
    mean = float(np.nanmean(arr))
    return mean if np.isfinite(mean) and mean > 1e-12 else 1.0


def _size_factor_normalizers(sample_arrays, sample_names):
    """Fit simple library-size-style factors and divide profiles by them."""

    means = [_mean_positive_signal(arr) for arr in sample_arrays]
    target = float(np.nanmean(means)) if means else 1.0
    if not np.isfinite(target) or target <= 1e-12:
        target = 1.0
    return {name: AggregateSizeFactorNorm(mean / target) for name, mean in zip(sample_names, means)}


def _aggregate_fp_score(profile):
    """Simple flank-minus-center score used only for drawing order."""

    arr = np.asarray(profile, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size < 6:
        return 0.0
    center_width = max(2, int(round(arr.size * 0.12)))
    flank_width = max(center_width, int(round(arr.size * 0.20)))
    mid = arr.size // 2
    half = center_width // 2
    center = arr[max(0, mid - half):min(arr.size, mid + half + center_width % 2)]
    flank = np.concatenate([arr[:flank_width], arr[-flank_width:]])
    if center.size == 0 or flank.size == 0:
        return 0.0
    return float(np.nanmean(flank) - np.nanmean(center))


def _aggregate_bed_paths(outdir, prefix, comparison, site_set):
    bed_dir = os.path.join(outdir, prefix, "beds")
    site_set = (site_set or "all").replace("_", "-")
    if site_set == "bound":
        return {cond: os.path.join(bed_dir, f"{prefix}_{cond}_bound.bed") for cond in comparison}
    return {cond: os.path.join(bed_dir, prefix + "_all.bed") for cond in comparison}


def _limit_aggregate_centers(centers, max_centers):
    if max_centers is None:
        return centers
    max_centers = int(max_centers)
    if max_centers <= 0 or len(centers) <= max_centers:
        return centers
    indices = np.linspace(0, len(centers) - 1, max_centers, dtype=int)
    return [centers[idx] for idx in indices]


def _cached_aggregate_centers_for_row(prefix, comparison, site_set, aggregate_site_maps, cond_groups, max_centers=None):
    site_set = (site_set or "all").replace("_", "-")

    def add_unique(centers, path):
        seen = set(centers)
        for center in _read_bed_centers(path):
            if center not in seen:
                seen.add(center)
                centers.append(center)
        return centers

    centers_by_condition = {}
    if site_set == "bound":
        for cond in comparison:
            centers = []
            for idx in cond_groups.get(cond, []):
                if idx >= len(aggregate_site_maps):
                    continue
                path = aggregate_site_maps[idx].get(prefix, {}).get("bound")
                if path:
                    centers = add_unique(centers, path)
            centers_by_condition[cond] = _limit_aggregate_centers(centers, max_centers)
    else:
        shared = []
        for sample_map in aggregate_site_maps:
            path = sample_map.get(prefix, {}).get("all")
            if path:
                shared = add_unique(shared, path)
                break
        shared = _limit_aggregate_centers(shared, max_centers)
        centers_by_condition = {cond: list(shared) for cond in comparison}

    unique = []
    seen = set()
    for centers in centers_by_condition.values():
        for center in centers:
            if center not in seen:
                seen.add(center)
                unique.append(center)
    return centers_by_condition, unique


def _aggregate_centers_for_row(
    outdir,
    prefix,
    comparison,
    site_set,
    max_centers=None,
    aggregate_site_maps=None,
    cond_groups=None,
):
    if aggregate_site_maps:
        centers_by_condition, _ = _cached_aggregate_centers_for_row(
            prefix,
            comparison,
            site_set,
            aggregate_site_maps,
            cond_groups or {},
            max_centers=max_centers,
        )
        # Project-layout summary mode intentionally leaves the reusable sample
        # folders without materialized per-motif BEDs. process_tfbs writes the
        # selected motif sites to a temporary comparison tree, so use those
        # files for any motif/condition absent from the reusable directories.
        fallback_paths = _aggregate_bed_paths(outdir, prefix, comparison, site_set)
        for cond in comparison:
            if centers_by_condition.get(cond):
                continue
            centers_by_condition[cond] = _limit_aggregate_centers(
                _read_bed_centers(fallback_paths[cond]),
                max_centers,
            )
        unique = []
        seen = set()
        for centers in centers_by_condition.values():
            for center in centers:
                if center not in seen:
                    seen.add(center)
                    unique.append(center)
        return centers_by_condition, unique
    paths = _aggregate_bed_paths(outdir, prefix, comparison, site_set)
    centers_by_condition = {cond: _limit_aggregate_centers(_read_bed_centers(path), max_centers) for cond, path in paths.items()}
    unique = []
    seen = set()
    for centers in centers_by_condition.values():
        for center in centers:
            if center not in seen:
                seen.add(center)
                unique.append(center)
    return centers_by_condition, unique


def _robust_affine_normalizers(sample_arrays, sample_names):
    """Fit robust linear scalers from sampled aggregate-track windows."""

    centers = []
    widths = []
    cleaned = []
    for arr in sample_arrays:
        values = np.asarray(arr, dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            values = np.array([0.0], dtype=float)
        q05, q50, q95 = np.nanquantile(values, [0.05, 0.5, 0.95])
        width = float(q95 - q05)
        if not np.isfinite(width) or width <= 1e-12:
            width = 1.0
        centers.append(float(q50) if np.isfinite(q50) else 0.0)
        widths.append(width)
        cleaned.append(values)
    if not cleaned:
        return {}
    target_center = float(np.nanmedian(centers)) if centers else 0.0
    target_width = float(np.nanmedian(widths)) if widths else 1.0
    if not np.isfinite(target_width) or target_width <= 1e-12:
        target_width = 1.0
    out = {}
    for name, center, width in zip(sample_names, centers, widths):
        out[name] = AggregateAffineNorm(center, target_width / width, target_center)
    return out


def _sample_bigwig_window_values(bigwig_path, centers, flank, max_values=500000):
    """Read a deterministic sample of cut-site values for report-level normalization."""

    window = flank * 2
    if not centers or window <= 0:
        return np.array([0.0], dtype=float)
    max_windows = max(1, int(max_values // window))
    if len(centers) > max_windows:
        indices = np.linspace(0, len(centers) - 1, max_windows, dtype=int)
        centers = [centers[idx] for idx in indices]
    values = []
    with pyBigWig.open(bigwig_path) as bw:
        chroms = bw.chroms()
        for chrom, center in centers:
            if chrom not in chroms:
                continue
            start = center - flank
            end = center + flank
            if start < 0 or end > chroms[chrom] or end <= start:
                continue
            arr = np.asarray(bw.values(chrom, start, end, numpy=True), dtype=float)
            if arr.size != window:
                continue
            values.append(np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0))
    if not values:
        return np.array([0.0], dtype=float)
    return np.concatenate(values)


def _fit_aggregate_normalizers(selected, outdir, aggregate_signals, cond_groups, comparison, flank, mode, site_set="all", logger=None, max_centers=None, aggregate_site_maps=None):
    """Fit one report-level normalizer for aggregate cut-site profiles.

    This uses pooled windows from all displayed motifs. Fitting a separate quantile
    curve for each motif aggregate changes that motif's within-profile rank
    structure and can create artificial footprint shapes.
    """

    mode = (mode or "none").replace("_", "-")
    sample_names = [f"sample_{idx + 1}" for idx in range(len(aggregate_signals))]
    if mode == "none" or len(sample_names) <= 1:
        return {}

    all_centers = []
    seen = set()
    for _, row in selected.iterrows():
        prefix = str(row["output_prefix"])
        _, row_centers = _aggregate_centers_for_row(
            outdir,
            prefix,
            comparison,
            site_set,
            max_centers=max_centers,
            aggregate_site_maps=aggregate_site_maps,
            cond_groups=cond_groups,
        )
        for key in row_centers:
            if key not in seen:
                seen.add(key)
                all_centers.append(key)
    if not all_centers:
        return {}

    sample_arrays = [_sample_bigwig_window_values(path, all_centers, flank) for path in aggregate_signals]
    if mode == "sample-quantile":
        return {"mode": mode, "sample": _robust_affine_normalizers(sample_arrays, sample_names)}

    if mode == "size-factor":
        return {"mode": mode, "sample": _size_factor_normalizers(sample_arrays, sample_names)}

    if mode == "condition-quantile":
        condition_arrays = []
        valid_conditions = []
        for cond in comparison:
            indices = [idx for idx in cond_groups.get(cond, []) if idx < len(sample_arrays)]
            if not indices:
                continue
            arrays = [sample_arrays[idx] for idx in indices]
            condition_arrays.append(np.concatenate(arrays) if arrays else np.array([0.0], dtype=float))
            valid_conditions.append(cond)
        if len(condition_arrays) <= 1:
            return {}
        return {"mode": mode, "condition": _robust_affine_normalizers(condition_arrays, valid_conditions)}

    raise ValueError(f"Unsupported aggregate normalization mode: {mode}")


def _normalize_aggregate_profiles(sample_profiles, sample_names, condition_names, cond_groups, mode, norm_spec=None):
    """Apply report-level aggregate signal normalizers to motif profiles."""

    mode = (mode or "none").replace("_", "-")
    sample_profiles = {name: np.asarray(profile, dtype=float) for name, profile in zip(sample_names, sample_profiles)}
    if mode == "none" or len(sample_profiles) <= 1:
        return sample_profiles
    norm_spec = norm_spec or {}

    if mode in {"sample-quantile", "size-factor"}:
        norm_objects = norm_spec.get("sample", {})
        return {
            name: norm_objects[name].normalize(profile) if name in norm_objects else profile
            for name, profile in sample_profiles.items()
        }

    if mode == "condition-quantile":
        norm_objects = norm_spec.get("condition", {})
        out = dict(sample_profiles)
        for cond in condition_names:
            norm = norm_objects.get(cond)
            if norm is None:
                continue
            for idx in cond_groups.get(cond, []):
                if idx >= len(sample_names):
                    continue
                name = sample_names[idx]
                if name in sample_profiles:
                    out[name] = norm.normalize(sample_profiles[name])
        return out

    raise ValueError(f"Unsupported aggregate normalization mode: {mode}")


def _aggregate_payload_for_row(task):
    if len(task) == 11:
        row, comparison, outdir, aggregate_signals, cond_groups, flank, x_len, base, normalization, aggregate_norm_spec, sample_names = task
        site_set = "all"
        max_centers = None
        aggregate_site_maps = None
    elif len(task) == 12:
        row, comparison, outdir, aggregate_signals, cond_groups, flank, x_len, base, normalization, aggregate_norm_spec, sample_names, site_set = task
        max_centers = None
        aggregate_site_maps = None
    elif len(task) == 13:
        row, comparison, outdir, aggregate_signals, cond_groups, flank, x_len, base, normalization, aggregate_norm_spec, sample_names, site_set, max_centers = task
        aggregate_site_maps = None
    else:
        row, comparison, outdir, aggregate_signals, cond_groups, flank, x_len, base, normalization, aggregate_norm_spec, sample_names, site_set, max_centers, aggregate_site_maps = task
    c1, c2 = comparison
    prefix = str(row["output_prefix"])
    centers_by_condition, all_centers = _aggregate_centers_for_row(
        outdir,
        prefix,
        comparison,
        site_set,
        max_centers=max_centers,
        aggregate_site_maps=aggregate_site_maps,
        cond_groups=cond_groups,
    )
    if not all_centers:
        return None

    if not sample_names or len(sample_names) != len(aggregate_signals):
        sample_names = [f"sample_{idx + 1}" for idx in range(len(aggregate_signals))]
    sample_norms = (aggregate_norm_spec or {}).get("sample", {})
    condition_norms = (aggregate_norm_spec or {}).get("condition", {})
    sample_to_condition = {idx: cond for cond, indices in cond_groups.items() for idx in indices}

    conditions = []
    for cond in (c1, c2):
        centers = centers_by_condition.get(cond, all_centers)
        sample_profiles = []
        samples = []
        for signal_idx in cond_groups.get(cond, []):
            sample_name = sample_names[signal_idx]
            norm = None
            if normalization in {"sample-quantile", "size-factor"}:
                norm = sample_norms.get(sample_name) or sample_norms.get(f"sample_{signal_idx + 1}")
            elif normalization == "condition-quantile":
                norm = condition_norms.get(sample_to_condition.get(signal_idx))
            sample_profile = np.asarray(_mean_profile(aggregate_signals[signal_idx], centers, flank, norm=norm), dtype=float)
            sample_profiles.append(sample_profile)
            samples.append({
                "name": sample_name,
                "profile": [round(float(v), 6) for v in sample_profile],
                "fp_score": round(float(_aggregate_fp_score(sample_profile)), 6),
            })
        if sample_profiles:
            mean_profile = np.nanmean(np.asarray(sample_profiles, dtype=float), axis=0)
            profile = [round(float(v), 6) for v in mean_profile]
            fp_score = round(float(_aggregate_fp_score(mean_profile)), 6)
        else:
            profile = [0.0] * x_len
            fp_score = 0.0
        conditions.append({"name": cond, "profile": profile, "samples": samples, "n_sites": len(centers), "fp_score": fp_score})
    return {
        "prefix": prefix,
        "name": str(row.get("name", prefix)),
        "motif_id": str(row.get("motif_id", "")),
        "change": float(row.get(base + "_change", 0.0)),
        "pvalue": float(row.get(base + "_pvalue_numeric", 1.0)),
        "n_sites": len(all_centers),
        "site_set": site_set,
        "max_sites_per_motif": max_centers,
        "conditions": conditions,
    }


def build_diff_footprint_aggregate_payload(motifs, info_table, comparison, args):
    """Build compact aggregate profiles for embedding in comparison HTML."""

    if not getattr(args, "aggregate_signals", None):
        return None
    if len(args.aggregate_signals) != len(args.signals):
        raise ValueError("--aggregate-signals must have the same length as --signals")
    if not 1 <= int(getattr(args, "default_aggregate_plots", 4)) <= 12:
        raise ValueError("--default-aggregate-plots must be between 1 and 12")

    c1, c2 = comparison
    base = f"{c1}_{c2}"
    rows = info_table.copy()
    rows[base + "_pvalue_numeric"] = pd.to_numeric(rows[base + "_pvalue"], errors="coerce").fillna(1.0)
    rows[base + "_abs_change"] = pd.to_numeric(rows[base + "_change"], errors="coerce").fillna(0.0).abs()
    mode = getattr(args, "plot_aggregate", "sig")
    top_n = max(1, int(getattr(args, "plot_aggregate_top_n", 20)))
    sig_only = bool(getattr(args, "aggregate_sig_only", False))
    no_fallback = bool(getattr(args, "aggregate_sig_no_fallback", False))
    select_rows = rows
    if sig_only:
        highlighted_col = base + "_highlighted"
        if highlighted_col in select_rows.columns:
            highlighted = select_rows[highlighted_col].astype(str).str.lower().isin({"true", "1", "yes"})
            select_rows = select_rows.loc[highlighted].copy()
        else:
            select_rows = select_rows.iloc[0:0].copy()
    requested_motifs = list(getattr(args, "plot_aggregate_motifs", None) or [])
    if requested_motifs:
        selected = select_aggregate_rows(rows, requested_motifs)
    elif mode == "all":
        selected = select_rows.sort_values([base + "_pvalue_numeric", base + "_abs_change"], ascending=[True, False])
    elif mode == "top":
        selected = select_rows.sort_values([base + "_pvalue_numeric", base + "_abs_change"], ascending=[True, False]).head(top_n)
    else:
        threshold = float(getattr(args, "aggregate_pvalue_threshold", 0.05))
        selected = select_rows[select_rows[base + "_pvalue_numeric"] <= threshold].sort_values([base + "_pvalue_numeric", base + "_abs_change"], ascending=[True, False])
        if not selected.empty:
            selected = selected.head(top_n)
        if selected.empty and not no_fallback:
            selected = select_rows.sort_values([base + "_pvalue_numeric", base + "_abs_change"], ascending=[True, False]).head(top_n)

    flank = max(1, int(getattr(args, "aggregate_flank", 100)))
    x = list(range(-flank, flank))
    requested_aggregate_norm = (getattr(args, "aggregate_normalization", "match") or "match").replace("_", "-")
    normalization = (getattr(args, "normalization", "none") or "none").replace("_", "-") if requested_aggregate_norm == "match" else requested_aggregate_norm
    site_set = (getattr(args, "aggregate_site_set", "all") or "all").replace("_", "-")
    max_centers = getattr(args, "aggregate_max_sites", None)
    cond_groups = {cond: list(indices) for cond, indices in getattr(args, "cond_groups", {}).items()}
    aggregate_site_maps = getattr(args, "aggregate_site_maps", None)
    aggregate_bed_root = getattr(args, "tmp_tfbs_root", None) or args.outdir
    aggregate_norm_spec = _fit_aggregate_normalizers(
        selected,
        aggregate_bed_root,
        list(args.aggregate_signals),
        cond_groups,
        (c1, c2),
        flank,
        normalization,
        site_set=site_set,
        logger=getattr(args, "logger", None),
        max_centers=max_centers,
        aggregate_site_maps=aggregate_site_maps,
    )
    sample_names = list(getattr(args, "sample_names", []) or [f"sample_{idx + 1}" for idx in range(len(args.aggregate_signals))])
    tasks = [(row.to_dict(), (c1, c2), aggregate_bed_root, list(args.aggregate_signals), cond_groups, flank, len(x), base, normalization, aggregate_norm_spec, sample_names, site_set, max_centers, aggregate_site_maps) for _, row in selected.iterrows()]

    cores = max(1, int(getattr(args, "cores", 1) or 1))
    if cores > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=min(cores, len(tasks))) as executor:
            payloads = list(executor.map(_aggregate_payload_for_row, tasks))
    else:
        payloads = [_aggregate_payload_for_row(task) for task in tasks]
    motifs_payload = [payload for payload in payloads if payload is not None]
    y_label = "Corrected cut-site signal"
    if normalization == "sample-quantile":
        y_label = "Quantile-scaled corrected cut-site signal"
    elif normalization == "condition-quantile":
        y_label = "Condition-quantile-scaled corrected cut-site signal"
    elif normalization == "size-factor":
        y_label = "Size-factor-scaled corrected cut-site signal"
    default_plot_count = max(1, min(12, int(getattr(args, "default_aggregate_plots", 4))))
    default_motifs = [payload["prefix"] for payload in motifs_payload]
    return {"x": x, "motifs": motifs_payload, "comparison": f"{c1} / {c2}", "normalization": normalization, "site_set": site_set, "max_sites_per_motif": max_centers, "default_motifs": default_motifs, "default_plot_count": min(default_plot_count, len(default_motifs)), "x_label": "Distance from motif center (bp)", "y_label": y_label}


def select_aggregate_rows(rows, selectors):
    """Resolve an ordered list of exact motif selectors against a result table."""

    columns = [column for column in ("output_prefix", "prefix", "motif_id", "name") if column in rows.columns]
    if not columns:
        raise ValueError("Aggregate-motif selection requires motif identifiers in the result table")
    selected_indices = []
    for selector in selectors:
        token = str(selector).strip().casefold()
        matches = set()
        for column in columns:
            values = rows[column].fillna("").astype(str).str.strip().str.casefold()
            matches.update(rows.index[values == token].tolist())
        if not matches:
            raise ValueError(f"Unknown aggregate motif: {selector}")
        if len(matches) > 1:
            choices = ", ".join(
                str(rows.loc[index].get("output_prefix", rows.loc[index].get("name", index)))
                for index in sorted(matches)
            )
            raise ValueError(f"Ambiguous aggregate motif {selector!r}; use one of: {choices}")
        index = next(iter(matches))
        if index not in selected_indices:
            selected_indices.append(index)
    return rows.loc[selected_indices].copy()


def _compressed_json_b64(payload):
    text = json.dumps(payload, separators=(",", ":"), ensure_ascii=False)
    return base64.b64encode(gzip.compress(text.encode("utf-8"), compresslevel=9)).decode("ascii")


def _benjamini_hochberg_values(pvalues):
    pvals = np.asarray(pvalues, dtype=float)
    qvals = np.ones(pvals.shape, dtype=float)
    finite = np.isfinite(pvals)
    if not finite.any():
        return qvals
    clipped = np.clip(pvals[finite], 0.0, 1.0)
    order = np.argsort(clipped)
    ranked = clipped[order]
    adjusted = ranked * float(len(ranked)) / np.arange(1, len(ranked) + 1, dtype=float)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    unsorted = np.empty_like(adjusted)
    unsorted[order] = adjusted
    qvals[finite] = unsorted
    return qvals


def _motif_logo_svg(motif, width=420, height=150):
    counts = getattr(motif, "counts", None)
    if counts is None:
        return getattr(motif, "logo_svg", "") or ""
    try:
        counts = np.asarray(counts, dtype=float)
        if counts.shape[0] != 4 or counts.shape[1] == 0:
            return ""
        col_sums = np.sum(counts, axis=0)
        col_sums = np.where(np.isclose(col_sums, 0.0), 1.0, col_sums)
        pfm = counts / col_sums
        entropy = -np.sum(np.where(pfm > 0, pfm * np.log2(np.maximum(pfm, 1e-12)), 0.0), axis=0)
        bits = pfm * np.maximum(0.0, 2.0 - entropy)
    except Exception:
        return ""
    bases = ["A", "C", "G", "T"]
    colors = {"A": "#198754", "C": "#0d6efd", "G": "#f59f00", "T": "#dc3545"}
    left, right, top, bottom = 46, 14, 16, 32
    plot_w, plot_h = width - left - right, height - top - bottom
    npos = bits.shape[1]
    col_w = plot_w / max(1, npos)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img">', '<rect width="100%" height="100%" fill="#ffffff"/>', f'<line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#3b4552" stroke-width="1.2"/>', f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#3b4552" stroke-width="1.2"/>', f'<text x="18" y="{top + plot_h / 2}" transform="rotate(-90 18 {top + plot_h / 2})" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="12" font-weight="700" fill="#152133">bits</text>', f'<text x="{left + plot_w / 2}" y="{height - 7}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="12" font-weight="700" fill="#152133">position</text>']
    for tick in [0, 1, 2]:
        y = top + plot_h - (tick / 2.0) * plot_h
        parts.append(f'<line x1="{left - 4}" y1="{y:.2f}" x2="{left}" y2="{y:.2f}" stroke="#3b4552" stroke-width="1"/>')
        parts.append(f'<text x="{left - 8}" y="{y + 4:.2f}" text-anchor="end" font-family="Arial,Helvetica,sans-serif" font-size="11" font-weight="700" fill="#56616f">{tick}</text>')
    for pos in range(npos):
        y_cursor = top + plot_h
        order = np.argsort(bits[:, pos])
        x_center = left + pos * col_w + col_w / 2.0
        if npos <= 18 or pos in {0, npos - 1} or (pos + 1) % 5 == 0:
            parts.append(f'<text x="{x_center:.2f}" y="{top + plot_h + 13}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="9" font-weight="700" fill="#56616f">{pos + 1}</text>')
        for base_idx in order:
            value = float(bits[base_idx, pos])
            if value <= 0.015:
                continue
            letter_h = max(3.0, value / 2.0 * plot_h)
            y_cursor -= letter_h
            base = bases[base_idx]
            font_size = max(8.0, min(40.0, letter_h * 1.25))
            parts.append(f'<text x="{x_center:.2f}" y="{y_cursor + letter_h * 0.88:.2f}" text-anchor="middle" font-family="Arial Black,Arial,Helvetica,sans-serif" font-size="{font_size:.2f}" font-weight="900" fill="{colors[base]}">{base}</text>')
    parts.append('</svg>')
    return "".join(parts)


def _motif_matrix_map(motifs):
    matrices = {}
    for motif in motifs:
        prefix = str(getattr(motif, "prefix", getattr(motif, "name", "")))
        counts = getattr(motif, "counts", None)
        if counts is None:
            continue
        try:
            counts = np.asarray(counts, dtype=float)
            if counts.shape[0] != 4 or counts.shape[1] == 0:
                continue
            matrices[prefix] = [
                [round(float(value), 4) for value in row]
                for row in counts.tolist()
            ]
        except Exception:
            continue
    return matrices


def _motif_logo_map(motifs):
    logos = {}
    for motif in motifs:
        prefix = str(getattr(motif, "prefix", getattr(motif, "name", "")))
        if getattr(motif, "counts", None) is not None:
            continue
        svg = _motif_logo_svg(motif)
        png = getattr(motif, "base", "") or ""
        entry = {}
        if svg:
            entry["svg"] = svg
        elif png:
            entry["png"] = "data:image/png;base64," + png
        if entry:
            logos[prefix] = entry
    return logos


def plot_interactive_diff_footprints(
    motifs,
    comparison,
    html_out,
    aggregate_data=None,
    title="Differential footprint report",
    report_label=None,
    change_label="Differential footprint score",
    results_table=None,
):
    cond1, cond2 = comparison
    display_title = f"{title} ({cond1} vs {cond2})" if title == "Differential footprint report" else title
    groups = [cond1 + "_up", cond2 + "_up", "n.s."]
    colors = {cond1 + "_up": "#dc2626", cond2 + "_up": "#2563eb", "n.s.": "#8a94a6"}
    points = []
    for motif in motifs:
        group = getattr(motif, "group", "n.s.")
        if group not in colors:
            group = "n.s."
        pvalue = max(float(getattr(motif, "pvalue", 1.0)), 1e-308)
        qvalue = getattr(motif, "qvalue", None)
        qvalue = float(qvalue) if qvalue is not None and np.isfinite(float(qvalue)) else np.nan
        point = {
            "prefix": str(getattr(motif, "prefix", getattr(motif, "name", ""))),
            "name": str(getattr(motif, "name", "")),
            "motif_id": str(getattr(motif, "id", "")),
            "group": group,
            "change": round(float(getattr(motif, "change", 0.0)), 6),
            "pvalue": pvalue,
            "fdr": qvalue,
            "neglog10p": round(float(-np.log10(pvalue)), 6),
        }
        region_stats = getattr(motif, "region_stats", None)
        if isinstance(region_stats, dict):
            for key in (
                "ci_lower", "ci_upper", "n_regions_set_1", "n_regions_set_2",
                "n_motif_regions_set_1", "n_motif_regions_set_2",
                "motif_prevalence_set_1", "motif_prevalence_set_2", "status",
                "statistical_method",
            ):
                if key in region_stats:
                    value = region_stats[key]
                    point[key] = None if isinstance(value, (float, np.floating)) and not np.isfinite(value) else value
        points.append(point)
    fallback_fdr = _benjamini_hochberg_values([point["pvalue"] for point in points])
    for point, qvalue in zip(points, fallback_fdr):
        if not np.isfinite(point["fdr"]):
            point["fdr"] = float(qvalue)
    payload = {
        "title": display_title,
        "report_label": report_label or "",
        "conditions": [cond1, cond2],
        "groups": groups,
        "colors": colors,
        "points": points,
        "motif_matrices": _motif_matrix_map(motifs),
        "logos": _motif_logo_map(motifs),
        "aggregate": aggregate_data or {"x": [], "motifs": []},
        "change_label": change_label,
        "results_tsv": results_table.to_csv(sep="\t", index=False, na_rep="NA") if results_table is not None else "",
    }
    from fp_tools.tools.static_comparison_browser import write_embedded_static_browser

    review = {"schema": "fp-tools.review-multi-comparisons.v1", "title": display_title,
              "comparisons": [{"label": display_title, "payload": payload}]}
    write_embedded_static_browser(review, html_out, source_payload=payload)
