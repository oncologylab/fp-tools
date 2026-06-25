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
import pyBigWig
import pysam

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
    fasta_obj = pysam.FastaFile(fasta)
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

    fasta_obj = pysam.FastaFile(args.genome)
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


def process_tfbs(TF_name, args, log2fc_params):
    """Split into bound/unbound, write per-TF BED/overview, return TF summary row."""
    logger = FpToolsLogger("", args.verbosity, args.log_q)

    bed_outdir = os.path.join(args.outdir, TF_name, "beds")
    filename = os.path.join(bed_outdir, TF_name + ".tmp")
    tmp_files = [filename]
    no_cond = len(args.cond_names)
    comparisons = args.comparisons
    diff_dist = scipy.stats.norm

    if args.output_peaks is not None:
        # Import lazily so diff-footprints --help and parser-only paths do not touch
        # pybedtools/genomepy cache initialization.
        from pybedtools import BedTool

        output_bt = BedTool(args.output_peaks)
        sites_bt = BedTool(filename)
        intersection = sites_bt.intersect(output_bt, u=True)
        filename = intersection.fn
        tmp_files.append(intersection.fn)

    stime = datetime.now()
    header = ["TFBS_chr", "TFBS_start", "TFBS_end", "TFBS_name", "TFBS_score", "TFBS_strand"] \
             + args.peak_header_list \
             + [f"{sample}_score" for sample in args.sample_names]
    with open(filename) as f:
        bedlines = [dict(zip(header, line.rstrip().split("\t"))) for line in f.readlines()]
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
    # write *_all.bed
    outfile = os.path.join(bed_outdir, TF_name + "_all.bed")
    dict_to_tab(bedlines, outfile, header + condition_columns + condition_sd_columns)

    # write bound/unbound per condition
    for condition in args.cond_names:
        chosen_columns = header[:-len(args.sample_names)] + [condition + "_score"]
        for state in ["bound", "unbound"]:
            chosen_bool = 1 if state == "bound" else 0
            subset = [bl for bl in bedlines if bl[condition + "_bound"] == chosen_bool]
            outfile = os.path.join(bed_outdir, f"{TF_name}_{condition}_{state}.bed")
            dict_to_tab(subset, outfile, chosen_columns)

    # overview (txt + optional xlsx)
    overview_columns = header + condition_columns + condition_sd_columns + [c + "_bound" for c in args.cond_names] \
                       + [f"{c1}_{c2}_delta_fp" for (c1, c2) in comparisons] \
                       + [f"{c1}_{c2}_log2fc" for (c1, c2) in comparisons]
    overview_txt = os.path.join(args.outdir, TF_name, TF_name + "_overview.txt")
    dict_to_tab(bedlines, overview_txt, overview_columns, header=True)

    bed_table = pd.DataFrame(bedlines, columns=overview_columns)
    logger.spam(f"Read table {bed_table.shape} for TF {TF_name}")

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
    info_columns += [f"{cond}_{metric}" for cond, metric in itertools.product(args.cond_names, ["mean_score", "score_sd", "n_replicates", "bound"])]
    info_columns += [f"{c1}_{c2}_{metric}" for (c1, c2), metric in itertools.product(comparisons, ["change", "pvalue", "mean_delta_fp", "mean_log2fc", "delta_fp_se", "log2fc_se"])]
    info_table = pd.DataFrame(np.nan, columns=info_columns, index=[TF_name])

    info_table.at[TF_name, "total_tfbs"] = n_rows
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
            info_table.at[TF_name, base + "_mean_delta_fp"] = np.round(float(np.mean(observed_deltas)), 5) if len(observed_deltas) else np.nan
            info_table.at[TF_name, base + "_mean_log2fc"] = np.round(float(np.mean(observed_log2fcs)), 5) if len(observed_log2fcs) else np.nan
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
                change = (obs_mean - bg_mean) / np.mean([obs_std, bg_std])
                info_table.at[TF_name, base + "_change"] = np.round(change, 5)
            else:
                info_table.at[TF_name, base + "_change"] = 0
                info_table.at[TF_name, base + "_pvalue"] = 1

            np.random.seed(n_obs)
            sample_changes = []
            for _ in range(100):
                sample = scipy.stats.norm.rvs(bg_mean, bg_std, size=n_obs)
                sm, ss = float(np.mean(sample)), float(np.std(sample))
                sample_changes.append((sm - bg_mean) / np.mean([ss, bg_std]))
            ttest = scipy.stats.ttest_1samp(sample_changes, float(info_table.at[TF_name, base + "_change"]))
            info_table.at[TF_name, base + "_pvalue"] = ttest[1]

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

    # cleanup tmp
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

    xvalues = np.array([v["change"] for v in diff_scores.values()])
    yvalues = np.array([v["log10pvalue"] for v in diff_scores.values()])

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


def _aggregate_centers_for_row(outdir, prefix, comparison, site_set, max_centers=None):
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


def _fit_aggregate_normalizers(selected, outdir, aggregate_signals, cond_groups, comparison, flank, mode, site_set="all", logger=None, max_centers=None):
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
        _, row_centers = _aggregate_centers_for_row(outdir, prefix, comparison, site_set, max_centers=max_centers)
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
    elif len(task) == 12:
        row, comparison, outdir, aggregate_signals, cond_groups, flank, x_len, base, normalization, aggregate_norm_spec, sample_names, site_set = task
        max_centers = None
    else:
        row, comparison, outdir, aggregate_signals, cond_groups, flank, x_len, base, normalization, aggregate_norm_spec, sample_names, site_set, max_centers = task
    c1, c2 = comparison
    prefix = str(row["output_prefix"])
    centers_by_condition, all_centers = _aggregate_centers_for_row(outdir, prefix, comparison, site_set, max_centers=max_centers)
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
    if mode == "all":
        selected = select_rows.sort_values([base + "_pvalue_numeric", base + "_abs_change"], ascending=[True, False])
    elif mode == "top":
        selected = select_rows.sort_values([base + "_pvalue_numeric", base + "_abs_change"], ascending=[True, False]).head(top_n)
    else:
        threshold = float(getattr(args, "aggregate_pvalue_threshold", 0.05))
        selected = select_rows[select_rows[base + "_pvalue_numeric"] <= threshold].sort_values([base + "_pvalue_numeric", base + "_abs_change"], ascending=[True, False])
        if selected.empty and not no_fallback:
            selected = select_rows.sort_values([base + "_pvalue_numeric", base + "_abs_change"], ascending=[True, False]).head(top_n)

    flank = max(1, int(getattr(args, "aggregate_flank", 100)))
    x = list(range(-flank, flank))
    requested_aggregate_norm = (getattr(args, "aggregate_normalization", "match") or "match").replace("_", "-")
    normalization = (getattr(args, "normalization", "none") or "none").replace("_", "-") if requested_aggregate_norm == "match" else requested_aggregate_norm
    site_set = (getattr(args, "aggregate_site_set", "all") or "all").replace("_", "-")
    max_centers = getattr(args, "aggregate_max_sites", None)
    cond_groups = {cond: list(indices) for cond, indices in getattr(args, "cond_groups", {}).items()}
    aggregate_norm_spec = _fit_aggregate_normalizers(
        selected,
        args.outdir,
        list(args.aggregate_signals),
        cond_groups,
        (c1, c2),
        flank,
        normalization,
        site_set=site_set,
        logger=getattr(args, "logger", None),
        max_centers=max_centers,
    )
    sample_names = list(getattr(args, "sample_names", []) or [f"sample_{idx + 1}" for idx in range(len(args.aggregate_signals))])
    tasks = [(row.to_dict(), (c1, c2), args.outdir, list(args.aggregate_signals), cond_groups, flank, len(x), base, normalization, aggregate_norm_spec, sample_names, site_set, max_centers) for _, row in selected.iterrows()]

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
    return {"x": x, "motifs": motifs_payload, "comparison": f"{c1} / {c2}", "normalization": normalization, "site_set": site_set, "max_sites_per_motif": max_centers, "x_label": "Distance from motif center (bp)", "y_label": y_label}


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
        points.append({
            "prefix": str(getattr(motif, "prefix", getattr(motif, "name", ""))),
            "name": str(getattr(motif, "name", "")),
            "motif_id": str(getattr(motif, "id", "")),
            "group": group,
            "change": round(float(getattr(motif, "change", 0.0)), 6),
            "pvalue": pvalue,
            "fdr": qvalue,
            "neglog10p": round(float(-np.log10(pvalue)), 6),
        })
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
    }
    payload_b64 = _compressed_json_b64(payload)

    html_template = '''<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>__TITLE_ATTR__</title><meta name="viewport" content="width=device-width,initial-scale=1"><style>
:root{color-scheme:light;--border:#d8e2ef;--ink:#172033;--muted:#64748b;--panel:#fff;--bg:#f6f9fc;--plot-row-height:640px;--options-height:300px}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font-family:Arial,Helvetica,sans-serif;font-size:13px}.wrap{padding:8px}.panel{max-width:1900px;margin:0 auto}.head{display:block;margin-bottom:7px}.head h1{margin:0;font-size:34px;line-height:1.05;font-weight:900}.sub{margin:3px 0 0;color:var(--muted);font-weight:800}.options{border:1px solid var(--border);background:var(--panel);border-radius:8px;margin-bottom:8px}.options>summary{cursor:pointer;padding:7px 10px;font-weight:900}.options-grid{display:grid;grid-template-columns:300px 520px minmax(0,1fr);grid-template-areas:"actions samples selected";gap:8px;padding:0 8px 8px;align-items:stretch}.options-actions{grid-area:actions}.options-samples{grid-area:samples}.selected-card{grid-area:selected}.option-col{min-width:0;height:var(--options-height);overflow:auto}.section-title{margin:0 0 5px;font-size:11px;text-transform:uppercase;letter-spacing:.04em;color:#46566c;font-weight:900}.card,.rank-card,.volcano-card,.selected-card,.aggregate-card{background:var(--panel);border:1px solid var(--border);border-radius:8px}.controls{display:flex;flex-wrap:wrap;gap:7px;padding:7px}.control-stack,.export-stack{display:grid;gap:7px;padding:7px}.controls label,.color-row,.sample-style-row{display:flex;align-items:center;gap:6px}#color-controls{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:6px}.color-row{justify-content:space-between;min-width:0}.color-row span{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.controls input,.controls select,.sample-style-row input,.sample-style-row select,.panel-tf,#rank-rows,#rank-rows-slider,#plot-count{height:24px;border:1px solid #c6d3e1;border-radius:5px;background:white;font-size:12px}.sample-visible{width:16px;height:16px}.rows-control,.selected-toolbar{display:flex;align-items:center;justify-content:space-between;gap:8px}.rows-control input[type=number]{width:58px}.rows-control input[type=range]{width:100%;min-width:90px}.selected-toolbar{margin-bottom:5px}.selected-toolbar label,.sort-note{display:flex;align-items:center;gap:6px;font-weight:800;color:#526176}.sort-note{font-size:11px;text-transform:uppercase;letter-spacing:.03em}button{height:26px;border:1px solid #b9c8da;border-radius:5px;background:#f8fbff;font-weight:800;color:#26364d;cursor:pointer}.export-stack button{width:100%;text-align:left}.sample-style-panel{display:grid;grid-template-columns:1fr;gap:8px}.sample-style-group{border:1px solid #e2e8f0;border-radius:6px;padding:6px}.sample-style-group-title{display:flex;gap:6px;align-items:center;font-weight:900;margin-bottom:5px}.sample-style-dot{width:10px;height:10px;border-radius:50%;display:inline-block}.sample-style-row{display:grid;grid-template-columns:30px minmax(92px,1fr) 42px 54px 54px 66px;gap:5px;min-height:28px}.sample-style-head{font-size:11px;font-weight:900;color:#526176}.sample-style-name{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.sample-style-row input[type=color]{width:42px;min-width:42px;padding:1px}.sample-style-row input[type=number],.sample-style-row select{width:100%;min-width:0;padding:0 3px}.dashboard{display:grid;grid-template-columns:430px 560px minmax(0,1fr);grid-template-areas:"rank volcano aggregate";gap:8px;align-items:start}.dashboard>*{min-width:0}.rank-card{grid-area:rank;padding:8px;height:var(--plot-row-height);overflow:auto}.volcano-card{grid-area:volcano;padding:6px;height:var(--plot-row-height);overflow:hidden;display:flex;align-items:flex-start;justify-content:center}.aggregate-card{grid-area:aggregate;padding:8px;min-width:0;height:var(--plot-row-height);overflow:hidden;position:relative}.selected-card{padding:8px;min-width:0}.rank-head{display:flex;align-items:center;justify-content:space-between;gap:8px;margin-bottom:5px}.selected-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:6px}.aggregate-grid{display:grid;grid-template-columns:repeat(var(--aggregate-cols,2),minmax(0,1fr));grid-template-rows:repeat(var(--aggregate-rows,2),minmax(0,1fr));gap:7px;height:100%;place-content:center}.selected-motif,.aggregate-tile{border:1px solid #e1e8f0;border-radius:7px;background:#fff;padding:6px;min-width:0}.selected-motif.active,.aggregate-tile.active{outline:2px solid #93c5fd}.aggregate-tile{display:flex;align-items:stretch;justify-content:center}.selected-head{display:block;margin-bottom:5px}.selected-head .panel-tf{width:100%;min-width:0}.detail-grid{display:grid;grid-template-columns:1fr;gap:4px}.detail-grid p{margin:0;font-size:12px;color:#526176}.motif-group{font-size:15px;font-weight:900}.metric-line{font-weight:800;color:#334155}.motif-logo{height:116px;display:flex;align-items:center;justify-content:center;overflow:hidden;margin-bottom:6px}.motif-logo svg,.motif-logo img{max-width:100%;max-height:112px}.logo-empty{color:#94a3b8}.panel-label{font-weight:900}.sample-picker summary{cursor:pointer;font-weight:800;color:#40516a}.sample-menu{display:grid;grid-template-columns:1fr 1fr;gap:3px;padding-top:4px}.sample-menu label{font-size:11px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}#chart,#rank-chart,.aggregate-panel{width:100%;height:auto;display:block}.volcano-card #chart{width:100%;max-width:548px;aspect-ratio:1/1}.aggregate-panel{height:100%;aspect-ratio:1/1}.aggregate-legend{position:absolute;top:12px;right:12px;z-index:4;display:grid;gap:3px;max-width:220px;padding:6px 7px;background:rgba(255,255,255,.9);border:1px solid #d8e2ef;border-radius:6px;box-shadow:0 1px 3px rgba(15,23,42,.08);font-size:10px;font-weight:800;color:#334155}.legend-row{display:grid;grid-template-columns:34px minmax(0,1fr);gap:5px;align-items:center}.legend-row span{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.legend-line{width:32px;height:0;border-top-style:solid}.plot-title{font-size:15px;font-weight:900;fill:#172033}.axis{stroke:#344256;stroke-width:1.2}.zero{stroke:#7c8798;stroke-width:1.1;stroke-dasharray:4 4}.grid{stroke:#e3eaf3;stroke-width:1}.tick{font-size:11px;fill:#526176;font-weight:700}.axis-label{font-size:12px;fill:#243247;font-weight:900}.summary-label{font-size:10px}.rank-bar{cursor:pointer}.rank-bar.active{stroke:#111827;stroke-width:1.5}.pt{cursor:pointer}.pt.selected{filter:drop-shadow(0 1px 2px rgba(15,23,42,.28))}@media(max-width:1620px){.dashboard{grid-template-columns:400px 520px minmax(0,1fr)}}@media(max-width:1500px){.options-grid{grid-template-columns:280px minmax(0,1fr);grid-template-areas:"actions samples" "selected selected"}.dashboard{grid-template-columns:400px minmax(0,1fr);grid-template-areas:"rank volcano" "aggregate aggregate"}.sample-style-panel{grid-template-columns:1fr}.aggregate-card{height:auto;overflow:auto}.aggregate-grid{height:auto}.aggregate-panel{height:auto}}@media(max-width:980px){.head,.dashboard{display:block}.head h1{font-size:28px}.option-col{height:auto;max-height:none}.rank-card,.volcano-card,.selected-card,.aggregate-card{margin-bottom:8px;height:auto}.aggregate-grid{grid-template-columns:1fr}.sample-style-group{margin-bottom:7px}.sample-style-row{grid-template-columns:28px minmax(86px,1fr) 40px 50px 50px 62px;gap:4px}.sample-style-row input[type=color]{width:40px;min-width:40px}.aggregate-legend{position:static;margin-bottom:8px;max-width:none}}@media(max-width:700px){.options-grid{grid-template-columns:1fr;grid-template-areas:"actions" "samples" "selected";gap:8px}.selected-grid{grid-template-columns:1fr}}
.color-row{font-size:10px;gap:3px;justify-content:flex-start}.color-row input[type=color]{width:28px;min-width:28px;padding:1px}
.dashboard{grid-template-columns:500px 640px minmax(0,1fr)}.aggregate-card{display:grid;grid-template-columns:minmax(0,1fr) 178px;grid-template-rows:minmax(0,1fr);column-gap:6px;position:relative}.aggregate-grid{grid-column:1;grid-row:1;gap:5px;min-height:0}.volcano-card #chart{max-width:628px}.aggregate-legend{grid-column:2;grid-row:1;align-self:start;justify-self:stretch;position:static;z-index:auto;margin-bottom:0;background:rgba(255,255,255,.94)}

</style></head><body><div class="wrap"><div class="panel"><div class="head"><h1>Differential footprint report (<span id="title-cond1">__COND1__</span> vs <span id="title-cond2">__COND2__</span>)</h1><p class="sub" id="report-method" style="display:__REPORT_LABEL_DISPLAY__">__REPORT_LABEL__</p></div><details class="options" open><summary>Show/Hide</summary><div class="options-grid"><div class="option-col options-actions"><p class="section-title">Export editable SVG</p><div class="card export-stack"><button id="download-logo">Download motif logo panel</button><button id="download-rank">Download bar plot panel</button><button id="download-volcano">Download volcano plot</button><button id="download-aggregate">Download motif aggregate panel</button><button id="download-panel">Download combined panel</button></div><p class="section-title" style="margin-top:8px">Groups</p><div class="card control-stack"><div class="controls" id="color-controls"></div><label class="rows-control">Rows <input id="rank-rows-slider" type="range" min="2" max="200" step="1" value="20"><input id="rank-rows" type="number" min="2" max="200" step="1" value="20"></label></div></div><div class="option-col options-samples"><p class="section-title">Sample line styles</p><div class="sample-style-panel" id="aggregate-sample-styles"></div></div><div class="selected-card option-col"><div class="selected-toolbar"><p class="section-title">Selected motifs</p><label>Plots <select id="plot-count"><option value="1">1</option><option value="2">2</option><option value="3">3</option><option value="4" selected>4</option><option value="5">5</option><option value="6">6</option><option value="7">7</option><option value="8">8</option><option value="9">9</option><option value="10">10</option><option value="11">11</option><option value="12">12</option></select></label></div><div id="selected-grid" class="selected-grid"></div></div></div></details><div class="dashboard" id="dashboard"><aside class="rank-card"><svg id="rank-chart" viewBox="0 0 330 680" aria-label="Top differential motifs"></svg></aside><main class="volcano-card"><svg id="chart" viewBox="0 0 760 760" aria-label="Differential footprint volcano plot"></svg></main><section class="aggregate-card"><div id="aggregate-legend" class="aggregate-legend"></div><div class="aggregate-grid" id="aggregate-grid"></div></section></div></div></div><script>
const reportPayloadB64="__PAYLOAD__",aggregateDisplayBp=60,plotSvgStyle='svg,text{font-family:Arial,Helvetica,sans-serif}.plot-title{font-size:15px;font-weight:900;fill:#172033}.axis{stroke:#344256;stroke-width:1.2}.zero{stroke:#7c8798;stroke-width:1.1;stroke-dasharray:4 4}.grid{stroke:#e3eaf3;stroke-width:1}.tick{font-size:11px;fill:#526176;font-weight:700}.axis-label{font-size:12px;fill:#243247;font-weight:900}.summary-label{font-size:10px}.rank-bar.active{stroke:#111827;stroke-width:1.5}.pt.selected{filter:drop-shadow(0 1px 2px rgba(15,23,42,.28))}';let payload=null,panelPrefixes=[],visibleSamples=null,activePanel=0,rankRows=null,sampleLineStyles={},groupColors={};const chart=document.getElementById('chart'),rankChart=document.getElementById('rank-chart'),selectedGrid=document.getElementById('selected-grid'),aggregateGrid=document.getElementById('aggregate-grid'),aggregateLegend=document.getElementById('aggregate-legend'),colorControls=document.getElementById('color-controls'),aggregateSampleStyles=document.getElementById('aggregate-sample-styles'),rankRowsSel=document.getElementById('rank-rows'),rankRowsSlider=document.getElementById('rank-rows-slider'),plotCountSel=document.getElementById('plot-count'),titleCond1=document.getElementById('title-cond1'),titleCond2=document.getElementById('title-cond2'),reportMethod=document.getElementById('report-method');
function escText(value){return String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]))}function b64ToBytes(b64){return Uint8Array.from(atob(b64),c=>c.charCodeAt(0))}async function decodePayload(){if(!('DecompressionStream'in window))throw new Error('This standalone report needs a modern browser with gzip DecompressionStream support.');const ds=new DecompressionStream('gzip');const stream=new Blob([b64ToBytes(reportPayloadB64)]).stream().pipeThrough(ds);return JSON.parse(await new Response(stream).text())}function motifLabel(item){if(!item)return'';const id=item.motif_id||item.id||'';return id?`${item.name} (${id})`:item.name}function pointByPrefix(prefix){return (payload.points||[]).find(p=>p.prefix===prefix)}function aggregateByPrefix(prefix){return ((payload.aggregate||{}).motifs||[]).find(m=>m.prefix===prefix)}function allAggregateMotifs(){return ((payload.aggregate||{}).motifs||[])}function allSampleLabels(){const seen=new Set();allAggregateMotifs().forEach(m=>(m.conditions||[]).forEach(c=>(c.samples||[]).forEach(s=>seen.add(s.name))));return [...seen]}function allSelectableMotifs(){const out=[],seen=new Set();[...allAggregateMotifs(),...(payload.points||[])].forEach(item=>{if(item&&item.prefix&&!seen.has(item.prefix)){seen.add(item.prefix);out.push(item)}});return out.sort((a,b)=>motifLabel(a).localeCompare(motifLabel(b),undefined,{sensitivity:'base'}))}function currentGroupColors(){return{...payload.colors,...groupColors}}function currentConditionColors(){const colors=currentGroupColors();return{[payload.conditions[0]]:colors[payload.conditions[0]+'_up'],[payload.conditions[1]]:colors[payload.conditions[1]+'_up']}}function reportSummary(){const label=String(payload.report_label||'').trim();if(label)return label;const norm=(payload.aggregate||{}).normalization||'none',bed=payload.input_beds_label||'all.bed';return `Aggregate normalization: ${norm}; input beds: ${bed}`}function renderHeader(){const colors=currentGroupColors();titleCond1.style.color=colors[payload.conditions[0]+'_up'];titleCond2.style.color=colors[payload.conditions[1]+'_up'];reportMethod.textContent=reportSummary();reportMethod.style.display='block'}function renderColorControls(){colorControls.innerHTML=payload.groups.map(group=>`<label class="color-row"><span>${escText(group)}</span><input type="color" data-color-group="${escText(group)}" value="${currentGroupColors()[group]}"></label>`).join('');colorControls.querySelectorAll('input').forEach(inp=>inp.addEventListener('input',()=>{groupColors[inp.dataset.colorGroup]=inp.value;renderAll()}))}
function motifLogoSvg(prefix,attrs=''){const counts=(payload.motif_matrices||{})[prefix];if(!Array.isArray(counts)||counts.length!==4||!counts[0]||!counts[0].length)return'';const width=420,height=150,bases=['A','C','G','T'],colors={A:'#198754',C:'#0d6efd',G:'#f59f00',T:'#dc3545'},left=46,right=14,top=16,bottom=32,plotW=width-left-right,plotH=height-top-bottom,npos=counts[0].length,colW=plotW/Math.max(1,npos),bits=[[],[],[],[]];for(let pos=0;pos<npos;pos++){let colSum=counts.reduce((acc,row)=>acc+(Number(row[pos])||0),0);if(Math.abs(colSum)<1e-12)colSum=1;const pfm=counts.map(row=>(Number(row[pos])||0)/colSum),entropy=-pfm.reduce((acc,p)=>acc+(p>0?p*Math.log2(Math.max(p,1e-12)):0),0),scale=Math.max(0,2-entropy);pfm.forEach((p,i)=>bits[i][pos]=p*scale)}const attr=attrs?` ${attrs}`:'',parts=[`<svg${attr} xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img">`,`<rect width="100%" height="100%" fill="#ffffff"/>`,`<line x1="${left}" y1="${top+plotH}" x2="${left+plotW}" y2="${top+plotH}" stroke="#3b4552" stroke-width="1.2"/>`,`<line x1="${left}" y1="${top}" x2="${left}" y2="${top+plotH}" stroke="#3b4552" stroke-width="1.2"/>`,`<text x="18" y="${top+plotH/2}" transform="rotate(-90 18 ${top+plotH/2})" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="12" font-weight="700" fill="#152133">bits</text>`,`<text x="${left+plotW/2}" y="${height-7}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="12" font-weight="700" fill="#152133">position</text>`];[0,1,2].forEach(tick=>{const y=top+plotH-(tick/2)*plotH;parts.push(`<line x1="${left-4}" y1="${y.toFixed(2)}" x2="${left}" y2="${y.toFixed(2)}" stroke="#3b4552" stroke-width="1"/>`,`<text x="${left-8}" y="${(y+4).toFixed(2)}" text-anchor="end" font-family="Arial,Helvetica,sans-serif" font-size="11" font-weight="700" fill="#56616f">${tick}</text>`)});for(let pos=0;pos<npos;pos++){let yCursor=top+plotH,order=[0,1,2,3].sort((a,b)=>bits[a][pos]-bits[b][pos]),xCenter=left+pos*colW+colW/2;if(npos<=18||pos===0||pos===npos-1||(pos+1)%5===0)parts.push(`<text x="${xCenter.toFixed(2)}" y="${top+plotH+13}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="9" font-weight="700" fill="#56616f">${pos+1}</text>`);order.forEach(baseIdx=>{const value=bits[baseIdx][pos];if(value<=0.015)return;const letterH=Math.max(3,value/2*plotH);yCursor-=letterH;const base=bases[baseIdx],fontSize=Math.max(8,Math.min(40,letterH*1.25));parts.push(`<text x="${xCenter.toFixed(2)}" y="${(yCursor+letterH*.88).toFixed(2)}" text-anchor="middle" font-family="Arial Black,Arial,Helvetica,sans-serif" font-size="${fontSize.toFixed(2)}" font-weight="900" fill="${colors[base]}">${base}</text>`)})}parts.push('</svg>');return parts.join('')}
function motifLogoHtml(prefix){const svg=motifLogoSvg(prefix);if(svg)return svg;const logo=(payload.logos||{})[prefix]||{};return logo.svg|| (logo.png?`<img alt="Motif logo" src="${logo.png}">`:'<span class="logo-empty">Motif logo unavailable</span>')}
function targetPlotCount(){return Math.max(1,Math.min(12,Number(plotCountSel.value)||4))}function visiblePanelPrefixes(){return panelPrefixes.slice(0,targetPlotCount())}function gridShape(n){const count=Math.max(1,Math.min(12,Number(n)||4)),cols=Math.ceil(Math.sqrt(count)),rows=Math.ceil(count/cols);return{cols,rows}}function setAggregateGridShape(){const shape=gridShape(visiblePanelPrefixes().length||Number(plotCountSel.value)||4);aggregateGrid.style.setProperty('--aggregate-cols',shape.cols);aggregateGrid.style.setProperty('--aggregate-rows',shape.rows)}
function defaultPanelPrefixes(target){const hasAgg=new Set(allAggregateMotifs().map(m=>m.prefix)),points=(payload.points||[]).filter(p=>p.prefix&&hasAgg.has(p.prefix)),negN=Math.floor(target/2),posN=target-negN,positive=points.filter(p=>p.change>0).sort((a,b)=>b.change-a.change||a.pvalue-b.pvalue).slice(0,posN),negative=points.filter(p=>p.change<0).sort((a,b)=>a.change-b.change||a.pvalue-b.pvalue).slice(0,negN),out=[];[...positive,...negative].forEach(p=>{if(!out.includes(p.prefix))out.push(p.prefix)});points.sort((a,b)=>Math.abs(b.change)-Math.abs(a.change)||a.pvalue-b.pvalue).forEach(p=>{if(out.length<target&&!out.includes(p.prefix))out.push(p.prefix)});allAggregateMotifs().forEach(m=>{if(out.length<target&&!out.includes(m.prefix))out.push(m.prefix)});return out.slice(0,target)}
function ensurePanels(){const motifs=allSelectableMotifs(),labels=allSampleLabels(),target=targetPlotCount();panelPrefixes=panelPrefixes.filter(Boolean).slice(0,12);if(!panelPrefixes.length)panelPrefixes=defaultPanelPrefixes(Math.max(target,4));defaultPanelPrefixes(12).forEach(prefix=>{if(panelPrefixes.length<target&&!panelPrefixes.includes(prefix))panelPrefixes.push(prefix)});for(const motif of motifs){if(panelPrefixes.length>=target)break;if(!panelPrefixes.includes(motif.prefix))panelPrefixes.push(motif.prefix)}while(panelPrefixes.length<target)panelPrefixes.push((motifs[panelPrefixes.length%Math.max(1,motifs.length)]||{}).prefix||'');activePanel=Math.max(0,Math.min(activePanel,target-1));visibleSamples=new Set([...(visibleSamples||new Set(labels))].filter(label=>labels.includes(label)));if(!visibleSamples.size)visibleSamples=new Set(labels)}
function setPanelMotif(idx,prefix){panelPrefixes[idx]=prefix;activePanel=idx;renderAll()}function setSelectedMotif(prefix,opts={}){setPanelMotif(Number.isInteger(opts.panel)?opts.panel:activePanel,prefix)}
function lineDash(type){return type==='dash'?'6 4':(type==='dot'?'1.2 3':'')}function dashAttr(type){const dash=lineDash(type);return dash?` stroke-dasharray="${dash}"`:''}function lineWidthValue(input,fallback){const v=Number(input&&input.value!==undefined?input.value:input);return Number.isFinite(v)&&v>0?Math.min(8,Math.max(.1,v)):fallback}function alphaValue(input,fallback){const v=Number(input&&input.value!==undefined?input.value:input);return Number.isFinite(v)?Math.min(1,Math.max(0.05,v)):fallback}function sampleStyleKey(name){return String(name||'sample').replace(/[^A-Za-z0-9_.-]+/g,'_')}function sampleLineStyle(name,defaults={}){const key=sampleStyleKey(name),stored=sampleLineStyles[key]||{};return{color:stored.color||defaults.color||'#2563eb',alpha:stored.alpha??0.9,width:stored.width||.7,type:stored.type||'solid'}}function setSampleLineStyle(name,patch){const key=sampleStyleKey(name);sampleLineStyles[key]={...sampleLineStyle(name),...patch}}
function renderSampleStyleControls(){const rows=[],colors=currentConditionColors(),by={};allAggregateMotifs().forEach(m=>(m.conditions||[]).forEach(c=>{(by[c.name]||(by[c.name]=new Set()));(c.samples||[]).forEach(s=>by[c.name].add(s.name))}));Object.entries(by).forEach(([cond,names],idx)=>{const defaultColor=colors[cond]||['#2563eb','#dc2626','#16a34a','#9333ea'][idx%4];rows.push(`<div class="sample-style-group" data-sample-group="${escText(cond)}"><div class="sample-style-group-title"><span class="sample-style-dot" style="background:${defaultColor}"></span><span>${escText(cond)} samples</span></div><div class="sample-style-row sample-style-head"><span>Show</span><span>Sample</span><span>Color</span><span>Alpha</span><span>Width</span><span>Line</span></div>`);[...names].forEach(name=>{const style=sampleLineStyle(name,{color:defaultColor}),checked=visibleSamples.has(name)?' checked':'';rows.push(`<label class="sample-style-row" title="Adjust ${escText(cond)} sample ${escText(name)}"><input class="sample-visible" data-sample-visible="${escText(name)}" type="checkbox" aria-label="Show ${escText(cond)} sample ${escText(name)}"${checked}><span class="sample-style-name">${escText(name)}</span><input data-sample-color="${escText(name)}" type="color" aria-label="Color for ${escText(cond)} sample ${escText(name)}" value="${style.color}"><input data-sample-alpha="${escText(name)}" type="number" aria-label="Alpha for ${escText(cond)} sample ${escText(name)}" min="0.05" max="1" step="0.05" value="${style.alpha}"><input data-sample-width="${escText(name)}" type="number" aria-label="Line width for ${escText(cond)} sample ${escText(name)}" min="0.2" max="5" step="0.1" value="${style.width}"><select data-sample-type="${escText(name)}"><option value="solid"${style.type==='solid'?' selected':''}>Solid</option><option value="dash"${style.type==='dash'?' selected':''}>Dash</option><option value="dot"${style.type==='dot'?' selected':''}>Dot</option></select></label>`)});rows.push('</div>')});aggregateSampleStyles.innerHTML=rows.join('');aggregateSampleStyles.querySelectorAll('[data-sample-visible]').forEach(el=>el.addEventListener('change',()=>{if(el.checked)visibleSamples.add(el.dataset.sampleVisible);else visibleSamples.delete(el.dataset.sampleVisible);if(!visibleSamples.size)allSampleLabels().forEach(label=>visibleSamples.add(label));renderAll()}));aggregateSampleStyles.querySelectorAll('[data-sample-color]').forEach(el=>el.addEventListener('input',()=>{setSampleLineStyle(el.dataset.sampleColor,{color:el.value});renderAll(false)}));aggregateSampleStyles.querySelectorAll('[data-sample-alpha]').forEach(el=>el.addEventListener('input',()=>{setSampleLineStyle(el.dataset.sampleAlpha,{alpha:alphaValue(el,.9)});renderAll(false)}));aggregateSampleStyles.querySelectorAll('[data-sample-width]').forEach(el=>el.addEventListener('input',()=>{setSampleLineStyle(el.dataset.sampleWidth,{width:lineWidthValue(el,.7)});renderAll(false)}));aggregateSampleStyles.querySelectorAll('[data-sample-type]').forEach(el=>el.addEventListener('change',()=>{setSampleLineStyle(el.dataset.sampleType,{type:el.value});renderAll(false)}))}
function niceStep(raw){if(!Number.isFinite(raw)||raw<=0)return 1;const pow=Math.pow(10,Math.floor(Math.log10(raw))),f=raw/pow;return (f<=1?1:f<=1.5?1.5:f<=2.5?2.5:f<=5?5:10)*pow}function niceLimit(value){const step=niceStep(Math.abs(value)/5);return Math.max(step,Math.ceil(Math.abs(value)/step)*step)}function niceTicks(min,max,n){const step=niceStep((max-min)/Math.max(1,n-1)),start=Math.ceil(min/step)*step,end=Math.floor(max/step)*step,out=[];for(let v=start;v<=end+step/2;v+=step)out.push(Number(v.toPrecision(12)));return out.length?out:[0]}function fmtTick(value){return Math.abs(value)>=1?value.toFixed(1).replace('-0.0','0.0'):value.toFixed(2).replace('-0.00','0.00')}function fmtShort(v){return Math.abs(v)>=1000?`${Math.round(v/100)/10}k`:String(Math.round(v*100)/100)}
function renderVolcano(showHighlights=true){const colors=currentGroupColors(),width=760,height=760,margin={top:34,right:48,bottom:60,left:84},innerW=width-margin.left-margin.right,innerH=height-margin.top-margin.bottom,plotX0=margin.left,plotY0=margin.top,plotX1=plotX0+innerW,plotY1=plotY0+innerH,changeLabel=payload.change_label||'Differential footprint score',points=payload.points||[];const xs=points.map(p=>p.change),ys=points.map(p=>p.neglog10p),xabs=niceLimit(Math.max(1e-9,Math.abs(Math.min(...xs,0)),Math.abs(Math.max(...xs,0)))*1.05),ymin=0,ymax=niceLimit(Math.max(1,Math.max(...ys,1)*1.03));const sx=x=>plotX0+((x+xabs)/(2*xabs))*innerW,sy=y=>plotY1-((y-ymin)/(ymax-ymin||1))*innerH,xTicks=niceTicks(-xabs,xabs,7),yTicks=niceTicks(ymin,ymax,7),selected=showHighlights?new Set(visiblePanelPrefixes()):new Set(),tickStyle='font-size:15px;font-weight:900;font-family:Arial,Helvetica,sans-serif',axisStyle='font-size:17px;font-weight:900;font-family:Arial,Helvetica,sans-serif';const parts=[`<rect width="${width}" height="${height}" fill="#ffffff"/>`,`<rect x="${plotX0}" y="${plotY0}" width="${innerW}" height="${innerH}" fill="#fbfdff" stroke="#d9e2ec"/>`];yTicks.forEach(v=>parts.push(`<line x1="${plotX0}" y1="${sy(v)}" x2="${plotX1}" y2="${sy(v)}" class="grid"/>`,`<text x="${plotX0-12}" y="${sy(v)+5}" class="tick" style="${tickStyle}" text-anchor="end">${fmtTick(v)}</text>`));xTicks.forEach(v=>parts.push(`<line x1="${sx(v)}" y1="${plotY0}" x2="${sx(v)}" y2="${plotY1}" class="grid"/>`,`<text x="${sx(v)}" y="${plotY1+24}" class="tick" style="${tickStyle}" text-anchor="middle">${fmtTick(v)}</text>`));parts.push(`<line x1="${sx(0)}" y1="${plotY0}" x2="${sx(0)}" y2="${plotY1}" class="zero"/>`,`<line x1="${plotX0}" y1="${plotY1}" x2="${plotX1}" y2="${plotY1}" class="axis"/>`,`<line x1="${plotX0}" y1="${plotY0}" x2="${plotX0}" y2="${plotY1}" class="axis"/>`,`<text x="${(plotX0+plotX1)/2}" y="${height-10}" class="axis-label" style="${axisStyle}" text-anchor="middle">${escText(changeLabel)}</text>`,`<text x="16" y="${plotY0+innerH/2}" class="axis-label" style="${axisStyle}" text-anchor="middle" transform="rotate(-90 16 ${plotY0+innerH/2})">-log10(p-value)</text>`,`<text x="${plotX0+18}" y="${plotY1-14}" font-size="24" font-weight="900" fill="${colors[payload.conditions[1]+'_up']}">${escText(payload.conditions[1]+'_up')}</text>`,`<text x="${plotX1-18}" y="${plotY1-14}" text-anchor="end" font-size="24" font-weight="900" fill="${colors[payload.conditions[0]+'_up']}">${escText(payload.conditions[0]+'_up')}</text>`);points.map((p,idx)=>({p,idx,selected:selected.has(p.prefix)})).sort((a,b)=>Number(a.selected)-Number(b.selected)).forEach(item=>{const p=item.p,selected=item.selected;parts.push(`<circle class="pt${selected?' selected':''}" data-prefix="${escText(p.prefix)}" cx="${sx(p.change).toFixed(2)}" cy="${sy(p.neglog10p).toFixed(2)}" r="${selected?7.2:4.2}" fill="${colors[p.group]||colors['n.s.']}" fill-opacity="${selected?0.98:0.76}" stroke="${selected?'#111827':'#ffffff'}" stroke-width="${selected?2.7:.9}"><title>${escText(motifLabel(p))}</title></circle>`)});chart.innerHTML=parts.join('');if(showHighlights)chart.querySelectorAll('.pt').forEach(el=>el.addEventListener('click',()=>setSelectedMotif(el.dataset.prefix)))}
function drawTopMotifs(showHighlights=true){const points=(payload.points||[]),limit=Math.max(2,Math.floor(Number(rankRowsSel.value||20))),perDir=Math.max(1,Math.floor(limit/2)),positive=points.filter(p=>p.change>0).sort((a,b)=>b.change-a.change||a.pvalue-b.pvalue).slice(0,perDir),negative=points.filter(p=>p.change<0).sort((a,b)=>a.change-b.change||a.pvalue-b.pvalue).slice(0,perDir),shown=[...negative,...positive],width=380,rowH=14,gap=3,sectionGap=8,margin={top:64,bottom:68,left:128,right:14},height=Math.max(430,margin.top+shown.length*(rowH+gap)+sectionGap+margin.bottom),xMid=246,xW=112,maxAbs=niceLimit(Math.max(...shown.map(p=>Math.abs(p.change)),1e-9)),colors=currentGroupColors(),sx=v=>xMid+(v/maxAbs)*xW,axisY=height-60,ticks=niceTicks(-maxAbs,maxAbs,5),selectedSet=showHighlights?new Set(visiblePanelPrefixes()):new Set();rankChart.setAttribute('viewBox',`0 0 ${width} ${height}`);let parts=[`<rect width="${width}" height="${height}" fill="#fff"/><text x="${width/2}" y="18" class="plot-title" text-anchor="middle">Top differential motifs</text><line x1="${xMid}" y1="${margin.top-36}" x2="${xMid}" y2="${axisY}" stroke="#172033" stroke-width="2.2"/>`,`<text x="${xMid-6}" y="${margin.top-22}" text-anchor="end" font-size="14" font-weight="900" fill="${colors[payload.conditions[1]+'_up']}">${escText(payload.conditions[1]+'_up')}</text>`,`<text x="${xMid+6}" y="${margin.top-22}" text-anchor="start" font-size="14" font-weight="900" fill="${colors[payload.conditions[0]+'_up']}">${escText(payload.conditions[0]+'_up')}</text>`];ticks.forEach(v=>parts.push(`<line x1="${sx(v)}" y1="${axisY-4}" x2="${sx(v)}" y2="${axisY+4}" class="axis"/>`,`<text x="${sx(v)}" y="${axisY+17}" class="tick" text-anchor="middle">${fmtShort(v)}</text>`));parts.push(`<line x1="${sx(-maxAbs)}" y1="${axisY}" x2="${sx(maxAbs)}" y2="${axisY}" class="axis"/><text x="${xMid}" y="${height-8}" class="axis-label" text-anchor="middle">${escText(payload.change_label||'Differential footprint score')}</text>`);let y=margin.top;function drawRows(rows){rows.forEach(p=>{const barW=Math.abs(p.change)/maxAbs*xW,x=p.change>=0?xMid:xMid-barW,color=colors[p.group]||colors['n.s.'],active=selectedSet.has(p.prefix),name=escText(motifLabel(p)).slice(0,20),nameY=y+rowH-2;parts.push(`<text class="rank-name${active?' active':''}" data-prefix="${escText(p.prefix)}" x="6" y="${nameY}" font-size="10" font-weight="${active?'900':'700'}" fill="${active?color:'#526176'}" style="cursor:pointer">${name}</text><rect class="rank-bar${active?' active':''}" data-prefix="${escText(p.prefix)}" x="${x}" y="${y}" width="${barW}" height="${rowH}" fill="${color}" fill-opacity="${active?0.95:0.72}"><title>${escText(motifLabel(p))}: ${p.change}</title></rect><text x="${p.change>=0?x-3:x+barW+3}" y="${nameY}" class="tick" text-anchor="${p.change>=0?'end':'start'}">${fmtShort(p.change)}</text>`);y+=rowH+gap})}drawRows(negative);y+=sectionGap;drawRows(positive);rankChart.innerHTML=parts.join('');if(showHighlights)rankChart.querySelectorAll('.rank-bar,.rank-name').forEach(el=>el.addEventListener('click',()=>setSelectedMotif(el.dataset.prefix)))}
function fmtSci(v){const n=Number(v);return Number.isFinite(n)?n.toExponential(0).replace('e-0','e-').replace('e+0','e+'):'NA'}function fmtDelta(v){const n=Number(v);return Number.isFinite(n)?n.toFixed(1):'NA'}function renderSelectedCards(){const motifs=allSelectableMotifs(),colors=currentGroupColors();selectedGrid.innerHTML=visiblePanelPrefixes().map((prefix,idx)=>{const point=pointByPrefix(prefix)||{},aggregate=aggregateByPrefix(prefix),motif=aggregate||point||{prefix},change=Number(point.change||motif.change||0),fdr=Number(point.fdr??point.pvalue??motif.pvalue??1),group=point.group||'n.s.',groupColor=colors[group]||colors['n.s.'];const options=motifs.map(m=>{const hasAgg=!!aggregateByPrefix(m.prefix),suffix=hasAgg?'':' - no aggregate';return `<option value="${escText(m.prefix)}" ${m.prefix===prefix?'selected':''}>${escText(motifLabel(m)+suffix)}</option>`}).join('');return `<div class="selected-motif${idx===activePanel?' active':''}" data-selected-panel="${idx}"><div class="selected-head"><select class="panel-tf" data-panel-tf="${idx}" aria-label="Motif for aggregate plot ${idx+1}">${options}</select></div><div class="motif-logo">${motifLogoHtml(prefix)}</div><div class="detail-grid"><p class="motif-group" style="color:${groupColor}">${escText(group)}</p><p class="metric-line">&#916;FP = ${fmtDelta(change)}</p><p class="metric-line">FDR = ${fmtSci(fdr)}</p></div></div>`}).join('');selectedGrid.querySelectorAll('[data-selected-panel]').forEach(el=>el.addEventListener('click',ev=>{if(ev.target.closest('select'))return;activePanel=Number(el.dataset.selectedPanel);renderAll(false)}));selectedGrid.querySelectorAll('[data-panel-tf]').forEach(sel=>sel.addEventListener('change',()=>setPanelMotif(Number(sel.dataset.panelTf),sel.value)))}
function legendSamples(){const colors=currentConditionColors(),seen=new Set(),rows=[];allAggregateMotifs().forEach(m=>(m.conditions||[]).forEach(c=>(c.samples||[]).forEach(s=>{if((visibleSamples||new Set(allSampleLabels())).has(s.name)&&!seen.has(s.name)){seen.add(s.name);rows.push({...s,condition:c.name,style:sampleLineStyle(s.name,{color:colors[c.name]||'#64748b'})})}})));return rows}function borderStyleForLine(type){return type==='dash'?'dashed':(type==='dot'?'dotted':'solid')}function renderAggregateLegend(){const rows=legendSamples();aggregateLegend.style.display=rows.length?'grid':'none';aggregateLegend.innerHTML=rows.map(row=>`<div class="legend-row"><i class="legend-line" style="border-top-color:${row.style.color};border-top-width:${lineWidthValue(row.style.width,.7)}px;border-top-style:${borderStyleForLine(row.style.type)};opacity:${alphaValue(row.style.alpha,.9)}"></i><span title="${escText(row.name)}">${escText(row.name)}</span></div>`).join('')}function samplesForPanel(motif){const allowed=visibleSamples||new Set(allSampleLabels()),out=[];(motif.conditions||[]).forEach(c=>(c.samples||[]).forEach(s=>{if(allowed.has(s.name))out.push({...s,condition:c.name})}));return out}function meansForPanel(motif,samples){const by={};samples.forEach(s=>(by[s.condition]||(by[s.condition]=[])).push(s));return Object.entries(by).map(([condition,rows])=>{const cond=(motif.conditions||[]).find(c=>c.name===condition)||{},len=Math.max(...rows.map(r=>r.profile.length),0),profile=[];for(let i=0;i<len;i++)profile.push(rows.reduce((acc,r)=>acc+(Number(r.profile[i])||0),0)/rows.length);return{name:condition,condition,profile,n_sites:cond.n_sites||motif.n_sites||0}})}
function pathD(profile,x,sx,sy){return profile.map((y,i)=>`${i?'L':'M'}${sx(x[i]).toFixed(2)},${sy(y).toFixed(2)}`).join(' ')}function visibleConditionNames(samples){return [...new Set(samples.map(s=>s.condition).filter(Boolean))]}function visibleSiteCount(motif,samples){const names=visibleConditionNames(samples);if(names.length===1){const cond=(motif.conditions||[]).find(c=>c.name===names[0]);if(cond&&cond.n_sites!==undefined)return cond.n_sites}return motif.n_sites||0}
function drawAggregatePanel(motif,idx){const rawX=(payload.aggregate||{}).x||[],keep=rawX.map((v,i)=>({v,i})).filter(p=>p.v>=-aggregateDisplayBp&&p.v<=aggregateDisplayBp),x=keep.length?keep.map(p=>p.v):rawX,samples=samplesForPanel(motif).map(s=>({...s,profile:keep.length?keep.map(p=>s.profile[p.i]):s.profile})),siteCount=visibleSiteCount(motif,samples),width=300,height=300,margin={top:30,right:8,bottom:34,left:36},innerW=width-margin.left-margin.right,innerH=height-margin.top-margin.bottom,colors=currentConditionColors(),allY=samples.flatMap(s=>s.profile||[]).filter(Number.isFinite);let ymin=Math.min(...allY,0),ymax=Math.max(...allY,1e-9);const pad=Math.max((ymax-ymin||1)*.18,1e-6);ymin-=pad;ymax+=pad;const sx=v=>margin.left+((v-x[0])/(x[x.length-1]-x[0]||1))*innerW,sy=v=>margin.top+innerH-((v-ymin)/(ymax-ymin||1))*innerH,yTicks=niceTicks(ymin,ymax,4),xTicks=[-aggregateDisplayBp,0,aggregateDisplayBp],sampleSeries=samples.slice().sort((a,b)=>(Number(b.fp_score||0)-Number(a.fp_score||0)));let parts=[`<svg class="aggregate-panel" data-panel="${idx}" viewBox="0 0 ${width} ${height}"><style>${plotSvgStyle}</style><rect width="${width}" height="${height}" fill="#fff"/><text x="${width/2}" y="18" class="plot-title" text-anchor="middle">${escText(motifLabel(motif))}</text>`];yTicks.forEach(v=>parts.push(`<line x1="${margin.left}" y1="${sy(v)}" x2="${margin.left+innerW}" y2="${sy(v)}" class="grid"/><text x="${margin.left-6}" y="${sy(v)+3}" class="tick" text-anchor="end">${fmtTick(v)}</text>`));xTicks.forEach(v=>parts.push(`<line x1="${sx(v)}" y1="${margin.top}" x2="${sx(v)}" y2="${margin.top+innerH}" class="grid"/><text x="${sx(v)}" y="${margin.top+innerH+17}" class="tick" text-anchor="middle">${v}</text>`));parts.push(`<line x1="${sx(0)}" y1="${margin.top}" x2="${sx(0)}" y2="${margin.top+innerH}" class="zero"/><line x1="${margin.left}" y1="${margin.top+innerH}" x2="${margin.left+innerW}" y2="${margin.top+innerH}" class="axis"/><line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${margin.top+innerH}" class="axis"/><text x="${margin.left+7}" y="${margin.top+innerH-8}" class="tick" fill="#94a3b8">${siteCount}</text>`);sampleSeries.forEach((s,i)=>{const style=sampleLineStyle(s.name,{color:colors[s.condition]||'#64748b'}),dash=dashAttr(style.type);parts.push(`<path d="${pathD(s.profile,x,sx,sy)}" fill="none" stroke="${style.color}" stroke-width="${lineWidthValue(style.width,.7)}"${dash} stroke-opacity="${alphaValue(style.alpha,.9)}"><title>${escText(s.name)} - ${escText(s.condition)}</title></path>`)});parts.push(`<text x="${margin.left+innerW/2}" y="${height-6}" class="axis-label" text-anchor="middle">${escText((payload.aggregate||{}).x_label||'Distance from motif center (bp)')}</text><text x="10" y="${margin.top+innerH/2}" class="axis-label" text-anchor="middle" transform="rotate(-90 10 ${margin.top+innerH/2})">${escText((payload.aggregate||{}).y_label||'Corrected cut-site signal')}</text></svg>`);return parts.join('')}
function drawMissingAggregatePanel(item,idx){const width=300,height=300,label=motifLabel(item)||item.prefix||'Selected motif';return `<svg class="aggregate-panel" data-panel="${idx}" viewBox="0 0 ${width} ${height}"><style>${plotSvgStyle}</style><rect width="${width}" height="${height}" fill="#fff"/><rect x="18" y="18" width="${width-36}" height="${height-36}" fill="#fbfdff" stroke="#d9e2ec"/><text x="${width/2}" y="86" class="plot-title" text-anchor="middle">${escText(label)}</text><text x="${width/2}" y="138" class="axis-label" text-anchor="middle">No aggregate profile</text><text x="${width/2}" y="158" class="tick" text-anchor="middle">in this HTML</text><text x="${width/2}" y="190" class="tick" text-anchor="middle">Use --plot-aggregate all</text><text x="${width/2}" y="208" class="tick" text-anchor="middle">or increase --plot-aggregate-top-n</text></svg>`}
function renderAggregateGrid(){setAggregateGridShape();aggregateGrid.innerHTML=visiblePanelPrefixes().map((prefix,idx)=>{const motif=aggregateByPrefix(prefix),point=pointByPrefix(prefix)||{prefix};return `<div class="aggregate-tile${idx===activePanel?' active':''}" data-tile="${idx}">${motif?drawAggregatePanel(motif,idx):drawMissingAggregatePanel(point,idx)}</div>`}).join('');aggregateGrid.querySelectorAll('[data-tile]').forEach(el=>el.addEventListener('click',()=>{activePanel=Number(el.dataset.tile);renderAll(false)}))}
function renderAll(refreshStyles=true){ensurePanels();renderHeader();if(refreshStyles)renderSampleStyleControls();drawTopMotifs();renderVolcano();renderSelectedCards();renderAggregateGrid();renderAggregateLegend()}function styledSvgClone(svgNode){const clone=svgNode.cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');clone.setAttribute('font-family','Arial,Helvetica,sans-serif');clone.style.fontFamily='Arial,Helvetica,sans-serif';if(!clone.querySelector('style'))clone.insertAdjacentHTML('afterbegin',`<style>${plotSvgStyle}</style>`);return clone}function svgBlob(svgNode){return new Blob([new XMLSerializer().serializeToString(styledSvgClone(svgNode))],{type:'image/svg+xml;charset=utf-8'})}function downloadBlob(blob,filename){const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=filename;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000)}function exportLegendSvg(x=0,y=8){const rows=legendSamples(),legendW=170,rowH=15,pad=7,h=rows.length*rowH+pad*2,parts=[];if(!rows.length)return'';parts.push(`<g class="aggregate-export-legend" transform="translate(${x},${y})"><rect width="${legendW}" height="${h}" rx="5" fill="#ffffff" fill-opacity="0.94" stroke="#d8e2ef"/>`);rows.forEach((row,i)=>{const yy=pad+i*rowH+9,dash=dashAttr(row.style.type);parts.push(`<line x1="8" y1="${yy-3}" x2="38" y2="${yy-3}" stroke="${row.style.color}" stroke-width="${lineWidthValue(row.style.width,.7)}"${dash} stroke-opacity="${alphaValue(row.style.alpha,.9)}"/><text x="44" y="${yy}" class="tick">${escText(row.name).slice(0,20)}</text>`)});parts.push('</g>');return parts.join('')}function motifLogoPanelSvg(){const cards=[...selectedGrid.querySelectorAll('.selected-motif')],cardW=240,cardH=220,gap=10,cols=Math.max(1,Math.min(4,cards.length)),rows=Math.ceil(cards.length/cols),colors=currentGroupColors(),parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${cols*cardW+(cols-1)*gap} ${rows*cardH+(rows-1)*gap}" font-family="Arial,Helvetica,sans-serif"><style>${plotSvgStyle}.motif-card-title{font-size:13px;font-weight:900;fill:#172033}.motif-card-metric{font-size:12px;font-weight:800;fill:#334155}</style><rect width="100%" height="100%" fill="#ffffff"/>`];visiblePanelPrefixes().forEach((prefix,idx)=>{const point=pointByPrefix(prefix)||{},aggregate=aggregateByPrefix(prefix),motif=aggregate||point||{prefix},logo=(payload.logos||{})[prefix]||{},group=point.group||'n.s.',groupColor=colors[group]||colors['n.s.'],change=Number(point.change||motif.change||0),fdr=Number(point.fdr??point.pvalue??motif.pvalue??1),x=(idx%cols)*(cardW+gap),y=Math.floor(idx/cols)*(cardH+gap),active=idx===activePanel;parts.push(`<g transform="translate(${x},${y})"><rect width="${cardW}" height="${cardH}" rx="7" fill="#ffffff" stroke="${active?'#93c5fd':'#d8e2ef'}" stroke-width="${active?3:1}"/><text x="10" y="20" class="motif-card-title">${escText(motifLabel(motif)).slice(0,32)}</text>`);const matrixSvg=motifLogoSvg(prefix,`x="24" y="34" width="${cardW-48}" height="96"`);if(matrixSvg){parts.push(matrixSvg)}else if(logo.svg){parts.push(logo.svg.replace(/<svg\b/i,`<svg x="24" y="34" width="${cardW-48}" height="96"`))}else if(logo.png){parts.push(`<image x="24" y="34" width="${cardW-48}" height="96" preserveAspectRatio="xMidYMid meet" href="${logo.png}"/>`)}else{parts.push(`<text x="${cardW/2}" y="84" class="tick" text-anchor="middle">Motif logo unavailable</text>`)}parts.push(`<text x="10" y="150" font-size="13" font-weight="900" fill="${groupColor}">${escText(group)}</text><text x="10" y="171" class="motif-card-metric">&#916;FP = ${fmtDelta(change)}</text><text x="10" y="192" class="motif-card-metric">FDR = ${fmtSci(fdr)}</text></g>`)});parts.push('</svg>');return parts.join('')}function downloadMotifLogoPanel(){downloadBlob(new Blob([motifLogoPanelSvg()],{type:'image/svg+xml;charset=utf-8'}),'diff_footprints_motif_logo_panel.svg')}function downloadAggregateGrid(){const svgs=[...document.querySelectorAll('.aggregate-panel')],w=300,h=300,shape=gridShape(svgs.length),cols=shape.cols,gridW=cols*w,gridH=shape.rows*h,legendW=legendSamples().length?182:0,gap=legendW?12:0,totalW=gridW+gap+legendW,totalH=gridH;let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${totalW} ${totalH}" font-family="Arial,Helvetica,sans-serif"><style>${plotSvgStyle}</style><rect width="${totalW}" height="${totalH}" fill="#ffffff"/>`];svgs.forEach((svg,i)=>{const clone=styledSvgClone(svg);clone.querySelector('style')?.remove();parts.push(`<g transform="translate(${(i%cols)*w},${Math.floor(i/cols)*h})">${clone.innerHTML}</g>`)});parts.push(exportLegendSvg(gridW+gap,8),'</svg>');downloadBlob(new Blob(parts,{type:'image/svg+xml;charset=utf-8'}),'diff_footprints_aggregate_grid.svg')}function downloadStandalonePlot(drawFn,node,filename){drawFn(false);downloadBlob(svgBlob(node),filename);renderAll(false)}function aggregateGridInnerSvg(){const svgs=[...document.querySelectorAll('.aggregate-panel')],w=300,h=300,shape=gridShape(svgs.length),gridW=shape.cols*w,gridH=shape.rows*h,legendW=legendSamples().length?182:0,gap=legendW?12:0,parts=[];svgs.forEach((svg,i)=>{const clone=styledSvgClone(svg);clone.querySelector('style')?.remove();parts.push(`<g transform="translate(${(i%shape.cols)*w},${Math.floor(i/shape.cols)*h})">${clone.innerHTML}</g>`)});parts.push(exportLegendSvg(gridW+gap,8));return{inner:parts.join(''),width:gridW+gap+legendW,height:gridH}}function downloadDashboardPanel(){const rank=styledSvgClone(rankChart),volcano=styledSvgClone(chart),agg=aggregateGridInnerSvg(),rankBox=rank.getAttribute('viewBox').split(/\\s+/).map(Number),volcanoBox=volcano.getAttribute('viewBox').split(/\\s+/).map(Number),rankW=rankBox[2]||330,rankH=rankBox[3]||760,volcanoW=volcanoBox[2]||760,volcanoH=volcanoBox[3]||760,panelH=Math.max(rankH,volcanoH,760),aggScale=panelH/Math.max(1,agg.height),aggW=agg.width*aggScale,gap=24,totalW=rankW+volcanoW+aggW+gap*2;rank.querySelector('style')?.remove();volcano.querySelector('style')?.remove();const svg=`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${totalW} ${panelH}" font-family="Arial,Helvetica,sans-serif"><style>${plotSvgStyle}</style><rect width="${totalW}" height="${panelH}" fill="#ffffff"/><g>${rank.innerHTML}</g><g transform="translate(${rankW+gap},0)">${volcano.innerHTML}</g><g transform="translate(${rankW+volcanoW+gap*2},0) scale(${aggScale})">${agg.inner}</g></svg>`;downloadBlob(new Blob([svg],{type:'image/svg+xml;charset=utf-8'}),'diff_footprints_panel.svg')}function syncRankRows(source){const max=Math.max(2,Number(rankRowsSel.max)||200),value=Math.max(2,Math.min(max,Math.floor(Number(source.value)||20)));rankRowsSel.value=value;rankRowsSlider.value=value;renderAll(false)}
document.getElementById('download-volcano').addEventListener('click',()=>downloadStandalonePlot(renderVolcano,chart,'diff_footprints_volcano.svg'));document.getElementById('download-rank').addEventListener('click',()=>downloadStandalonePlot(drawTopMotifs,rankChart,'diff_footprints_barplot.svg'));document.getElementById('download-aggregate').addEventListener('click',downloadAggregateGrid);document.getElementById('download-panel').addEventListener('click',downloadDashboardPanel);document.getElementById('download-logo').addEventListener('click',downloadMotifLogoPanel);rankRowsSel.addEventListener('input',()=>syncRankRows(rankRowsSel));rankRowsSlider.addEventListener('input',()=>syncRankRows(rankRowsSlider));plotCountSel.addEventListener('change',()=>renderAll(false));decodePayload().then(data=>{payload=data;rankRows=rankRowsSel;groupColors={...payload.colors};const maxRows=Math.max(20,(payload.points||[]).length);rankRowsSel.max=maxRows;rankRowsSlider.max=maxRows;renderColorControls();renderAll()}).catch(err=>{selectedGrid.innerHTML=`<div class="selected-motif"><h2>Could not open report payload</h2><p>${escText(err.message)}</p></div>`});
</script></body></html>'''
    html_str = (html_template
        .replace('__PAYLOAD__', payload_b64)
        .replace('__TITLE_ATTR__', html.escape(display_title, quote=True))
        .replace('__TITLE__', html.escape(display_title))
        .replace('__COND1__', html.escape(cond1))
        .replace('__COND2__', html.escape(cond2))
        .replace('__REPORT_LABEL__', html.escape(report_label or ''))
        .replace('__REPORT_LABEL_DISPLAY__', 'block' if report_label else 'none'))
    with open(html_out, 'w') as f:
        f.write(html_str)
