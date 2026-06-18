#!/usr/bin/env python
"""
Helper functions for BINDetect scoring, summaries, and output generation.

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
            sample_name = f"{condition}_rep{rep_no}"
            sample_bigwigs[sample_name] = pyBigWig.open(args.signals[signal_idx], "rb")
            signal_to_sample[signal_idx] = sample_name

    fasta_obj = pysam.FastaFile(args.genome)
    chrom_boundaries = dict(zip(fasta_obj.references, fasta_obj.lengths))

    rand_window = 200
    background_signal = {
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
        # Import lazily so BINDetect --help and parser-only paths do not touch
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

    log2fc_pdf.close()

    # cleanup tmp
    for fn in tmp_files:
        try:
            os.remove(fn)
        except Exception:
            logger.error(f"Could not remove temporary file {fn} (harmless).")

    return info_table


# ------------------------------ plotting utils ------------------------------ #
def plot_bindetect(motifs, cluster_obj, conditions, args):
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

    # Volcano plot lives in the main BINDetect PDF
    volcano_fig, ax1 = plt.subplots(figsize=(4.0, 4.0))
    ax1.set_title("BINDetect volcano plot", fontsize=PDF_FONT_SIZE, fontweight="bold", pad=12)
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


def build_bindetect_aggregate_payload(motifs, info_table, comparison, args):
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
        else:
            selected = selected.head(top_n)

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


def _motif_logo_map(motifs):
    logos = {}
    for motif in motifs:
        prefix = str(getattr(motif, "prefix", getattr(motif, "name", "")))
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


def plot_interactive_bindetect(
    motifs,
    comparison,
    html_out,
    aggregate_data=None,
    title="Differential footprint report",
    report_label=None,
    change_label="Differential footprint score",
):
    cond1, cond2 = comparison
    groups = [cond1 + "_up", cond2 + "_up", "n.s."]
    colors = {cond1 + "_up": "#2563eb", cond2 + "_up": "#dc2626", "n.s.": "#8a94a6"}
    points = []
    for motif in motifs:
        group = getattr(motif, "group", "n.s.")
        if group not in colors:
            group = "n.s."
        pvalue = max(float(getattr(motif, "pvalue", 1.0)), 1e-308)
        points.append({
            "prefix": str(getattr(motif, "prefix", getattr(motif, "name", ""))),
            "name": str(getattr(motif, "name", "")),
            "motif_id": str(getattr(motif, "id", "")),
            "group": group,
            "change": round(float(getattr(motif, "change", 0.0)), 6),
            "pvalue": pvalue,
            "neglog10p": round(float(-np.log10(pvalue)), 6),
        })
    payload = {
        "title": title,
        "report_label": report_label or "",
        "conditions": [cond1, cond2],
        "groups": groups,
        "colors": colors,
        "points": points,
        "logos": _motif_logo_map(motifs),
        "aggregate": aggregate_data or {"x": [], "motifs": []},
        "change_label": change_label,
    }
    payload_b64 = _compressed_json_b64(payload)

    html_template = '''<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>__TITLE_ATTR__</title><meta name="viewport" content="width=device-width,initial-scale=1"><style>
:root{color-scheme:light;--border:#d8e2ef;--ink:#172033;--muted:#64748b;--panel:#fff;--bg:#f6f9fc}*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font-family:Arial,Helvetica,sans-serif;font-size:13px}.wrap{padding:8px}.panel{max-width:1900px;margin:0 auto}.head{display:grid;grid-template-columns:minmax(280px,1fr) auto auto;gap:10px;align-items:end;margin-bottom:6px}.head h1{margin:0;font-size:24px;line-height:1.1}.sub{margin:0;color:var(--muted);font-weight:700}.options{border:1px solid var(--border);background:var(--panel);border-radius:8px;margin-bottom:8px}.options>summary{cursor:pointer;padding:7px 10px;font-weight:900}.options-grid{display:grid;grid-template-columns:1fr 1.1fr 1fr;gap:8px;padding:0 8px 8px}.section-title{margin:0 0 5px;font-size:11px;text-transform:uppercase;letter-spacing:.04em;color:#46566c;font-weight:900}.card,.rank-card,.volcano-card,.selected-card,.aggregate-card{background:var(--panel);border:1px solid var(--border);border-radius:8px}.controls{display:flex;flex-wrap:wrap;gap:7px;padding:7px}.controls label,.color-row,.sample-style-row{display:flex;align-items:center;gap:6px}.controls input,.controls select,.sample-style-row input,.sample-style-row select,.panel-tools select,#rank-rows{height:24px;border:1px solid #c6d3e1;border-radius:5px;background:white;font-size:12px}.button-row{display:flex;gap:7px;padding:7px;flex-wrap:wrap}button{height:26px;border:1px solid #b9c8da;border-radius:5px;background:#f8fbff;font-weight:800;color:#26364d;cursor:pointer}.sample-style-panel{display:grid;grid-template-columns:repeat(2,minmax(280px,1fr));gap:8px;padding:0 8px 8px}.sample-style-group{border:1px solid #e2e8f0;border-radius:6px;padding:6px}.sample-style-group-title{display:flex;gap:6px;align-items:center;font-weight:900;margin-bottom:5px}.sample-style-dot{width:10px;height:10px;border-radius:50%;display:inline-block}.sample-style-row{display:grid;grid-template-columns:minmax(88px,1fr) 44px 58px 58px 74px;gap:7px;min-height:28px}.sample-style-head{font-size:11px;font-weight:900;color:#526176}.sample-style-name{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.dashboard{display:grid;grid-template-columns:340px minmax(0,1fr) 520px;grid-template-areas:"rank volcano right";gap:8px;align-items:start}.dashboard>*{min-width:0}.rank-card{grid-area:rank;padding:8px}.volcano-card{grid-area:volcano;padding:6px}.right-rail{grid-area:right;display:grid;gap:8px;align-content:start;min-width:0}.selected-card,.aggregate-card{padding:8px;min-width:0}.rank-head{display:flex;align-items:center;justify-content:space-between;gap:8px;margin-bottom:5px}.selected-grid,.aggregate-grid{display:grid;grid-template-columns:1fr 1fr;gap:7px}.selected-motif,.aggregate-tile{border:1px solid #e1e8f0;border-radius:7px;background:#fff;padding:6px;min-width:0}.selected-motif.active,.aggregate-tile.active{outline:2px solid #93c5fd}.selected-motif h2{font-size:14px;line-height:1.15;margin:0 0 5px}.detail-grid{display:grid;grid-template-columns:1fr 1fr 1fr;gap:5px}.detail-grid p{margin:0;font-size:11px;color:#526176}.motif-logo{height:116px;display:flex;align-items:center;justify-content:center;overflow:hidden}.motif-logo svg,.motif-logo img{max-width:100%;max-height:112px}.logo-empty{color:#94a3b8}.panel-tools{display:grid;grid-template-columns:50px minmax(0,1fr) 96px;gap:5px;align-items:center;margin-bottom:5px}.panel-label{font-weight:900}.sample-picker{grid-column:1/-1}.sample-picker summary{cursor:pointer;font-weight:800;color:#40516a}.sample-menu{display:grid;grid-template-columns:1fr 1fr;gap:3px;padding-top:4px}.sample-menu label{font-size:11px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}#chart,#rank-chart,.aggregate-panel{width:100%;height:auto;display:block}.aggregate-panel{aspect-ratio:1/1}.plot-title{font-size:15px;font-weight:900;fill:#172033}.axis{stroke:#344256;stroke-width:1.2}.zero{stroke:#7c8798;stroke-width:1.1;stroke-dasharray:4 4}.grid{stroke:#e3eaf3;stroke-width:1}.tick{font-size:11px;fill:#526176;font-weight:700}.axis-label{font-size:12px;fill:#243247;font-weight:900}.summary-label{font-size:10px}.rank-bar{cursor:pointer}.rank-bar.active{stroke:#111827;stroke-width:1.5}.pt{cursor:pointer}.pt.selected{filter:drop-shadow(0 1px 2px rgba(15,23,42,.28))}@media(max-width:1500px){.dashboard{grid-template-columns:320px minmax(0,1fr);grid-template-areas:"rank volcano" "right right"}.right-rail{grid-template-columns:1fr 1fr}.sample-style-panel{grid-template-columns:1fr}}@media(max-width:980px){.head,.options-grid,.dashboard,.right-rail{display:block}.rank-card,.volcano-card,.selected-card,.aggregate-card{margin-bottom:8px}.selected-grid,.aggregate-grid{grid-template-columns:1fr}.sample-style-panel{display:block}.sample-style-group{margin-bottom:7px}}


</style></head><body><div class="wrap"><div class="panel"><div class="head"><h1>__TITLE__</h1><p class="sub">__COND1__ / __COND2__</p><p class="sub" id="report-method" style="display:__REPORT_LABEL_DISPLAY__">__REPORT_LABEL__</p></div><details class="options" open><summary>Plot options</summary><div class="options-grid"><div><p class="section-title">Groups</p><div class="card controls" id="color-controls"></div></div><div><p class="section-title">Aggregate profile</p><div class="card controls"><label class="mean-toggle"><input id="aggregate-show-mean" type="checkbox">Show mean</label><label>Mean <input id="aggregate-mean-width" type="number" min="0.2" max="6" step="0.1" value="1.05"></label><label>Mean type <select id="aggregate-mean-type"><option value="solid">Solid</option><option value="dash">Dash</option><option value="dot">Dot</option></select></label><label>Width <select id="aggregate-width"><option value="normal">Normal</option><option value="wide">Wide</option><option value="full">Full width</option></select></label></div></div><div><p class="section-title">Export editable SVG</p><div class="card button-row"><button id="download-volcano">Download volcano SVG</button><button id="download-aggregate">Download aggregate SVG</button><button id="download-logo">Download motif logo SVG</button></div></div></div><div><p class="section-title">Sample line styles</p><div class="sample-style-panel" id="aggregate-sample-styles"></div></div></details><div class="dashboard" id="dashboard"><aside class="rank-card"><div class="rank-head"><p class="section-title">Top differential motifs</p><label>Rows <select id="rank-rows"><option value="20" selected>20</option><option value="50">50</option><option value="100">100</option><option value="all">All</option></select></label></div><svg id="rank-chart" viewBox="0 0 330 680" aria-label="Top differential motifs"></svg></aside><main class="volcano-card"><svg id="chart" viewBox="0 0 1100 760" aria-label="Differential footprint volcano plot"></svg></main><aside class="right-rail"><div class="selected-card"><p class="section-title">Selected motifs</p><div id="selected-grid" class="selected-grid"></div></div><section class="aggregate-card"><div class="aggregate-grid" id="aggregate-grid"></div></section></aside></div></div></div><script>
const reportPayloadB64="__PAYLOAD__";let payload=null,panelPrefixes=[],panelSamples=[],activePanel=0,rankRows=null,sampleLineStyles={};const chart=document.getElementById('chart'),rankChart=document.getElementById('rank-chart'),selectedGrid=document.getElementById('selected-grid'),aggregateGrid=document.getElementById('aggregate-grid'),colorControls=document.getElementById('color-controls'),aggregateWidth=document.getElementById('aggregate-width'),aggregateShowMean=document.getElementById('aggregate-show-mean'),aggregateMeanWidth=document.getElementById('aggregate-mean-width'),aggregateMeanType=document.getElementById('aggregate-mean-type'),aggregateSampleStyles=document.getElementById('aggregate-sample-styles'),rankRowsSel=document.getElementById('rank-rows');
function escText(value){return String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]))}function b64ToBytes(b64){return Uint8Array.from(atob(b64),c=>c.charCodeAt(0))}async function decodePayload(){if(!('DecompressionStream'in window))throw new Error('This standalone report needs a modern browser with gzip DecompressionStream support.');const ds=new DecompressionStream('gzip');const stream=new Blob([b64ToBytes(reportPayloadB64)]).stream().pipeThrough(ds);return JSON.parse(await new Response(stream).text())}function motifLabel(item){if(!item)return'';const id=item.motif_id||item.id||'';return id?`${item.name} (${id})`:item.name}function pointByPrefix(prefix){return (payload.points||[]).find(p=>p.prefix===prefix)}function aggregateByPrefix(prefix){return ((payload.aggregate||{}).motifs||[]).find(m=>m.prefix===prefix)}function allAggregateMotifs(){return ((payload.aggregate||{}).motifs||[])}function allSampleLabels(){const seen=new Set();allAggregateMotifs().forEach(m=>(m.conditions||[]).forEach(c=>(c.samples||[]).forEach(s=>seen.add(s.name))));return [...seen]}function currentGroupColors(){const out={...payload.colors};document.querySelectorAll('[data-color-group]').forEach(inp=>out[inp.dataset.colorGroup]=inp.value);return out}function currentConditionColors(){const groupColors=currentGroupColors();return{[payload.conditions[0]]:groupColors[payload.conditions[0]+'_up'],[payload.conditions[1]]:groupColors[payload.conditions[1]+'_up']}}function renderColorControls(){colorControls.innerHTML=payload.groups.map(group=>`<label class="color-row"><span>${escText(group)}</span><input type="color" data-color-group="${escText(group)}" value="${payload.colors[group]}"></label>`).join('');colorControls.querySelectorAll('input').forEach(inp=>inp.addEventListener('input',renderAll))}
function ensurePanels(){const motifs=allAggregateMotifs(),labels=allSampleLabels();if(!panelPrefixes.length)panelPrefixes=motifs.slice(0,4).map(m=>m.prefix);while(panelPrefixes.length<4)panelPrefixes.push((motifs[panelPrefixes.length%Math.max(1,motifs.length)]||{}).prefix||'');while(panelSamples.length<4)panelSamples.push(new Set(labels));panelSamples=panelSamples.slice(0,4).map(set=>{const clean=new Set([...set].filter(label=>labels.includes(label)));return clean.size?clean:new Set(labels)})}
function setPanelMotif(idx,prefix){panelPrefixes[idx]=prefix;activePanel=idx;renderAll()}function setSelectedMotif(prefix,opts={}){setPanelMotif(Number.isInteger(opts.panel)?opts.panel:activePanel,prefix)}
function lineDash(type){return type==='dash'?'6 4':(type==='dot'?'1.2 3':'')}function dashAttr(type){const dash=lineDash(type);return dash?` stroke-dasharray="${dash}"`:''}function lineWidthValue(input,fallback){const v=Number(input&&input.value!==undefined?input.value:input);return Number.isFinite(v)&&v>0?Math.min(8,Math.max(.1,v)):fallback}function alphaValue(input,fallback){const v=Number(input&&input.value!==undefined?input.value:input);return Number.isFinite(v)?Math.min(1,Math.max(0.05,v)):fallback}function sampleStyleKey(name){return String(name||'sample').replace(/[^A-Za-z0-9_.-]+/g,'_')}function sampleLineStyle(name,defaults={}){const key=sampleStyleKey(name),stored=sampleLineStyles[key]||{};return{color:stored.color||defaults.color||'#2563eb',alpha:stored.alpha??0.9,width:stored.width||.7,type:stored.type||'solid'}}function setSampleLineStyle(name,patch){const key=sampleStyleKey(name);sampleLineStyles[key]={...sampleLineStyle(name),...patch}}
function renderSampleStyleControls(){const rows=[],colors=currentConditionColors(),by={};allAggregateMotifs().forEach(m=>(m.conditions||[]).forEach(c=>{(by[c.name]||(by[c.name]=new Set()));(c.samples||[]).forEach(s=>by[c.name].add(s.name))}));Object.entries(by).forEach(([cond,names],idx)=>{const defaultColor=colors[cond]||['#2563eb','#dc2626','#16a34a','#9333ea'][idx%4];rows.push(`<div class="sample-style-group" data-sample-group="${escText(cond)}"><div class="sample-style-group-title"><span class="sample-style-dot" style="background:${defaultColor}"></span><span>${escText(cond)} samples</span></div><div class="sample-style-row sample-style-head"><span>Sample</span><span>Color</span><span>Alpha</span><span>Width</span><span>Line</span></div>`);[...names].forEach(name=>{const style=sampleLineStyle(name,{color:defaultColor});rows.push(`<label class="sample-style-row" title="Adjust ${escText(cond)} sample ${escText(name)}"><span class="sample-style-name">${escText(name)}</span><input data-sample-color="${escText(name)}" type="color" aria-label="Color for ${escText(cond)} sample ${escText(name)}" value="${style.color}"><input data-sample-alpha="${escText(name)}" type="number" aria-label="Alpha for ${escText(cond)} sample ${escText(name)}" min="0.05" max="1" step="0.05" value="${style.alpha}"><input data-sample-width="${escText(name)}" type="number" aria-label="Line width for ${escText(cond)} sample ${escText(name)}" min="0.2" max="5" step="0.1" value="${style.width}"><select data-sample-type="${escText(name)}"><option value="solid"${style.type==='solid'?' selected':''}>Solid</option><option value="dash"${style.type==='dash'?' selected':''}>Dash</option><option value="dot"${style.type==='dot'?' selected':''}>Dot</option></select></label>`)});rows.push('</div>')});aggregateSampleStyles.innerHTML=rows.join('');aggregateSampleStyles.querySelectorAll('[data-sample-color]').forEach(el=>el.addEventListener('input',()=>{setSampleLineStyle(el.dataset.sampleColor,{color:el.value});renderAll(false)}));aggregateSampleStyles.querySelectorAll('[data-sample-alpha]').forEach(el=>el.addEventListener('input',()=>{setSampleLineStyle(el.dataset.sampleAlpha,{alpha:alphaValue(el,.9)});renderAll(false)}));aggregateSampleStyles.querySelectorAll('[data-sample-width]').forEach(el=>el.addEventListener('input',()=>{setSampleLineStyle(el.dataset.sampleWidth,{width:lineWidthValue(el,.7)});renderAll(false)}));aggregateSampleStyles.querySelectorAll('[data-sample-type]').forEach(el=>el.addEventListener('change',()=>{setSampleLineStyle(el.dataset.sampleType,{type:el.value});renderAll(false)}))}
function niceTicks(min,max,n){const out=[];for(let i=0;i<n;i++)out.push(min+(max-min)*(i/Math.max(1,n-1)));return out}function fmtTick(value){return Math.abs(value)>=1?value.toFixed(1).replace('-0.0','0.0'):value.toFixed(2).replace('-0.00','0.00')}function fmtShort(v){return Math.abs(v)>=1000?`${Math.round(v/100)/10}k`:String(Math.round(v*100)/100)}
function renderVolcano(){const colors=currentGroupColors(),width=1100,height=760,margin={top:52,right:54,bottom:78,left:92},innerW=width-margin.left-margin.right,innerH=height-margin.top-margin.bottom-58,plotX0=margin.left,plotY0=margin.top,plotX1=plotX0+innerW,plotY1=plotY0+innerH,changeLabel=payload.change_label||'Differential footprint score',points=payload.points||[];const xs=points.map(p=>p.change),ys=points.map(p=>p.neglog10p),xabs=Math.max(1,Math.abs(Math.min(...xs,0)),Math.abs(Math.max(...xs,0)))*1.1,ymin=0,ymax=Math.max(1,Math.max(...ys,1)*1.08);const sx=x=>plotX0+((x+xabs)/(2*xabs))*innerW,sy=y=>plotY1-((y-ymin)/(ymax-ymin||1))*innerH,xTicks=niceTicks(-xabs,xabs,7),yTicks=niceTicks(ymin,ymax,7),selected=new Set(panelPrefixes);const parts=[`<rect width="${width}" height="${height}" fill="#ffffff"/>`,`<text x="${(plotX0+plotX1)/2}" y="26" class="plot-title" text-anchor="middle">Differential footprint evidence</text>`,`<rect x="${plotX0}" y="${plotY0}" width="${innerW}" height="${innerH}" fill="#fbfdff" stroke="#d9e2ec"/>`];yTicks.forEach(v=>parts.push(`<line x1="${plotX0}" y1="${sy(v)}" x2="${plotX1}" y2="${sy(v)}" class="grid"/>`,`<text x="${plotX0-12}" y="${sy(v)+4}" class="tick" text-anchor="end">${v.toFixed(1)}</text>`));xTicks.forEach(v=>parts.push(`<line x1="${sx(v)}" y1="${plotY0}" x2="${sx(v)}" y2="${plotY1}" class="grid"/>`,`<text x="${sx(v)}" y="${plotY1+28}" class="tick" text-anchor="middle">${fmtTick(v)}</text>`));parts.push(`<line x1="${sx(0)}" y1="${plotY0}" x2="${sx(0)}" y2="${plotY1}" class="zero"/>`,`<line x1="${plotX0}" y1="${plotY1}" x2="${plotX1}" y2="${plotY1}" class="axis"/>`,`<line x1="${plotX0}" y1="${plotY0}" x2="${plotX0}" y2="${plotY1}" class="axis"/>`,`<text x="${(plotX0+plotX1)/2}" y="${height-18}" class="axis-label" text-anchor="middle">${escText(changeLabel)}</text>`,`<text x="26" y="${plotY0+innerH/2}" class="axis-label" text-anchor="middle" transform="rotate(-90 26 ${plotY0+innerH/2})">-log10(p-value)</text>`,`<text x="${plotX0+20}" y="${plotY0+24}" font-size="14" font-weight="900" fill="${colors[payload.conditions[1]+'_up']}">${escText(payload.conditions[1]+'_up')}</text>`,`<text x="${plotX1-20}" y="${plotY0+24}" text-anchor="end" font-size="14" font-weight="900" fill="${colors[payload.conditions[0]+'_up']}">${escText(payload.conditions[0]+'_up')}</text>`);points.map((p,idx)=>({p,idx,selected:selected.has(p.prefix)})).sort((a,b)=>Number(a.selected)-Number(b.selected)).forEach(item=>{const p=item.p,selected=item.selected;parts.push(`<circle class="pt${selected?' selected':''}" data-prefix="${escText(p.prefix)}" cx="${sx(p.change).toFixed(2)}" cy="${sy(p.neglog10p).toFixed(2)}" r="${selected?7.2:4.2}" fill="${colors[p.group]||colors['n.s.']}" fill-opacity="${selected?0.98:0.76}" stroke="${selected?'#111827':'#ffffff'}" stroke-width="${selected?2.7:.9}"><title>${escText(motifLabel(p))}</title></circle>`)});chart.innerHTML=parts.join('');chart.querySelectorAll('.pt').forEach(el=>el.addEventListener('click',()=>setSelectedMotif(el.dataset.prefix)))}
function drawTopMotifs(){const points=(payload.points||[]),limit=rankRowsSel.value==='all'?Infinity:Number(rankRowsSel.value||20),perDir=rankRowsSel.value==='all'?Infinity:Math.max(1,Math.floor(limit/2)),positive=points.filter(p=>p.change>0).sort((a,b)=>b.change-a.change||a.pvalue-b.pvalue).slice(0,perDir),negative=points.filter(p=>p.change<0).sort((a,b)=>a.change-b.change||a.pvalue-b.pvalue).slice(0,perDir),shown=[...positive,...negative],width=330,rowH=14,gap=3,headerH=17,margin={top:50,bottom:68,left:116,right:10},height=Math.max(430,margin.top+shown.length*(rowH+gap)+headerH*2+margin.bottom),xMid=218,xW=94,maxAbs=Math.max(...shown.map(p=>Math.abs(p.change)),1e-9),colors=currentGroupColors(),sx=v=>xMid+(v/maxAbs)*xW,axisY=height-48,ticks=[-maxAbs,-maxAbs/2,0,maxAbs/2,maxAbs];rankChart.setAttribute('viewBox',`0 0 ${width} ${height}`);let parts=[`<rect width="${width}" height="${height}" fill="#fff"/><text x="${width/2}" y="18" class="plot-title" text-anchor="middle">Top differential motifs</text><line x1="${xMid}" y1="${margin.top-8}" x2="${xMid}" y2="${axisY}" class="zero"/>`];ticks.forEach(v=>parts.push(`<line x1="${sx(v)}" y1="${axisY-4}" x2="${sx(v)}" y2="${axisY+4}" class="axis"/>`,`<text x="${sx(v)}" y="${axisY+17}" class="tick" text-anchor="middle">${fmtShort(v)}</text>`));parts.push(`<line x1="${sx(-maxAbs)}" y1="${axisY}" x2="${sx(maxAbs)}" y2="${axisY}" class="axis"/><text x="${xMid}" y="${height-8}" class="axis-label" text-anchor="middle">${escText(payload.change_label||'Differential footprint score')}</text>`);let y=margin.top;function drawSection(label,rows){parts.push(`<text x="6" y="${y+10}" class="tick" font-weight="900" fill="${colors[label]||'#526176'}">${escText(label)}</text>`);y+=headerH;rows.forEach(p=>{const barW=Math.abs(p.change)/maxAbs*xW,x=p.change>=0?xMid:xMid-barW,color=colors[p.group]||colors['n.s.'],active=panelPrefixes.includes(p.prefix);parts.push(`<text x="6" y="${y+rowH-2}" class="tick summary-label">${escText(motifLabel(p)).slice(0,18)}</text><rect class="rank-bar${active?' active':''}" data-prefix="${escText(p.prefix)}" x="${x}" y="${y}" width="${barW}" height="${rowH}" fill="${color}" fill-opacity="${active?0.95:0.72}"><title>${escText(motifLabel(p))}: ${p.change}</title></rect><text x="${p.change>=0?x+barW+3:x-3}" y="${y+rowH-2}" class="tick" text-anchor="${p.change>=0?'start':'end'}">${fmtShort(p.change)}</text>`);y+=rowH+gap})}drawSection(payload.conditions[0]+'_up',positive);y+=4;drawSection(payload.conditions[1]+'_up',negative);rankChart.innerHTML=parts.join('');rankChart.querySelectorAll('.rank-bar').forEach(el=>el.addEventListener('click',()=>setSelectedMotif(el.dataset.prefix)))}
function renderSelectedCards(){selectedGrid.innerHTML=panelPrefixes.map((prefix,idx)=>{const point=pointByPrefix(prefix)||{},motif=aggregateByPrefix(prefix)||point,logo=payload.logos[prefix]||{},label=motifLabel(motif)||prefix,change=Number(point.change||motif.change||0),pvalue=Number(point.pvalue||motif.pvalue||1);return `<div class="selected-motif${idx===activePanel?' active':''}" data-selected-panel="${idx}"><p class="section-title">Panel ${idx+1}</p><h2>${escText(label)}</h2><div class="detail-grid"><p><strong>Group:</strong><br>${escText(point.group||'')}</p><p><strong>${escText(payload.change_label||'Differential footprint score')}:</strong><br>${change.toFixed(4)}</p><p><strong>P-value:</strong><br>${pvalue.toExponential(3)}</p></div><div class="motif-logo">${logo.svg|| (logo.png?`<img alt="Motif logo" src="${logo.png}">`:'<span class="logo-empty">Motif logo unavailable</span>')}</div></div>`}).join('');selectedGrid.querySelectorAll('[data-selected-panel]').forEach(el=>el.addEventListener('click',()=>{activePanel=Number(el.dataset.selectedPanel);renderAll(false)}))}
function samplesForPanel(motif,idx){const allowed=panelSamples[idx]||new Set(allSampleLabels()),out=[];(motif.conditions||[]).forEach(c=>(c.samples||[]).forEach(s=>{if(allowed.has(s.name))out.push({...s,condition:c.name})}));return out}function meansForPanel(motif,samples){const by={};samples.forEach(s=>(by[s.condition]||(by[s.condition]=[])).push(s));return Object.entries(by).map(([condition,rows])=>{const cond=(motif.conditions||[]).find(c=>c.name===condition)||{},len=Math.max(...rows.map(r=>r.profile.length),0),profile=[];for(let i=0;i<len;i++)profile.push(rows.reduce((acc,r)=>acc+(Number(r.profile[i])||0),0)/rows.length);return{name:condition,condition,profile,n_sites:cond.n_sites||motif.n_sites||0}})}
function pathD(profile,x,sx,sy){return profile.map((y,i)=>`${i?'L':'M'}${sx(x[i]).toFixed(2)},${sy(y).toFixed(2)}`).join(' ')}function bedLabel(motif){const siteSet=motif.site_set||(payload.aggregate||{}).site_set||'motif-site set';return siteSet.endsWith('bed')?siteSet:`${siteSet} motif-site set`}function panelSubtitle(motif,samples){return `${samples.length} sample${samples.length===1?'':'s'} - ${motif.n_sites||0} sites - ${bedLabel(motif)}`}
function drawAggregatePanel(motif,idx){const samples=samplesForPanel(motif,idx),means=meansForPanel(motif,samples),x=(payload.aggregate||{}).x||[],width=360,height=360,margin={top:42,right:18,bottom:42,left:58},innerW=width-margin.left-margin.right,innerH=height-margin.top-margin.bottom,colors=currentConditionColors(),series=aggregateShowMean.checked?samples.concat(means):samples,allY=series.flatMap(s=>s.profile||[]).filter(Number.isFinite);let ymin=Math.min(...allY,0),ymax=Math.max(...allY,1e-9);const pad=Math.max((ymax-ymin||1)*.18,1e-6);ymin-=pad;ymax+=pad;const sx=v=>margin.left+((v-x[0])/(x[x.length-1]-x[0]||1))*innerW,sy=v=>margin.top+innerH-((v-ymin)/(ymax-ymin||1))*innerH,yTicks=niceTicks(ymin,ymax,4),xTicks=[x[0],0,x[x.length-1]],sampleSeries=samples.slice().sort((a,b)=>(Number(b.fp_score||0)-Number(a.fp_score||0)));let parts=[`<svg class="aggregate-panel" data-panel="${idx}" viewBox="0 0 ${width} ${height}"><rect width="${width}" height="${height}" fill="#fff"/><text x="${width/2}" y="17" class="plot-title" text-anchor="middle">${escText(motifLabel(motif))}</text><text x="${width/2}" y="32" class="tick" text-anchor="middle">${escText(panelSubtitle(motif,samples))}</text>`];yTicks.forEach(v=>parts.push(`<line x1="${margin.left}" y1="${sy(v)}" x2="${margin.left+innerW}" y2="${sy(v)}" class="grid"/><text x="${margin.left-8}" y="${sy(v)+3}" class="tick" text-anchor="end">${fmtTick(v)}</text>`));xTicks.forEach(v=>parts.push(`<line x1="${sx(v)}" y1="${margin.top}" x2="${sx(v)}" y2="${margin.top+innerH}" class="grid"/><text x="${sx(v)}" y="${margin.top+innerH+18}" class="tick" text-anchor="middle">${v}</text>`));parts.push(`<line x1="${sx(0)}" y1="${margin.top}" x2="${sx(0)}" y2="${margin.top+innerH}" class="zero"/><line x1="${margin.left}" y1="${margin.top+innerH}" x2="${margin.left+innerW}" y2="${margin.top+innerH}" class="axis"/><line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${margin.top+innerH}" class="axis"/>`);sampleSeries.forEach((s,i)=>{const style=sampleLineStyle(s.name,{color:colors[s.condition]||'#64748b'}),dash=dashAttr(style.type);parts.push(`<path d="${pathD(s.profile,x,sx,sy)}" fill="none" stroke="${style.color}" stroke-width="${lineWidthValue(style.width,.7)}"${dash} stroke-opacity="${alphaValue(style.alpha,.9)}"><title>${escText(s.name)} - ${escText(s.condition)}</title></path>`)});if(aggregateShowMean.checked)means.forEach((s,i)=>{const color=colors[s.condition]||'#64748b';parts.push(`<path d="${pathD(s.profile,x,sx,sy)}" fill="none" stroke="${color}" stroke-width="${lineWidthValue(aggregateMeanWidth,1.05)}"${dashAttr(aggregateMeanType.value)} stroke-opacity="0.95" stroke-linecap="round"/><text x="${margin.left+8}" y="${margin.top+12+i*12}" font-size="9" font-weight="900" fill="${color}">${escText(s.condition)} mean</text>`)});parts.push(`<text x="${margin.left+innerW/2}" y="${height-8}" class="axis-label" text-anchor="middle">${escText((payload.aggregate||{}).x_label||'Distance from motif center (bp)')}</text><text x="15" y="${margin.top+innerH/2}" class="axis-label" text-anchor="middle" transform="rotate(-90 15 ${margin.top+innerH/2})">${escText((payload.aggregate||{}).y_label||'Corrected cut-site signal')}</text></svg>`);return parts.join('')}
function samplePickerHtml(idx){const selected=panelSamples[idx]||new Set(allSampleLabels());return `<details class="sample-picker"><summary>Samples: ${selected.size}</summary><div class="sample-menu">${allSampleLabels().map(label=>`<label><input type="checkbox" data-panel-sample="${idx}" data-sample="${escText(label)}" ${selected.has(label)?'checked':''}> ${escText(label)}</label>`).join('')}</div></details>`}function renderAggregateGrid(){const motifs=allAggregateMotifs();aggregateGrid.innerHTML=panelPrefixes.map((prefix,idx)=>{const motif=aggregateByPrefix(prefix)||motifs[0]||{conditions:[]};return `<div class="aggregate-tile${idx===activePanel?' active':''}" data-tile="${idx}"><div class="panel-tools"><span class="panel-label">Panel ${idx+1}</span><select class="panel-tf" data-panel-tf="${idx}">${motifs.map(m=>`<option value="${escText(m.prefix)}" ${m.prefix===motif.prefix?'selected':''}>${escText(motifLabel(m))}</option>`).join('')}</select><button data-download-panel="${idx}">Download SVG</button>${samplePickerHtml(idx)}</div>${drawAggregatePanel(motif,idx)}</div>`}).join('');aggregateGrid.querySelectorAll('[data-tile]').forEach(el=>el.addEventListener('click',ev=>{if(ev.target.closest('button,select,details,input,label,summary'))return;activePanel=Number(el.dataset.tile);renderAll(false)}));aggregateGrid.querySelectorAll('[data-panel-tf]').forEach(sel=>sel.addEventListener('change',()=>setPanelMotif(Number(sel.dataset.panelTf),sel.value)));aggregateGrid.querySelectorAll('[data-panel-sample]').forEach(inp=>inp.addEventListener('change',()=>{const idx=Number(inp.dataset.panelSample),set=panelSamples[idx];if(inp.checked)set.add(inp.dataset.sample);else set.delete(inp.dataset.sample);if(!set.size)allSampleLabels().forEach(label=>set.add(label));renderAll(false)}));aggregateGrid.querySelectorAll('[data-download-panel]').forEach(btn=>btn.addEventListener('click',ev=>{ev.stopPropagation();const svg=document.querySelector(`.aggregate-panel[data-panel="${btn.dataset.downloadPanel}"]`);if(svg)downloadBlob(svgBlob(svg),`diff_footprints_panel_${Number(btn.dataset.downloadPanel)+1}.svg`)}))}
function renderAll(refreshStyles=true){ensurePanels();renderColorControls();if(refreshStyles)renderSampleStyleControls();drawTopMotifs();renderVolcano();renderSelectedCards();renderAggregateGrid()}function svgBlob(svgNode){const clone=svgNode.cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');return new Blob([new XMLSerializer().serializeToString(clone)],{type:'image/svg+xml;charset=utf-8'})}function downloadBlob(blob,filename){const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=filename;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000)}function downloadAggregateGrid(){const svgs=[...document.querySelectorAll('.aggregate-panel')],w=360,h=360,cols=2;let parts=[`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${cols*w} ${Math.ceil(svgs.length/cols)*h}">`];svgs.forEach((svg,i)=>parts.push(`<g transform="translate(${(i%cols)*w},${Math.floor(i/cols)*h})">${svg.innerHTML}</g>`));parts.push('</svg>');downloadBlob(new Blob(parts,{type:'image/svg+xml;charset=utf-8'}),'diff_footprints_aggregate_grid.svg')}
document.getElementById('download-volcano').addEventListener('click',()=>downloadBlob(svgBlob(chart),'diff_footprints_volcano.svg'));document.getElementById('download-aggregate').addEventListener('click',downloadAggregateGrid);document.getElementById('download-logo').addEventListener('click',()=>{const svg=selectedGrid.querySelector(`.selected-motif[data-selected-panel="${activePanel}"] svg`);if(svg)downloadBlob(svgBlob(svg),'diff_footprints_motif_logo.svg')});[aggregateShowMean,aggregateMeanWidth,aggregateMeanType,aggregateWidth,rankRowsSel].forEach(el=>el.addEventListener('change',()=>renderAll(false)));aggregateMeanWidth.addEventListener('input',()=>renderAll(false));decodePayload().then(data=>{payload=data;rankRows=rankRowsSel;renderAll()}).catch(err=>{selectedGrid.innerHTML=`<div class="selected-motif"><h2>Could not open report payload</h2><p>${escText(err.message)}</p></div>`});
</script></body></html>'''
    html_str = (html_template
        .replace('__PAYLOAD__', payload_b64)
        .replace('__TITLE_ATTR__', html.escape(f'{title} {cond1} / {cond2}', quote=True))
        .replace('__TITLE__', html.escape(title))
        .replace('__COND1__', html.escape(cond1))
        .replace('__COND2__', html.escape(cond2))
        .replace('__REPORT_LABEL__', html.escape(report_label or ''))
        .replace('__REPORT_LABEL_DISPLAY__', 'block' if report_label else 'none'))
    with open(html_out, 'w') as f:
        f.write(html_str)
