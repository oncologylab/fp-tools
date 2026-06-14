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


def _fit_aggregate_normalizers(selected, outdir, aggregate_signals, cond_groups, comparison, flank, mode, logger=None):
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
        bed_path = os.path.join(outdir, prefix, "beds", prefix + "_all.bed")
        for chrom, center in _read_bed_centers(bed_path):
            key = (chrom, center)
            if key not in seen:
                seen.add(key)
                all_centers.append(key)
    if not all_centers:
        return {}

    sample_arrays = [_sample_bigwig_window_values(path, all_centers, flank) for path in aggregate_signals]
    if mode == "sample-quantile":
        return {"mode": mode, "sample": _robust_affine_normalizers(sample_arrays, sample_names)}

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

    if mode == "sample-quantile":
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
    row, comparison, outdir, aggregate_signals, cond_groups, flank, x_len, base, normalization, aggregate_norm_spec, sample_names = task
    c1, c2 = comparison
    prefix = str(row["output_prefix"])
    bed_path = os.path.join(outdir, prefix, "beds", prefix + "_all.bed")
    centers = _read_bed_centers(bed_path)
    if not centers:
        return None

    if not sample_names or len(sample_names) != len(aggregate_signals):
        sample_names = [f"sample_{idx + 1}" for idx in range(len(aggregate_signals))]
    normalized_profiles = {}
    sample_norms = (aggregate_norm_spec or {}).get("sample", {})
    condition_norms = (aggregate_norm_spec or {}).get("condition", {})
    sample_to_condition = {idx: cond for cond, indices in cond_groups.items() for idx in indices}
    for signal_idx, signal_path in enumerate(aggregate_signals):
        sample_name = sample_names[signal_idx]
        norm = None
        if normalization == "sample-quantile":
            norm = sample_norms.get(sample_name) or sample_norms.get(f"sample_{signal_idx + 1}")
        elif normalization == "condition-quantile":
            norm = condition_norms.get(sample_to_condition.get(signal_idx))
        normalized_profiles[sample_name] = np.asarray(_mean_profile(signal_path, centers, flank, norm=norm), dtype=float)

    conditions = []
    for cond in (c1, c2):
        sample_profiles = []
        samples = []
        for signal_idx in cond_groups.get(cond, []):
            sample_name = sample_names[signal_idx]
            sample_profile = normalized_profiles.get(sample_name, np.zeros(x_len, dtype=float))
            sample_profiles.append(sample_profile)
            samples.append({"name": sample_name, "profile": [round(float(v), 6) for v in sample_profile]})
        if sample_profiles:
            profile = [round(float(v), 6) for v in np.nanmean(np.asarray(sample_profiles, dtype=float), axis=0)]
        else:
            profile = [0.0] * x_len
        conditions.append({"name": cond, "profile": profile, "samples": samples})
    return {
        "prefix": prefix,
        "name": str(row.get("name", prefix)),
        "motif_id": str(row.get("motif_id", "")),
        "change": float(row.get(base + "_change", 0.0)),
        "pvalue": float(row.get(base + "_pvalue_numeric", 1.0)),
        "n_sites": len(centers),
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
    if mode == "all":
        selected = rows.sort_values([base + "_pvalue_numeric", base + "_abs_change"], ascending=[True, False])
    elif mode == "top":
        selected = rows.sort_values([base + "_pvalue_numeric", base + "_abs_change"], ascending=[True, False]).head(top_n)
    else:
        threshold = float(getattr(args, "aggregate_pvalue_threshold", 0.05))
        selected = rows[rows[base + "_pvalue_numeric"] <= threshold].sort_values([base + "_pvalue_numeric", base + "_abs_change"], ascending=[True, False])
        if selected.empty:
            selected = rows.sort_values([base + "_pvalue_numeric", base + "_abs_change"], ascending=[True, False]).head(top_n)

    flank = max(1, int(getattr(args, "aggregate_flank", 100)))
    x = list(range(-flank, flank))
    normalization = (getattr(args, "normalization", "none") or "none").replace("_", "-")
    cond_groups = {cond: list(indices) for cond, indices in getattr(args, "cond_groups", {}).items()}
    aggregate_norm_spec = _fit_aggregate_normalizers(
        selected,
        args.outdir,
        list(args.aggregate_signals),
        cond_groups,
        (c1, c2),
        flank,
        normalization,
        logger=getattr(args, "logger", None),
    )
    sample_names = list(getattr(args, "sample_names", []) or [f"sample_{idx + 1}" for idx in range(len(args.aggregate_signals))])
    tasks = [(row.to_dict(), (c1, c2), args.outdir, list(args.aggregate_signals), cond_groups, flank, len(x), base, normalization, aggregate_norm_spec, sample_names) for _, row in selected.iterrows()]

    cores = max(1, int(getattr(args, "cores", 1) or 1))
    if cores > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=min(cores, len(tasks))) as executor:
            payloads = list(executor.map(_aggregate_payload_for_row, tasks))
    else:
        payloads = [_aggregate_payload_for_row(task) for task in tasks]
    motifs_payload = [payload for payload in payloads if payload is not None]
    y_label = "Corrected cut-site signal (a.u.)"
    if normalization == "sample-quantile":
        y_label = "Quantile-scaled corrected cut-site signal (a.u.)"
    elif normalization == "condition-quantile":
        y_label = "Condition-quantile-scaled corrected cut-site signal (a.u.)"
    return {"x": x, "motifs": motifs_payload, "comparison": f"{c1} / {c2}", "normalization": normalization, "x_label": "Distance from motif center (bp)", "y_label": y_label}


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




def plot_interactive_bindetect(motifs, comparison, html_out, aggregate_data=None, title='Differential footprint report'):
    cond1, cond2 = comparison
    groups = [cond1 + '_up', cond2 + '_up', 'n.s.']
    colors = {cond1 + '_up': '#2563eb', cond2 + '_up': '#dc2626', 'n.s.': '#8a94a6'}
    points = []
    for motif in motifs:
        group = getattr(motif, 'group', 'n.s.')
        if group not in colors:
            group = 'n.s.'
        pvalue = max(float(getattr(motif, 'pvalue', 1.0)), 1e-308)
        points.append({
            'prefix': str(getattr(motif, 'prefix', getattr(motif, 'name', ''))),
            'name': str(getattr(motif, 'name', '')),
            'motif_id': str(getattr(motif, 'id', '')),
            'group': group,
            'change': round(float(getattr(motif, 'change', 0.0)), 6),
            'pvalue': pvalue,
            'neglog10p': round(float(-np.log10(pvalue)), 6),
        })
    payload = {
        'title': title,
        'comparison': f'{cond1} / {cond2}',
        'conditions': [cond1, cond2],
        'groups': groups,
        'colors': colors,
        'points': points,
        'logos': _motif_logo_map(motifs),
        'aggregate': aggregate_data or {'motifs': [], 'x': [], 'normalization': 'none'},
    }
    payload_b64 = _compressed_json_b64(payload)
    html_template = '''<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>__TITLE_ATTR__</title><style>
:root{--ink:#152133;--muted:#596579;--line:#d9e2ec;--grid:#e8eef5;--panel:#fff;--bg:#eef3f8;--accent:#173b73;--soft:#f7fafc}*{box-sizing:border-box}body{margin:0;font-family:Arial,Helvetica,sans-serif;background:var(--bg);color:var(--ink);font-weight:700}.wrap{max-width:min(1840px,calc(100vw - 28px));margin:10px auto;padding:0 10px}.panel{background:var(--panel);border:1px solid var(--line);border-radius:8px;box-shadow:0 14px 34px rgba(21,33,51,.10);overflow:hidden}.head{padding:10px 18px 8px;border-bottom:1px solid var(--line);background:linear-gradient(180deg,#fff 0%,#f7fafc 100%)}h1{margin:0;font-size:21px;line-height:1.12;font-weight:900}.sub{margin:3px 0 0;color:var(--muted);font-size:12px;font-weight:700}.top-row{display:grid;grid-template-columns:360px minmax(360px,1fr) 360px 180px;gap:10px;align-items:stretch;padding:9px 12px;border-bottom:1px solid var(--line);background:#fbfdff}.main-row{display:grid;grid-template-columns:minmax(720px,1fr) 520px;gap:0;align-items:stretch}.main-row.aggregate-wide{grid-template-columns:minmax(720px,1fr) 760px}.main-row.aggregate-full{grid-template-columns:1fr}.main-row.aggregate-full .chart-box{border-right:0;border-bottom:1px solid var(--line)}.chart-box{padding:8px 12px 10px;border-right:1px solid var(--line)}#chart{width:100%;height:min(640px,calc(100vh - 260px));min-height:560px;display:block;background:#fff}.aggregate-pane{padding:8px 12px 10px;background:#fff}#aggregate-chart{width:100%;height:315px;display:block;background:#fff}.section-title{font-size:11px;line-height:1.1;text-transform:uppercase;letter-spacing:.08em;color:#728197;margin:0 0 5px;font-weight:900}#logo-title{text-transform:none;letter-spacing:0}.card{border:1px solid var(--line);border-radius:7px;background:#fff;padding:7px;min-height:86px}.controls{display:flex;flex-wrap:wrap;gap:6px;align-content:flex-start}.color-row{display:flex;align-items:center;justify-content:space-between;gap:6px;font-size:12px;color:#334e68;font-weight:800;border:1px solid #e6edf5;border-radius:999px;padding:4px 6px;background:#fbfdff}.color-row input{width:28px;height:20px;border:1px solid var(--line);background:#fff;padding:0;border-radius:4px}.detail h2{margin:0 0 4px;font-size:16px;line-height:1.05;font-weight:900}.detail-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:6px}.detail p{margin:0;font-size:11px;color:var(--muted);font-weight:800;line-height:1.18}.logo{height:86px;min-height:86px;display:flex;align-items:center;justify-content:center;overflow:hidden}.logo svg,.logo img{max-width:100%;height:auto;display:block}.logo-empty{color:#728197;font-size:13px;font-weight:800}.button-row{display:grid;gap:5px}button{border:1px solid #b8c5d6;background:#fff;color:var(--accent);border-radius:6px;padding:5px 8px;font-family:Arial,Helvetica,sans-serif;font-size:11px;font-weight:900;cursor:pointer}button:hover{background:#f2f6fb}.agg-toolbar{display:flex;align-items:center;justify-content:space-between;gap:8px;margin-bottom:8px}.agg-toolbar label{font-size:11px;line-height:1;color:var(--muted);font-weight:900;text-transform:uppercase;letter-spacing:.06em}.agg-toolbar select{border:1px solid var(--line);border-radius:6px;background:#fff;color:var(--ink);font-family:Arial,Helvetica,sans-serif;font-size:12px;font-weight:800;padding:5px 8px}.combo{position:relative;margin-bottom:10px}.combo input{width:100%;border:1px solid var(--line);border-radius:6px;padding:6px 9px;font-family:Arial,Helvetica,sans-serif;font-size:12px;font-weight:800;color:var(--ink);background:#fff}.combo-list{display:none;position:absolute;z-index:20;left:0;right:0;top:32px;max-height:240px;overflow:auto;border:1px solid var(--line);border-radius:6px;background:#fff;box-shadow:0 10px 24px rgba(21,33,51,.16)}.combo-list.open{display:block}.combo-option{padding:7px 10px;font-size:13px;font-weight:800;cursor:pointer;border-bottom:1px solid #eef3f8}.combo-option small{display:block;color:#728197;font-size:11px;font-weight:700}.combo-option:hover,.combo-option.active{background:#edf4ff}.axis{stroke:#3b4552;stroke-width:1.35}.grid{stroke:var(--grid);stroke-width:1}.zero{stroke:#677386;stroke-width:1.4;stroke-dasharray:4 4}.tick{font-family:Arial,Helvetica,sans-serif;font-size:12px;fill:var(--muted);font-weight:800}.axis-label{font-family:Arial,Helvetica,sans-serif;font-size:14px;fill:var(--ink);font-weight:900}.plot-title{font-family:Arial,Helvetica,sans-serif;font-size:15px;fill:var(--ink);font-weight:900}.pt{cursor:pointer}.pt:hover{stroke:#111827;stroke-width:1.5}.pt.selected{stroke:#111827;stroke-width:2.6;fill-opacity:.95}@media(max-width:1260px){.top-row{grid-template-columns:1fr 1fr}.main-row{grid-template-columns:1fr}.chart-box{border-right:0;border-bottom:1px solid var(--line)}#chart{height:auto;min-height:0}}@media(max-width:760px){.top-row{grid-template-columns:1fr}.detail-grid{grid-template-columns:1fr}}
</style></head><body><div class="wrap"><div class="panel"><div class="head"><h1>__TITLE__</h1><p class="sub">__COND1__ / __COND2__</p></div><div class="top-row"><div><p class="section-title">Groups</p><div class="card controls" id="color-controls"></div></div><div><p class="section-title">Selected motif</p><div class="card detail" id="detail"><h2>Loading report</h2><div class="detail-grid"><p>Data are embedded in this standalone HTML file.</p></div></div></div><div><p class="section-title" id="logo-title">Motif logo</p><div class="card logo" id="logo-box"><span class="logo-empty">Motif logo</span></div></div><div><p class="section-title">Export editable SVG</p><div class="card button-row"><button id="download-volcano">Download volcano SVG</button><button id="download-aggregate">Download aggregate SVG</button><button id="download-logo">Download motif logo SVG</button></div></div></div><div class="main-row" id="main-row"><div class="chart-box"><svg id="chart" viewBox="0 0 980 620" aria-label="Differential footprint volcano plot"></svg></div><div class="aggregate-pane"><div class="agg-toolbar"><p class="section-title">Aggregate profile</p><label>Width <select id="aggregate-width"><option value="normal">Normal</option><option value="wide">Wide</option><option value="full">Full width</option></select></label></div><div class="combo" id="aggregate-combo" style="display:none"><input id="aggregate-search" type="text" autocomplete="off" placeholder="Search motif"><div class="combo-list" id="aggregate-options"></div></div><svg id="aggregate-chart" viewBox="0 0 520 315" aria-label="Aggregate footprint profile"></svg></div></div></div></div><script>
const reportPayloadB64="__PAYLOAD__";let payload=null,selectedPrefix=null,activeOptionIndex=0,sortedAggregateMotifs=[];const chart=document.getElementById('chart'),aggregateChart=document.getElementById('aggregate-chart'),detail=document.getElementById('detail'),logoBox=document.getElementById('logo-box'),logoTitle=document.getElementById('logo-title'),colorControls=document.getElementById('color-controls'),mainRow=document.getElementById('main-row'),aggregateWidth=document.getElementById('aggregate-width'),aggregateCombo=document.getElementById('aggregate-combo'),aggregateSearch=document.getElementById('aggregate-search'),aggregateOptions=document.getElementById('aggregate-options');
function escText(value){return String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]))}function b64ToBytes(b64){return Uint8Array.from(atob(b64),c=>c.charCodeAt(0))}async function decodePayload(){if(!('DecompressionStream'in window))throw new Error('This standalone report needs a modern browser with gzip DecompressionStream support.');const ds=new DecompressionStream('gzip');const stream=new Blob([b64ToBytes(reportPayloadB64)]).stream().pipeThrough(ds);return JSON.parse(await new Response(stream).text())}function fmtTick(value){return Math.abs(value)>=1?value.toFixed(1).replace('-0.0','0.0'):value.toFixed(2).replace('-0.00','0.00')}function niceTicks(min,max,n){const out=[];for(let i=0;i<n;i++)out.push(min+(max-min)*(i/Math.max(1,n-1)));return out}function currentGroupColors(){const out={...payload.colors};document.querySelectorAll('[data-color-group]').forEach(inp=>out[inp.dataset.colorGroup]=inp.value);return out}function currentConditionColors(){const groupColors=currentGroupColors();return{[payload.conditions[0]]:groupColors[payload.conditions[0]+'_up'],[payload.conditions[1]]:groupColors[payload.conditions[1]+'_up']}}function pointByPrefix(prefix){return payload.points.find(p=>p.prefix===prefix)}function aggregateByPrefix(prefix){return(payload.aggregate.motifs||[]).find(m=>m.prefix===prefix)}
function renderColorControls(){colorControls.innerHTML=payload.groups.map(group=>`<label class="color-row"><span>${escText(group)}</span><input type="color" data-color-group="${escText(group)}" value="${payload.colors[group]}"></label>`).join('');colorControls.querySelectorAll('input').forEach(inp=>inp.addEventListener('input',()=>{renderVolcano();renderAggregate(selectedPrefix)}))}function setAggregateLayout(mode){mainRow.classList.toggle('aggregate-wide',mode==='wide');mainRow.classList.toggle('aggregate-full',mode==='full');renderAggregate(selectedPrefix)}aggregateWidth.addEventListener('change',()=>setAggregateLayout(aggregateWidth.value));
function renderVolcano(){const colors=currentGroupColors(),width=980,height=620,margin={top:58,right:54,bottom:72,left:92},innerW=width-margin.left-margin.right,innerH=390,plotX0=margin.left,plotY0=margin.top,plotX1=plotX0+innerW,plotY1=plotY0+innerH;const xs=payload.points.map(p=>p.change),ys=payload.points.map(p=>p.neglog10p),xabs=Math.max(1,Math.abs(Math.min(...xs,0)),Math.abs(Math.max(...xs,0)))*1.1,ymin=0,ymax=Math.max(1,Math.max(...ys,1)*1.08);const sx=x=>plotX0+((x+xabs)/(2*xabs))*innerW,sy=y=>plotY1-((y-ymin)/(ymax-ymin||1))*innerH,xTicks=niceTicks(-xabs,xabs,7),yTicks=niceTicks(ymin,ymax,7);const parts=[`<rect width="${width}" height="${height}" fill="#ffffff"/>`,`<text x="${(plotX0+plotX1)/2}" y="32" class="plot-title" text-anchor="middle">Differential footprint evidence</text>`,`<rect x="${plotX0}" y="${plotY0}" width="${innerW}" height="${innerH}" fill="#fbfdff" stroke="#d9e2ec"/>`];yTicks.forEach(v=>parts.push(`<line x1="${plotX0}" y1="${sy(v)}" x2="${plotX1}" y2="${sy(v)}" class="grid"/>`,`<text x="${plotX0-12}" y="${sy(v)+4}" class="tick" text-anchor="end">${v.toFixed(1)}</text>`));xTicks.forEach(v=>parts.push(`<line x1="${sx(v)}" y1="${plotY0}" x2="${sx(v)}" y2="${plotY1}" class="grid"/>`,`<text x="${sx(v)}" y="${plotY1+25}" class="tick" text-anchor="middle">${fmtTick(v)}</text>`));parts.push(`<line x1="${sx(0)}" y1="${plotY0}" x2="${sx(0)}" y2="${plotY1}" class="zero"/>`,`<line x1="${plotX0}" y1="${plotY1}" x2="${plotX1}" y2="${plotY1}" class="axis"/>`,`<line x1="${plotX0}" y1="${plotY0}" x2="${plotX0}" y2="${plotY1}" class="axis"/>`,`<text x="${(plotX0+plotX1)/2}" y="${plotY1+62}" class="axis-label" text-anchor="middle">Differential footprint score</text>`,`<text x="28" y="${plotY0+innerH/2}" class="axis-label" text-anchor="middle" transform="rotate(-90 28 ${plotY0+innerH/2})">-log10(p-value)</text>`);payload.points.forEach((p,idx)=>{const selected=p.prefix===selectedPrefix;parts.push(`<circle class="pt${selected?' selected':''}" data-prefix="${escText(p.prefix)}" data-index="${idx}" cx="${sx(p.change).toFixed(2)}" cy="${sy(p.neglog10p).toFixed(2)}" r="${selected?6.0:4.3}" fill="${colors[p.group]||colors['n.s.']}" fill-opacity="${selected?.95:.76}" stroke="#ffffff" stroke-width="0.9"/>`)});chart.innerHTML=parts.join('');chart.querySelectorAll('.pt').forEach(el=>el.addEventListener('click',()=>setSelectedMotif(el.dataset.prefix,{from:'volcano'})))}
function motifLabel(item){if(!item)return'';const id=item.motif_id||item.id||'';return id?`${item.name} (${id})`:item.name}function renderDetail(point){if(!point){detail.innerHTML='<h2>No motif selected</h2><div class="detail-grid"><p>Select a motif from the volcano or aggregate search.</p></div>';return}detail.innerHTML=`<h2>${escText(motifLabel(point))}</h2><div class="detail-grid"><p><strong>Group:</strong><br>${escText(point.group)}</p><p><strong>Change:</strong><br>${Number(point.change).toFixed(4)}</p><p><strong>P-value:</strong><br>${Number(point.pvalue).toExponential(3)}</p></div>`}function renderLogo(prefix){const logo=payload.logos[prefix],agg=aggregateByPrefix(prefix),point=pointByPrefix(prefix),label=motifLabel(point)||motifLabel(agg)||'';logoTitle.textContent=label?`${label} Motif logo`:'Motif logo';if(!logo){logoBox.innerHTML='<span class="logo-empty">Motif logo unavailable</span>';return}if(logo.svg)logoBox.innerHTML=logo.svg;else logoBox.innerHTML=`<img alt="Motif logo" src="${logo.png}">`}
function renderAggregate(prefix){const agg=payload.aggregate||{motifs:[],x:[]};if(!agg.motifs||agg.motifs.length===0){aggregateChart.innerHTML='<text x="260" y="180" text-anchor="middle" class="tick">Aggregate profiles unavailable</text>';return}const motif=aggregateByPrefix(prefix)||agg.motifs[0];const mode=aggregateWidth.value,width=mode==='full'?980:(mode==='wide'?760:520),height=315;aggregateChart.setAttribute('viewBox',`0 0 ${width} ${height}`);const x=agg.x,margin={top:44,right:26,bottom:62,left:92},innerW=width-margin.left-margin.right,innerH=height-margin.top-margin.bottom;const sampleProfiles=motif.conditions.flatMap(c=>(c.samples||[]).flatMap(s=>s.profile));const meanProfiles=motif.conditions.flatMap(c=>c.profile);const allY=[...meanProfiles,...sampleProfiles].filter(Number.isFinite);let ymin=Math.min(...allY,0),ymax=Math.max(...allY,1e-9);const minPad=Math.max(1e-4,Math.abs(ymax)*.05,Math.abs(ymin)*.05);const pad=Math.max((ymax-ymin||1)*.22,minPad);ymin-=pad;ymax+=pad;const sx=v=>margin.left+((v-x[0])/(x[x.length-1]-x[0]||1))*innerW,sy=v=>margin.top+innerH-((v-ymin)/(ymax-ymin||1))*innerH,colors=currentConditionColors(),xTicks=[x[0],Math.round(x[0]/2),0,Math.round(x[x.length-1]/2),x[x.length-1]],yTicks=niceTicks(ymin,ymax,5),lineD=profile=>profile.map((y,i)=>`${i===0?'M':'L'}${sx(x[i]).toFixed(2)},${sy(y).toFixed(2)}`).join(' ');const parts=[`<rect width="${width}" height="${height}" fill="#ffffff"/>`,`<text x="${width/2}" y="22" class="plot-title" text-anchor="middle">${escText(motifLabel(motif))} (${motif.n_sites} sites)</text>`,`<text x="${width/2}" y="37" class="tick" text-anchor="middle">${escText(agg.normalization||'none')} normalization</text>`];yTicks.forEach(v=>parts.push(`<line x1="${margin.left}" y1="${sy(v)}" x2="${margin.left+innerW}" y2="${sy(v)}" class="grid"/>`,`<text x="${margin.left-10}" y="${sy(v)+4}" class="tick" text-anchor="end">${v.toPrecision(2)}</text>`));xTicks.forEach(v=>parts.push(`<line x1="${sx(v)}" y1="${margin.top}" x2="${sx(v)}" y2="${margin.top+innerH}" class="grid"/>`,`<text x="${sx(v)}" y="${margin.top+innerH+25}" class="tick" text-anchor="middle">${v}</text>`));parts.push(`<line x1="${sx(0)}" y1="${margin.top}" x2="${sx(0)}" y2="${margin.top+innerH}" class="zero"/>`,`<line x1="${margin.left}" y1="${margin.top+innerH}" x2="${margin.left+innerW}" y2="${margin.top+innerH}" class="axis"/>`,`<line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${margin.top+innerH}" class="axis"/>`);motif.conditions.forEach((cond,idx)=>{const color=colors[cond.name]||['#2563eb','#dc2626','#16a34a','#9333ea'][idx%4];(cond.samples||[]).forEach(sample=>{parts.push(`<path d="${lineD(sample.profile)}" fill="none" stroke="${color}" stroke-width="0.9" stroke-opacity="0.30"><title>${escText(sample.name)}</title></path>`)});});motif.conditions.forEach((cond,idx)=>{const color=colors[cond.name]||['#2563eb','#dc2626','#16a34a','#9333ea'][idx%4];parts.push(`<path d="${lineD(cond.profile)}" fill="none" stroke="${color}" stroke-width="2.2"/>`,`<text x="${margin.left+8}" y="${margin.top+16+idx*16}" font-family="Arial,Helvetica,sans-serif" font-size="12" font-weight="900" fill="${color}">${escText(cond.name)} mean</text>`)});parts.push(`<text x="${margin.left+innerW/2}" y="${height-16}" class="axis-label" text-anchor="middle">${escText(agg.x_label||'Distance from motif center (bp)')}</text>`,`<text x="22" y="${margin.top+innerH/2}" class="axis-label" text-anchor="middle" transform="rotate(-90 22 ${margin.top+innerH/2})">${escText(agg.y_label||'Corrected cut-site signal (a.u.)')}</text>`);aggregateChart.innerHTML=parts.join('')}
function filteredAggregateMotifs(){const q=aggregateSearch.value.trim().toLowerCase();const list=sortedAggregateMotifs.filter(m=>!q||(`${m.name} ${m.motif_id||''} ${m.prefix}`).toLowerCase().includes(q));return list.slice(0,80)}function renderAggregateOptions(){const list=filteredAggregateMotifs();activeOptionIndex=Math.min(activeOptionIndex,Math.max(0,list.length-1));aggregateOptions.innerHTML=list.map((m,idx)=>`<div class="combo-option${idx===activeOptionIndex?' active':''}" data-prefix="${escText(m.prefix)}"><span>${escText(m.name)}</span><small>${escText(m.prefix)} · ${m.n_sites||0} sites</small></div>`).join('');aggregateOptions.classList.add('open');aggregateOptions.querySelectorAll('.combo-option').forEach(el=>el.addEventListener('mousedown',ev=>{ev.preventDefault();setSelectedMotif(el.dataset.prefix,{from:'search'});aggregateOptions.classList.remove('open')}))}function setSearchValue(prefix){const motif=aggregateByPrefix(prefix);if(motif)aggregateSearch.value=motif.name}function setupAggregateSearch(){sortedAggregateMotifs=[...(payload.aggregate.motifs||[])].sort((a,b)=>(a.name||'').localeCompare(b.name||'')||(a.prefix||'').localeCompare(b.prefix||''));if(!sortedAggregateMotifs.length)return;aggregateCombo.style.display='block';aggregateSearch.addEventListener('focus',()=>{activeOptionIndex=0;renderAggregateOptions()});aggregateSearch.addEventListener('input',()=>{activeOptionIndex=0;renderAggregateOptions()});aggregateSearch.addEventListener('keydown',ev=>{const list=filteredAggregateMotifs();if(ev.key==='ArrowDown'){ev.preventDefault();activeOptionIndex=Math.min(activeOptionIndex+1,Math.max(0,list.length-1));renderAggregateOptions()}else if(ev.key==='ArrowUp'){ev.preventDefault();activeOptionIndex=Math.max(activeOptionIndex-1,0);renderAggregateOptions()}else if(ev.key==='Enter'){ev.preventDefault();if(list[activeOptionIndex]){setSelectedMotif(list[activeOptionIndex].prefix,{from:'search'});aggregateOptions.classList.remove('open')}}else if(ev.key==='Escape'){aggregateOptions.classList.remove('open')}});document.addEventListener('click',ev=>{if(!aggregateCombo.contains(ev.target))aggregateOptions.classList.remove('open')})}
function setSelectedMotif(prefix,opts={}){const agg=aggregateByPrefix(prefix);const point=pointByPrefix(prefix)||payload.points.find(p=>p.name===(agg&&agg.name))||payload.points[0];const selected=(agg&&agg.prefix)||(point&&point.prefix);if(!selected)return;selectedPrefix=selected;renderDetail(point);renderLogo(selected);renderVolcano();renderAggregate(selected);if(opts.from!=='search')setSearchValue(selected)}function svgBlob(svgNode){const clone=svgNode.cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');const text=new XMLSerializer().serializeToString(clone);return new Blob([text],{type:'image/svg+xml;charset=utf-8'})}function downloadBlob(blob,filename){const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=filename;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000)}document.getElementById('download-volcano').addEventListener('click',()=>downloadBlob(svgBlob(chart),'diff_footprints_volcano.svg'));document.getElementById('download-aggregate').addEventListener('click',()=>downloadBlob(svgBlob(aggregateChart),'diff_footprints_aggregate.svg'));document.getElementById('download-logo').addEventListener('click',()=>{const svg=logoBox.querySelector('svg');if(svg)downloadBlob(svgBlob(svg),'diff_footprints_motif_logo.svg')});
decodePayload().then(data=>{payload=data;renderColorControls();setupAggregateSearch();const first=(payload.aggregate.motifs&&payload.aggregate.motifs[0]&&payload.aggregate.motifs[0].prefix)||(payload.points[0]&&payload.points[0].prefix);if(first)setSelectedMotif(first);else renderVolcano()}).catch(err=>{detail.innerHTML=`<h2>Could not open report payload</h2><div class="detail-grid"><p>${escText(err.message)}</p></div>`});
</script></body></html>'''
    html_str = (html_template
        .replace('__PAYLOAD__', payload_b64)
        .replace('__TITLE_ATTR__', html.escape(f'{title} {cond1} / {cond2}', quote=True))
        .replace('__TITLE__', html.escape(title))
        .replace('__COND1__', html.escape(cond1))
        .replace('__COND2__', html.escape(cond2)))
    with open(html_out, 'w') as f:
        f.write(html_str)
