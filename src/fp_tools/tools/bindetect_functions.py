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
import html
import json
import random
import itertools
from datetime import datetime
import os

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


def plot_interactive_bindetect(motifs, comparison, html_out):
    cond1, cond2 = comparison
    groups = [cond1 + "_up", cond2 + "_up", "n.s."]
    colors = {
        cond1 + "_up": "#2ca25f",
        cond2 + "_up": "#de2d26",
        "n.s.": "#9e9e9e",
    }
    points_by_group = {group: [] for group in groups}
    for group in groups:
        for motif in motifs:
            if getattr(motif, "group", "n.s.") != group:
                continue
            points_by_group[group].append({
                "group": group,
                "color": colors[group],
                "x": float(getattr(motif, "change", 0)),
                "y": float(-np.log10(max(float(getattr(motif, "pvalue", 1)), 1e-308))),
                "name": motif.name,
                "base": ("data:image/png;base64," + getattr(motif, "base", "")) if getattr(motif, "base", "") else "",
            })

    points = [point for group in groups for point in points_by_group[group]]
    xvals = [p["x"] for p in points] or [0.0]
    yvals = [p["y"] for p in points] or [1.0]
    xabs = max(abs(min(xvals)), abs(max(xvals)), 1.0)
    xlim = xabs * 1.1
    ymin = 0.0
    ymax = max(yvals) * 1.08 if max(yvals) > 0 else 1.0

    width = 980
    height = 720
    margin = {"top": 72, "right": 70, "bottom": 90, "left": 90}
    inner_w = width - margin["left"] - margin["right"]
    inner_h = 430
    plot_x0 = margin["left"]
    plot_y0 = margin["top"]
    plot_x1 = plot_x0 + inner_w
    plot_y1 = plot_y0 + inner_h

    def sx(x):
        return plot_x0 + ((x + xlim) / (2 * xlim)) * inner_w

    def sy(y):
        return plot_y1 - ((y - ymin) / (ymax - ymin or 1.0)) * inner_h

    def fmt_tick(value):
        if abs(value) >= 1:
            out = f"{value:.1f}"
        else:
            out = f"{value:.2f}"
        return out.replace("-0.00", "0.00").replace("-0.0", "0.0")

    def esc(text):
        return html.escape(str(text), quote=True)

    y_ticks = []
    for i in range(7):
        yval = ymin + (ymax - ymin) * (i / 6.0)
        y_ticks.append((yval, sy(yval)))

    x_ticks = []
    for i in range(7):
        xval = -xlim + (2 * xlim) * (i / 6.0)
        x_ticks.append((xval, sx(xval)))

    grid_parts = []
    for _, ypos in y_ticks:
        grid_parts.append(
            f'<line x1="{plot_x0:.2f}" y1="{ypos:.2f}" x2="{plot_x1:.2f}" y2="{ypos:.2f}" class="grid" />'
        )
    for _, xpos in x_ticks:
        grid_parts.append(
            f'<line x1="{xpos:.2f}" y1="{plot_y0:.2f}" x2="{xpos:.2f}" y2="{plot_y1:.2f}" class="grid" />'
        )
    zero_x = sx(0.0)
    grid_parts.append(
        f'<line x1="{zero_x:.2f}" y1="{plot_y0:.2f}" x2="{zero_x:.2f}" y2="{plot_y1:.2f}" class="zero" />'
    )

    label_parts = []
    for yval, ypos in y_ticks:
        label_parts.append(
            f'<text x="{plot_x0 - 12}" y="{ypos + 4:.2f}" class="tick" text-anchor="end">{esc(f"{yval:.1f}")}</text>'
        )
    for xval, xpos in x_ticks:
        label_parts.append(
            f'<text x="{xpos:.2f}" y="{plot_y1 + 24}" class="tick" text-anchor="middle">{esc(fmt_tick(xval))}</text>'
        )

    point_parts = []
    for group in groups:
        for point in points_by_group[group]:
            point_parts.append(
                "<circle "
                f'cx="{sx(point["x"]):.2f}" '
                f'cy="{sy(point["y"]):.2f}" '
                'r="4.2" '
                f'fill="{point["color"]}" '
                'fill-opacity="0.72" '
                'stroke="#ffffff" stroke-width="0.8" class="pt" '
                f'data-group="{esc(point["group"])}" '
                f'data-name="{esc(point["name"])}" '
                f'data-change="{point["x"]:.5f}" '
                f'data-pvalue="{10 ** (-point["y"]):.12g}" '
                f'data-logo="{esc(point["base"])}"'
                " />"
            )

    legend_items = []
    for label in groups:
        legend_items.append(
            '<div class="legend-item">'
            f'<span class="swatch" style="background:{colors[label]}"></span>'
            f'<span>{esc(label)}</span>'
            '</div>'
        )

    html_str = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>BINDetect """ + f"{cond1} / {cond2}" + """</title>
  <style>
    body { margin: 0; font-family: Arial, Helvetica, sans-serif; background: #edf2f7; color: #152133; font-weight: 700; }
    .wrap { max-width: 1080px; margin: 24px auto; padding: 0 18px; }
    .panel {
      background: #ffffff;
      border: 1px solid #d9e2ec;
      border-radius: 14px;
      box-shadow: 0 10px 30px rgba(21, 33, 51, 0.08);
      overflow: hidden;
    }
    .head {
      padding: 18px 24px 12px;
      border-bottom: 1px solid #e7edf4;
      background: linear-gradient(180deg, #ffffff 0%, #f7fafc 100%);
    }
    h1 { margin: 0; font-size: 24px; font-weight: 700; }
    .sub { margin: 6px 0 0; color: #52606d; font-size: 14px; }
    .body {
      display: grid;
      grid-template-columns: 1fr 260px;
      gap: 0;
      align-items: start;
    }
    .chart-box { padding: 10px 16px 18px; }
    #chart { width: 100%; height: auto; display: block; }
    .grid { stroke: #e7edf4; stroke-width: 1; }
    .axis { stroke: #5b6875; stroke-width: 1.2; }
    .zero { stroke: #7b8794; stroke-width: 1.4; stroke-dasharray: 4 4; }
    .tick { font-size: 12px; fill: #52606d; font-weight: 700; }
    .axis-label { font-size: 14px; fill: #152133; font-weight: 700; }
    .pt { cursor: pointer; }
    .pt:hover { stroke: #152133; stroke-width: 1.4; }
    .side {
      border-left: 1px solid #e7edf4;
      padding: 18px 18px 20px;
      background: #fbfdff;
    }
    .legend {
      display: flex;
      flex-direction: column;
      gap: 8px;
      margin-bottom: 16px;
    }
    .legend-item {
      display: flex;
      align-items: center;
      gap: 10px;
      font-size: 13px;
      color: #334e68;
      font-weight: 700;
    }
    .swatch {
      width: 12px;
      height: 12px;
      border-radius: 999px;
      display: inline-block;
      border: 1px solid rgba(0,0,0,0.12);
    }
    .meta-title {
      font-size: 12px;
      text-transform: uppercase;
      letter-spacing: 0.08em;
      color: #7b8794;
      margin: 0 0 8px;
      font-weight: 700;
    }
    .detail {
      min-height: 146px;
      border: 1px solid #d9e2ec;
      border-radius: 12px;
      background: #ffffff;
      padding: 12px;
      box-sizing: border-box;
    }
    .detail h2 {
      margin: 0 0 6px;
      font-size: 16px;
      line-height: 1.2;
    }
    .detail p {
      margin: 4px 0;
      font-size: 13px;
      color: #52606d;
      font-weight: 700;
    }
    .logo {
      margin-top: 14px;
      border: 1px solid #d9e2ec;
      border-radius: 12px;
      background: #ffffff;
      min-height: 140px;
      display: flex;
      align-items: center;
      justify-content: center;
      overflow: hidden;
      padding: 10px;
      box-sizing: border-box;
    }
    .logo img {
      max-width: 100%;
      max-height: 120px;
      display: none;
    }
    .logo span {
      color: #7b8794;
      font-size: 13px;
      font-weight: 700;
    }
    @media (max-width: 920px) {
      .body { grid-template-columns: 1fr; }
      .side { border-left: 0; border-top: 1px solid #e7edf4; }
    }
  </style>
</head>
<body>
  <div class="wrap">
    <div class="panel">
      <div class="head">
        <h1>BINDetect Volcano Plot</h1>
        <p class="sub">""" + f"{cond1} / {cond2}" + """</p>
      </div>
      <div class="body">
        <div class="chart-box">
          <svg id="chart" viewBox="0 0 980 720" aria-label="BINDetect volcano plot">
            <rect x=\"0\" y=\"0\" width=\"980\" height=\"720\" fill=\"#ffffff\" />
            <rect x=\"""" + f"{plot_x0:.2f}" + """\" y=\"""" + f"{plot_y0:.2f}" + """\" width=\"""" + f"{inner_w:.2f}" + """\" height=\"""" + f"{inner_h:.2f}" + """\" fill=\"#fbfdff\" stroke=\"#d9e2ec\" />
            """ + "".join(grid_parts) + """
            <line x1=\"""" + f"{plot_x0:.2f}" + """\" y1=\"""" + f"{plot_y1:.2f}" + """\" x2=\"""" + f"{plot_x1:.2f}" + """\" y2=\"""" + f"{plot_y1:.2f}" + """\" class=\"axis\" />
            <line x1=\"""" + f"{plot_x0:.2f}" + """\" y1=\"""" + f"{plot_y0:.2f}" + """\" x2=\"""" + f"{plot_x0:.2f}" + """\" y2=\"""" + f"{plot_y1:.2f}" + """\" class=\"axis\" />
            """ + "".join(label_parts) + """
            <text x=\"""" + f"{(plot_x0 + plot_x1) / 2:.2f}" + """\" y=\"""" + f"{plot_y1 + 56:.2f}" + """\" class=\"axis-label\" text-anchor=\"middle\">Differential binding score</text>
            <text x=\"26\" y=\"""" + f"{plot_y0 + inner_h / 2:.2f}" + """\" class=\"axis-label\" text-anchor=\"middle\" transform=\"rotate(-90 26 """ + f"{plot_y0 + inner_h / 2:.2f}" + """)\">-log10(pvalue)</text>
            """ + "".join(point_parts) + """
          </svg>
        </div>
        <aside class="side">
          <p class="meta-title">Groups</p>
          <div class="legend">""" + "".join(legend_items) + """</div>
          <p class="meta-title">Selected motif</p>
          <div class="detail" id="detail">
            <h2>Hover over a point</h2>
            <p>Motif name, group, change, and p-value appear here.</p>
          </div>
          <div class="logo" id="logo-box">
            <img id="logo-img" alt="Motif logo">
            <span id="logo-empty">Motif logo</span>
          </div>
        </aside>
      </div>
    </div>
  </div>
  <script>
    const detail = document.getElementById('detail');
    const logoImg = document.getElementById('logo-img');
    const logoEmpty = document.getElementById('logo-empty');
    document.querySelectorAll('.pt').forEach((el) => {
      el.addEventListener('mouseenter', () => {
        const name = el.dataset.name;
        const group = el.dataset.group;
        const change = Number(el.dataset.change);
        const pvalue = Number(el.dataset.pvalue);
        detail.innerHTML =
          `<h2>${name}</h2>` +
          `<p><strong>Group:</strong> ${group}</p>` +
          `<p><strong>Change:</strong> ${change.toFixed(4)}</p>` +
          `<p><strong>P-value:</strong> ${pvalue.toExponential(3)}</p>`;
        if (el.dataset.logo) {
          logoImg.src = el.dataset.logo;
          logoImg.style.display = 'block';
          logoEmpty.style.display = 'none';
        } else {
          logoImg.removeAttribute('src');
          logoImg.style.display = 'none';
          logoEmpty.style.display = 'block';
        }
      });
    });
  </script>
</body>
</html>"""
    with open(html_out, "w") as f:
        f.write(html_str)
