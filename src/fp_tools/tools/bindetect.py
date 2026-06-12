#!/usr/bin/env python
"""
BINDetect command driver for TF binding and differential binding analysis.

This module handles:
- motif/scored-signal integration
- bound versus unbound site calling
- differential binding statistics between conditions
- summary tables, PDFs, and interactive HTML outputs

It also includes replicate grouping support and skewness report integration.
"""

import os
import sys
import argparse
import time
import numpy as np
import multiprocessing as mp
import itertools
import pandas as pd
import seaborn as sns
from collections import Counter
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
from fp_tools.parsers import add_bindetect_arguments
from fp_tools.tools.bindetect_functions import *
from fp_tools.tools import bindetect_skew_report as skewrep
from fp_tools.tools.bindetect_replicate_report import build_replicate_report

from fp_tools.utils.utilities import (
    check_required, check_files, make_directory, merge_dicts,
    monitor_progress, expand_dirs, check_cores, file_writer
)
from fp_tools.utils.regions import *
from fp_tools.utils.motifs import *
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.plotting_style import PDF_FONT_SIZE, apply_pdf_style, apply_ascii_minus_to_figure

# tame some noisy warnings during curve fitting
from scipy.optimize import OptimizeWarning

warnings.simplefilter("ignore", OptimizeWarning)
warnings.simplefilter("ignore", RuntimeWarning)
apply_pdf_style()


def norm_fit(x, mean, std, scale):
    return scale * scipy.stats.norm.pdf(x, mean, std)


# ----------------------------------------------------------------------------- #
def run_bindetect(args):
    """Run the BINDetect pipeline from parsed CLI arguments."""
    check_required(args, ["signals", "motifs", "genome", "peaks"])

    # derive condition names if not given (one per input bigwig)
    args.cond_names = (
        [os.path.basename(os.path.splitext(bw)[0]) for bw in args.signals]
        if args.cond_names is None else args.cond_names
    )
    args.outdir = os.path.abspath(args.outdir)

    # ----------------------------- replicate grouping ------------------------- #
    # Allow passing identical names in --cond-names to denote replicates.
    # Build: { condition -> [indices into args.signals] }
    orig = list(args.cond_names)
    if len(orig) != len(args.signals):
        raise ValueError("--cond-names must have the same length as --signals")
    if getattr(args, "norm_off", False):
        args.normalization = "none"
    idxs = {}
    for i, nm in enumerate(orig):
        idxs.setdefault(nm, []).append(i)
    args.cond_groups = idxs                   # condition -> signal indices
    args.cond_names = list(idxs.keys())       # unique condition list
    args.condition_replicates = {cond: len(indices) for cond, indices in idxs.items()}
    args.sample_names = []
    args.sample_to_condition = {}
    args.condition_samples = {cond: [] for cond in args.cond_names}
    for cond, indices in idxs.items():
        for rep_no, signal_idx in enumerate(indices, start=1):
            sample_name = f"{cond}_rep{rep_no}"
            args.sample_names.append(sample_name)
            args.sample_to_condition[sample_name] = cond
            args.condition_samples[cond].append(sample_name)
    # ------------------------------------------------------------------------- #

    # outputs we’ll create
    states = ["bound", "unbound"]
    outfiles = [os.path.abspath(os.path.join(
        args.outdir, "*", "beds", f"*_{cond}_{state}.bed"))
        for (cond, state) in itertools.product(args.cond_names, states)]
    outfiles += [
        os.path.abspath(os.path.join(args.outdir, "*", "beds", "*_all.bed")),
        os.path.abspath(os.path.join(args.outdir, "*", "plots", "*_log2fcs.pdf")),
        os.path.abspath(os.path.join(args.outdir, "*", "*_overview.txt")),
        os.path.abspath(os.path.join(args.outdir, "*", "*_overview.xlsx")),
        os.path.abspath(os.path.join(args.outdir, args.prefix + "_distances.txt")),
        os.path.abspath(os.path.join(args.outdir, args.prefix + "_results.txt")),
        os.path.abspath(os.path.join(args.outdir, args.prefix + "_results.xlsx")),
        os.path.abspath(os.path.join(args.outdir, args.prefix + "_figures.pdf")),
        os.path.abspath(os.path.join(args.outdir, args.prefix + "_clusters.pdf")),
    ]

    # ------------------------------ logger/pools ------------------------------ #
    logger = FpToolsLogger("BINDetect", args.verbosity)
    logger.begin()
    parser = add_bindetect_arguments(argparse.ArgumentParser())
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
    check_files([args.signals, args.motifs, args.genome, args.peaks], action="r")
    check_files(outfiles[-4:], action="w")
    make_directory(args.outdir)

    # condition comparisons
    no_conditions = len(args.cond_names)  # NOTE: use unique conditions (not #signals)
    if args.time_series:
        comparisons = list(zip(args.cond_names[:-1], args.cond_names[1:]))
    else:
        comparisons = list(itertools.combinations(args.cond_names, 2))
    args.comparisons = comparisons

    # debug/fig PDFs
    if args.debug:
        debug_out = os.path.join(args.outdir, args.prefix + "_debug.pdf")
        debug_pdf = PdfPages(debug_out, keep_empty=True)

    fig_out = os.path.join(args.outdir, args.prefix + "_figures.pdf")
    figure_pdf = PdfPages(fig_out, keep_empty=True)
    cluster_out = os.path.join(args.outdir, args.prefix + "_clusters.pdf")
    cluster_pdf = PdfPages(cluster_out, keep_empty=True)

    plt.figure()
    plt.axis('off')
    plt.text(0.5, 0.8, "BINDETECT FIGURES", ha="center", va="center", fontsize=PDF_FONT_SIZE, fontweight="bold")
    titles = ["Raw score distributions"]
    if no_conditions > 1 and not args.norm_off:
        titles.append("Normalized score distributions")
    if args.debug:
        for (c1, c2) in comparisons:
            titles.append(f"Background log2FCs ({c1} / {c2})")
    for (c1, c2) in comparisons:
        titles.append(f"BINDetect volcano plot ({c1} / {c2})")
    plt.text(0.1, 0.6, "\n".join([f"Page {i+2}) {t}" for i, t in enumerate(titles)]) + "\n\n", va="top", fontsize=PDF_FONT_SIZE, fontweight="bold")
    apply_ascii_minus_to_figure(plt.gcf())
    figure_pdf.savefig(bbox_inches='tight')
    plt.close()

    plt.figure()
    plt.axis('off')
    plt.text(0.5, 0.8, "BINDETECT CLUSTERS", ha="center", va="center", fontsize=PDF_FONT_SIZE, fontweight="bold")
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
    logger.info("Reading motifs from file")
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

    logger.info("Creating folder structure for each TF")
    for TF in motif_names:
        make_directory(os.path.join(args.outdir, TF))
        make_directory(os.path.join(args.outdir, TF, "beds"))
        make_directory(os.path.join(args.outdir, TF, "plots"))

    # logos
    logo_filenames = {m.prefix: os.path.join(args.outdir, m.prefix, m.prefix + ".png") for m in motif_list}
    logger.info("Plotting sequence logos for each motif")
    task_list = [pool.apply_async(OneMotif.logo_to_file, (m, logo_filenames[m.prefix],)) for m in motif_list]
    monitor_progress(task_list, logger)
    _ = [t.get() for t in task_list]
    logger.comment("")

    logger.debug("Getting base64 strings per motif")
    for m in motif_list:
        with open(logo_filenames[m.prefix], "rb") as png:
            m.base = base64.b64encode(png.read()).decode("utf-8")

    # --------------------- scan motifs + match to signals --------------------- #
    logger.comment("")
    logger.start_logger_queue()
    args.log_q = logger.queue
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
    apply_ascii_minus_to_figure(fig)
    figure_pdf.savefig(fig, bbox_inches='tight'); plt.close()

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
        apply_ascii_minus_to_figure(fig)
        figure_pdf.savefig(fig, bbox_inches='tight'); plt.close()
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
        apply_ascii_minus_to_figure(fig)
        figure_pdf.savefig(fig, bbox_inches='tight'); plt.close()

    # ---------------------- threshold (bound/unbound) ------------------------ #
    logger.info("Estimating bound/unbound threshold")
    bg_values = np.array([background["signal"][c] for c in args.cond_names]).flatten()
    logger.debug(f"Size of background array collected: {bg_values.size}")
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
    if args.cores == 1:
        for name in motif_names:
            logger.info(f"- {name}")
            results.append(process_tfbs(name, args, log2fc_params))
    else:
        tlist = [pool.apply_async(process_tfbs, (name, args, log2fc_params)) for name in motif_names]
        monitor_progress(tlist, logger)
        results = [t.get() for t in tlist]

    logger.info("Concatenating results from subsets")
    info_table = pd.concat(results)

    pool.terminate()
    pool.join()
    logger.stop_logger_queue()

    # ---------------------- cluster TF overlaps & outputs -------------------- #
    clustering = RegionCluster(TF_overlaps)
    clustering.cluster(threshold=args.cluster_threshold)

    convert = {m.prefix: m.name for m in motif_list}
    for cluster in clustering.clusters:
        for name in convert:
            clustering.clusters[cluster]["cluster_name"] = clustering.clusters[cluster]["cluster_name"].replace(name, convert[name])

    matrix_out = os.path.join(args.outdir, args.prefix + "_distances.txt")
    clustering.write_distance_mat(matrix_out)

    logger.comment("")
    logger.info("Writing all_bindetect files")

    names, ids = [], []
    for prefix in info_table.index:
        m = [m for m in motif_list if m.prefix == prefix]
        names.append(m[0].name); ids.append(m[0].id)
    info_table.insert(0, "output_prefix", info_table.index)
    info_table.insert(1, "name", names)
    info_table.insert(2, "motif_id", ids)

    cluster_names = []
    for name in info_table.index:
        for cluster in clustering.clusters:
            if name in clustering.clusters[cluster]["member_names"]:
                cluster_names.append(clustering.clusters[cluster]["cluster_name"])
    info_table.insert(3, "cluster", cluster_names)

    info_table_clustered = info_table.groupby("cluster").mean(numeric_only=True).reset_index()

    info_table["total_tfbs"] = info_table["total_tfbs"].map(int)
    for cond in args.cond_names:
        info_table[f"{cond}_bound"] = info_table[f"{cond}_bound"].map(int)
        if f"{cond}_n_replicates" in info_table.columns:
            info_table[f"{cond}_n_replicates"] = info_table[f"{cond}_n_replicates"].map(int)

    for (c1, c2) in comparisons:
        base = f"{c1}_{c2}"
        info_table[base + "_change"] = info_table[base + "_change"].astype(float).round(5)
        info_table[base + "_pvalue"] = info_table[base + "_pvalue"].map("{:.5E}".format, na_action="ignore")

        names_series = info_table["output_prefix"]
        changes = info_table[base + "_change"].astype(float)
        pvals = info_table[base + "_pvalue"].astype(float)
        filtered_p = pvals[pvals > 0]
        pval_min = np.percentile(filtered_p, 5) if len(filtered_p) >= 1 else 1.0
        change_min, change_max = np.percentile(changes, [5, 95])

        for i, (chg, p) in enumerate(zip(changes, pvals)):
            # info_table.at[names_series[i], base + "_highlighted"] = (chg < change_min) or (chg > change_max) or (p < pval_min)
            name_key = names_series.iloc[i] if hasattr(names_series, "iloc") else names_series[i]
            info_table.at[name_key, f"{base}_highlighted"] = (chg < change_min) or (chg > change_max) or (p < pval_min)

    bindetect_out = os.path.join(args.outdir, args.prefix + "_results.txt")
    info_table.to_csv(bindetect_out, sep="\t", index=False, header=True, na_rep="NA")

    repeated_conditions = any(count > 1 for count in args.condition_replicates.values())
    write_replicate_report = args.replicate_report == "on" or (
        args.replicate_report == "auto" and (repeated_conditions or args.replicate_map is not None)
    )
    if write_replicate_report and len(args.cond_names) > 1:
        report_out = args.replicate_report_out or os.path.join(args.outdir, args.prefix + "_replicate_report.tsv")
        summary_out = args.replicate_summary_out or os.path.join(args.outdir, args.prefix + "_replicate_summary.tsv")
        figure_out = args.replicate_figure_out or os.path.join(args.outdir, args.prefix + "_replicate_report.png")
        try:
            build_replicate_report(
                bindetect_out,
                report_out,
                summary_output=summary_out,
                figure_output=figure_out,
                replicate_map=args.replicate_map,
            )
            logger.info(f"Wrote replicate-aware BINDetect report to {report_out}")
        except Exception as exc:
            logger.warning(f"Could not write replicate-aware BINDetect report: {exc}")

    if not args.skip_excel:
        bindetect_excel = os.path.join(args.outdir, args.prefix + "_results.xlsx")
        with pd.ExcelWriter(bindetect_excel, engine='xlsxwriter') as writer:
            info_table.to_excel(writer, index=False, sheet_name="Individual motifs")
            info_table_clustered.to_excel(writer, index=False, sheet_name="Motif clusters")
            for sheet in writer.sheets:
                ws = writer.sheets[sheet]
                n_rows = ws.dim_rowmax
                n_cols = ws.dim_colmax
                ws.autofilter(0, 0, n_rows, n_cols)

    # BEGIN EDIT: emit skew/shift PDF right next to *_results.txt
    skew_pdf = os.path.join(args.outdir, args.prefix + "_results_skewness_report.pdf")
    try:
        # Prefer a programmatic API if available
        if hasattr(skewrep, "generate_skew_report"):
            skewrep.generate_skew_report(
                results_tsv=bindetect_out,
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
                    results_tsv=bindetect_out,
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
                    "bindetect_skew_report has no generate_skew_report() or main_from_kwargs(); skipping PDF.")
    except Exception as e:
        logger.warning(f"Could not generate skew/shift report: {e}")
    # END EDIT

    # ------------------------------ plots ------------------------------------ #
    if no_conditions > 1:
        logger.info("Creating BINDetect plot(s)")
        change_cols = [c for c in info_table.columns if "_change" in c]
        pvalue_cols = [c for c in info_table.columns if "_pvalue" in c]
        info_table[change_cols] = info_table[change_cols].fillna(0)
        info_table[pvalue_cols] = info_table[pvalue_cols].fillna(1)

        for (c1, c2) in comparisons:
            logger.info(f"- {c1} / {c2} (static plot)")
            base = f"{c1}_{c2}"
            for m in motif_list:
                name = m.prefix
                m.change = float(info_table.at[name, base + "_change"])
                m.pvalue = float(info_table.at[name, base + "_pvalue"])
                m.logpvalue = -np.log10(m.pvalue) if m.pvalue > 0 else -np.log10(1e-308)
                m.highlighted = info_table.at[name, base + "_highlighted"]
                if m.highlighted:
                    m.group = f"{c2}_up" if m.change < 0 else f"{c1}_up"
                else:
                    m.group = "n.s."
            volcano_fig, cluster_fig = plot_bindetect(motif_list, clustering, [c1, c2], args)
            apply_ascii_minus_to_figure(volcano_fig)
            apply_ascii_minus_to_figure(cluster_fig)
            figure_pdf.savefig(volcano_fig, bbox_inches='tight'); plt.close(volcano_fig)
            cluster_pdf.savefig(cluster_fig, bbox_inches='tight'); plt.close(cluster_fig)

            logger.info(f"- {c1} / {c2} (interactive plot)")
            html_out = os.path.join(args.outdir, args.prefix + "_" + base + ".html")
            aggregate_data = None
            if getattr(args, "aggregate_signals", None) and getattr(args, "plot_aggregate", "off") != "off":
                try:
                    aggregate_data = build_bindetect_aggregate_payload(motif_list, info_table, [c1, c2], args)
                except Exception as exc:
                    logger.warning(f"Could not build aggregate payload for interactive HTML: {exc}")
            plot_interactive_bindetect(motif_list, [c1, c2], html_out, aggregate_data=aggregate_data, title="Differential footprint report")

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
    figure_pdf.close()
    cluster_pdf.close()
    logger.end()


# ----------------------------------------------------------------------------- #
def run_cli():
    parser = argparse.ArgumentParser()
    parser = add_bindetect_arguments(parser)
    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()
    run_bindetect(args)


def match_motifs_cli():
    parser = argparse.ArgumentParser(prog="match-motifs")
    parser = add_bindetect_arguments(parser)
    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()
    if not args.signals or len(args.signals) != 1:
        parser.error("match-motifs expects exactly one --signals bigWig for single-sample motif matching")
    if args.cond_names is not None and len(args.cond_names) != 1:
        parser.error("match-motifs expects exactly one --cond-names value when provided")
    if args.prefix == "bindetect":
        args.prefix = "motif_matches"
    args.replicate_report = "off"
    run_bindetect(args)


def diff_footprints_cli():
    parser = argparse.ArgumentParser(prog="diff-footprints")
    parser = add_bindetect_arguments(parser)
    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()
    if not args.signals or len(args.signals) < 2:
        parser.error("diff-footprints expects at least two --signals bigWigs")
    if args.prefix == "bindetect":
        args.prefix = "diff_footprints"
    run_bindetect(args)


if __name__ == '__main__':
    run_cli()
