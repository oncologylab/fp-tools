#!/usr/bin/env python
"""
FootprintScores command driver for scoring cutsite bigWig inputs.

This module is responsible for:
- reading corrected cutsite bigWigs
- computing footprint, sum, or mean scores across regions
- writing the scored output bigWig
"""
import os
import signal
import sys
import argparse
import numpy as np
import pyBigWig
import multiprocessing as mp
from contextlib import closing

# Internal functions and classes (fp_tools namespace)
from fp_tools.parsers import add_scorebigwig_arguments
from fp_tools.utils.utilities import (
    check_required, check_files, bigwig_writer, check_cores, monitor_progress
)
from fp_tools.utils.regions import OneRegion, RegionList
from fp_tools.utils.sequences import *          # kept for parity (even if not used here directly)
from fp_tools.utils.signals import *            # fast_rolling_math, footprint_score_array, FOS_score
from fp_tools.utils.multiscale import (
    multiscale_depletion, parse_scales, summarize_multiscale,
    trim_multiscale_features, write_multiscale_npz,
)
from fp_tools.utils.logger import FpToolsLogger


def _normalize_paths(args):
    """
    Make file paths absolute and ensure the output directory exists
    before any subprocess touches them.
    """
    # input files
    if getattr(args, "signal", None):
        args.signal = os.path.abspath(args.signal)
    if getattr(args, "regions", None):
        args.regions = os.path.abspath(args.regions)

    # output files
    if getattr(args, "output", None):
        out_dir = os.path.dirname(args.output) or "."
        os.makedirs(out_dir, exist_ok=True)   # <-- key bit: create parent dir
        args.output = os.path.abspath(args.output)

    if getattr(args, "output_multiscale_npz", None):
        out_dir = os.path.dirname(args.output_multiscale_npz) or "."
        os.makedirs(out_dir, exist_ok=True)
        args.output_multiscale_npz = os.path.abspath(args.output_multiscale_npz)

    if getattr(args, "output_bed", None):
        out_dir = os.path.dirname(args.output_bed) or "."
        os.makedirs(out_dir, exist_ok=True)
        args.output_bed = os.path.abspath(args.output_bed)

    return args


# ----------------------------------------------------------------------------- #
def calculate_scores(regions, args):

    logger = FpToolsLogger("", args.verbosity, args.log_q)

    pybw_signal = pyBigWig.open(args.signal)        # cutsites signal
    pybw_header = pybw_signal.chroms()
    chrom_lengths = {chrom: int(pybw_header[chrom]) for chrom in pybw_header}

    # Set flank to enable scoring in ends of regions
    flank = args.region_flank

    multiscale_records = []

    # Go through each region
    for i, region in enumerate(regions):

        logger.debug(f"Calculating scores for region: {region}")

        # Extend region with necessary flank
        region.extend_reg(flank)
        reg_key = (region.chrom, region.start + flank, region.end - flank)   # output region

        # Get bigwig signal in region
        signal = region.get_signal(pybw_signal, logger=logger)
        signal = np.nan_to_num(signal).astype("float64")

        # -------- Prepare signal for score calculation ------- #
        if getattr(args, "absolute", False):
            signal = np.abs(signal)

        if args.min_limit is not None:
            signal[signal < args.min_limit] = args.min_limit
        if args.max_limit is not None:
            signal[signal > args.max_limit] = args.max_limit

        # ------------------ Calculate scores ----------------- #
        if args.score == "sum":
            scores = fast_rolling_math(signal, args.window, "sum")

        elif args.score == "mean":
            scores = fast_rolling_math(signal, args.window, "mean")

        elif args.score == "footprint":
            scores = footprint_score_array(signal, args.flank_min, args.flank_max, args.fp_min, args.fp_max)

        elif args.score == "multiscale":
            features = multiscale_depletion(signal, args.scales)
            scores = summarize_multiscale(features, args.multiscale_summary)

        elif args.score == "FOS":
            scores = FOS_score(signal, args.flank_min, args.flank_max, args.fp_min, args.fp_max)

        elif args.score == "none":
            scores = signal

        else:
            sys.exit(f"Scoring {args.score} not found")

        # ----------------- Post-process scores --------------- #

        # Smooth signal with args.smooth bp
        if args.smooth and args.smooth > 1:
            scores = fast_rolling_math(scores, args.smooth, "mean")

        # Remove ends to prevent overlap with other regions
        if flank > 0:
            scores = scores[flank:-flank]

        if args.score == "multiscale" and getattr(args, "output_multiscale_npz", None):
            multiscale_records.append((reg_key, trim_multiscale_features(features, flank)))

        args.writer_qs["scores"].put(("scores", reg_key, scores))

    if getattr(args, "output_multiscale_npz", None):
        return multiscale_records
    return 1




def _local_maxima(values):
    """Yield 0-based offsets that are local maxima in a one-dimensional score array."""

    if len(values) == 0:
        return
    if len(values) == 1:
        if np.isfinite(values[0]):
            yield 0
        return
    for idx, value in enumerate(values):
        if not np.isfinite(value):
            continue
        left = values[idx - 1] if idx > 0 else -np.inf
        right = values[idx + 1] if idx < len(values) - 1 else -np.inf
        if value >= left and value >= right and (value > left or value > right):
            yield idx


def _write_candidate_bed(score_bigwig, regions, output_bed, args, chrom_info, logger):
    """Call ranked local footprint candidates from the scored bigWig."""

    if args.score not in {"footprint", "FOS", "multiscale"}:
        logger.warning("--output-bed is intended for footprint-like scores; writing calls from the selected score anyway.")

    min_score = args.min_score
    call_width = max(1, int(args.call_width))
    half_width = max(1, call_width // 2)
    min_distance = max(0, int(args.min_distance))
    raw_calls = []

    with closing(pyBigWig.open(score_bigwig)) as bw:
        for region_idx, region in enumerate(regions, start=1):
            values = np.asarray(bw.values(region.chrom, region.start, region.end, numpy=True), dtype=float)
            values = np.nan_to_num(values, nan=-np.inf, posinf=-np.inf, neginf=-np.inf)
            candidates = []
            for offset in _local_maxima(values):
                score = float(values[offset])
                if min_score is not None and score < min_score:
                    continue
                center = int(region.start + offset)
                candidates.append((score, center, region_idx, region.chrom, region.start, region.end))
            candidates.sort(key=lambda item: (-item[0], item[1]))
            kept_centers = []
            for score, center, reg_i, chrom, reg_start, reg_end in candidates:
                if any(abs(center - other) < min_distance for other in kept_centers):
                    continue
                kept_centers.append(center)
                start = max(0, center - half_width)
                end = min(int(chrom_info[chrom]), start + call_width)
                start = max(0, end - call_width)
                raw_calls.append({
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "score": score,
                    "center": center,
                    "source_region": f"{chrom}:{reg_start}-{reg_end}",
                })

    raw_calls.sort(key=lambda row: (-row["score"], row["chrom"], row["start"]))
    if args.top_n is not None:
        raw_calls = raw_calls[:max(0, int(args.top_n))]

    with open(output_bed, "w", encoding="utf-8") as handle:
        handle.write("#chrom\tstart\tend\tname\tscore\tstrand\tsource_region\tcenter\traw_score\n")
        for idx, row in enumerate(raw_calls, start=1):
            name = f"footprint_{idx}"
            score = f"{row['score']:.6g}"
            handle.write(
                f"{row['chrom']}\t{row['start']}\t{row['end']}\t{name}\t{score}\t.\t"
                f"{row['source_region']}\t{row['center']}\t{score}\n"
            )
    logger.info(f"Wrote {len(raw_calls)} footprint candidate calls to {output_bed}")

def _validate_bigwig_output(path, logger):
    """Fail loudly if the writer did not produce a readable bigWig."""

    if not path or not os.path.exists(path):
        raise RuntimeError(f"Expected score bigWig was not created: {path}")
    if os.path.getsize(path) == 0:
        raise RuntimeError(f"Expected score bigWig is empty: {path}")

    try:
        with closing(pyBigWig.open(path)) as bw:
            if not bw or not bw.chroms():
                raise RuntimeError(f"Expected score bigWig has no chromosome header: {path}")
    except Exception as exc:
        raise RuntimeError(f"Expected score bigWig is not readable: {path}") from exc

# ----------------------------------------------------------------------------- #
def run_scorebigwig(args):
    # Ensure paths are sane and the output directory exists
    args = _normalize_paths(args)

    check_required(args, ["signal", "output", "regions"])
    check_files([args.signal, args.regions], "r")
    check_files([args.output, getattr(args, "output_multiscale_npz", None), getattr(args, "output_bed", None)], "w")

    if getattr(args, "output_multiscale_npz", None) and args.score != "multiscale":
        sys.exit("--output-multiscale-npz requires --score multiscale")

    # ------------------------------------------------------------------------- #
    # Logger
    # ------------------------------------------------------------------------- #
    logger = FpToolsLogger("FootprintScores", args.verbosity)
    logger.begin()
    parser = add_scorebigwig_arguments(argparse.ArgumentParser())
    logger.arguments_overview(parser, args)
    logger.output_files([args.output, getattr(args, "output_multiscale_npz", None), getattr(args, "output_bed", None)])

    logger.debug("Setting up listener for log")
    logger.start_logger_queue()
    args.log_q = logger.queue

    # ------------------------------------------------------------------------- #
    # I/O
    # ------------------------------------------------------------------------- #
    logger.info("Processing input files")

    logger.info("- Opening input cutsite bigwig")
    pybw_signal = pyBigWig.open(args.signal)
    pybw_header = pybw_signal.chroms()
    chrom_info = {chrom: int(pybw_header[chrom]) for chrom in pybw_header}
    logger.debug(f"Chromosome lengths from input bigwig: {chrom_info}")

    # Decide regions
    logger.info("- Getting output regions ready")
    if args.regions:
        regions = RegionList().from_bed(args.regions)

        # Exclude regions not present in bigwig
        not_in_bigwig = list(set(regions.get_chroms()) - set(chrom_info.keys()))
        if len(not_in_bigwig) > 0:
            logger.warning(
                f"Contigs {not_in_bigwig} were found in input --regions, but were not found in input --signal. "
                f"These regions cannot be scored and will therefore be excluded from output."
            )
            regions = regions.remove_chroms(not_in_bigwig)

        regions.apply_method(OneRegion.extend_reg, args.extend)
        regions.merge()
        regions.apply_method(OneRegion.check_boundary, chrom_info, "cut")

    else:
        regions = RegionList().from_chrom_lengths(chrom_info)

    # Set flank to enable scoring in ends of regions
    if args.score == "sum":
        args.region_flank = int(args.window / 2.0)
    elif args.score in ("footprint", "FOS"):
        args.region_flank = int(args.flank_max)
    elif args.score == "multiscale":
        args.scales = list(parse_scales(args.scales))
        args.region_flank = int(max(args.scales) * 2)
    else:
        args.region_flank = 0

    # Double-check boundaries with flank
    for i, region in enumerate(regions):
        region.extend_reg(args.region_flank)
        region = region.check_boundary(chrom_info, "cut")
        region.extend_reg(-args.region_flank)

    # Output bigwig header
    reference_chroms = sorted(list(chrom_info.keys()))
    header = [(chrom, chrom_info[chrom]) for chrom in reference_chroms]
    regions.loc_sort(reference_chroms)

    # ------------------------------------------------------------------------- #
    # Calculate & write
    # ------------------------------------------------------------------------- #
    logger.info("Calculating footprints in regions...")
    regions_chunks = regions.chunks(args.split)

    args.cores = check_cores(args.cores, logger)
    logger.debug(f"Worker cores: {args.cores}")
    logger.debug("Writer cores: 1")

    manager = mp.Manager()
    pool = None
    writer_pool = None

    # Start bigwig writer
    q = manager.Queue()
    writer_pool = mp.Pool(processes=1)
    writer_result = writer_pool.apply_async(bigwig_writer, args=(q, {"scores": args.output}, header, regions, args))
    writer_pool.close()  # no more jobs to writer_pool
    writer_qs = {"scores": q}
    args.writer_qs = writer_qs

    try:
        # Start workers
        pool = mp.Pool(processes=args.cores)
        task_list = [pool.apply_async(calculate_scores, args=[chunk, args]) for chunk in regions_chunks]
        pool.close()
        monitor_progress(task_list, logger)
        results = [task.get() for task in task_list]
        pool.join()
        if getattr(args, "output_multiscale_npz", None):
            records = [record for chunk_records in results for record in chunk_records]
            records.sort(key=lambda item: (reference_chroms.index(item[0][0]), item[0][1]))
            write_multiscale_npz(args.output_multiscale_npz, records, args.scales, args.multiscale_summary)
            logger.info(f"Wrote multiscale tensor sidecar: {args.output_multiscale_npz}")
    except KeyboardInterrupt:
        logger.warning("Interrupted by user; shutting down workers/writer...")
        if pool is not None:
            pool.terminate()
            pool.join()
        raise
    finally:
        # Tell writer to finish and clean up either way
        for q in writer_qs.values():
            q.put((None, None, None))
        if writer_pool is not None:
            writer_pool.join()

    writer_result.get()
    _validate_bigwig_output(args.output, logger)
    if getattr(args, "output_bed", None):
        _write_candidate_bed(args.output, regions, args.output_bed, args, chrom_info, logger)

    logger.stop_logger_queue()
    logger.end()


# ----------------------------------------------------------------------------- #
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = add_scorebigwig_arguments(parser)
    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()
    run_scorebigwig(args)
