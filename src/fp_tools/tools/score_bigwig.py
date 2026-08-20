#!/usr/bin/env python
"""
call-footprints command driver for scoring cutsite bigWig inputs.

This module is responsible for:
- reading corrected cutsite bigWigs
- computing footprint, sum, or mean scores across regions
- writing the scored output bigWig
"""
import os
import signal
import sys
import argparse
import copy
import bisect
import queue
import numpy as np
from fp_tools.utils import bigwig as pyBigWig
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing

# Internal functions and classes (fp_tools namespace)
from fp_tools.parsers import add_scorebigwig_arguments
from fp_tools.utils.utilities import (
    check_required, check_files, bigwig_writer, check_cores, monitor_progress,
)
from fp_tools.utils.regions import OneRegion, RegionList
from fp_tools.utils.sequences import *          # kept for parity (even if not used here directly)
from fp_tools.utils.signals import *            # fast_rolling_math, footprint_score_array, FOS_score
from fp_tools.utils.signals import local_maxima_indices
from fp_tools.utils.multiscale import (
    multiscale_depletion, parse_scales, summarize_multiscale,
    trim_multiscale_features, write_multiscale_npz,
)
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.project_layout import (
    corrected_bigwig_path,
    is_project_layout,
    normalized_bigwig_path,
    project_analysis_peaks,
    project_root,
    read_sample_table,
    samples_root,
)


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


def _stem(path):
    stem = os.path.basename(os.path.splitext(str(path))[0])
    for suffix in ("_corrected", ".corrected", "_cutsites", ".cutsites"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


def _scorebigwig_batch_items(args):
    if is_project_layout(getattr(args, "layout", None)) and getattr(args, "sample_table", None):
        if not getattr(args, "outdir", None):
            sys.exit("--layout project requires --outdir")
        project = project_root(getattr(args, "outdir", None))
        samples = read_sample_table(args.sample_table)
        args.sample_names = [row.sample for row in samples]
        args.sample_output_root = str(samples_root(project))
        args.signals = [
            str(normalized_bigwig_path(project, row.sample))
            if normalized_bigwig_path(project, row.sample).exists()
            else str(corrected_bigwig_path(project, row.sample))
            for row in samples
        ]
        args.regions = str(project_analysis_peaks(project, getattr(args, "regions", None)))

    signals = _as_list(getattr(args, "signals", None))
    if not signals and getattr(args, "signal", None):
        signals = [args.signal]

    if not signals:
        sys.exit("call-footprints requires --signal or --signals")
    sample_output_root = getattr(args, "sample_output_root", None)
    sample_names = _as_list(getattr(args, "sample_names", None))
    if sample_output_root:
        if not sample_names:
            sys.exit("--sample-output-root requires --sample-names")
        if len(sample_names) != len(signals):
            sys.exit("--sample-names must have the same number of values as --signals")
        outputs = [
            os.path.join(sample_output_root, sample, "footprints", f"{sample}_footprints.bw")
            for sample in sample_names
        ]
    elif len(signals) == 1 and getattr(args, "output", None) and not getattr(args, "outputs", None):
        outputs = [args.output]
    else:
        outputs = _as_list(getattr(args, "outputs", None))
        if not outputs:
            outdir = getattr(args, "outdir", None)
            if not outdir:
                sys.exit("call-footprints with --signals requires --outputs or --outdir")
            outputs = [os.path.join(outdir, f"{_stem(path)}_footprints.bw") for path in signals]

    if len(outputs) != len(signals):
        sys.exit("--outputs must have the same number of files as --signals")

    output_beds = _as_list(getattr(args, "output_beds", None))
    if output_beds and len(output_beds) != len(signals):
        sys.exit("--output-beds must have the same number of files as --signals")
    if not output_beds and sample_output_root and getattr(args, "call_candidates", False):
        output_beds = [
            os.path.join(sample_output_root, sample, "footprints", f"{sample}_candidate_footprints.bed")
            for sample in sample_names
        ]
    elif not output_beds and getattr(args, "output_bed", None):
        if len(signals) > 1:
            sys.exit("use --output-beds or --output-bed-dir when scoring multiple --signals")
        output_beds = [args.output_bed]
    if not output_beds and getattr(args, "output_bed_dir", None):
        output_beds = [os.path.join(args.output_bed_dir, f"{_stem(path)}_candidates.bed") for path in signals]
    if not output_beds:
        output_beds = [None] * len(signals)

    output_npzs = _as_list(getattr(args, "output_multiscale_npzs", None))
    if output_npzs and len(output_npzs) != len(signals):
        sys.exit("--output-multiscale-npzs must have the same number of files as --signals")
    if not output_npzs and getattr(args, "output_multiscale_npz", None):
        if len(signals) > 1:
            sys.exit("use --output-multiscale-npzs when scoring multiple --signals with --score multiscale")
        output_npzs = [args.output_multiscale_npz]
    if not output_npzs:
        output_npzs = [None] * len(signals)

    return list(zip(signals, outputs, output_beds, output_npzs))


def _is_batch_request(args):
    signals = _as_list(getattr(args, "signals", None))
    return bool(signals) or bool(getattr(args, "outputs", None)) or bool(getattr(args, "outdir", None))


def _sample_worker_plan(n_items, cores, requested=None):
    """Return (sample_workers, cores_per_sample) for multi-signal dispatch."""

    if n_items <= 1:
        return 1, cores
    if requested is not None:
        workers = max(1, min(int(requested), n_items))
    else:
        if cores is None:
            return 1, cores
        workers = min(n_items, max(1, int(cores) // 8))
    cores_per_sample = cores
    if cores is not None and workers > 1:
        cores_per_sample = max(1, int(cores) // workers)
    return workers, cores_per_sample


def _run_scorebigwig_batch_item(single_args):
    _run_scorebigwig_single(single_args)
    return single_args.output


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
            if getattr(args, "footprint_kernel", "fast") in {"reference", "legacy"}:
                scores = footprint_score_array(signal, args.flank_min, args.flank_max, args.fp_min, args.fp_max)
            else:
                scores = footprint_score_array_fast(signal, args.flank_min, args.flank_max, args.fp_min, args.fp_max)

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

    values = np.asarray(values, dtype=float)
    return iter(local_maxima_indices(values.astype("float64", copy=False)))


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
                insert_at = bisect.bisect_left(kept_centers, center)
                left_too_close = insert_at > 0 and center - kept_centers[insert_at - 1] < min_distance
                right_too_close = insert_at < len(kept_centers) and kept_centers[insert_at] - center < min_distance
                if left_too_close or right_too_close:
                    continue
                kept_centers.insert(insert_at, center)
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
def _run_scorebigwig_single(args):
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
    logger = FpToolsLogger("call-footprints", args.verbosity)
    logger.begin()
    parser = add_scorebigwig_arguments(argparse.ArgumentParser())
    logger.arguments_overview(parser, args)
    logger.output_files([args.output, getattr(args, "output_multiscale_npz", None), getattr(args, "output_bed", None)])

    # A one-core run stays entirely in-process.  This is materially faster for
    # frozen desktop applications on Windows because spawning each worker
    # would otherwise unpack the full one-file executable again.
    args.cores = check_cores(args.cores, logger)
    logger.debug("Setting up listener for log")
    if args.cores > 1:
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

    logger.debug(f"Worker cores: {args.cores}")
    logger.debug("Writer cores: 1")

    if args.cores == 1:
        writer_queue = queue.Queue()
        writer_args = copy.copy(args)
        worker_args = copy.copy(args)
        worker_args.writer_qs = {"scores": writer_queue}
        results = []
        with ThreadPoolExecutor(max_workers=1) as writer_executor:
            writer_future = writer_executor.submit(
                bigwig_writer,
                writer_queue,
                {"scores": args.output},
                header,
                copy.deepcopy(regions),
                writer_args,
            )
            try:
                for chunk in regions_chunks:
                    results.append(calculate_scores(chunk, copy.copy(worker_args)))
            finally:
                writer_queue.put((None, None, None))
            writer_future.result()

        if getattr(args, "output_multiscale_npz", None):
            records = [record for chunk_records in results for record in chunk_records]
            records.sort(key=lambda item: (reference_chroms.index(item[0][0]), item[0][1]))
            write_multiscale_npz(args.output_multiscale_npz, records, args.scales, args.multiscale_summary)

        _validate_bigwig_output(args.output, logger)
        if getattr(args, "output_bed", None):
            _write_candidate_bed(args.output, regions, args.output_bed, args, chrom_info, logger)
        logger.end()
        return

    manager = mp.Manager()
    pool = None
    writer_pool = None

    # Start bigwig writer
    q = manager.Queue()
    writer_pool = mp.Pool(processes=1)
    writer_args = copy.copy(args)
    writer_result = writer_pool.apply_async(bigwig_writer, args=(q, {"scores": args.output}, header, regions, writer_args))
    writer_pool.close()  # no more jobs to writer_pool
    writer_qs = {"scores": q}

    try:
        # Start workers. Each task receives its own shallow argument snapshot;
        # reusing one mutable Namespace here can race with multiprocessing
        # pickling in nested project-mode runs.
        pool = mp.Pool(processes=args.cores)
        worker_args = copy.copy(args)
        worker_args.writer_qs = dict(writer_qs)
        task_list = [
            pool.apply_async(calculate_scores, args=[chunk, copy.copy(worker_args)])
            for chunk in regions_chunks
        ]
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


def run_scorebigwig(args):
    if _is_batch_request(args):
        items = _scorebigwig_batch_items(args)
        single_args_list = []
        for signal_path, output_path, output_bed, output_npz in items:
            single_args = copy.copy(args)
            single_args.signals = None
            single_args.outputs = None
            single_args.output_beds = None
            single_args.output_bed_dir = None
            single_args.output_multiscale_npzs = None
            single_args.outdir = None
            single_args.signal = signal_path
            single_args.output = output_path
            single_args.output_bed = output_bed
            single_args.output_multiscale_npz = output_npz
            single_args_list.append(single_args)
        sample_workers, sample_cores = _sample_worker_plan(
            len(single_args_list),
            getattr(args, "cores", None),
            getattr(args, "sample_workers", None),
        )
        for single_args in single_args_list:
            single_args.cores = sample_cores
        if sample_workers == 1:
            for single_args in single_args_list:
                _run_scorebigwig_batch_item(single_args)
        else:
            with ProcessPoolExecutor(max_workers=sample_workers) as executor:
                futures = [executor.submit(_run_scorebigwig_batch_item, single_args) for single_args in single_args_list]
                for future in as_completed(futures):
                    future.result()
        return

    _run_scorebigwig_single(args)


# ----------------------------------------------------------------------------- #
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser = add_scorebigwig_arguments(parser)
    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit()
    run_scorebigwig(args)
