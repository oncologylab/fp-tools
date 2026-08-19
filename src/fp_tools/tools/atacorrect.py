#!/usr/bin/env python

"""
atac-correct command driver for bias estimation and cutsite correction.

This module keeps the main command-line execution path for:
- reading BAM/FASTA/peak inputs
- estimating Tn5 sequence bias
- writing corrected and auxiliary bigWig tracks
- generating the atac-correct summary PDF
"""

#--------------------------------------------------------------------------------------------------------#
#----------------------------------------- Import libraries ---------------------------------------------# 
#--------------------------------------------------------------------------------------------------------#

import os
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 

import sys
import argparse
import numpy as np
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
import re

from collections import OrderedDict
import itertools
import matplotlib
matplotlib.use("Agg")  #non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#Bio-specific packages
from fp_tools.utils import bigwig as pyBigWig
from fp_tools.utils.alignment import index_alignment, open_alignment
from fp_tools.utils.fasta import open_fasta

#Internal functions and classes
from fp_tools.parsers import add_atacorrect_arguments
from fp_tools.tools.atacorrect_functions import *
from fp_tools.utils.utilities import *
from fp_tools.utils.regions import OneRegion, RegionList
from fp_tools.utils.ngs import OneRead, ReadList
from fp_tools.utils.sequences import *
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.project_layout import (
	analysis_peaks_path,
	is_project_layout,
	merged_peaks_path,
	project_root,
	read_sample_table,
	samples_root,
	write_analysis_peaks,
)
from fp_tools.tools.normalize_bigwig import corrected_scaled_output_path, normalize_bigwigs

#np.seterr(divide='raise', invalid='raise')


def _maybe_scale_corrected_bigwigs(args, corrected_bigwigs, logger):
	"""Scale corrected bigWigs with q95 background matching when requested."""
	mode = getattr(args, "scale_corrected", "auto")
	if mode == "none":
		return []
	cohort = list(getattr(args, "scale_corrected_bigwigs", None) or [])
	if not cohort:
		cohort = list(corrected_bigwigs)
	else:
		seen = {os.path.abspath(path) for path in cohort}
		for path in corrected_bigwigs:
			abs_path = os.path.abspath(path)
			if abs_path not in seen:
				cohort.append(path)
				seen.add(abs_path)
	if mode == "auto" and len(cohort) < 2:
		logger.info("Skipping q95 corrected-track scaling in auto mode because fewer than two corrected bigWigs were provided")
		return []
	if not getattr(args, "scale_background", None):
		message = "--scale-background is required when q95 corrected-track scaling is requested"
		if mode == "auto":
			logger.info(message + "; skipping scaling")
			return []
		raise ValueError(message)
	logger.info("Q95-scaling corrected bigWigs")
	output_paths = [corrected_scaled_output_path(path) for path in cohort]
	norm_workers = getattr(args, "sample_workers", None)
	if norm_workers is None and getattr(args, "cores", None):
		norm_workers = min(len(cohort), max(1, int(args.cores) // 4))
	norm_workers = max(1, min(int(norm_workers or 1), len(cohort)))
	rows = normalize_bigwigs(
		cohort,
		args.scale_background,
		args.outdir,
		method="background-scale",
		stat="q95",
		target=getattr(args, "scale_target", "median"),
		chrom_sizes=getattr(args, "scale_chrom_sizes", None),
		output_paths=output_paths,
		workers=norm_workers,
	)
	for row in rows:
		logger.stats("Q95_SCALE\t{0}\t{1:.5f}\t{2}".format(row.sample, row.scale_factor, row.output_bigwig))
	return rows


def _compact_list(values):
	if values is None:
		return []
	if isinstance(values, (str, os.PathLike)):
		values = [values]
	return [str(value) for value in values if str(value).strip()]


def _sample_names_from_bams(bams, requested=None):
	names = _compact_list(requested)
	if names and len(names) != len(bams):
		sys.exit("Error: --sample-names must contain one value per --bams input")
	if not names:
		names = [os.path.splitext(os.path.basename(bam))[0] for bam in bams]
	clean = []
	for name in names:
		if not name or os.path.basename(name) != name or name in (".", ".."):
			sys.exit("Error: --sample-names values must be non-empty names, not paths: {0}".format(name))
		clean.append(name)
	if len(set(clean)) != len(clean):
		sys.exit("Error: --sample-names values must be unique")
	return clean


def _input_name_key(path_or_name):
	name = os.path.basename(str(path_or_name)).lower()
	for suffix in (".narrowpeak", ".broadpeak", ".bed", ".bam", ".sam", ".cram"):
		if name.endswith(suffix):
			name = name[:-len(suffix)]
	name = re.sub(r"([._-]?(merged|all|peaks?|regions?|replicate|rep))+$", "", name)
	name = re.sub(r"[^a-z0-9]+", "", name)
	return name


def _warn_peak_sample_mismatches(bams, sample_names, peak_files, logger):
	if len(bams) <= 1 or len(peak_files) <= 1:
		return
	if len(peak_files) != len(bams):
		logger.warning(
			"Multiple BAMs and multiple peak BEDs were supplied, but their counts differ "
			"({0} BAMs vs {1} peak BEDs). fp-tools will merge all peak BEDs into merged_peaks.bed; "
			"please confirm these peaks are intended for the same analysis.".format(len(bams), len(peak_files))
		)
		return
	mismatches = []
	for sample_name, peak_file in zip(sample_names, peak_files):
		sample_key = _input_name_key(sample_name)
		peak_key = _input_name_key(peak_file)
		if sample_key and peak_key and sample_key not in peak_key and peak_key not in sample_key:
			mismatches.append("{0} vs {1}".format(sample_name, os.path.basename(peak_file)))
	if mismatches:
		logger.warning(
			"Sample names and peak BED filenames do not appear to match by position: {0}. "
			"fp-tools will still merge all peak BEDs into merged_peaks.bed.".format("; ".join(mismatches))
		)


def _merge_peak_files(peak_files, outdir, merged_peaks_out=None):
	peak_files = _compact_list(peak_files)
	if not peak_files:
		sys.exit("Error: No .peaks-file given")
	if len(peak_files) == 1:
		return peak_files[0]
	if merged_peaks_out:
		merged_path = os.path.abspath(merged_peaks_out)
		make_directory(os.path.dirname(merged_path) or ".")
	else:
		make_directory(outdir)
		merged_path = os.path.join(outdir, "merged_peaks.bed")
	merged = RegionList()
	for peak_file in peak_files:
		merged.extend(RegionList().from_bed(peak_file))
	merged.merge()
	merged.write_bed(merged_path)
	return merged_path


def _selected_output_tracks(args):
	tracks = ["uncorrected", "bias", "expected", "corrected"]
	requested = getattr(args, "write_tracks", None)
	if requested:
		requested = list(requested)
		if "all" in requested:
			selected = list(tracks)
		else:
			selected = [track for track in tracks if track in requested]
	else:
		selected = list(tracks)
	selected = [track for track in selected if track not in getattr(args, "track_off", [])]
	return selected


def _corrected_output_paths_for_args(args):
	tracks = _selected_output_tracks(args)
	if "corrected" not in tracks:
		return []
	strands = ["forward", "reverse"] if args.split_strands else ["both"]
	paths = []
	for strand in strands:
		elements = [args.prefix, "corrected"] if strand == "both" else [args.prefix, "corrected", strand]
		paths.append(os.path.join(args.outdir, "{0}.bw".format("_".join(elements))))
	return paths


def _sample_worker_plan(n_items, cores, requested=None):
	"""Return (sample_workers, cores_per_sample) for multi-sample dispatch."""
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


def _run_atacorrect_sample(sample_args):
	_run_atacorrect_single(sample_args)
	return sample_args.prefix

#--------------------------------------------------------------------------------------------------------#
#-------------------------------------- Main pipeline function ------------------------------------------#
#--------------------------------------------------------------------------------------------------------#

def run_atacorrect(args):

	"""
	Batch-aware atac-correct entry point.

	``--bams`` accepts one or more BAM files. Single-BAM runs keep the legacy
	output layout; multi-BAM runs write one subdirectory per sample.
	"""

	if is_project_layout(getattr(args, "layout", None)) and getattr(args, "sample_table", None):
		if not getattr(args, "outdir", None):
			sys.exit("Error: --layout project requires --outdir")
		project = project_root(getattr(args, "outdir", None))
		samples = read_sample_table(args.sample_table)
		args.bams = [row.bam for row in samples]
		args.peaks = [row.peaks for row in samples]
		args.sample_names = [row.sample for row in samples]
		args.sample_output_root = str(samples_root(project))
		args.merged_peaks_out = str(merged_peaks_path(project))
		args.outdir = str(project)

	bams = _compact_list(getattr(args, "bams", None))
	fragments = _compact_list(getattr(args, "fragments", None))
	if bams and fragments:
		sys.exit("Error: use either --bams or --fragments, not both")
	inputs = fragments or bams
	if not inputs:
		sys.exit("Error: provide --bams <reads.bam> or --fragments <fragments.tsv.gz>")
	if args.genome == None:
		sys.exit("Error: No .fasta-file given")
	if not _compact_list(getattr(args, "peaks", None)):
		sys.exit("Error: No .peaks-file given")
	if len(inputs) > 1 and getattr(args, "prefix", None):
		sys.exit("Error: --prefix can only be used with a single input. Use --sample-names for multi-sample runs.")

	base_outdir = os.path.abspath(args.outdir) if args.outdir != None else os.path.abspath(os.getcwd())
	sample_output_root = getattr(args, "sample_output_root", None)
	if sample_output_root:
		sample_output_root = os.path.abspath(sample_output_root)
		if not getattr(args, "outdir", None):
			base_outdir = sample_output_root
	sample_names = _sample_names_from_bams(inputs, getattr(args, "sample_names", None))
	peak_files = _compact_list(args.peaks)
	preflight_logger = FpToolsLogger("atac-correct", getattr(args, "verbosity", 3))
	_warn_peak_sample_mismatches(inputs, sample_names, peak_files, preflight_logger)
	peaks_for_run = _merge_peak_files(args.peaks, base_outdir, getattr(args, "merged_peaks_out", None))
	if is_project_layout(getattr(args, "layout", None)) and getattr(args, "sample_table", None):
		project = project_root(getattr(args, "outdir", None))
		write_analysis_peaks(peaks_for_run, analysis_peaks_path(project), tuple(getattr(args, "drop_chroms", []) or []))
	corrected_bigwigs = []

	sample_args_list = []
	for bam, sample_name in zip(inputs, sample_names):
		sample_args = deepcopy(args)
		sample_args.bam = bam
		sample_args.bams = [] if fragments else [bam]
		sample_args.peaks = peaks_for_run
		sample_args.input_type = "fragments" if fragments else "bam"
		if fragments:
			sample_args.read_shift = [0, 0]
		sample_args.prefix = args.prefix if len(inputs) == 1 and args.prefix else sample_name
		if sample_output_root:
			sample_args.outdir = os.path.join(sample_output_root, sample_name, "atac_correct")
		else:
			sample_args.outdir = base_outdir if len(inputs) == 1 else os.path.join(base_outdir, sample_name)
		sample_args.scale_corrected = "none"
		sample_args._scale_after_single = False
		corrected_bigwigs.extend(_corrected_output_paths_for_args(sample_args))
		sample_args_list.append(sample_args)

	sample_workers, sample_cores = _sample_worker_plan(
		len(sample_args_list),
		getattr(args, "cores", None),
		getattr(args, "sample_workers", None),
	)
	for sample_args in sample_args_list:
		sample_args.cores = sample_cores

	if sample_workers == 1:
		for sample_args in sample_args_list:
			_run_atacorrect_sample(sample_args)
	else:
		preflight_logger.info(
			"Processing {0} BAMs with {1} concurrent sample workers and {2} cores per sample".format(
				len(sample_args_list), sample_workers, sample_cores
			)
		)
		with ProcessPoolExecutor(max_workers=sample_workers) as executor:
			futures = [executor.submit(_run_atacorrect_sample, sample_args) for sample_args in sample_args_list]
			for future in as_completed(futures):
				future.result()

	args.outdir = base_outdir
	if len(inputs) == 1:
		args.prefix = args.prefix if args.prefix else sample_names[0]
	mode = getattr(args, "scale_corrected", "auto")
	if mode != "none":
		logger = FpToolsLogger("atac-correct", args.verbosity)
		logger.begin()
		try:
			_maybe_scale_corrected_bigwigs(args, corrected_bigwigs, logger)
		except Exception as exc:
			logger.error("Corrected-track q95 scaling failed: {0}".format(exc))
			raise
		finally:
			logger.end()


def _run_atacorrect_single(args):

	"""
	Function for bias correction of input .bam files
	Calls functions in atac-correct_functions and several internal classes
	"""

	#Test if required arguments were given:
	if args.bam == None:
		sys.exit("Error: No .bam-file given")
	if args.genome == None:
		sys.exit("Error: No .fasta-file given")
	if args.peaks == None:
		sys.exit("Error: No .peaks-file given")

	#Adjust some parameters depending on input
	args.prefix = os.path.splitext(os.path.basename(args.bam))[0] if args.prefix == None else args.prefix
	args.outdir = os.path.abspath(args.outdir) if args.outdir != None else os.path.abspath(os.getcwd())

	#Set output bigwigs based on input
	tracks = _selected_output_tracks(args)

	if args.split_strands == True:
		strands = ["forward", "reverse"]
	else:
		strands = ["both"]

	output_bws = {}
	for track in tracks:
		output_bws[track] = {}
		for strand in strands:
			elements = [args.prefix, track] if strand == "both" else [args.prefix, track, strand]
			output_bws[track][strand] = {"fn": os.path.join(args.outdir, "{0}.bw".format("_".join(elements)))}

	#Set all output files
	bam_out = os.path.join(args.outdir, args.prefix + "_atacorrect.bam")
	bigwigs = [output_bws[track][strand]["fn"] for (track, strand) in itertools.product(tracks, strands)]
	figures_f = os.path.join(args.outdir, "{0}_atacorrect.pdf".format(args.prefix))

	output_files = list(bigwigs)
	if not getattr(args, "skip_qc", False):
		output_files.append(figures_f)
	output_files = list(OrderedDict.fromkeys(output_files)) 	#remove duplicates due to "both" option

	strands = ["forward", "reverse"]

	#----------------------------------------------------------------------------------------------------#
	# Print info on run
	#----------------------------------------------------------------------------------------------------#

	logger = FpToolsLogger("atac-correct", args.verbosity)
	logger.begin()

	parser = add_atacorrect_arguments(argparse.ArgumentParser())
	logger.arguments_overview(parser, args)
	logger.output_files(output_files)

	args.cores = check_cores(args.cores, logger)

	#----------------------------------------------------------------------------------------------------#
	# Test input file availability for reading 
	#----------------------------------------------------------------------------------------------------#

	logger.info("----- Processing input data -----")

	logger.debug("Testing input file availability")
	check_files([args.bam, args.genome, args.peaks], "r")

	logger.debug("Testing output directory/file writeability")
	make_directory(args.outdir)
	check_files(output_files, "w")

	#Open pdf for figures
	figure_pdf = None if getattr(args, "skip_qc", False) else PdfPages(figures_f, keep_empty=False)

	#----------------------------------------------------------------------------------------------------#
	# Read information in bam/fasta
	#----------------------------------------------------------------------------------------------------#

	logger.info("Reading info from .bam file")
	bamfile = open_alignment(args.bam, "rb")
	if bamfile.has_index() == False:
		if index_alignment(args.bam):
			logger.warning("No index found for bamfile; created one for faster access.")
		else:
			logger.warning("No BAM index found; using the portable sequential-scan cache.")

	bam_references = bamfile.references 	#chromosomes in correct order
	bam_chrom_info = dict(zip(bamfile.references, bamfile.lengths))
	logger.debug("bam_chrom_info: {0}".format(bam_chrom_info))
	bamfile.close()

	logger.info("Reading info from .fasta file")
	fastafile = open_fasta(args.genome)
	fasta_chrom_info = dict(zip(fastafile.references, fastafile.lengths))
	logger.debug("fasta_chrom_info: {0}".format(fasta_chrom_info))
	fastafile.close()

	# Compare chrom lengths for BAM input. Fragment files do not carry chromosome
	# sizes, so the FASTA header is authoritative for their synthetic read view.
	if getattr(args, "input_type", "bam") == "fragments":
		bam_references = list(fasta_chrom_info)
		bam_chrom_info = dict(fasta_chrom_info)

	#Compare chrom lengths
	chrom_in_common = set(bam_chrom_info.keys()).intersection(fasta_chrom_info.keys())
	for chrom in chrom_in_common:
		bamlen = bam_chrom_info[chrom]
		fastalen = fasta_chrom_info[chrom]
		if bamlen != fastalen:
			logger.warning("(Fastafile)\t{0} has length {1}".format(chrom, fasta_chrom_info[chrom]))
			logger.warning("(Bamfile)\t{0} has length {1}".format(chrom, bam_chrom_info[chrom]))
			sys.exit("Error: .bam and .fasta have different chromosome lengths. Please make sure the genome file is similar to the one used in mapping.")

	#Subset bam_references to those for which there are sequences in fasta
	chrom_not_in_fasta = set(bam_references) - set(fasta_chrom_info.keys())
	if len(chrom_not_in_fasta) > 1:
		logger.warning("The following contigs in the input BAM did not have sequences in --fasta: {0}. NOTE: These contigs will be skipped in calculation and output.".format(chrom_not_in_fasta))

	bam_references = [ref for ref in bam_references if ref in fasta_chrom_info]
	chrom_in_common = [ref for ref in chrom_in_common if ref in bam_references]

	#Drop mitochrondrial (or other chroms) from list
	chrom_in_common_orig = chrom_in_common
	chrom_in_common = [chrom for chrom in chrom_in_common if chrom not in args.drop_chroms]
	bam_references = [chrom for chrom in bam_references if chrom in chrom_in_common]

	dropped = set(chrom_in_common_orig) - set(chrom_in_common)
	if len(dropped) > 0:
		logger.info("The following contigs were dropped from analysis because they were found in '--drop-chroms': {0}".format(list(dropped)))
	else:
		logger.warning("No additional chromosomes were removed. Consider using '--drop-chroms' to remove mitochondrial and/or other unwanted contigs.")

	#Check if any contigs were left; else exit
	if len(chrom_in_common) == 0:
		logger.error("No common contigs left to run atac-correct on. Please check that the input BAM and '--fasta' are matching.")
		sys.exit()

	#----------------------------------------------------------------------------------------------------#
	# Read regions from bedfiles
	#----------------------------------------------------------------------------------------------------#

	logger.info("Processing input/output regions")

	#Chromosomes included in analysis
	genome_regions = RegionList().from_list([OneRegion([chrom, 0, bam_chrom_info[chrom]]) for chrom in chrom_in_common]) #full genome length
	logger.debug("CHROMS\t{0}".format("; ".join(["{0} ({1})".format(reg.chrom, reg.end) for reg in genome_regions])))
	genome_bp = sum([region.get_length() for region in genome_regions])

	# Process peaks
	peak_regions = RegionList().from_bed(args.peaks)
	peak_regions.merge()
	for i in range(len(peak_regions)-1, -1, -1):
		region = peak_regions[i]

		peak_regions[i] = region.check_boundary(bam_chrom_info, "cut")	#regions are cut/removed from list
		if peak_regions[i] is None:
			logger.warning("Peak region {0} was removed at it is either out of bounds or not in the chromosomes given in genome/bam.".format(region.tup(), i+1))
			del peak_regions[i]

	nonpeak_regions = deepcopy(genome_regions).subtract(peak_regions)

	# Process specific input regions if given
	if args.regions_in != None:
		input_regions = RegionList().from_bed(args.regions_in)
		input_regions.merge()
		input_regions.apply_method(OneRegion.check_boundary, bam_chrom_info, "cut")
	else:
		input_regions = nonpeak_regions

	# Process specific output regions
	if args.regions_out != None:
		output_regions = RegionList().from_bed(args.regions_out)
	else:
		output_regions = deepcopy(peak_regions)

	#Extend regions to make sure extend + flanking for window/flank are within boundaries
	flank_extend = args.k_flank + int(args.window/2.0)
	output_regions.apply_method(OneRegion.extend_reg, args.extend + flank_extend)
	output_regions.merge()	
	output_regions.apply_method(OneRegion.check_boundary, bam_chrom_info, "cut")
	output_regions.apply_method(OneRegion.extend_reg, -flank_extend)	#Cut to needed size knowing that the region will be extended in function

	#Remove blacklisted regions and chromosomes not in common
	blacklist_regions = RegionList().from_bed(args.blacklist) if args.blacklist != None else RegionList([])	 #fill in with regions from args.blacklist
	regions_dict = {"genome": genome_regions, "input_regions":input_regions, "output_regions":output_regions, "peak_regions":peak_regions, "nonpeak_regions":nonpeak_regions, "blacklist_regions": blacklist_regions}
	for sub in ["input_regions", "output_regions", "peak_regions", "nonpeak_regions"]:
		regions_sub = regions_dict[sub]
		regions_sub.subtract(blacklist_regions)
		regions_sub = regions_sub.apply_method(OneRegion.split_region, 50000)

		regions_sub.keep_chroms(chrom_in_common)
		regions_dict[sub] = regions_sub
	
	#write beds to look at in igv
	#input_regions.write_bed(os.path.join(args.outdir, "input_regions.bed"))
	#output_regions.write_bed(os.path.join(args.outdir, "output_regions.bed"))
	#peak_regions.write_bed(os.path.join(args.outdir, "peak_regions.bed"))
	#nonpeak_regions.write_bed(os.path.join(args.outdir, "nonpeak_regions.bed"))

	#Sort according to order in bam_references:
	output_regions.loc_sort(bam_references)
	chrom_order = {bam_references[i]:i for i in range(len(bam_references))}	 #for use later when sorting output 

	#### Statistics about regions ####
	genome_bp = sum([region.get_length() for region in regions_dict["genome"]])
	for key in regions_dict:
		total_bp = sum([region.get_length() for region in regions_dict[key]])
		logger.stats("{0}: {1} regions | {2} bp | {3:.2f}% coverage".format(key, len(regions_dict[key]), total_bp, total_bp/genome_bp*100))

	#Estallish variables for regions to be used
	input_regions = regions_dict["input_regions"]
	output_regions = regions_dict["output_regions"]
	peak_regions = regions_dict["peak_regions"]
	nonpeak_regions = regions_dict["nonpeak_regions"]

	#Exit if no input/output regions were found
	if len(input_regions) == 0 or len(output_regions) == 0 or len(peak_regions) == 0 or len(nonpeak_regions) == 0:
		logger.error("No regions found - exiting!")
		sys.exit()

	#----------------------------------------------------------------------------------------------------#
	# Estimate normalization factors
	#----------------------------------------------------------------------------------------------------#

	#Setup logger queue
	logger.debug("Setting up listener for log")
	logger.start_logger_queue()
	args.log_q = logger.queue

	#----------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.info("----- Estimating normalization factors -----")

	#If normalization is to be calculated
	if not args.norm_off:

		#Reads in peaks/nonpeaks
		logger.info("Counting reads in peak regions")
		peak_region_chunks = peak_regions.chunks(args.split)
		reads_peaks = sum(run_parallel(count_reads, peak_region_chunks, [args], args.cores, logger, "Counting (peaks):"))
		logger.comment("")

		logger.info("Counting reads in nonpeak regions")
		nonpeak_region_chunks = nonpeak_regions.chunks(args.split)
		reads_nonpeaks = sum(run_parallel(count_reads, nonpeak_region_chunks, [args], args.cores, logger, "Counting (nonpeaks):"))

		reads_total = reads_peaks + reads_nonpeaks

		logger.stats("TOTAL_READS\t{0}".format(reads_total))
		logger.stats("PEAK_READS\t{0}".format(reads_peaks))
		logger.stats("NONPEAK_READS\t{0}".format(reads_nonpeaks))

		lib_norm = 10000000/reads_total
		frip = reads_peaks/reads_total
		correct_factor = lib_norm*(1/frip)

		logger.stats("LIB_NORM\t{0:.5f}".format(lib_norm))
		logger.stats("FRiP\t{0:.5f}".format(frip))
	else:
		logger.info("Normalization was switched off")
		correct_factor = 1.0

	logger.stats("CORRECTION_FACTOR:\t{0:.5f}".format(correct_factor))

	#----------------------------------------------------------------------------------------------------#
	# Estimate sequence bias
	#----------------------------------------------------------------------------------------------------#

	logger.comment("")

	if args.bias_pkl is None:

		logger.info("Started estimation of sequence bias...")

		input_region_chunks = input_regions.chunks(args.split)										# split to 100 chunks (also decides the step of output)
		out_lst = run_parallel(bias_estimation, input_region_chunks, [args], args.cores, logger, "Estimating bias:")	# Output is list of AtacBias objects

		#Join objects
		estimated_bias = out_lst[0]		#initialize object with first output
		for output in out_lst[1:]:
			estimated_bias.join(output)		#bias object contains bias/background SequenceMatrix objects

		logger.debug("Bias estimated\tno_reads: {0}".format(estimated_bias.no_reads))

	else:
		logger.info("Loading sequence bias from '--bias-pkl' file...")
		estimated_bias = AtacBias().from_pickle(args.bias_pkl)

	#----------------------------------------------------------------------------------------------------#
	# Join estimations from all chunks of regions
	#----------------------------------------------------------------------------------------------------#

	bias_obj = estimated_bias
	bias_obj.correction_factor = correct_factor

	### Bias motif ###
	logger.info("Finalizing bias motif for scoring")
	for strand in strands:
		bias_obj.bias[strand].prepare_mat()

		if figure_pdf is not None:
			logger.debug("Saving pssm to figure pdf")
			fig = plot_pssm(bias_obj.bias[strand].pssm, "Tn5 insertion bias of reads ({0})".format(strand))
			figure_pdf.savefig(fig)

	
	#Write bias motif to pickle
	out_f = os.path.join(args.outdir, args.prefix + "_AtacBias.pickle")
	logger.debug("Saving bias object to pickle ({0})".format(out_f))
	bias_obj.to_pickle(out_f)
	
	#----------------------------------------------------------------------------------------------------#
	# Correct read bias and write to bigwig
	#----------------------------------------------------------------------------------------------------#

	logger.comment("")
	logger.info("----- Correcting reads from .bam within output regions -----")

	output_regions.loc_sort(bam_references)		#sort in order of references
	output_regions_chunks = output_regions.chunks(args.split)
	no_tasks = float(len(output_regions_chunks))
	chunk_sizes = [len(chunk) for chunk in output_regions_chunks]
	logger.debug("All regions chunked: {0} ({1})".format(len(output_regions), chunk_sizes))

	### Create key-file linking for bigwigs 
	key2file = {}
	for track in output_bws:
		for strand in output_bws[track]:
			filename = output_bws[track][strand]["fn"]
			key = "{}:{}".format(track, strand)
			key2file[key] = filename

	#Start correction/write cores
	n_bigwig = len(key2file.values())
	writer_cores = min(n_bigwig, max(1,int(args.cores*0.1)))	#at most one core per bigwig or 10% of cores (or 1)
	worker_cores = max(1, args.cores - writer_cores) 				
	logger.debug("Worker cores: {0}".format(worker_cores))
	logger.debug("Writer cores: {0}".format(writer_cores))

	worker_pool = mp.Pool(processes=worker_cores)
	writer_pool = mp.Pool(processes=writer_cores)
	manager = mp.Manager()

	#Start bigwig file writers
	writer_tasks = []
	header = [(chrom, bam_chrom_info[chrom]) for chrom in bam_references]
	key_chunks = [list(key2file.keys())[i::writer_cores] for i in range(writer_cores)]
	qs_list = []
	qs = {}
	for chunk in key_chunks:
		logger.debug("Creating writer queue for {0}".format(chunk))

		q = manager.Queue()
		qs_list.append(q)

		files = [key2file[key] for key in chunk]
		writer_tasks.append(writer_pool.apply_async(bigwig_writer, args=(q, dict(zip(chunk, files)), header, output_regions, args)))	 #, callback = lambda x: finished.append(x) print("Writing time: {0}".format(x)))
		for key in chunk:
			qs[key] = q

	args.qs = qs
	writer_pool.close() #no more jobs applied to writer_pool

	#Start correction
	logger.debug("Starting correction")
	task_list = [worker_pool.apply_async(bias_correction, args=[chunk, args, bias_obj]) for chunk in output_regions_chunks]
	worker_pool.close()
	monitor_progress(task_list, logger, "Correction progress:")	#does not exit until tasks in task_list finished
	results = [task.get() for task in task_list]

	#Get all results
	if getattr(args, "skip_qc", False):
		pre_bias = None
		post_bias = None
	else:
		pre_bias = results[0][0]	#initialize with first result
		post_bias = results[0][1]	#initialize with first result
		for result in results[1:]:
			pre_bias_chunk = result[0]
			post_bias_chunk = result[1]

			for direction in strands:
				pre_bias[direction].add_counts(pre_bias_chunk[direction])
				post_bias[direction].add_counts(post_bias_chunk[direction])

	#Stop all queues for writing
	logger.debug("Stop all queues by inserting None")
	for q in qs_list:
		q.put((None, None, None))

	#Fetch error codes from bigwig writers
	logger.debug("Fetching possible errors from bigwig_writer tasks")
	results = [task.get() for task in writer_tasks]	#blocks until writers are finished

	logger.debug("Joining bigwig_writer queues")
	
	qsum = sum([q.qsize() for q in qs_list])
	while qsum != 0:
		qsum = sum([q.qsize() for q in qs_list])
		logger.spam("- Queue sizes {0}".format([(key, qs[key].qsize()) for key in qs]))
		time.sleep(0.5)

	#Waits until all queues are closed
	writer_pool.join() 
	worker_pool.terminate()
	worker_pool.join()

	#Stop multiprocessing logger	
	logger.stop_logger_queue()

	#----------------------------------------------------------------------------------------------------#
	# Information and verification of corrected read frequencies
	#----------------------------------------------------------------------------------------------------#		

	if getattr(args, "skip_qc", False):
		logger.info("Skipped atac-correct diagnostic PDF and bias-correction verification counts (--skip-qc)")
	else:
		logger.comment("")
		logger.info("Verifying bias correction")

		#Calculating variance per base
		for strand in strands:

			#Invert negative counts
			abssum = np.abs(np.sum(post_bias[strand].neg_counts, axis=0))
			post_bias[strand].neg_counts = post_bias[strand].neg_counts + abssum

			#Join negative/positive counts
			post_bias[strand].counts += post_bias[strand].neg_counts	#now pos

			pre_bias[strand].prepare_mat()
			post_bias[strand].prepare_mat()

			pre_var = np.mean(np.var(pre_bias[strand].bias_pwm, axis=1)[:4])   #mean of variance per nucleotide
			post_var = np.mean(np.var(post_bias[strand].bias_pwm, axis=1)[:4])
			logger.stats("BIAS\tpre-bias variance {0}:\t{1:.7f}".format(strand, pre_var))
			logger.stats("BIAS\tpost-bias variance {0}:\t{1:.7f}".format(strand, post_var))

			#Plot figure
			fig_title = "Nucleotide frequencies in corrected reads\n({0} strand)".format(strand)
			figure_pdf.savefig(plot_correction(pre_bias[strand].bias_pwm, post_bias[strand].bias_pwm, fig_title))

	
	#----------------------------------------------------------------------------------------------------#
	# Finish up
	#----------------------------------------------------------------------------------------------------#

	plt.close('all')
	if figure_pdf is not None:
		figure_pdf.close()
	corrected_bigwigs = []
	if "corrected" in output_bws:
		corrected_bigwigs = [output_bws["corrected"][strand]["fn"] for strand in output_bws["corrected"]]
	if getattr(args, "_scale_after_single", True):
		try:
			_maybe_scale_corrected_bigwigs(args, corrected_bigwigs, logger)
		except Exception as exc:
			logger.error("Corrected-track q95 scaling failed: {0}".format(exc))
			raise
	logger.end()

#--------------------------------------------------------------------------------------------------------#
if __name__ == '__main__':

	parser = argparse.ArgumentParser()
	parser = add_atacorrect_arguments(parser)
	args = parser.parse_args()

	if len(sys.argv[1:]) == 0:
		parser.print_help()
		sys.exit()

	run_atacorrect(args)
