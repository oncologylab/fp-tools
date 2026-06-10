#!/usr/bin/env python

"""
ATACorrect command driver for bias estimation and cutsite correction.

This module keeps the main command-line execution path for:
- reading BAM/FASTA/peak inputs
- estimating Tn5 sequence bias
- writing corrected and auxiliary bigWig tracks
- generating the ATACorrect summary PDF
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
from copy import deepcopy

from collections import OrderedDict
import itertools
import matplotlib
matplotlib.use("Agg")  #non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#Bio-specific packages
import pyBigWig
import pysam

#Internal functions and classes
from fp_tools.parsers import add_atacorrect_arguments
from fp_tools.tools.atacorrect_functions import *
from fp_tools.utils.utilities import *
from fp_tools.utils.regions import OneRegion, RegionList
from fp_tools.utils.ngs import OneRead, ReadList
from fp_tools.utils.sequences import *
from fp_tools.utils.logger import FpToolsLogger

#np.seterr(divide='raise', invalid='raise')

#--------------------------------------------------------------------------------------------------------#
#-------------------------------------- Main pipeline function ------------------------------------------# 
#--------------------------------------------------------------------------------------------------------#

def run_atacorrect(args):

	"""
	Function for bias correction of input .bam files
	Calls functions in ATACorrect_functions and several internal classes
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
	tracks = ["uncorrected", "bias", "expected", "corrected"]
	tracks = [track for track in tracks if track not in args.track_off] 	# switch off printing

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
	
	output_files = bigwigs + [figures_f]
	output_files = list(OrderedDict.fromkeys(output_files)) 	#remove duplicates due to "both" option

	strands = ["forward", "reverse"]

	#----------------------------------------------------------------------------------------------------#
	# Print info on run
	#----------------------------------------------------------------------------------------------------#

	logger = FpToolsLogger("ATACorrect", args.verbosity)
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
	figure_pdf = PdfPages(figures_f, keep_empty=False)

	#----------------------------------------------------------------------------------------------------#
	# Read information in bam/fasta
	#----------------------------------------------------------------------------------------------------#

	logger.info("Reading info from .bam file")
	bamfile = pysam.AlignmentFile(args.bam, "rb")
	if bamfile.has_index() == False:
		logger.warning("No index found for bamfile - creating one via pysam.")
		pysam.index(args.bam)

	bam_references = bamfile.references 	#chromosomes in correct order
	bam_chrom_info = dict(zip(bamfile.references, bamfile.lengths))
	logger.debug("bam_chrom_info: {0}".format(bam_chrom_info))
	bamfile.close()

	logger.info("Reading info from .fasta file")
	fastafile = pysam.FastaFile(args.genome)
	fasta_chrom_info = dict(zip(fastafile.references, fastafile.lengths))
	logger.debug("fasta_chrom_info: {0}".format(fasta_chrom_info))
	fastafile.close()

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
		logger.warning("The following contigs in --bam did not have sequences in --fasta: {0}. NOTE: These contigs will be skipped in calculation and output.".format(chrom_not_in_fasta))

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
		logger.error("No common contigs left to run ATACorrect on. Please check that '--bam' and '--fasta' are matching.")
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
	figure_pdf.close()
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
