#!/usr/bin/env python

"""
Helper classes and routines used by the ATACorrect workflow.

This module contains the lower-level ATACorrect logic for:
- storing learned bias models
- estimating expected cutsite signal from sequence bias
- correcting observed signal
- generating ATACorrect diagnostic plots
"""

import os
import sys
import gc
import numpy as np
import multiprocessing as mp
import time
from datetime import datetime
import matplotlib
matplotlib.use('Agg')	#non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
import pickle

#Bio-specific packages
import pysam

#Internal functions and classes
from fp_tools.utils.sequences import SequenceMatrix, GenomicSequence
from fp_tools.utils.signals import fast_rolling_math
from fp_tools.utils.utilities import * 
from fp_tools.utils.regions import OneRegion, RegionList
from fp_tools.utils.ngs import OneRead, ReadList
from fp_tools.utils.logger import FpToolsLogger
from fp_tools.utils.plotting_style import PDF_FONT_SIZE, apply_pdf_style

#Catch warnings from curve_fit
import warnings
from scipy.optimize import OptimizeWarning
warnings.simplefilter("error", OptimizeWarning)
apply_pdf_style()

#--------------------------------------------------------------------------------------------------#
class AtacBias:
	""" Class for storing information about estimated bias """

	def __init__(self, L=10, stype="DWM"):

		self.stype = stype	#type of score matrix
		self.bias = {"forward": SequenceMatrix.create(L, self.stype),
					 "reverse": SequenceMatrix.create(L, self.stype),
					 "both": SequenceMatrix.create(L, self.stype)}
		self.no_reads = 0

	def join(self, obj):
		""" Join counts from AtacBias obj with other AtacBias obj """

		self.bias["forward"].add_counts(obj.bias["forward"])
		self.bias["reverse"].add_counts(obj.bias["reverse"])
		self.bias["both"].add_counts(obj.bias["both"])
		self.no_reads += obj.no_reads

	def to_pickle(self, f):
		""" Pickle an AtacBias object to a .pickle file """

		handle = open(f, "wb")
		pickle.dump(self, handle)
		handle.close()

		return(self)

	def from_pickle(self, f):
		""" Read an AtacBias object from a .pickle file """

		handle = open(f, "rb")
		self = pickle.load(handle)
		return(self)

#--------------------------------------------------------------------------------------------------#

def count_reads(regions_list, params):
	""" Count reads from bam within regions (counts position of cutsite to prevent double-counting) """

	bam_f = params.bam
	read_shift = params.read_shift
	bam_obj = pysam.AlignmentFile(bam_f, "rb")

	log_q = params.log_q
	logger = FpToolsLogger("", params.verbosity, log_q) #sending all logger calls to log_q

	#Count per region
	read_count = 0
	logger.spam("Started counting region_chunk ({0} -> {1})".format("_".join([str(element) for element in regions_list[0]]), "_".join([str(element) for element in regions_list[-1]])))
	
	for region in regions_list:
		read_lst = ReadList().from_bam(bam_obj, region)
		logger.spam("- {0} ({1} reads)".format(region, len(read_lst)))
		for read in read_lst:  
			read.get_cutsite(read_shift)
			if read.cutsite > region.start and read.cutsite <= region.end:  #only reads within borders
				read_count += 1
				
	logger.spam("Finished counting region_chunk ({0} -> {1})".format("_".join([str(element) for element in regions_list[0]]), "_".join([str(element) for element in regions_list[-1]])))
	bam_obj.close()

	return(read_count)


#--------------------------------------------------------------------------------------------------#
def bias_estimation(regions_list, params):
	""" Estimates bias of insertions within regions """

	#Info on run
	bam_f = params.bam
	fasta_f = params.genome
	k_flank = params.k_flank
	bg_shift = params.bg_shift
	read_shift = params.read_shift
	L = 2 * k_flank + 1

	logger = FpToolsLogger("", params.verbosity, params.log_q) 	#sending all logger calls to log_q

	#Open objects for reading
	bam_obj = pysam.AlignmentFile(bam_f, "rb")
	fasta_obj = pysam.FastaFile(fasta_f)
	chrom_lengths = dict(zip(bam_obj.references, bam_obj.lengths))  #Chromosome boundaries from bam_obj

	bias_obj = AtacBias(L, params.score_mat)

	strands = ["forward", "reverse"]

	#Estimate bias at each region
	for region in regions_list:

		read_lst = ReadList().from_bam(bam_obj, region)  
		for read in read_lst:
			read.get_cutsite(read_shift)

		## Kmer cutting bias ##
		if len(read_lst) > 0:

			#Extract sequence
			extended_region = region.extend_reg(k_flank + bg_shift) 				#Extend to allow full kmers
			extended_region.check_boundary(chrom_lengths, "cut")
			sequence_obj = GenomicSequence(extended_region).from_fasta(fasta_obj)
			
			#Split reads forward/reverse
			for_lst, rev_lst = read_lst.split_strands()
			read_lst_strand = {"forward": for_lst, "reverse": rev_lst} 
			logger.spam("Region: {0}. Forward reads: {1}. Reverse reads: {2}".format(region, len(for_lst), len(rev_lst)))

			for strand in strands:

				#Map reads to positions
				read_per_pos = {}
				for read in read_lst_strand[strand]:
					if read.cigartuples is not None:
						first_tuple = read.cigartuples[-1] if read.is_reverse else read.cigartuples[0]
						if first_tuple[0] == 0:	#Only include non-clipped reads
							read_per_pos[read.cutsite] = read_per_pos.get(read.cutsite, []) + [read]

				n_reads_kept = sum([len(l) for l in read_per_pos.values()]) #total number of reads across all cutsites in dict
				logger.spam("Region: {0} ({1}) | Number of positions cut: {2} | Sum of read_per_pos: {3}".format(region.tup(), strand, len(read_per_pos), n_reads_kept))

				#Get kmer for each position
				for cutsite in read_per_pos:
					if cutsite > region.start and cutsite < region.end:  	#only reads within borders
						read = read_per_pos[cutsite][0] 					#use first read in list to establish kmer
						no_cut = min(len(read_per_pos[cutsite]), 10)		#put cap on number of cuts to limit influence of outliers

						read.get_kmer(sequence_obj, k_flank)

						bias_obj.bias[strand].add_sequence(read.kmer, no_cut)
						read.shift_cutsite(-bg_shift)			#upstream of read; ensures that bg is not within fragment
						read.get_kmer(sequence_obj, k_flank)	#kmer updated to kmer for shifted read
						bias_obj.bias[strand].add_background(read.kmer, no_cut)

						bias_obj.no_reads += no_cut
	bam_obj.close()
	fasta_obj.close()

	return(bias_obj) 	#object containing information collected on bias 	


#--------------------------------------------------------------------------------------------------#
def relu(x, a, b):
	""" a and b are components of a linear curve (y=a*x+b) """
	y = np.maximum(0.0, a*x + b)
	return(y)

#--------------------------------------------------------------------------------------------------#
def bias_correction(regions_list, params, bias_obj):
	""" Corrects bias in cutsites (from bamfile) using estimated bias """

	logger = FpToolsLogger("", params.verbosity, params.log_q)

	bam_f = params.bam
	fasta_f = params.genome
	k_flank = params.k_flank
	read_shift = params.read_shift
	L = 2 * k_flank + 1
	w = params.window
	f = int(w/2.0)
	qs = params.qs

	f_extend = k_flank + f

	strands = ["forward", "reverse"]
	pre_bias = {strand: SequenceMatrix.create(L, "PWM") for strand in strands}
	post_bias = {strand: SequenceMatrix.create(L, "PWM") for strand in strands}

	#Open bamfile and fasta
	bam_obj = pysam.AlignmentFile(bam_f, "rb")
	fasta_obj = pysam.FastaFile(fasta_f)

	out_signals = {}
	
	#Go through each region
	for region_obj in regions_list:

		region_obj.extend_reg(f_extend)
		reg_len = region_obj.get_length()	#length including flanking
		reg_key = (region_obj.chrom, region_obj.start+f_extend, region_obj.end-f_extend)	#output region
		out_signals[reg_key] = {"uncorrected":{}, "bias":{}, "expected":{}, "corrected":{}}


		################################
		####### Uncorrected reads ######
		################################

		#Get cutsite positions for each read
		read_lst = ReadList().from_bam(bam_obj, region_obj) 
		for read in read_lst:
			read.get_cutsite(read_shift)
		logger.spam("Read {0} reads from region {1}".format(len(read_lst), region_obj))

		#Exclude reads with cutsites outside region
		read_lst = ReadList([read for read in read_lst if read.cutsite > region_obj.start and read.cutsite < region_obj.end])
		for_lst, rev_lst = read_lst.split_strands()
		read_lst_strand = {"forward": for_lst, "reverse": rev_lst}

		for strand in strands:
			out_signals[reg_key]["uncorrected"][strand] = read_lst_strand[strand].signal(region_obj) 
			out_signals[reg_key]["uncorrected"][strand] = np.round(out_signals[reg_key]["uncorrected"][strand], 5)


		################################
		###### Estimation of bias ######
		################################

		#Get sequence in this region
		sequence_obj = GenomicSequence(region_obj).from_fasta(fasta_obj)
		
		#Score sequence using forward/reverse motifs
		for strand in strands:
			if strand == "forward": 
				seq = sequence_obj.sequence
				bias = bias_obj.bias[strand].score_sequence(seq)
			elif strand == "reverse":
				seq = sequence_obj.revcomp
				bias = bias_obj.bias[strand].score_sequence(seq)[::-1]  #3'-5'

			out_signals[reg_key]["bias"][strand] = np.nan_to_num(bias) 		#convert any nans to 0


		#################################
		###### Correction of reads ######
		#################################
		
		reg_end = reg_len - k_flank
		step = 10
		overlaps = int(params.window / step)
		window_starts = list(range(k_flank, reg_end-params.window, step))
		window_ends = list(range(k_flank+params.window, reg_end, step))
		windows = list(zip(window_starts, window_ends))

		for strand in strands:

			########### Estimate bias threshold ###########
			bias_predictions = np.zeros((overlaps,reg_len))
			bias_predictions[k_flank:reg_end] = np.nan #flanks stay 0 as no windows overlap
			row = 0

			for window in windows:

				signal_w = out_signals[reg_key]["uncorrected"][strand][window[0]:window[1]]
				bias_w = out_signals[reg_key]["bias"][strand][window[0]:window[1]]

				signalmax = np.max(signal_w)

				if signalmax > 0:
					try:
						popt, pcov = curve_fit(relu, bias_w, signal_w)
						bias_predict = relu(bias_w, *popt)

					except (OptimizeWarning, RuntimeError):
						cut_positions = np.logical_not(np.isclose(signal_w, 0))
						bias_min = np.min(bias_w[cut_positions])
						bias_predict = bias_w - bias_min
						bias_predict[bias_predict < 0] = 0

					if np.max(bias_predict) > 0:
						bias_predict = bias_predict / np.max(bias_predict)
				else:
					bias_predict = np.zeros(window[1]-window[0])	

				bias_predictions[row, window[0]:window[1]] = bias_predict
				row = row + 1 if row < overlaps - 1 else 0

				bias_prediction = np.nanmean(bias_predictions, axis=0) #nanmean because ends of array contain nan (before windows are completely overlapping)

			bias = bias_prediction 

			######## Calculate expected signal ######
			signal_sum = fast_rolling_math(out_signals[reg_key]["uncorrected"][strand], w, "sum")
			signal_sum[np.isnan(signal_sum)] = 0 	#f-width ends of region

			bias_sum = fast_rolling_math(bias, w, "sum")	#ends of arr are nan
			nulls = np.logical_or(np.isclose(bias_sum, 0), np.isnan(bias_sum))
			bias_sum[nulls] = 1 		# N-regions will give stretches of 0-bias
			bias_probas = bias / bias_sum
			bias_probas[nulls] = 0 		#nan to 0

			out_signals[reg_key]["expected"][strand] = signal_sum * bias_probas 

			######## Correct signal ########
			out_signals[reg_key]["uncorrected"][strand] *= bias_obj.correction_factor
			out_signals[reg_key]["expected"][strand] *= bias_obj.correction_factor
			out_signals[reg_key]["corrected"][strand] = out_signals[reg_key]["uncorrected"][strand] - out_signals[reg_key]["expected"][strand]

			######## Rescale signal to fit uncorrected sum ########
			uncorrected_sum = fast_rolling_math(out_signals[reg_key]["uncorrected"][strand], w, "sum")
			uncorrected_sum[np.isnan(uncorrected_sum)] = 0
			corrected_sum = fast_rolling_math(np.abs(out_signals[reg_key]["corrected"][strand]), w, "sum")	#negative values count as positive
			corrected_sum[np.isnan(corrected_sum)] = 0

			#Positive signal left after correction
			corrected_pos = np.copy(out_signals[reg_key]["corrected"][strand])
			corrected_pos[corrected_pos < 0] = 0 
			corrected_pos_sum = fast_rolling_math(corrected_pos, w, "sum")
			corrected_pos_sum[np.isnan(corrected_pos_sum)] = 0 
			corrected_neg_sum = corrected_sum - corrected_pos_sum

			#The corrected sum is less than the signal sum, so scale up positive cuts
			zero_sum = corrected_pos_sum == 0
			corrected_pos_sum[zero_sum] = np.nan 	#allow for zero division
			scale_factor = (uncorrected_sum - corrected_neg_sum) / corrected_pos_sum
			scale_factor[zero_sum] = 1			#Scale factor is 1 (which will be multiplied to the 0 values)
			scale_factor[scale_factor < 1] = 1	#Only scale up if needed
			pos_bool = out_signals[reg_key]["corrected"][strand] > 0
			out_signals[reg_key]["corrected"][strand][pos_bool] *= scale_factor[pos_bool]

			
		#######################################
		########   Verify correction   ########
		#######################################
		
		#Verify correction across all reads
		for strand in strands:
			for idx in range(k_flank,reg_len - k_flank -1): 
				if idx > k_flank and idx < reg_len-k_flank:

					orig = out_signals[reg_key]["uncorrected"][strand][idx]
					correct = out_signals[reg_key]["corrected"][strand][idx]

					if orig != 0 or correct != 0:	#if both are 0, don't add to pre/post bias
						if strand == "forward":
							kmer = sequence_obj.sequence[idx-k_flank:idx+k_flank+1]
						else:
							kmer = sequence_obj.revcomp[reg_len-idx-k_flank-1:reg_len-idx+k_flank]

						#Save kmer for bias correction verification
						pre_bias[strand].add_sequence(kmer, orig)
						post_bias[strand].add_sequence(kmer, correct)


		#######################################
		########    Write to queue    #########
		#######################################

		#Set size back to original
		for track in out_signals[reg_key]:
			for strand in out_signals[reg_key][track]:
				out_signals[reg_key][track][strand] = out_signals[reg_key][track][strand][f_extend:-f_extend]

		#Calculate "both" if split_strands == False
		if params.split_strands == False:
			for track in out_signals[reg_key]:
				out_signals[reg_key][track]["both"] = out_signals[reg_key][track]["forward"] + out_signals[reg_key][track]["reverse"]

		#Send to queue
		strands_to_write = ["forward", "reverse"] if params.split_strands == True else ["both"]
		for track in out_signals[reg_key]:

			#Send to writer per strand
			for strand in strands_to_write:
				key = "{0}:{1}".format(track, strand)

				if key in qs: #only write the signals where the files were initialized
					logger.spam("Sending {0} signal from region {1} to writer queue".format(key, reg_key))
					qs[key].put((key, reg_key, out_signals[reg_key][track][strand]))

		#Sent to qs - delete from this process
		out_signals[reg_key] = None

	bam_obj.close()
	fasta_obj.close()
	
	gc.collect()

	return([pre_bias, post_bias])



####################################################################################################
######################################## Plot functions ############################################
####################################################################################################

colors = {0:"green", 1:"red", 2:"blue", 3:"darkkhaki"}
names = {0:"A", 1:"T", 2:"C", 3:"G"}

def plot_pssm(matrix, title):
	""" Plot pssm in matrix """

	#Make figure
	fig, ax = plt.subplots()
	fig.suptitle(title, fontsize=PDF_FONT_SIZE, weight="bold")
	
	#Formatting of x axis
	length = matrix.shape[1]
	flank = int(length/2.0)		
	xvals = np.arange(length)  # each position corresponds to i in mat
	
	#Customize minor tick labels
	xtick_pos = xvals[:-1] + 0.5
	xtick_labels = list(range(-flank, flank))
	ax.xaxis.set_major_locator(ticker.FixedLocator(xvals))
	ax.xaxis.set_major_formatter(ticker.FixedFormatter(xtick_labels))	
	ax.xaxis.set_minor_locator(ticker.FixedLocator(xtick_pos))			#locate minor ticks between major ones (cutsites)
	ax.xaxis.set_minor_formatter(ticker.NullFormatter())

	#Make background grid on major ticks
	plt.grid(color='0.8', which="minor", ls="--", axis="x")

	plt.xlim([0, length-1])
	plt.xlabel('Position from cutsite')
	plt.ylabel('PSSM score')

	######## Plot data #######
	#Plot PSSM / bias motif
	for nuc in range(4):
		plt.plot(xvals, matrix[nuc,:], color=colors[nuc], label=names[nuc])

	#Cutsite-line
	plt.axvline(flank-0.5, linewidth=2, color="black", zorder=100)

	#Finish up
	plt.legend(loc="lower right", prop={'size': PDF_FONT_SIZE, 'weight': 'bold'})
	plt.tight_layout()
	fig.subplots_adjust(top=0.88, hspace=0.5)

	return(fig)

#----------------------------------------------------------------------------------------------------#
def plot_correction(pre_mat, post_mat, title):
	""" Plot comparison of pre-correction and post-correction matrices """

	#Make figure
	fig, ax = plt.subplots()
	fig.suptitle(title, fontsize=PDF_FONT_SIZE, weight="bold")

	L = pre_mat.shape[1]
	flank = int(L/2.0)		
	xvals = np.arange(L)  # each position corresponds to i in mat
	
	#Customize minor tick labels
	xtick_pos = xvals[:-1] + 0.5
	xtick_labels = list(range(-flank, flank))		#-flank - flank without 0
	ax.xaxis.set_major_locator(ticker.FixedLocator(xvals))
	ax.xaxis.set_major_formatter(ticker.FixedFormatter(xtick_labels))	
	ax.xaxis.set_minor_locator(ticker.FixedLocator(xtick_pos))			#locate minor ticks between major ones (cutsites)
	ax.xaxis.set_minor_formatter(ticker.NullFormatter())
	
	#PWMs for all mats
	pre_pwm = pre_mat
	post_pwm = post_mat

	#Pre correction
	for nuc in range(4):
		yvals = [pre_pwm[nuc, m] for m in range(L)]
		plt.plot(xvals, yvals, linestyle="--", color=colors[nuc], linewidth=1, alpha=0.5)
	
	#Post correction
	for nuc in range(4):
		yvals = [post_pwm[nuc, m] for m in range(L)]
		plt.plot(xvals, yvals, color=colors[nuc], linewidth=2, label=names[nuc])
	
	plt.xlim([0, L-1])
	plt.xlabel('Position from cutsite')
	plt.ylabel('Nucleotide frequency')

	#Set legend
	plt.plot([0],[0], linestyle="--", linewidth=1, color="black", label="pre-correction")
	plt.plot([0],[0], color="black", label="post-correction")
	plt.legend(loc="lower right", prop={'size': PDF_FONT_SIZE, 'weight': 'bold'})

	plt.tight_layout()
	fig.subplots_adjust(top=0.88, hspace=0.4)

	return(fig)
