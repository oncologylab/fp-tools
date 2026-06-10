# cython: language_level=3

"""Sequence utilities for counting and scoring with PWMs and DWMs."""

import numpy as np
cimport numpy as np
import cython
from libc.math cimport log2
import sys

import pysam
complement = {0:1, 1:0, 2:3, 3:2}


#--------------------------------------------------------------------------------------------------#
class SequenceMatrix:
	""" Class-builder """

	@staticmethod
	def create(length, stype):
		if stype == "PWM":
			return(NucleotideMatrix(length))
		elif stype == "DWM":
			return(DiNucleotideMatrix(length))
		else:
			sys.exit("Unkown score type {0}".format(stype))

	def add_counts(self, obj):
		
		self.counts += obj.counts
		self.bg_counts += obj.bg_counts
		self.neg_counts += obj.neg_counts
		
		self.no_bias += obj.no_bias
		self.no_bg += obj.no_bg


#-----------------------------------------------------#
class NucleotideMatrix(SequenceMatrix):

	def __init__(self, length):

		self.length = length
		self.counts = np.zeros((5, length))			#bias counts
		self.neg_counts = np.zeros((5, length))		#bias counts for negative cuts
		self.bg_counts = np.zeros((5, length))		#background counts
		
		self.no_bias = 0	#number of sequences in bias mat
		self.no_bg = 0		#number of sequences in bg mat

		#Added by prepare_mat
		self.pssm = np.zeros((5, length)) 	#ATCGN

	@cython.boundscheck(False)
	@cython.wraparound(False)
	def add_sequence(self, np.ndarray[np.int64_t, ndim=1] sequence, double amount = 1.0):
		""" Adds sequence to SequenceMatrix the number of times specified in amount """

		cdef int Sm, m
		cdef int length = self.length
		cdef np.ndarray[np.float64_t, ndim=2] counts = self.counts
		cdef np.ndarray[np.float64_t, ndim=2] neg_counts = self.neg_counts
		cdef int So 	#S_other

		#Check if seq contains N
		for m in range(length):
			Sm = sequence[m]	#ATCG(N)
			if Sm == 4:
				return

		for m in range(length):
			Sm = sequence[m]

			if amount > 0:
				counts[Sm,m] += amount
			else:
				neg_counts[Sm,m] += amount	#addition of minus = subtraction

		self.counts = counts
		self.neg_counts = neg_counts
		self.no_bias += amount


	@cython.boundscheck(False)
	@cython.wraparound(False)		
	def add_background(self, np.ndarray[np.int64_t, ndim=1] sequence, double amount = 1.0):
		""" Adds sequence to count of background nucleotides """

		cdef int Sm, m
		cdef int length = self.length
		cdef np.ndarray[np.float64_t, ndim=2] bg_counts = self.bg_counts

		for m in range(length):
			Sm = sequence[m]
			bg_counts[Sm, m] += amount

		self.bg_counts = bg_counts
		self.no_bg += amount

	def prepare_mat(self):
		""" Prepare matrices for scoring with score_sequence """

		self.no_bias = np.mean(np.sum(self.counts[:4,:], axis=0))
		self.no_bg = np.mean(np.sum(self.bg_counts[:4,:], axis=0))

		self.bias_pwm = np.true_divide(self.counts + 1, self.no_bias + 4)
		self.bg_pwm = np.true_divide(self.bg_counts + 1, self.no_bg + 4)

		#Prepare PSSM for scoring
		self.pssm = np.log2(np.true_divide(self.bias_pwm, self.bg_pwm))


	@cython.boundscheck(False)	#dont check boundaries
	@cython.cdivision(True)		#no check for zero division
	@cython.wraparound(False) 	#dont deal with negative indices
	def score_sequence(self, np.ndarray[np.int64_t, ndim=1] sequence):
		""" Score nucleotide sequence against motif """ 

		cdef np.ndarray[np.float64_t, ndim=2] pssm = self.pssm
		cdef int reglen = len(sequence)
		cdef int L = pssm.shape[1]
		cdef int F = int(L/2.0) 
		cdef np.ndarray[np.float64_t, ndim=1] scores = np.zeros(len(sequence))		#array to fill in
		scores[:] = np.nan
		cdef int i, m, Sm, contain_N
		cdef double score

		#PSSM matrix scoring
		for i in range(reglen - L + 1):
			contain_N = 0
			score = 0.0
			for m in range(L): 	#m loops from start to end of motif
				Sm = sequence[i + m]
				if Sm > 3:
					contain_N = 1

				score += pssm[Sm, m]

			#Assign score
			if contain_N == 0:
				scores[i+F] = score 	#logscale

		#Return score array
		return(scores)	


#-----------------------------------------------------#
class DiNucleotideMatrix(SequenceMatrix):
	""" Dinucleotide matrix """

	def __init__(self, length):
		self.length = length
		self.counts = np.zeros((5, 5, length, length))
		self.neg_counts = np.zeros((5, length))		#not used, only aded to enable joining of counts in add_counts
		self.bg_counts = np.zeros((5, 5, length, length))

		self.no_bias = 0	#number of sequences in bias mat
		self.no_bg = 0		#number of sequences in bg mat

		#Added later by prepare_mat
		self.pssm = np.zeros((5,length))

		self.bias_pwm = np.zeros((5,5,length,length))
		self.bg_pwm = np.zeros((5,5,length,length))
		self.bias_dwm = np.zeros((5,5,length,length))
		self.bg_dwm = np.zeros((5,5,length,length))


	@cython.boundscheck(False)
	@cython.wraparound(False) 	#dont deal with negative indices
	def add_sequence(self, np.ndarray[np.int64_t, ndim=1] sequence, double amount = 1.0):

		cdef np.ndarray[np.float64_t, ndim=4] bias_counts = self.counts
		cdef int L = self.length
		cdef int Sm, Sn, m, n

		#Check if seq contains N
		for m in range(L):
			Sm = sequence[m]	#ATCG(N)
			if Sm == 4:
				return

		#Add counts
		for m in range(L):
			Sm = sequence[m]	#ATCG(N)
			for n in range(L):
				Sn = sequence[n]
				bias_counts[Sm,Sn,m,n] += amount

		self.no_bias += amount
		self.counts = bias_counts


	@cython.boundscheck(False)
	@cython.wraparound(False)		
	def add_background(self, np.ndarray[np.int64_t, ndim=1] sequence, double amount = 1.0):
		""" Adds sequence to count of background nucleotides """

		cdef np.ndarray[np.float64_t, ndim=4] bg_counts = self.bg_counts
		cdef int L = self.length
		cdef int Sm, Sn, m, n
		
		#Check if seq contains N
		for m in range(L):
			Sm = sequence[m]	#ATCG(N)
			if Sm == 4:
				return

		for m in range(L):
			Sm = sequence[m]	#ATCG(N)
			for n in range(L):
				Sn = sequence[n]
				bg_counts[Sm,Sn,m,n] += amount

		self.no_bg += amount
		self.bg_counts = bg_counts


	def prepare_mat(self):

		L = self.length

		#PWM bias
		bias_cm = np.sum(self.counts, axis=(0,2)) / float(L)	    #divide by length because every position was added *length times
		self.no_bias = np.mean(np.sum(bias_cm[:4,:], axis=0))
		self.bias_pwm = np.true_divide(bias_cm + 1, self.no_bias + 4)

		#PWM bg
		bg_cm = np.sum(self.bg_counts, axis=(0,2)) / float(L)
		self.no_bg = np.mean(np.sum(bg_cm[:4,:], axis=0))
		self.bg_pwm = np.true_divide(bg_cm + 1, self.no_bg + 4)
		
		#Pssm for visualization
		self.pssm = np.log2(np.true_divide(self.bias_pwm, self.bg_pwm))	

		#Create dwms from counts
		pseudo_bias = max([16,self.no_bias])
		pseudo_bg = max([16,self.no_bg])
		for m in range(L):
			for n in range(L):
				for Sm in range(5):
					for Sn in range(5):
						self.bias_dwm[Sm,Sn,m,n] = (self.counts[Sm,Sn,m,n] + pseudo_bias * self.bias_pwm[Sn,n] * self.bias_pwm[Sm,m]) / float(self.no_bias + pseudo_bias)
						self.bg_dwm[Sm,Sn,m,n] = (self.bg_counts[Sm,Sn,m,n] + pseudo_bg * self.bg_pwm[Sn,n] * self.bg_pwm[Sm,m]) / float(self.no_bg + pseudo_bg)

		#Log for scoring
		self.bias_pwm_log = np.log2(self.bias_pwm)
		self.bg_pwm_log = np.log2(self.bg_pwm)
		self.bias_dwm_log = np.log2(self.bias_dwm)
		self.bg_dwm_log = np.log2(self.bg_dwm)


	@cython.boundscheck(False)
	@cython.cdivision(True)		#no check for zero division
	@cython.wraparound(False) 	#dont deal with negative indices
	def score_sequence(self, np.ndarray[np.int64_t, ndim=1] sequence):
		#Score nucleotide sequence against dinucleotide motif

		cdef np.ndarray[np.float64_t, ndim=2] bias_PWM = self.bias_pwm_log
		cdef np.ndarray[np.float64_t, ndim=2] bg_PWM = self.bg_pwm_log
		cdef np.ndarray[np.float64_t, ndim=4] bias_DWM = self.bias_dwm_log
		cdef np.ndarray[np.float64_t, ndim=4] bg_DWM = self.bg_dwm_log

		cdef int L = self.length
		cdef int F = int(L/2.0) 
		cdef int reg_len = len(sequence)
		cdef np.ndarray[np.float64_t, ndim=1] scores = np.zeros(reg_len)		#array to fill in
		scores[:] = np.nan
		cdef int i, m, n, a, Sm, Sn, contain_N
		cdef double bg_score_log, bias_score_log, Cn, prod_log, cond_log, probas

		#Index along sequence
		for i in range(reg_len - L + 1): #i is 1 if reglen == L
			
			contain_N = 0
			bg_score_log = 0
			bias_score_log = 0
			for n in range(L):
				Sn = sequence[i + n]	#For nucleotide Sn on position n in kmer

				#If N in kmer:
				if Sn > 3:
					contain_N = 1

				else:
					
					###Score against background motif###
					#Sum of all conditional probabilities of possible nucleotides on position n (Sn))
					Cn = 0.0
					for a in range(4):
						prod_log = 0.0
						for m in range(L):
							if m != n:
								Sm = sequence[i + m]
								prod_log += bg_DWM[Sm,a,m,n] - bg_PWM[a,n]
						
						Cn += 2**(prod_log + bg_PWM[a,n])  #back to linear to add up

					#Conditional probability for true Sn:
					cond_log = 0.0
					for m in range(L):
						if m != n:
							Sm = sequence[i + m]
							cond_log += bg_DWM[Sm,Sn,m,n] - bg_PWM[Sn,n]

					proba = cond_log + bg_PWM[Sn,n] - log2(Cn)
					bg_score_log += proba
					
					###Score against bias motif###
					#Sum of all conditional probabilities of possible nucleotides on position n (Sn))
					Cn = 0.0
					for a in range(4):
						prod_log = 0.0
						for m in range(L):
							if m != n:
								Sm = sequence[i + m]
								prod_log += bias_DWM[Sm,a,m,n] - bias_PWM[a,n]

						Cn += 2**(prod_log + bias_PWM[a,n])  #back to linear to add up
					
					#Conditional probability for true Sn:
					cond_log = 0.0
					for m in range(L):
						if m != n:
							Sm = sequence[i + m]
							cond_log += bias_DWM[Sm,Sn,m,n] - bias_PWM[Sn,n]
				
					proba = cond_log + bias_PWM[Sn,n] - log2(Cn)
					bias_score_log += proba

			#Score for index i
			if contain_N == 0:	
				scores[i+F] = bias_score_log - bg_score_log		#log-scale

		return(scores)

#--------------------------------------------------------------------------------------------------#
def nuc_to_num(str sequence):
	""" Convert DNA sequence string to internal number format """

	cdef int length = len(sequence)
	cdef np.ndarray[np.int64_t, ndim=1] num_sequence = np.zeros(length, dtype=int)
	cdef int i, num
	cdef str nuc

	for i in range(length):
		nuc = sequence[i]

		#number format
		if nuc == "A" or nuc == "a":
			num = 0
		elif nuc == "T" or nuc == "t":
			num = 1
		elif nuc == "C" or nuc == "c":
			num = 2
		elif nuc == "G" or nuc == "g":
			num = 3
		else:
			num = 4

		num_sequence[i] = num

	return(num_sequence)

#--------------------------------------------------------------------------------------------------#
class GenomicSequence:
	""" Used in ATACorrect for fast scoring of PWM/DWM """

	def __init__(self, region):

		self.region = region
		self.length = region.end - region.start
		self.sequence = np.zeros(self.length, dtype=int) #in number format
		self.revcomp = np.zeros(self.length, dtype=int)	 #in number format

	@cython.boundscheck(False)	#dont check boundaries
	@cython.wraparound(False) 	#dont deal with negative indices
	def from_fasta(self, fasta_obj):
		""" Fill sequence from fasta object """

		cdef str fasta = fasta_obj.fetch(self.region.chrom, self.region.start, self.region.end)
		cdef int length = self.length
		cdef np.ndarray[np.int64_t, ndim=1] sequence = self.sequence
		cdef np.ndarray[np.int64_t, ndim=1] revcomp_sequence = self.revcomp
		cdef int i, num, comp
		cdef str nuc

		for i in range(length):
			nuc = fasta[i]

			#number format
			if nuc == "A" or nuc == "a":
				num = 0
			elif nuc == "T" or nuc == "t":
				num = 1
			elif nuc == "C" or nuc == "c":
				num = 2
			elif nuc == "G" or nuc == "g":
				num = 3
			else:
				num = 4

			#complement
			if num == 0:
				comp = 1
			elif num == 1:
				comp = 0
			elif num == 2:
				comp = 3
			elif num == 3:
				comp = 2
			else:
				comp = 4

			sequence[i] = num
			revcomp_sequence[length-i-1] = comp

		self.sequence = sequence
		self.revcomp = revcomp_sequence

		return(self)


def get_gc_content(regions, fasta):
	""" Get GC content from regions in fasta """
	gc = 0
	total = 0
	fasta_obj = pysam.FastaFile(fasta)
	for region in regions:
		seq = fasta_obj.fetch(region.chrom, region.start, region.end)
		gc += sum([1 for nuc in seq if (nuc == "G" or nuc == "g" or nuc == "C" or nuc == "c")])
		total += region.end - region.start
	fasta_obj.close()
	gc_content = gc / float(total)

	return(gc_content)
