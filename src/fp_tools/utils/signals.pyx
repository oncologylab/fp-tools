# cython: language_level=3

"""Signal-array utilities for bias, cutsite, and footprint-score workflows."""

import numpy as np
cimport numpy as np
import math
import cython
from libc.math cimport fabs, isfinite

#--------------------------------------------------------------------------------------------------#
class OneSignal(np.ndarray):
	""" Work in progress; placeholder for future development """

#--------------------------------------------------------------------------------------------------#
class SignalList(list):
	""" Work in progress; placeholder for future development """

	def __new__(cls):
		pass

	def __init__(self, matrix=None, name=""):

		self.aggregate = ""
		self.name = name
		self.mat = matrix
		self.n = 0 #n regions

	def from_regions(self, regions, bigwig):
		""" 
			Read from regions and bigwig 
			Assumes .signal in regions
		"""
		pass
	
	def filter_outliers(self, lower=0.0, upper=1.0):
		""" Filter rows based on outlier values """

		#Exclude outlier rows 
		max_values = np.max(self.mat, axis=1)
		upper_limit = np.percentile(max_values, [100*upper])[0]	#remove-outliers is a fraction
		logical = max_values <= upper_limit 
		#logger.debug("{0}:\tUpper limit: {1} (regions removed: {2})".format(self.name, upper_limit, self.n - sum(logical)))
		#signalmat = signalmat[logical]

	def aggregate(self, normalize=False, smooth=1):
		""" Makes aggregate across all rows """

		self.aggregate = ""
		return(self.aggregate)

	def correlate():
		""" """
		pass

	def footprint():
		""" """
		pass

#--------------------------------------------------------------------------------------------------#
def shuffle_array(np.ndarray[np.float64_t, ndim=1] arr, int no_rand, np.ndarray[np.int64_t, ndim=1] shift_options):
	""" Shuffles array of values within the boundaries given in shift """

	cdef int max_shift = max([abs(np.min(shift_options)), abs(np.max(shift_options))])
	cdef np.ndarray[np.float64_t, ndim=1] ext_arr = np.concatenate((np.zeros(max_shift), arr, np.zeros(max_shift)))	 #pad with max shift to allow for shuffling outside borders
	cdef int ext_arr_len = len(ext_arr)

	cdef np.ndarray[np.int64_t, ndim=1] nonzero_index = np.nonzero(ext_arr)[0]
	cdef int no_shift = len(nonzero_index)

	cdef np.ndarray[np.int64_t, ndim=2] rand_rel_positions = np.random.choice(shift_options, size=(no_shift, no_rand)) 		#positions of shuffled reads
	cdef np.ndarray[np.float64_t, ndim=2] rand_mat = np.zeros((no_rand, ext_arr_len))

	cdef double value
	cdef int i, j, pos, rand_pos

	#Get all relative placements
	for i in range(no_shift):

		pos = nonzero_index[i] 	#original index in ext_arr
		value = ext_arr[pos]	#value at original position

		#Get new random positions
		for j in range(no_rand):
			rand_pos = pos + rand_rel_positions[i,j]  #i for cut, j for each randomization
			rand_mat[j, rand_pos] = rand_mat[j, rand_pos] + value

	return(rand_mat[:,max_shift:-max_shift])



#--------------------------------------------------------------------------------------------------#
@cython.boundscheck(False)	#dont check boundaries
@cython.cdivision(True)		#no check for zero division
@cython.wraparound(False) 	#dont deal with negative indices
def fast_rolling_math(np.ndarray[np.float64_t, ndim=1] arr, int w, operation):

	"""
	Rolling operation of arr with window size w 
	Possible operations are: "max", "min", "mean", "sum"
	Returns array the same size as arr with "NaN" in flanking positions (for sum/mean)
	"""

	cdef int L = arr.shape[0]
	cdef np.ndarray[np.float64_t, ndim=1] roll_arr = np.zeros(L)
	roll_arr[:] = np.nan
	cdef int i, j, start_i
	cdef int lf = int(np.floor(w / 2.0))
	cdef int rf = int(np.ceil(w / 2.0))
	cdef float minval, maxval, valsum
	cdef np.ndarray[np.float64_t, ndim=1] prefix
	cdef np.ndarray[np.int64_t, ndim=1] nan_prefix
	cdef double value
	
	#Max in window
	if operation == "max":
		for i in range(L):		#centered on reg

			maxval = arr[i]		#initialize with first value
			start_i = i - lf 	#start_i + w gives full region
		
			for j in range(w):
				if start_i + j > 0 and start_i + j < L:
					if arr[start_i+j] > maxval:
						maxval = arr[start_i+j]
			
			#Assign maxval to index
			roll_arr[i] = maxval
	
	#Min in window
	elif operation == "min":
		for i in range(L):

			minval = arr[i]
			start_i = i - lf 	

			for j in range(w):
				if start_i + j > 0 and start_i + j < L:
					if arr[start_i+j] < minval:
						minval = arr[start_i+j]

			#Assign maxval to index
			roll_arr[i] = minval

	#Sum of window
	elif operation == "sum":
		prefix = np.zeros(L + 1)
		nan_prefix = np.zeros(L + 1, dtype=np.int64)
		for i in range(L):
			value = arr[i]
			prefix[i + 1] = prefix[i]
			nan_prefix[i + 1] = nan_prefix[i]
			if value != value:
				nan_prefix[i + 1] += 1
			else:
				prefix[i + 1] += value
		for i in range(L-w+1):
			if nan_prefix[i + w] - nan_prefix[i] > 0:
				roll_arr[i+lf] = np.nan
			else:
				roll_arr[i+lf] = prefix[i + w] - prefix[i]

	#Mean in window
	elif operation == "mean":
		prefix = np.zeros(L + 1)
		nan_prefix = np.zeros(L + 1, dtype=np.int64)
		for i in range(L):
			value = arr[i]
			prefix[i + 1] = prefix[i]
			nan_prefix[i + 1] = nan_prefix[i]
			if value != value:
				nan_prefix[i + 1] += 1
			else:
				prefix[i + 1] += value
		for i in range(L-w+1):
			if nan_prefix[i + w] - nan_prefix[i] > 0:
				roll_arr[i+lf] = np.nan
			else:
				roll_arr[i+lf] = (prefix[i + w] - prefix[i])/(w*1.0)

	#Product of values in window
	elif operation == "prod":
		for i in range(L):
			prod = 1.0
			start_i = i - lf
			for j in range(w):
				if start_i + j > 0 and start_i + j < L:
					prod *= arr[start_i+j]

			roll_arr[i] = prod

	return roll_arr


#--------------------------------------------------------------------------------------------------#
@cython.boundscheck(False)	#dont check boundaries
@cython.cdivision(True)		#no check for zero division
@cython.wraparound(False) 	#dont deal with negative indices
def footprint_score_array(np.ndarray[np.float64_t, ndim=1] arr, int flank_min, int flank_max, int fp_min, int fp_max):

	cdef int L = arr.shape[0]
	cdef np.ndarray[np.float64_t, ndim=1] footprint_scores = np.zeros(L)
	cdef int i, j, footprint_w, flank_w

	cdef float fp_sum, fp_mean, fp_score, flank_mean
	cdef float left_sum, right_sum, val

	#Each position in array, starting at i, going until last possible start of region
	for i in range(L-2*flank_max-fp_max):
		for flank_w in range(flank_min, flank_max+1):
			for footprint_w in range(fp_min, fp_max+1):

				#Sum of left flank
				left_sum = 0.0
				for j in range(flank_w):
					val = arr[i+j]
					if val > 0.0:
						left_sum += val

				#Sum of footprint (only negative counts)
				fp_sum = 0.0
				for j in range(footprint_w):
					val = arr[i+flank_w+j]	
					if val < 0.0:
						fp_sum += val 	#val is negative

				#Sum of right flank
				right_sum = 0.0
				for j in range(flank_w):
					val = arr[i+flank_w+footprint_w+j]
					if val > 0.0:
						right_sum += val

				#Calculate score
				fp_mean = fp_sum/(1.0*footprint_w)	#will be minus
				flank_mean = (right_sum + left_sum)/(2.0*flank_w)

				fp_score = flank_mean - fp_mean 	#-- = +

				#Save score across footprint
				for pos in range(i+flank_w, i+flank_w+footprint_w):
					if fp_score > footprint_scores[pos]:
						footprint_scores[pos] = fp_score

	return(footprint_scores)


#--------------------------------------------------------------------------------------------------#
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def footprint_score_array_fast(np.ndarray[np.float64_t, ndim=1] arr, int flank_min, int flank_max, int fp_min, int fp_max):

	cdef int L = arr.shape[0]
	cdef np.ndarray[np.float64_t, ndim=1] footprint_scores = np.zeros(L)
	cdef np.ndarray[np.float64_t, ndim=1] pos_prefix = np.zeros(L + 1)
	cdef np.ndarray[np.float64_t, ndim=1] neg_prefix = np.zeros(L + 1)
	cdef int i, pos, footprint_w, flank_w
	cdef double val, fp_sum, fp_mean, fp_score, flank_mean
	cdef double left_sum, right_sum

	for i in range(L):
		val = arr[i]
		pos_prefix[i + 1] = pos_prefix[i]
		neg_prefix[i + 1] = neg_prefix[i]
		if val > 0.0:
			pos_prefix[i + 1] += val
		elif val < 0.0:
			neg_prefix[i + 1] += val

	# Same window traversal and update semantics as footprint_score_array,
	# but left/center/right sums are O(1) prefix differences.
	for i in range(L - 2 * flank_max - fp_max):
		for flank_w in range(flank_min, flank_max + 1):
			left_sum = pos_prefix[i + flank_w] - pos_prefix[i]
			for footprint_w in range(fp_min, fp_max + 1):
				fp_sum = neg_prefix[i + flank_w + footprint_w] - neg_prefix[i + flank_w]
				right_sum = (
					pos_prefix[i + flank_w + footprint_w + flank_w]
					- pos_prefix[i + flank_w + footprint_w]
				)

				fp_mean = fp_sum / (1.0 * footprint_w)
				flank_mean = (right_sum + left_sum) / (2.0 * flank_w)
				fp_score = flank_mean - fp_mean

				for pos in range(i + flank_w, i + flank_w + footprint_w):
					if fp_score > footprint_scores[pos]:
						footprint_scores[pos] = fp_score

	return(footprint_scores)


#--------------------------------------------------------------------------------------------------#
@cython.cdivision(True)		#no check for zero division
@cython.boundscheck(False)	#dont check boundaries
@cython.wraparound(False) 	#dont deal with negative indices
def FOS_score(np.ndarray[np.float64_t, ndim=1] arr, int flank_min, int flank_max, int fp_min, int fp_max):
	
	cdef int L = arr.shape[0]
	cdef np.ndarray[np.float64_t, ndim=1] footprint_scores = np.zeros(L)
	footprint_scores[:] = 10000
	cdef int i, j, footprint_w, flank_w

	cdef float Lm, Cm, Rm, fos_score
	cdef float left_sum, right_sum, fp_sum

	#Each position in array, starting at i, going until last possible start of region
	for i in range(L-2*flank_max-fp_max):
		for flank_w in range(flank_min, flank_max):
			for footprint_w in range(fp_min, fp_max):
			
				#Sum of left flank
				left_sum = 0.0
				for j in range(flank_w):
					left_sum += arr[i+j]

				#Sum of footprint
				fp_sum = 0.0
				for j in range(footprint_w):
					fp_sum += arr[i+flank_w+j]

				#Sum of right flank
				right_sum = 0.0
				for j in range(flank_w):
					right_sum += arr[i+flank_w+footprint_w+j]

				#Calculate score
				Lm = left_sum / (flank_w*1.0)		#left mean
				Cm = fp_sum / (footprint_w*1.0)		#center mean
				Rm = right_sum / (flank_w*1.0)		#right mean

				if Cm < Rm and Cm < Lm and Lm > 0 and Rm > 0:
					fos_score = (Cm+1.0)/Lm + (Cm+1.0)/Rm
				else:
					fos_score = 10000

				#Save score to arr (smallest are best)
				for j in range(footprint_w):
					if fos_score < footprint_scores[i+flank_w+j]:
						footprint_scores[i+flank_w+j] = fos_score

	return(footprint_scores)


#--------------------------------------------------------------------------------------------------#
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def add_bias_prediction_window(
		np.ndarray[np.float64_t, ndim=1] sum_arr,
		np.ndarray[np.int64_t, ndim=1] count_arr,
		np.ndarray[np.float64_t, ndim=1] prediction,
		int start,
		int end):
	"""Accumulate one bias-prediction window for later nanmean-equivalent finalization."""

	cdef int i, j
	cdef int L = sum_arr.shape[0]
	cdef double value
	if start < 0:
		start = 0
	if end > L:
		end = L
	for i in range(start, end):
		j = i - start
		value = prediction[j]
		if value == value:
			sum_arr[i] += value
			count_arr[i] += 1
	return None


#--------------------------------------------------------------------------------------------------#
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def finalize_bias_prediction(
		np.ndarray[np.float64_t, ndim=1] sum_arr,
		np.ndarray[np.int64_t, ndim=1] count_arr,
		int k_flank,
		int reg_end):
	"""Finalize accumulated bias predictions with the previous nanmean semantics."""

	cdef int i
	cdef int L = sum_arr.shape[0]
	cdef np.ndarray[np.float64_t, ndim=1] out = np.zeros(L)
	for i in range(L):
		if count_arr[i] > 0:
			out[i] = sum_arr[i] / count_arr[i]
		elif i >= k_flank and i < reg_end:
			out[i] = np.nan
		else:
			out[i] = 0.0
	return out


cdef inline void _rolling_sum_nan_to_zero(
		double[:] arr,
		int w,
		double[:] out) noexcept:
	cdef Py_ssize_t L = arr.shape[0]
	cdef Py_ssize_t i
	cdef int lf = int(math.floor(w / 2.0))
	cdef double running = 0.0
	cdef long nan_count = 0
	cdef double value

	for i in range(L):
		out[i] = 0.0
	if w <= 0 or L <= 0 or w > L:
		return
	for i in range(w):
		value = arr[i]
		if value != value:
			nan_count += 1
		else:
			running += value
	out[lf] = 0.0 if nan_count > 0 else running
	for i in range(1, L - w + 1):
		value = arr[i - 1]
		if value != value:
			nan_count -= 1
		else:
			running -= value
		value = arr[i + w - 1]
		if value != value:
			nan_count += 1
		else:
			running += value
		out[i + lf] = 0.0 if nan_count > 0 else running


#--------------------------------------------------------------------------------------------------#
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def atac_correct_arrays(
		np.ndarray[np.float64_t, ndim=1] uncorrected,
		np.ndarray[np.float64_t, ndim=1] bias,
		int window,
		double correction_factor):
	"""Apply ATAC expected-signal correction and positive-signal rescaling."""

	cdef Py_ssize_t L = uncorrected.shape[0]
	cdef Py_ssize_t i
	cdef double bsum, cpos_sum, scale
	cdef np.ndarray[np.float64_t, ndim=1] uncorrected_scaled = np.empty(L, dtype=np.float64)
	cdef np.ndarray[np.float64_t, ndim=1] expected = np.empty(L, dtype=np.float64)
	cdef np.ndarray[np.float64_t, ndim=1] corrected = np.empty(L, dtype=np.float64)
	cdef np.ndarray[np.float64_t, ndim=1] signal_sum = np.empty(L, dtype=np.float64)
	cdef np.ndarray[np.float64_t, ndim=1] bias_sum = np.empty(L, dtype=np.float64)
	cdef np.ndarray[np.float64_t, ndim=1] corrected_abs = np.empty(L, dtype=np.float64)
	cdef np.ndarray[np.float64_t, ndim=1] corrected_pos = np.empty(L, dtype=np.float64)
	cdef np.ndarray[np.float64_t, ndim=1] uncorrected_sum = np.empty(L, dtype=np.float64)
	cdef np.ndarray[np.float64_t, ndim=1] corrected_sum = np.empty(L, dtype=np.float64)
	cdef np.ndarray[np.float64_t, ndim=1] corrected_pos_sum = np.empty(L, dtype=np.float64)

	if bias.shape[0] != L:
		raise ValueError("uncorrected and bias arrays must have the same length")

	_rolling_sum_nan_to_zero(uncorrected, window, signal_sum)
	_rolling_sum_nan_to_zero(bias, window, bias_sum)

	for i in range(L):
		bsum = bias_sum[i]
		if bsum != bsum or fabs(bsum) <= 1e-8:
			expected[i] = 0.0
		else:
			expected[i] = signal_sum[i] * (bias[i] / bsum)
		uncorrected_scaled[i] = uncorrected[i] * correction_factor
		expected[i] *= correction_factor
		corrected[i] = uncorrected_scaled[i] - expected[i]
		if corrected[i] != corrected[i]:
			corrected_abs[i] = np.nan
			corrected_pos[i] = np.nan
		else:
			corrected_abs[i] = fabs(corrected[i])
			corrected_pos[i] = corrected[i] if corrected[i] > 0.0 else 0.0

	_rolling_sum_nan_to_zero(uncorrected_scaled, window, uncorrected_sum)
	_rolling_sum_nan_to_zero(corrected_abs, window, corrected_sum)
	_rolling_sum_nan_to_zero(corrected_pos, window, corrected_pos_sum)

	for i in range(L):
		cpos_sum = corrected_pos_sum[i]
		if cpos_sum == 0.0:
			scale = 1.0
		else:
			scale = (uncorrected_sum[i] - (corrected_sum[i] - cpos_sum)) / cpos_sum
			if scale != scale or scale < 1.0:
				scale = 1.0
		if corrected[i] > 0.0:
			corrected[i] *= scale

	return uncorrected_scaled, expected, corrected


#--------------------------------------------------------------------------------------------------#
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def local_maxima_indices(np.ndarray[np.float64_t, ndim=1] values):
	"""Return local-maximum offsets using the same plateau semantics as call-footprints."""

	cdef Py_ssize_t L = values.shape[0]
	cdef Py_ssize_t i
	cdef double current, left, right
	cdef list out = []
	if L == 0:
		return out
	if L == 1:
		current = values[0]
		if current == current and isfinite(current):
			out.append(0)
		return out
	for i in range(L):
		current = values[i]
		if current != current or not isfinite(current):
			continue
		left = -np.inf if i == 0 else values[i - 1]
		right = -np.inf if i == L - 1 else values[i + 1]
		if current >= left and current >= right and (current > left or current > right):
			out.append(i)
	return out
