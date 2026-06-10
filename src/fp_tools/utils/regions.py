#!/usr/bin/env python

"""Classes for working with genomic regions."""

import numpy as np
import sys
import re
from copy import deepcopy
import pyBigWig
from collections import Counter
import logging
import traceback

from fp_tools.utils.logger import FpToolsLogger

#Clustering
import sklearn.preprocessing as preprocessing
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform

#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

class OneRegion(list):

	nuc_to_pos = {"A":0, "T":1, "C":2, "G":3}

	def __init__(self, lst=["",0,0]):

		super(OneRegion, self).__init__(iter(lst))
		no_fields = len(lst)

		#Required
		self.chrom = lst[0]
		self.start = int(lst[1])		#exclude start
		self.end = int(lst[2])			#include end

		#Optional
		self.name = lst[3] if no_fields > 3 else ""
		self.score = lst[4] if no_fields > 4 else ""
		self.strand = lst[5] if no_fields > 5 else "."

	def __str__(self):
		return("{0}".format("\t".join(str(x) for x in self)))
	
	def tup(self):
		return((self.chrom, self.start, self.end, self.name, self.strand))
	
	def pretty(self):
		""" Pretty print of coordinates in the format chr:start-stop(strand) """
		return("{0}:{1}-{2}({3})".format(self.chrom, self.start, self.end, self.strand))

	def update(self):

		self[0] = self.chrom
		self[1] = self.start
		self[2] = self.end

	def get_length(self):
		return(self.end - self.start)

	def get_width(self):
		return(self.end - self.start)

	def extend_reg(self, bp):
		""" Extend region with bp to either side """

		self.start -= bp
		self[1] = self.start

		self.end += bp
		self[2] = self.end

		return(self)

	def set_width(self, bp):
		""" Set with of region centered on original region """

		flank_5, flank_3 = int(np.floor(bp/2.0)), int(np.ceil(bp/2.0))  #flank 5', flank 3'

		if self.strand == "-":
			mid = int(np.ceil((self.start + self.end) / 2.0))
			self.start = mid - flank_5
			self.end = mid + flank_3
		else:
			mid = int(np.floor((self.start + self.end) / 2.0))
			self.start = mid - flank_3
			self.end = mid + flank_5

		self[1] = self.start
		self[2] = self.end

		return(self)


	def split_region(self, max_length):
		""" Split genomic region into smaller subsets. Returns a RegionList of OneRegion objects """

		regions = RegionList() 	#Empty regions object

		starts = list(range(self.start, self.end, max_length))
		ends = list(range(self.start + max_length, self.end, max_length)) + [self.end]

		for (start, end) in zip(starts, ends):
			regions.append(OneRegion([self.chrom, start, end]))

		return(regions)


	def check_boundary(self, boundaries_dict, action="cut", logger=FpToolsLogger()):
		""" Check if region is within chromosome boundaries. Actions:
				- "cut": cut region to bounds. If the chromosome is not in boundaries_dict, "cut" falls back on "remove"
				- "remove": remove region outside bounds (returns None)
				- "exit": exit the program with error message through logger
		"""
		
		#Establish if region is outside of bounds:
		outside = 0
		if self.chrom not in boundaries_dict:
			if action == "exit":
				logger.error("Chromosome for region \"{0}\" is not found in list of available chromosomes ({1})".format(self, list(boundaries_dict.keys())))
				sys.exit(1)

			self = None	#cannot cut to bounds when boundaries are not known; remove
			return(self)

		elif self.start < 0:
			outside = 1
		elif self.end > int(boundaries_dict[self.chrom]):
			outside = 1

		#Perform action if region is outside of bounds
		if outside == 1:
			if action == "cut":
				self.start = max([0, self.start])
				self.end = min([boundaries_dict[self.chrom], self.end])
				
				#Update positions in list
				self[1] = self.start
				self[2] = self.end

				#If the region has been cut to be 0 of less length; remove
				if self.get_length() <= 0:
					self = None

			elif action == "remove":
				self = None
				
			elif action == "exit":
				logger.error("Region \"{0}\" is outside of the chromosome boundaries ({1}: {2})".format(self, self.chrom, boundaries_dict[self.chrom]))
				sys.exit(1)

		return(self)

	def get_signal(self, pybw, numpy_bool = True, logger=FpToolsLogger(), key=None):
		""" Get signal from bigwig in region. key is a string which will be written in case of an error """

		try:
			#Define whether pybigwig was compiled with numpy
			if pyBigWig.numpy == 1:
				values = pybw.values(self.chrom, self.start, self.end, numpy=numpy_bool)
			else:
				values = np.array(pybw.values(self.chrom, self.start, self.end)) #fetch list of values and convert to numpy arr
			values = np.nan_to_num(values)	#nan to 0
			
			if self.strand == "-":
				signal = values[::-1]
			else:
				signal = values
				
		except Exception as e:
			if key is not None:
				logger.error("Error reading region: {0} from pybigwig object ({1}). Exception is: {2}".format(self.tup(), key, e))
			else:
				logger.error("Error reading region: {0} from pybigwig object. Exception is: {1}".format(self.tup(), e))

			traceback.print_tb(e.__traceback__)
			raise e
			
		return(signal)	


#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#

class RegionList(list):
	""" Class Regions is a list of OneRegion objects """

	def __str__(self):
		return("\n".join([str(oneregion) for oneregion in self]))
	
	def __init__(self, lst=[]):
		super(RegionList, self).__init__(iter(lst))

	# ---- I/O ---- #
	def from_list(self, lst):
		""" Initialize object from list of OneRegion objects """

		for obj in lst:
			self.append(obj)
		return(self)

	def from_bed(self, bedfile_f, logger=FpToolsLogger()):
		""" Initialize Object from bedfile """

		#Read all lines
		try:
			bedlines = open(bedfile_f).readlines()
		except UnicodeDecodeError as e:
			logger.error("Could not read input file '{0}'. It looks like the file might be compressed. fp-tools reads .bed files in raw text format. Exception is: '{1}'".format(bedfile_f, e))
			sys.exit(1)
		except Exception as e:
			logger.error("Could not read input file '{0}'. Exception is: '{1}'".format(bedfile_f, e))
			sys.exit(1)

		self = RegionList([None]*len(bedlines))	#intialize to prevent appending for large bedfiles
		for i, line in enumerate(bedlines):
		
			if line.startswith("#"): #comment lines are excluded
				continue

			#Test line format
			if re.match(r"[^\s]+\t\d+\t\d+\b.*", line) == None:
				logger.error("Line {0} in {1} is not proper bed format:\n{2}".format(i+1, bedfile_f, line))
				sys.exit(1)

			columns = line.split("\t") #split (and then .rstrip() afterwards)
			columns[-1] = columns[-1].rstrip() # remove line-ending from last col
			columns[1] = int(columns[1]) #start
			columns[2] = int(columns[2]) #end
			
			if columns[1] >= columns[2]:
				logger.error("Start position is larger than end position in line {0} in {1}:\n{2}".format(i+1, bedfile_f, line))
				sys.exit(1)

			region = OneRegion(columns)
			self[i] = region
		
		#Remove any None's in list
		for i in range(len(self)-1, -1, -1):	#counting down so that index in list does not change as None's are removed
			if self[i] is None:
				del self[i]

		return(self)

	def from_chrom_lengths(self, chrom_length_dict):
		""" Initialize regions with dict of chromosome lengths """

		for chrom in sorted(chrom_length_dict):
			self.append(OneRegion([chrom, 0, chrom_length_dict[chrom]]))

		return(self)

	def as_bed(self, additional=True):
		""" Return bedfile formatting as string """

		bed = ""
		for region in self:

			line = "{0}\n".format("\t".join([str(col) for col in region]))
			bed += line

		return(bed)

	def write_bed(self, bed_f):
		""" Write regions to bedfile """

		bed_content = self.as_bed()
		out = open(bed_f, "w")
		out.write(bed_content)
		out.close()


	# ---- Working with regions ---- #
	def count(self):
		""" Number of regions """
		return(len(self))

	def loc_sort(self, contig_list=[]):
		""" 
		Sorts list of region objects based on genomic location. 
		contig_list gives the correct order of contigs; else normal sort is applied
		"""

		#Get all chroms of list; join with contig_list to prevent key error during sort
		#all_chroms = self.get_chroms()

		#Sort
		if len(contig_list) > 0:
			order_dict = dict(zip(contig_list, range(len(contig_list))))
			self.sort(key=lambda region: (order_dict.get(region.chrom, -1), region.start, region.end, region.name))
		else:
			self.sort(key=lambda region: (region.chrom, region.start, region.end, region.name))

	def score_sort(self):
		""" Sorts list of region objects based on score in last column """
		self.sort(key=lambda region: float(region[-1]), reverse=True)

	def get_chroms(self):

		chroms = {}
		for region in self:
			chroms[region.chrom] = ""
			
		return(list(chroms.keys()))


	def split_chrom(self):
		""" Split list of regions into separate chromosomes (output is RegionLists in dict) """

		output_dict = {}
		for region in self:
			if region.chrom not in output_dict:
				output_dict[region.chrom] = RegionList()

			output_dict[region.chrom].append(region)

		return(output_dict)


	def chunks(self, n):
		""" Split list of regions into n chunks (e.g. to use in multiprocessing) """

		no_reg = self.count()

		if no_reg > 0:
			per_chunk = int(np.ceil(no_reg/float(n)))
			chunk_lst = [RegionList(self[i:i+per_chunk]) for i in range(0, no_reg, per_chunk)]
		else:
			chunk_lst = [self]
		return(chunk_lst)


	def merge(self, name=False): 	
		""" Merge overlapping genomic regions. If name == True, regions will only be merged if the name is the same (used in count_overlaps) """

		self.loc_sort()		#sort before overlapping
		no = self.count()	#number of regions

		prev_chrom, prev_start, prev_end = "", 0, 0
		i = 1
		while i < no:

			prev_chrom, prev_start, prev_end = self[i-1].chrom, self[i-1].start, self[i-1].end
			curr_chrom, curr_start, curr_end = self[i].chrom, self[i].start, self[i].end

			#Check naming
			merge_flag = True
			if name == True:
				merge_flag = True if self[i-1].name == self[i].name else False

			if curr_chrom == prev_chrom and curr_start < prev_end and merge_flag == True:		#if overlapping regions
					self[i].start = prev_start
					self[i][1] = prev_start

					self[i].end = curr_end
					self[i][2] = curr_end

					del self[i-1]
					no -= 1
			else:
				i += 1  #compare next

		return(self)

	def remove_chroms(self, chromlist):
		""" Removes regions within chromlist chromosomes """

		no = self.count()
		i = 0
		while i < no:
			if self[i].chrom in chromlist:
				del self[i]
				no -= 1
			else:
				i += 1
		
		return(self)

	def keep_chroms(self, chromlist):
		""" Keep only regions within chromlist chromosomes """
		no = self.count()
		i = 0
		while i < no:
			if self[i].chrom not in chromlist:
				del self[i]
				no -= 1
			else:
				i += 1
		
		return(self)

	def subset(self, no):
		""" Take no-size subset of regions from RegionList """

		if self.count() > no:
			del self[no:]
		return(self)


	def remove_duplicates(self):
		""" Returns a unique list of sites """

		unique = deepcopy(self)
		unique.loc_sort()
		prev_chrom, prev_start, prev_end, prev_strand = "", 0, 0, ""

		for i, region in enumerate(unique):

			curr_chrom, curr_start, curr_end, curr_strand = region.chrom, region.start, region.end, region.strand

			if curr_chrom == prev_chrom and curr_start == prev_start and curr_end == prev_end and curr_strand == prev_strand:  #not unique
				unique[i] = None

			#Save as previous for next comparison
			prev_chrom, prev_start, prev_end, prev_strand = curr_chrom, curr_start, curr_end, curr_strand

		#Remove all None
		unique = RegionList([region for region in unique if region is not None])

		return(unique)


	def subtract(self, b_regions):
		""" Subtract b_regions from self regions """

		#Sort before starting
		self.loc_sort()
		b_regions.loc_sort()

		#a_len = self.count()
		#b_len = b_regions.count()

		chromlist = sorted(list(set([region.chrom for region in self] + [region.chrom for region in b_regions])))
		chrom_pos = {chrom:chromlist.index(chrom) for chrom in chromlist}

		a = self
		#b = b_regions

		a_i = 0
		b_i = 0
		while a_i < self.count() and b_i < b_regions.count():

			a_chrom, a_start, a_end = a[a_i].chrom, a[a_i].start, a[a_i].end
			b_chrom, b_start, b_end = b_regions[b_i].chrom, b_regions[b_i].start, b_regions[b_i].end

			if a_chrom == b_chrom:

				if a_end <= b_start:	#current a is placed before current b
					a_i += 1	

				elif a_start >= b_end:	#current a is placed after current b 
					b_i += 1

				elif a_start >= b_start and a_end <= b_end: 	#a is completely contained within b (a is removed)
					del a[a_i]

				elif a_start < b_start and a_end > b_end: 		#b is completely contained within a (a is split into two)
					a[a_i] = OneRegion([a_chrom, a_start, b_start])				#first new region
					a.insert(a_i + 1, OneRegion([a_chrom, b_end, a_end]))			#second new region

				elif a_end > b_end: 		#delete left side
					a[a_i] = OneRegion([a_chrom, b_end, a_end])

				elif a_start < b_start: 	#delete right side
					a[a_i] = OneRegion([a_chrom, a_start, b_start])
			
			elif chrom_pos[a_chrom] > chrom_pos[b_chrom]: 	#if a_chrom is after current b_chrom
				b_i += 1

			elif chrom_pos[b_chrom] > chrom_pos[a_chrom]:	#if b_chrom is after current a_chrom
				a_i += 1

		return(self)

	def apply_method(self, method, *args):
		""" Applies a method to every OneRegion object in regions list """ 

		reglist = RegionList()
		for i in range(self.count()):
			out = method(self[i], *args) 	#apply method to region

			if out == None:			#no return, change OneRegion object itself
				pass
			else:
				if type(out) is OneRegion:
					reglist.append(out)
				elif type(out) is RegionList:
					reglist.extend(out)

		return(reglist)


	def resolve_overlaps(self, priority="higher"):
		""" Priority "higher" or "lower" """

		self = deepcopy(self) #ensures that self is not changed
		self.loc_sort()
		no_regions = len(self)

		i = 0
		j = 1

		while i + j < no_regions:

			reg_a = self[i]
			reg_b = self[i+j]

			if reg_a == None:
				i += 1
				j = 1
			elif reg_b == None:
				j += 1
			else:
				if reg_a.chrom == reg_b.chrom:
					if reg_b.start < reg_a.end:
						scores = [reg_a.score, reg_b.score]

						if priority == "higher":
							worst = scores.index(min(scores))
						else:
							worst = scores.index(max(scores))

						if worst == 0:
							self[i] = None
						else:
							self[i+j] = None
					else:
						i += 1
						j = 1
				else:
					i += 1
					j = 1

		#Remove all None
		non_overlapping = RegionList([reg for reg in self if reg is not None])

		return(non_overlapping)


	def get_signal_dict(self, bigwigs):
		""" Get dict of signal[region.tup][bigwig] = signal """

		signal_dict = {region.tup():{bigwig:[] for bigwig in bigwigs} for region in self}
		for bigwig in bigwigs:
			pybw = pyBigWig.open(bigwig, "rb")
			for region in self:
				signal_dict[region.tup()][bigwig] = region.get_signal(pybw)
			pybw.close()

		return(signal_dict)


	def get_width_distri(self):

		c = Counter()
		for region in self:
			reglen = region.get_width()
			c[reglen] += 1
		return(c)


	
	def count_overlaps(self):
		""" Returns a dictionary of strings and tuples - string-keys represent total bp per TF, and tuples represent the overlap in bp between two TFs """

		#Join all with similar name first
		self.loc_sort()
		self = self.merge(name=True) 

		#Count overlap
		i = 0  #index of current site
		j = 1  #index of next/compared site

		no_sites = len(self)
		
		overlap = {}

		while i < no_sites:

			s1_chrom, s1_start, s1_end, s1_id = self[i].chrom, self[i].start, self[i].end, self[i].name

			reglen = s1_end - s1_start
			overlap[s1_id] = overlap.get(s1_id,0) + reglen

			#Compare to next sites
			flag = 1	#overlapping
			j = 1
			while flag == 1 and i+j < no_sites:

				s2_chrom, s2_start, s2_end, s2_id = self[i+j][:4]

				#Do the regions overlap?
				if s1_chrom == s2_chrom:
					if s1_start < s2_end and s1_end > s2_start+1:	#if overlap
						bp_overlap = min([s1_end, s2_end]) - max([s1_start, s2_start]) 
						overlap[(s1_id, s2_id)] = overlap.get((s1_id, s2_id), 0) + bp_overlap
						overlap[(s2_id, s1_id)] = overlap.get((s2_id, s1_id), 0) + bp_overlap
						j += 1
					else:
						flag = 0
				else:
					flag = 0
			i += 1

		return(overlap)



#--------------------------------------------------------------------------------------------------#
#------------------------------ Additional functions related to regions ---------------------------#
#--------------------------------------------------------------------------------------------------#

class RegionCluster:

	def __init__(self, overlap_dict): 

		self.overlaps = overlap_dict	#overlap dict is from RegionList.count_overlaps
		
		#Added later
		self.clusters = {}	# ID:{"member_idx":[], "member_names":[], "cluster_name":"", "representative"}

	def cluster(self, threshold=0.5, method="average"):
		""" Main function to cluster the overlap dictionary into clusters"""

		self.overlap_to_distance()

		if len(self.names) > 1:
			self.linkage_mat = linkage(squareform(self.distance_mat), method)
			self.labels = fcluster(self.linkage_mat, threshold, criterion="distance")		#ordering of the dendrogram

			#Find clusters below threshold
			self.linkage_clusters = dict(zip(range(self.n), [[num] for num in range(self.n)]))
			for i, row in enumerate(self.linkage_mat):
				ID1 = int(row[0])
				ID2 = int(row[1])
				new = self.n + i
				dist = row[2]

				if dist <= threshold:
					self.linkage_clusters[new] = self.linkage_clusters[ID1] + self.linkage_clusters[ID2] + [new]
					del self.linkage_clusters[ID1]
					del self.linkage_clusters[ID2]

			#Add member-names to clusters
			for cluster in self.linkage_clusters:

				self.clusters[cluster] = {"member_idx": [idx for idx in self.linkage_clusters[cluster] if idx < self.n]}
				self.clusters[cluster]["member_names"] = [self.names[idx] for idx in self.clusters[cluster]["member_idx"]]
		
		else:	#only one TF
			self.linkage_clusters = {0:[0]}
			self.linkage_mat = np.array([[0]])
			self.clusters[0] = {"member_idx":[0]}
			self.clusters[0]["member_names"] = [self.names[idx] for idx in self.clusters[0]["member_idx"]]

		self.get_cluster_names()	#Set names of clusters
		self.assign_colors()
		
	def overlap_to_distance(self):
		""" Convert overlap-dict from count_overlaps to distance matrix """
		
		#Find all region names
		names = [key for key in self.overlaps.keys() if isinstance(key, str)]
		names = sorted(names)
		self.n = len(names)

		distance_mat = np.zeros((self.n, self.n))

		for i, id1 in enumerate(names): #rows
			for j, id2 in enumerate(names): #columns
				if i != j:
					tot_overlap = self.overlaps.get((id1,id2), 0)	#Find key-pair otherwise no overlap

					tot_id1 = self.overlaps[id1]
					tot_id2 = self.overlaps[id2]

					id1_dist = 1 - tot_overlap / float(tot_id1) if tot_id1 > 0 else 1 
					id2_dist = 1 - tot_overlap / float(tot_id2) if tot_id2 > 0 else 1

					dist = min([id1_dist, id2_dist])
					distance_mat[i,j] = dist

		self.distance_mat = distance_mat
		self.names = names


	def get_cluster_names(self):
		""" Names each cluster based on the members of the cluster """

		self.cluster_names = {}
		self.name2cluster = {}

		#Sort clusters by distance scores to other motifs in cluster
		#cluster_i = 1
		for cluster in self.clusters:

			members_idx = self.clusters[cluster]["member_idx"]

			if len(members_idx) > 1:

				### Code to assign a specific name to clusters
				members_idx = self.clusters[cluster]["member_idx"]
				distances = {}
				for member in members_idx:
					distances[member] = np.sum(self.distance_mat[member,:])	#0 if member is one-member cluster

				#Sort cluster by distance; smallest values=most representative first
				ordered_idx = sorted(distances.keys(), key=lambda member: distances[member])

				self.clusters[cluster]["representative"] = self.names[ordered_idx[0]]
				self.cluster_names[cluster] = "C_" + self.names[ordered_idx[0]] #cluster is the idx for cluster
			else:		
				self.cluster_names[cluster] = "C_" + self.clusters[cluster]["member_names"][0]
				self.clusters[cluster]["representative"] = self.cluster_names[cluster]
				
			self.clusters[cluster]["cluster_name"] = self.cluster_names[cluster]

			#Assign every member to its cluster
			for member in members_idx:
				self.name2cluster[self.names[member]] = self.cluster_names[cluster] 



	def write_distance_mat(self, out_f):
		""" Writes out distance matrix to out_f file """
		np.savetxt(out_f, self.distance_mat, delimiter="\t", header="\t".join(self.names), fmt="%.4f")


	def assign_colors(self):
		""" Assign colors for plotting the dendrogram """

		clusters = self.linkage_clusters
		no_IDS = self.n

		colorlist = ["blue", "green", "red", "orange"]
		node_color = ["black"] * (2*no_IDS-1)
		i = 0
		for cluster in sorted(list(clusters.keys())):
			if len(clusters[cluster]) > 1:
				color = colorlist[i]
				for node in clusters[cluster]:
					node_color[node] = color
				i += 1 

				if i == len(colorlist):
					i = 0

		self.node_color = node_color #list corresponding to each possible clustering in tree
