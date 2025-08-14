import sys
import os
import parasail
import numpy as np
import re

class alignment_tsd_tir_finder:
	def __init__(self, min_ok_length = 10, max_mismatch = 1, polyAT_ok = False, polyAT_threshold = 1, 
				check_inverts = False, gap_penalty = 1, extension_penalty = 0):
				
		self.revcmp_table = str.maketrans('ACGTacgt', 'TGCAtgca')
		
		self.np_encoding = {'-':0, 'A':1, 'C':2, 'G':3, 'T':4}
		self.np_decoding = {0:'-', 1:'A', 2:'C', 3:'G', 4:'T'}
		
		self.min_ok_length = min_ok_length
		self.max_mismatch = max_mismatch
		self.max_consecutive_mismatches = 1
		self.polyAT_ok = polyAT_ok
		self.polyAT_threshold = polyAT_threshold
		
		self.check_inverts = check_inverts
		self.gap_penalty = gap_penalty
		self.ext_penalty = extension_penalty

	#To find inverted repeats, revcmp the right string and realign
	def revcomp(self, sequence):
		sequence = sequence[::-1]
		sequence = sequence.translate(self.revcmp_table)
		return sequence

	#Convert a string to a numeric representation to make vector ops a little easier to work with in python
	def encode_numpy(self, sequence):
		sequence = re.sub('[^ATCGatcg]', '-', sequence)
		sequence = np.array([self.np_encoding[c] for c in sequence], dtype = np.int32)
		return sequence

	#Vectorized numpy run-length encoding function;
	#Credit to #https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encodi
	#Here, used to encode runs of matching and mismatching characters from encoded numeric numpy array
	def rle(self, ia):
		""" run length encoding. Partial credit to R rle function. 
			Multi datatype arrays catered for including non Numpy
			returns: tuple (runlengths, startpositions, values) """
		n = len(ia)
		if n == 0: 
			return None, None, None
		else:
			y = ia[1:] != ia[:-1]               # pairwise unequal (string safe)
			i = np.append(np.where(y), n - 1)   # must include last element posi
			z = np.diff(np.append(-1, i))       # run lengths
			p = np.cumsum(np.append(0, z))[:-1] # positions
			return z, p, ia[i]

	def purge_poly_AT_seq(self, lseq, rseq):
		left_counts = np.bincount(lseq, minlength = 5)
		right_counts  = np.bincount(rseq, minlength = 5)
		
		#Correspond to sum of counts of C and G in the sequence
		left_cg = left_counts[2] + left_counts[3]
		right_cg = right_counts[2] + right_counts[3]
		
		sequence_is_ok = left_cg > self.polyAT_threshold and right_cg > self.polyAT_threshold
		
		return sequence_is_ok
	
	#Sequence alignment based approach using parasail
	def find_similar_sequences(self, left, right):
		#Convert characters to integers and represent with numpy arrays
		#We ultimately convert back to strings at the end of this, but this makes RLE work, gap finding, etc. much easier
		left = self.encode_numpy(left)
		right = self.encode_numpy(right)
		
		gaps_left = np.cumsum(left == 0)  #counts of '-' characters for position adjustments later
		gaps_right = np.cumsum(right == 0)
		
		all_eq = left == right
		
		run_lengths, start_positions, values = self.rle(all_eq)

		winning_left = None
		winning_right = None
		winning_distance = 1_000_000
		
		winning_lstart = None
		winning_lend = None
		winning_rstart = None
		winning_rend = None
		
		if run_lengths.shape[0] < 3:
			'''
			Edge cases
			All matches or all mismatches, run_lengths.shape[0] = 1
			Run of mismatches - > matches or vice versa
			'''
			
			sufficient_length = False
			
			#There are either only matches or only mismatches
			if run_lengths.shape[0] == 1:
				#If it's all matches, return the strings unmodified and at full length;
				#if this is false, there is nothing to return
				if values[0]:
					if run_lengths[0] >= self.min_ok_length:
						sufficient_length = True
						winning_lstart = 0
						winning_rstart = 0
						winning_lend = left.shape[0]
						winning_rend = right.shape[0]
						
				
			#All matches followed by all mismatches or all mismatches followed by all matches
			if run_lengths.shape[0] == 2:
				#First run is matches
				if values[0]:
					if run_lengths[0] >= self.min_ok_length:
						sufficient_length = True
						winning_lstart = 0
						winning_rstart = 0
						winning_lend = int(run_lengths[0])
						winning_rend = int(run_lengths[0])

				#First run is mismatches
				else:
					if run_lengths[1] >= self.min_ok_length:
						sufficient_length = True
						winning_lstart = int(start_positions[1])
						winning_rstart = int(start_positions[1])
						winning_lend = int(winning_lstart + run_lengths[1])
						winning_rend = int(winning_rstart + run_lengths[1])
						
			if sufficient_length:
				left_indices = left[winning_lstart:winning_lend]
				#Skip polyAT check; can return TSDs which are polyAT
				if self.polyAT_ok:
					is_not_poly_AT_sequence = True
				#Check sequence for A/T percentage; 
				else:
					#Check to see if the recovered sequence is a polyA or polyT or poly AT repeat; these are not TSD candidates
					is_not_poly_AT_sequence = self.purge_poly_AT_seq(left_indices, left_indices)
				
				if is_not_poly_AT_sequence:
					winning_left = ''.join([self.np_decoding[c] for c in left_indices])
					winning_right = winning_left
				else:
					winning_lstart = None
					winning_lend = None
					winning_rstart = None
					winning_rend = None
							
		else:
			next_group = []					
			#Group runs of aligned sequences separated by no more than one mismatch
			for l, v, s in zip(run_lengths, values, start_positions):
				if v: #the strings have the same character over the next run_length positions
					next_group.append((v, l, s,))
				else: #The strings have a run of one or more non-matching chars
					if l == self.max_consecutive_mismatches and len(next_group) > 0: #The size of the no-match run is exactly 1
						next_group.append((v, l, s,))
					else: #The size of the no-match run is > 1 / this is a new run
						#Process the group to clean to actual putatitve TSDs
						if len(next_group) == 0:
							continue
						else:
							#Because of the way we build the above, these arrays always start and end with matches
							df = np.array(next_group)
							df_size = df.shape[0]
							#Sum the lengths of mismatch runs over all mismatches in the array
							num_mismatches = np.sum(df[df[:,0] == 0][:,1])

							#Find only the longest run of matched characters satisfying the rules
							if num_mismatches > self.max_mismatch:
								'''Rules: 
									(1) There cannot be more than max_mismatch mismatches in the string
									(2) The output must start and end with a match
									(3) Output must be at least min_ok_length characters
									(4) Return only the longest string satisfying these conditions, if any.
								'''
								#Find each consecutive sub-dataframe satisfying the above conditions
								win_start = 0
								win_length = 0
								#This logic will need updating if there's a different number of max acceptable mismatches
								check_size = 2*self.max_mismatch
								cut_size = check_size + 1
								
								for i in range(0, df_size - check_size, 2):
									#Sum of run_lengths for this dataframe = output string length; this a check to find which run of sequences is most satisfactory
									this_length = np.sum(df[i:i+cut_size, 1])
									if this_length > win_length:
										win_length = this_length
										win_start = i
								
								#Select out the winning dataframe
								df = df[win_start:win_start + cut_size, :]
							
							#The place in the input strings where the shared subsequence is found
							start = df[0, 2] #First start index of a group
							end = df[-1, 1] + df[-1, 2] #last start index of a group + run length
							
							#The longest shared subsequence passing rules still has a minimum OK size
							run_size = end - start
							if run_size >= self.min_ok_length:
								left_indices = left[start:end]
								right_indices = right[start:end]
								
								#Skip polyAT check; can return TSDs which are polyAT
								if self.polyAT_ok:
									is_not_poly_AT_sequence = True
								#Check sequence for A/T percentage; 
								else:
									#Check to see if the recovered sequence is a polyA or polyT or poly AT repeat; these are not TSD candidates
									is_not_poly_AT_sequence = self.purge_poly_AT_seq(left_indices, right_indices)
								
								if is_not_poly_AT_sequence:
									#Relative locations of start, end in UNGAPPED left/right input strings
									lstart = int(start - gaps_left[start])
									lend   = int(end - gaps_left[end-1]) #ends are 1-indexed in string slices, but the gap counts are still 0-indexed; left offset by 1
									rstart = int(start - gaps_right[start])
									rend   = int(end - gaps_right[end-1])
									
									#Find the longest TSD
									length_of_seq_between_tsds = rstart - lend
									if length_of_seq_between_tsds < winning_distance:
										winning_distance = length_of_seq_between_tsds
								
										winning_lstart = lstart
										winning_lend = lend
										winning_rstart = rstart
										winning_rend = rend
								
										winning_left = []
										winning_right = []
										#Convert numpy ints back to characters; format mismatches accordingly
										for c1, c2 in zip([self.np_decoding[c] for c in left_indices], [self.np_decoding[c] for c in right_indices]):
											if c1 != c2:
												c1 = c1.lower()
												c2 = c2.lower()
											winning_left.append(c1)
											winning_right.append(c2)
										
										#Update the result to return
										winning_left = ''.join(winning_left)
										winning_right = ''.join(winning_right)	
										
						#Reset to continue processing
						next_group = []
		
		return winning_left, winning_lstart, winning_lend, winning_right, winning_rstart, winning_rend
	
	#Sequence alignment based approach using parasail
	def tsd_by_sequence_alignment(self, left_seq, right_seq):
		forward = None
		reverse = None
		#Low penalty semi-global sequence alignment to find repeats within substrings
		res = parasail.sg_trace_striped_sat(left_seq, right_seq, self.gap_penalty, self.ext_penalty, parasail.blosum62)
		
		#Fortunately this never encodes a double '-', so we can ignore a match on that character - it will always be false
		left = res.traceback.query
		right = res.traceback.ref
				
		lseq, lst, lend, rseq, rst, rend = self.find_similar_sequences(left, right)
		
		if lseq is not None:
			forward = (lseq, lst, lend, rseq, rst, rend,)

		if self.check_inverts:
			rsl = len(right_seq)
			og_rseq = right_seq
			right_seq = self.revcomp(right_seq)
			#Low penalty semi-global sequence alignment to find repeats within substrings
			#res = parasail.sg_trace_striped_sat(left_seq, right_seq, self.gap_penalty, self.ext_penalty, parasail.blosum62)
			res = parasail.sw_trace_striped_sat(left_seq, right_seq, self.gap_penalty, self.ext_penalty, parasail.blosum62)
			
			#Fortunately this never encodes a double '-', so we can ignore a match on that character - it will always be false
			left = res.traceback.query
			right = res.traceback.ref
						
			lseq, lst, lend, rseq, rst, rend = self.find_similar_sequences(left, right)

			#Flip the right sequence, start and stop indices to forward orientation
			if rseq is not None:
				fst, fend = rsl-rend, rsl-rst
				rseq = self.revcomp(rseq)
				rst, end = fst, fend
					
				reverse = (lseq, lst, lend, rseq, rst, rend,)
			
		return forward, reverse
	

#This code doesn't actually need this, I'm using it for testing
import pyfastx
mn = alignment_tsd_tir_finder(check_inverts = False)		
fa = pyfastx.Fasta('output/Step1_extend_tsd_input.fa')
for seq in fa:
	l = seq.seq[0:50]
	r = seq.seq[-70:]
	#r = seq.seq
	f, r = mn.tsd_by_sequence_alignment(l, r)
	print(f)


	
	