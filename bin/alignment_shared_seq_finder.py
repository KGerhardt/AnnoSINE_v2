import sys
import os
import parasail
import numpy as np
import re

class alignment_tsd_tir_finder:
	def __init__(self, method = 'tsd_searcher', min_ok_length = 10, max_mismatch = 1, polyAT_ok = False, polyAT_threshold = 1, 
				check_inverts = False, gap_penalty = 1, extension_penalty = 0, sf_score_thresh = 10, 
				sf_mismatch_thresh = 2, sf_mismatch_penalty = 1, return_best_only = True):
				
				
		self.method = method
		self.revcmp_table = str.maketrans('ACGTacgt', 'TGCAtgca')
		
		self.np_encoding = {'-':0, 'A':1, 'C':2, 'G':3, 'T':4}
		self.np_decoding = {0:'-', 1:'A', 2:'C', 3:'G', 4:'T'}
		
		self.min_ok_length = min_ok_length
		self.max_mismatch = max_mismatch
		self.max_consecutive_mismatches = 1
		
		self.sf_mm = sf_mismatch_thresh
		self.sf_pen = sf_mismatch_penalty
		self.sf_score = sf_score_thresh
		
		self.polyAT_ok = polyAT_ok
		self.polyAT_threshold = polyAT_threshold
		
		self.check_inverts = check_inverts
		self.gap_penalty = gap_penalty
		self.ext_penalty = extension_penalty
		
		self.best = return_best_only
		
		self.candidates = []

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
	
	def extract_hit_from_df(self, df, l_enc, r_enc, gaps_l, gaps_r):
		#The place in the input strings where the shared subsequence is found
		start = df[0, 2] #First start index of a group
		end = df[-1, 2] + df[-1, 1] #last start index of a group + run length
		left_indices = l_enc[start:end]
		right_indices = r_enc[start:end]

		#Skip polyAT check; can return TSDs which are polyAT
		if self.polyAT_ok:
			is_not_poly_AT_sequence = True
		#Check sequence for A/T percentage; 
		else:
			#Check to see if the recovered sequence is a polyA or polyT or poly AT repeat; these are not TSD candidates
			is_not_poly_AT_sequence = self.purge_poly_AT_seq(left_indices, right_indices)
			
		if is_not_poly_AT_sequence:
			#Relative locations of start, end in UNGAPPED left/right input strings
			lstart = int(start - gaps_l[start])
			lend   = int(end - gaps_l[end-1]) #ends are 1-indexed in string slices, but the gap counts are still 0-indexed; left offset by 1
			rstart = int(start - gaps_r[start])
			rend   = int(end - gaps_r[end-1])
			
			tsd_length = lend - lstart
			tsd_mismatches = int(np.sum(df[df[:,0] == 0][:,1]))
			
			left = []
			right = []
			#Convert numpy ints back to characters; format mismatches accordingly
			for c1, c2 in zip([self.np_decoding[c] for c in left_indices], [self.np_decoding[c] for c in right_indices]):
				if c1 != c2:
					c1 = c1.lower()
					c2 = c2.lower()
				left.append(c1)
				right.append(c2)
			
			#Update the result to return
			left = ''.join(left)
			right = ''.join(right)
		else:
			#Default return case
			left, lstart, lend, right, rstart, rend, tsd_length, tsd_mismatches = None, None, None, None, None, None, None, None
			
		return left, lstart, lend, right, rstart, rend, tsd_length, tsd_mismatches 
	
	#Score similar sequences with sinefinder logic: 
	#TSD score = num_matches - num_mismatches; must exceed score threshold (def. 10) and must start and end with match
	def find_similar_sequences_sinefinder(self, left, right, is_forward = True):
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
	
		current_mismatch_score = 0
		current_match_score = 0
		
		possible_candidates = []
		current_group = []
		for v, s, l in zip(values, start_positions, run_lengths):
			#Next group of matches; this will always be added
			if v:
				current_match_score += l
				current_group.append((v, l, s,))
			#Next group of mismatches; this may be added or skipped
			else:
				#the next mismatch run would exceed acceptable mismatch count
				if current_mismatch_score + (l * self.sf_pen) > self.sf_mm:
					#The current candidate is acceptable; add it to a list
					if current_match_score - current_mismatch_score >= self.sf_score and len(current_group) > 0:
						#This must end in a true because of the way the current group is constructed
						possible_candidates.append(np.array(current_group))
					
					#Reset
					current_mismatch_score = 0
					current_match_score = 0
					
					#If the next mismatch sequence couldn't be added to any segment, just proceed
					if (self.sf_pen * l) > self.sf_mm:
						current_group = []
					
					#The next mismatch segment could possibly be added to a sequence of previously observed matches
					else:
						#Pop previous true + false pairs, check if the next set of mismatches could be added
						while len(current_group) > 2:
							current_group = current_group[2:]
							for i in current_group:
								running_match = 0
								running_mismatch = 0
								#If it's a match
								if i[0]:
									#Add the number of matches
									running_match += i[1]
								else:
									#add the number of mismatches
									running_mismatch += (self.sf_pen * i[1])
							
							#If there is a remaining set of matches, 
							if running_mismatch + (self.sf_pen * l) < self.sf_mm:
								current_group.append((v, l, s,))
								current_match_score = running_match
								current_mismatch_score = running_mismatch + (self.sf_pen * l)
								break
								
						if len(current_group) == 1:
							current_mismatch_score = (self.sf_pen * l)
							current_match_score = current_group[0][1]
							current_group.append((v, l, s,))
					
				else:
					if len(current_group) > 0:
						current_group.append((v, l, s,))
						current_mismatch_score += (l * self.sf_pen)
		#Leftover group not yet added
		if len(current_group) > 0:
			#If the last element was a mismatch, pop it
			if not current_group[-1][0]:
				current_mismatch_score -= current_group[-1][1]
				current_group = current_group[:-1]
			#Add a candidate	
			if current_match_score - current_mismatch_score >= self.sf_score and len(current_group) > 0:
				#This must end in a true because of the way the current group is constructed
				possible_candidates.append(np.array(current_group))

		for df in possible_candidates:
			l, ls, le, r, rs, rend, tsdl, tsd_mm = self.extract_hit_from_df(df, left, right, gaps_left, gaps_right)
			if l is not None:
				next_candidate = (l, ls, le, r, rs, rend, tsdl, tsd_mm, is_forward,)
				self.candidates.append(next_candidate)
	
	#Score similar sequences with TSD searcher logic
	def find_similar_sequences_tsd_searcher(self, left, right, is_forward = True):
		#Convert characters to integers and represent with numpy arrays
		#We ultimately convert back to strings at the end of this, but this makes RLE work, gap finding, etc. much easier
		left = self.encode_numpy(left)
		right = self.encode_numpy(right)
		
		gaps_left = np.cumsum(left == 0)  #counts of '-' characters for position adjustments later
		gaps_right = np.cumsum(right == 0)
		
		all_eq = left == right
		
		run_lengths, start_positions, values = self.rle(all_eq)
		
		lseq = None
		rseq = None
		lstart = None
		rstart = None
		lend = None
		rend = None
		tsd_length = None
		mismatch_length = None
		
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
						tsd_length = run_lengths[0]
						lstart = 0
						rstart = 0
						lend = left.shape[0]
						rend = right.shape[0]
						
				
			#All matches followed by all mismatches or all mismatches followed by all matches
			if run_lengths.shape[0] == 2:
				#First run is matches
				if values[0]:
					if run_lengths[0] >= self.min_ok_length:
						sufficient_length = True
						tsd_length = run_lengths[0]
						lstart = 0
						rstart = 0
						lend = int(run_lengths[0])
						rend = int(run_lengths[0])

				#First run is mismatches
				else:
					if run_lengths[1] >= self.min_ok_length:
						sufficient_length = True
						tsd_length = run_lengths[1]
						lstart = int(start_positions[1])
						rstart = int(start_positions[1])
						lend = int(lstart + run_lengths[1])
						rend = int(rstart + run_lengths[1])
						
			if sufficient_length:
				mismatch_length = 0 #it always will be in these cases
				left_indices = left[lstart:lend]
				#Skip polyAT check; can return TSDs which are polyAT
				if self.polyAT_ok:
					is_not_poly_AT_sequence = True
				#Check sequence for A/T percentage; 
				else:
					#Check to see if the recovered sequence is a polyA or polyT or poly AT repeat; these are not TSD candidates
					is_not_poly_AT_sequence = self.purge_poly_AT_seq(left_indices, left_indices)
				
				if is_not_poly_AT_sequence:
					lseq = ''.join([self.np_decoding[c] for c in left_indices])
					rseq = lseq
					
					next_candidate = (lseq, lstart, lend, rseq, rstart, rend, tsd_length, tsd_mismatches, is_forward,)
					self.candidates.append(next_candidate)
					
					
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
								l, ls, le, r, rs, rend, tsdl, tsd_mm = self.extract_hit_from_df(df, left, right, gaps_left, gaps_right)
								if l is not None:
									next_candidate = (l, ls, le, r, rs, rend, tsdl, tsd_mm, is_forward,)
									self.candidates.append(next_candidate)

						#Reset to continue processing
						next_group = []
		

	#If requested, return only the best matching TSD
	#Biologically, this must be the closest to the original candidate, irrespective of length
	def get_best_hit(self):
		pass
	
	#Sequence alignment based approach using parasail
	def tsd_by_sequence_alignment(self, left_seq, right_seq):
		self.candidates = []
		forward = None
		reverse = None
		#Low penalty semi-global sequence alignment to find repeats within substrings
		res = parasail.sg_trace_striped_sat(left_seq, right_seq, self.gap_penalty, self.ext_penalty, parasail.blosum62)
		
		#Fortunately this never encodes a double '-', so we can ignore a match on that character - it will always be false
		left = res.traceback.query
		right = res.traceback.ref
		
		if self.method == 'tsd_searcher':
			self.find_similar_sequences_tsd_searcher(left, right, forward = True)
		if self.method == 'sinefinder':
			self.find_similar_sequences_sinefinder(left, right, forward = True)

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
						
			if self.method == 'tsd_searcher':
				self.find_similar_sequences_tsd_searcher(left, right)
			if self.method == 'sinefinder':
				self.find_similar_sequences_sinefinder(left, right)
			
			#Flip the right sequence, start and stop indices to forward orientation
			if rseq is not None:
				fst, fend = rsl-rend, rsl-rst
				rseq = self.revcomp(rseq)
				rst, end = fst, fend
					
				reverse = (lseq, lst, lend, rseq, rst, rend,)
			
		return forward, reverse
	

#This code doesn't actually need this, I'm using it for testing
import pyfastx
mn = alignment_tsd_tir_finder(check_inverts = False, polyAT_ok = False, method = 'sinefinder')		
fa = pyfastx.Fasta('output/Step1_extend_tsd_input.fa')
for seq in fa:
	l = seq.seq[0:50]
	r = seq.seq[-70:]
	#r = seq.seq
	f, r = mn.tsd_by_sequence_alignment(l, r)
	if len(mn.candidates) > 0:
		for c in mn.candidates:
			print(c)
	


	
	