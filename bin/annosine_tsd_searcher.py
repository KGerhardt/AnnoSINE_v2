import sys
import os
import pyfastx
import multiprocessing
import parasail
import numpy as np
import re

#Fast numpy RLE	
#https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encodi
def rle(ia):
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

class tsd_finder:
	def __init__(self, search_sequences_file):
		self.sequences_file = search_sequences_file
		self.fa = pyfastx.Fasta(self.sequences_file)
		self.np_encoding = {'-':-1, 'A':0, 'C':1, 'G':2, 'T':3}
		self.np_decoding = {-1:'-', 0:'A', 1:'C', 2:'G', 3:'T'}
		self.ok_chars = set('ATCG')

	def encode_numpy(self, sequence):
		sequence = sequence.upper()
		sequence = re.sub('[^ATCG]', '-', sequence)
		sequence = np.array([self.np_encoding[c] for c in sequence], dtype = np.int32)
		return sequence
	
	#Sequence alignment based approach using parasail
	def find_similar_sequences(self, left, right, right_offset, min_ok_length = 10, max_mismatch = 1):
		#Convert characters to integers and represent with numpy arrays
		#We ultimately convert back to strings at the end of this, but this makes RLE work, gap finding, etc. much easier
		left = self.encode_numpy(left)
		right = self.encode_numpy(right)
	
		gaps_left = np.cumsum(left == -1)  #counts of '-' characters for position adjustments later
		gaps_right = np.cumsum(right == -1)
		
		all_eq = left == right
		run_lengths, start_positions, values = rle(all_eq)
		
		
		cleanup = []
		next_group = []
		#Group runs of aligned sequences separated by no more than one mismatch
		for l, v, s in zip(run_lengths, values, start_positions):
			if v: #the strings have the same character over the next run_length positions
				next_group.append((v, l, s,))
			else: #The strings have a run of one or more non-matching chars
				if l == 1: #The size of the no-match run is exactly 1
					next_group.append((v, l, s,))
				else: #The size of the no-match run is > 1 / this is a new run
					#Process the group to clean to actual putatitve TSDs
					if len(next_group) == 0:
						continue
					else:
						#Because of the way we build the above, these arrays always start and end with matches
						df = np.array(next_group)
						df_size = df.shape[0]
						#The dataframe always looks like True, (True, False) repeating, True, so this math finds the count of 1-length mismatches
						num_mismatches = ((df_size - 1) / 2)

						#Find only the longest run of matched characters satisfying the rules
						if num_mismatches > max_mismatch:
							'''Rules: 
								(1) There cannot be more than max_mismatch mismatches in the string
								(2) The output must start and end with a match
								(3) Output must be at least min_ok_length characters
								(4) Return only the longest string satisfying these conditions, if any.
							'''
							#Find each consecutive sub-dataframe satisfying the above conditions
							win_start = 0
							win_length = 0
							check_size = 2*max_mismatch
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
						if run_size >= min_ok_length:
							lslice = []
							rslice = []
							#Convert numpy ints back to characters; format mismatches accordingly
							for c1, c2 in zip([self.np_decoding[c] for c in left[start:end]], [self.np_decoding[c] for c in right[start:end]]):
								if c1 != c2:
									c1 = c1.lower()
									c2 = c2.lower()
								lslice.append(c1)
								rslice.append(c2)
							
							#Make char arrays into strings
							lslice = ''.join(lslice)
							rslice = ''.join(rslice)	
							
							#Relative locations of start, end in UNGAPPED input strings
							#+1 to switch back to 1-indexing
							lstart = int(start - gaps_left[start]) + 1
							lend   = int(end - gaps_left[end]) + 1
							rstart = right_offset + int(start - gaps_right[start]) + 1
							rend   = right_offset + int(end - gaps_right[end]) + 1
							
							#encode info
							next_run = (lslice, lstart, lend, rslice, rstart, rend)
							
							cleanup.append(next_run)
					
					#Reset to continue processing
					next_group = []
		
		return cleanup
	
	#Sequence alignment based approach using parasail
	def tsd_by_sequence_alignment(self, left_chunk = 50, right_chunk = 70, min_ok_length = 10, max_mismatch = 1):
		for record in self.fa:
			seqlen = len(record.seq)

			left_sequence = record.seq[0:(left_chunk + 1)]
			right_offset = seqlen - right_chunk
			right_sequence = record.seq[right_offset:]
			
			#Low penalty semi-global sequence alignment to find repeats within substrings
			res = parasail.sg_trace_striped_sat(left_sequence, right_sequence, 1, 0, parasail.blosum62)
			
			#Fortunately this never encodes a double '-', so we can ignore a match on that character - it will always be false
			left = res.traceback.query
			right = res.traceback.ref
			
			seqid = record.description
			
			my_repeats = self.find_similar_sequences(left, right, right_offset, min_ok_length = min_ok_length, max_mismatch = max_mismatch)
			print(f'>{seqid}')
			printout = []
			for r in my_repeats:
				lseq, lst, lend, rseq, rst, rend = r
				printout.append(f'{lseq} ({lst}-{lend}) {rseq} ({rst}-{rend})')
			printout = ' '.join(printout)
			print(printout)
			
				
			
			
mn = tsd_finder('Step1_extend_tsd_input.fa')
mn.tsd_by_sequence_alignment()

'''
l = 'ATCG'
r = 'ATCGTTT'

res = parasail.sg_trace_scan(l, r, 11, 1, parasail.blosum62)
print(dir(res))
#print(res.cigar.seq)
print(dir(res.traceback))
print(res.traceback.query)
print(res.traceback.ref)
'''