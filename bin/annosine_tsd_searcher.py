import sys
import os
import pyfastx
import multiprocessing
import parasail
import numpy as np
import re
import shutil

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
		self.np_encoding = {'-':0, 'A':1, 'C':2, 'G':3, 'T':4}
		self.np_decoding = {0:'-', 1:'A', 2:'C', 3:'G', 4:'T'}
		self.ok_chars = set('ATCG')

	def encode_numpy(self, sequence):
		sequence = sequence.upper()
		sequence = re.sub('[^ATCG]', '-', sequence)
		sequence = np.array([self.np_encoding[c] for c in sequence], dtype = np.int32)
		return sequence
	
	def purge_poly_AT_seq(self, lseq, rseq, threshold = 0):
		left_counts = np.bincount(lseq, minlength = 5)
		right_counts  = np.bincount(rseq, minlength = 5)
		
		#Correspond to sum of counts of C and G in the sequence
		left_cg = left_counts[2] + left_counts[3]
		right_cg = right_counts[2] + right_counts[3]
		
		sequence_is_ok = left_cg > threshold and right_cg > threshold
		
		return sequence_is_ok
	
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
		
		winning_left = None
		winning_right = None
		winning_distance = 1_000_000
		lstart = None
		lend = None
		rstart = None
		rend = None
		
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
							left_indices = left[start:end]
							right_indices = right[start:end]
							
							#Check to see if the recovered sequence is a polyA or polyT or poly AT repeat; these are not TSD candidates
							is_not_poly_AT_sequence = self.purge_poly_AT_seq(left_indices, right_indices, threshold = 1)
							
							if is_not_poly_AT_sequence:
								#Relative locations of start, end in UNGAPPED input strings
								#+1 to switch back to 1-indexing
								lstart = int(start - gaps_left[start]) + 1
								lend   = int(end - gaps_left[end]) + 1
								rstart = right_offset + int(start - gaps_right[start]) + 1
								rend   = right_offset + int(end - gaps_right[end]) + 1
								
								#Check if the sequence is more likely to be biologically correct;
								#i.e., TSDs are closer together / closer to the ends of the SINE
								length_of_seq_between_tsds = rstart - lend
								if length_of_seq_between_tsds < winning_distance:
									winning_distance = length_of_seq_between_tsds
							
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
		
		return winning_left, lstart, lend, winning_right, rstart, rend
	
	#Sequence alignment based approach using parasail
	def tsd_by_sequence_alignment(self, seq_low_id, seq_high_id, output_file, left_chunk = 50, right_chunk = 70, min_ok_length = 10, max_mismatch = 1):
		with open(output_file, 'w') as out:
			for seqid in range(seq_low_id, seq_high_id):
				record = self.fa[seqid]
				
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
				
				print(f'>{seqid}', file = out)
				lseq, lst, lend, rseq, rst, rend = self.find_similar_sequences(left, right, right_offset, min_ok_length = min_ok_length, max_mismatch = max_mismatch)
				if lseq is not None:
					print(f'{lseq} ({lst}-{lend}) {rseq} ({rst}-{rend})', file = out)
				else:
					print('', file = out)
		

def para_search_init(genome_file, output_prefix):
	global gf
	global op
	gf = genome_file
	op = output_prefix
	

def para_search(sequence_range_tuple):
	seqid_low, seqid_high = sequence_range_tuple
	output_file = f'{op}_{seqid_low}_{seqid_high}.txt'
	mn = tsd_finder(gf)
	mn.tsd_by_sequence_alignment(seqid_low, seqid_high, output_file)
	return output_file

#Given
def split_seq_indices(size, num_grps):
	splitsize = 1.0/num_grps*size
	newseq = [((int(round(i*splitsize)), int(round((i+1)*splitsize)),)) for i in range(num_grps)]
	return newseq

def tsd_search_parallel(sequence_file, output_prefix, lc = 50, rc = 70, minlen = 10, max_mis = 1, threads = 1):
	fa = pyfastx.Fasta(sequence_file, build_index = True)
	number_of_seqs = len(fa)
	print(number_of_seqs)
	
	output_file = 'Step2_tsd_new.txt'
	#Don't bother with parallelization if there's only a few sequences
	if number_of_seqs < 10000:
		mn = tsd_finder(sequence_file)
		mn.tsd_by_sequence_alignment(0, number_of_seqs, output_file)
	else:	
		#Sequence index ranges
		groups = split_seq_indices(number_of_seqs, threads)
		with open(output_file, 'wb') as out:
			with multiprocessing.Pool(threads, initializer = para_search_init, initargs = (sequence_file, output_file,)) as pool:
				for result in pool.map(para_search, groups):
					with open(result, 'rb') as inf:
						shutil.copyfileobj(inf, out, 128*1024)
					os.remove(result)

f = 'Step1_extend_tsd_input.fa'
tsd_search_parallel(f, 'test', threads = 10)
