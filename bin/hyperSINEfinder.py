import sys
import os
import multiprocessing
import sqlite3
import pyfastx
import hyperscan
import numpy as np

from itertools import product

#Set OMP threads for TSD searcher to 1
os.environ['OMP_NUM_THREADS'] = '1'
from tsd_searcher import alignment_tsd_tir_finder

def on_match(id: int, fr: int, to: int, flags: int, context:  None):
	context.append((id, fr, to,))

class hyperSINEfinder:
	def __init__(self, genome):
		self.fa = pyfastx.Fasta(genome, build_index = True)
		conn = sqlite3.connect(f'{genome}.fxi')
		curs = conn.cursor()
		self.seqlens = {seqid:seqlen for seqid, seqlen in curs.execute("SELECT chrom, blen FROM seq").fetchall()}
		curs.close()
		conn.close()
		
		#All possible SINE forward patterns exploded out to explicit regex with no ambiguity
		self.patterns = None
		self.pattern_reverser = None
		#Genome sequence to search
		self.sequence = None
		#Hyperscan database of regex patterns
		self.hsdb = None
		
		#self.tsd_checker = alignment_tsd_tir_finder(method = 'sinefinder', return_best_only = True)
		#self.tsd_checker = alignment_tsd_tir_finder(method = 'tsd_searcher', return_best_only = True)
		self.tsd_checker = 	alignment_tsd_tir_finder(min_ok_length = 10, max_mismatch = 2, return_best_only = True, best_hit_approach = 'closest')

	
	def revcomp(self, string):
		switcher = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
		revcmp = ''.join([switcher[c] if c in switcher else 'N' for c in string[::-1]])
		return revcmp
		
	def set_patterns(self):
		'''
		#Original pattern - I omit spacers because we do not need to explicitly search for them
		pattern = (
		('TSD_region_1', ".{,40}"),
		('a_box', "[GA][CGA]TGG"),
		('spacer_1', ".{25,50}"),
		('b_box', "GTTC[AG]A"),
		('spacer_2', ".{20,500}?"),
		('polyA', "A{6,}|T{6,}"),
		('TSD_region_2', ".{,40}"),
		'''
	
		a_box_fwd = ('AATGG',
					 'ACTGG',
					 'AGTGG',
					 'GATGG',
					 'GCTGG',
					 'GGTGG')
					
		b_box_fwd = ('GTTCAA',
					 'GTTCGA')
					
		poly_AT = ('A{6,}',
				   'T{6,}')
				   
		a_box_rev = [self.revcomp(rgx) for rgx in a_box_fwd]
		b_box_rev = [self.revcomp(rgx) for rgx in b_box_fwd]
		
		self.patterns = (a_box_fwd, b_box_fwd, poly_AT, b_box_rev, a_box_rev)
		
	def make_hsdb(self):
		patterns = []
		index = 0
		for group in self.patterns:
			for rgx in group:
				next_pat = (rgx.encode(encoding = 'ascii'), index, hyperscan.HS_FLAG_SOM_LEFTMOST,)
				patterns.append(next_pat)
			index += 1
		
		self.hsdb = hyperscan.Database(mode=hyperscan.HS_MODE_BLOCK)
		
		expressions, ids, flags = zip(*patterns)
								
		self.hsdb.compile(
			expressions=expressions, ids=ids, elements=len(patterns), flags=flags
		)
		
	def prep(self):
		self.set_patterns()
		self.make_hsdb()
	
	def remove_ungreedy(self, hits_array):
		#Remove ungreedy hits		
		#i row and i+1 row have the same regex pattern id (here, both are polyAT)
		same_id = hits_array[:,0][:-1] == hits_array[:,0][1:] 
		#i row and i+1 row have the same start location (hyperscan returned end-truncated hits for a polyAT)
		same_start = hits_array[:,1][:-1] == hits_array[:,1][1:] 
		#i row has a greater end index than i+1 row; i+1 row is an insufficiently greedy match if both previous conditiions
		greatest_end = hits_array[:,2][:-1] > hits_array[:,2][1:]
	
		#for this program, this only affects polyAT and keeps the longest polyAT for a given genomic region
		filterer = np.zeros(shape = hits_array.shape[0], dtype = np.bool_)
		#Find indices with greedy matches
		filterer[1:] = np.logical_and(np.logical_and(same_id, same_start), greatest_end)
		filterer = np.logical_not(filterer) #invert it so that true = most greedy indices only
		
		hits_array = hits_array[np.where(filterer)] #Remove all instances of insufficiently greedy hits

		del filterer
		
		return hits_array
	
	'''
	#Fast numpy run length encoding
	#https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encodi
	def rle(self, inarray):
			""" run length encoding. Partial credit to R rle function. 
				Multi datatype arrays catered for including non Numpy
				returns: tuple (runlengths, startpositions, values) """
			n = inarray.shape[0]
			y = inarray[1:] != inarray[:-1]               # pairwise unequal (string safe)
			i = np.append(np.where(y), n - 1)   # must include last element posi
			z = np.diff(np.append(-1, i))       # run lengths
			p = np.cumsum(np.append(0, z))[:-1] # positions
			indices = np.arange(z.shape[0], dtype = np.int32)
			
			rle_arr = np.column_stack([indices, inarray[i], p, z])
			
			return rle_arr
	
	'''	
	'''
	#Find zeros close enough to ones, should record which pairs pass this 
	
	In this, I refer to zero pattern, one pattern, and two pattern. This is a_box, b_box, polyAT for forward, polyAT, b_box, a_box for reverse
	'''
	def find_plausible_sine_patterns(self, hits, forward = True):
		
		if forward:
			#I encoded these patterns as 0, 1, 2 forward, 2 3, 4 reverse. PolyAT is identical both directions, so grouping 2 with both works
			first = 0
			second = 1
			third = 2
			#Original SINEfinder pattern searcher uses these distances to identify what is OK offset
			dist1_low = 25
			dist1_high = 50
			dist2_low = 20
			dist2_high = 500
		else:
			first = 2
			second = 3
			third = 4
			#In the reverse search, distance orders are inverted
			dist1_low = 20
			dist1_high = 500
			dist2_low = 25
			dist2_high = 50
		
		#Final return
		sine_candidates = None
		
		#First pass - identify first and second patterns close enough together; keep all twos
		zero_mask = hits[:, 0] == first
		zeros = hits[zero_mask]
		zero_orig_indices = np.flatnonzero(zero_mask)
		
		one_mask = hits[:, 0] == second
		ones = hits[one_mask]
		one_orig_indices = np.flatnonzero(one_mask)
		
		sort_order = zeros[:, 2].argsort()
		sorted_zeros = zeros[sort_order]
		sorted_zero_orig_indices = zero_orig_indices[sort_order]
		sorted_zeros_end = sorted_zeros[:, 2]
		
		#Find loci where the end of a zero pattern would be inserted into offset arrays of the start of a one pattern
		left = np.searchsorted(sorted_zeros_end, ones[: , 1]-dist1_high, side = 'left')
		right = np.searchsorted(sorted_zeros_end, ones[: , 1]-dist1_low, side = 'right')
		
		#Places where the indices differ indicate that a zero pattern is close enough to at least one one pattern
		worthwhile_indices = np.where(left < right)[0]
		#passing_ones = one_orig_indices[worthwhile_indices]
			
		acceptable_zero_ones = []
		
		for i in worthwhile_indices:
			l, r = left[i], right[i]
			for zidx in sorted_zero_orig_indices[l:r]:
				acceptable_zero_ones.append((zidx, one_orig_indices[i]))
		
		acceptable_zero_ones = np.array(acceptable_zero_ones, dtype = np.int32)
		
		#two_mask = np.flatnonzero(hits[:, 0] == third)
		
		#ok_first_pass = np.union1d(acceptable_zero_ones.flatten(), two_mask)
		
		#if ok_first_pass.shape[0] > 0:
		if acceptable_zero_ones.shape[0] > 0:
			#Second pass - of the remaining ones, identify all ones and twos close enough together; keep all zeros
			
			#No need to redo this
			#one_mask = hits[:, 0] == second
			#ones = hits[one_mask]
			#one_orig_indices = np.flatnonzero(one_mask)
			
			two_mask = hits[:, 0] == third
			twos = hits[two_mask]
			two_orig_indices = np.flatnonzero(two_mask)
			
			sort_order = ones[:, 2].argsort()
			sorted_ones = ones[sort_order]
			sorted_one_orig_indices = one_orig_indices[sort_order]
			sorted_ones_end = sorted_ones[:, 2]
			
			left = np.searchsorted(sorted_ones_end, twos[: , 1]- dist2_high, side = 'left')
			right = np.searchsorted(sorted_ones_end, twos[: , 1]- dist2_low, side = 'right')
			worthwhile_indices = np.where(left < right)[0]
			
			acceptable_one_twos = []
			
			for i in worthwhile_indices:
				l, r = left[i], right[i]
				for oidx in sorted_one_orig_indices[l:r]:
					acceptable_one_twos.append((oidx, two_orig_indices[i],))

			acceptable_one_twos = np.array(acceptable_one_twos, dtype = np.int32)
			
			#zero_mask = np.flatnonzero(hits[:, 0] == first)			
			
			#ok_second_pass = np.union1d(acceptable_one_twos.flatten(), zero_mask)
						
			#if ok_second_pass.shape[0] > 0:
			if acceptable_one_twos.shape[0] > 0:
				
				#Cross-check
				#which zeros have a one in the one_two index
				acceptable_zero_ones = acceptable_zero_ones[np.isin(acceptable_zero_ones[:, 1], acceptable_one_twos[:, 0])]
				#Which ones have a two in the zero_one index after filtering
				acceptable_one_twos = acceptable_one_twos[np.isin(acceptable_one_twos[:, 0], acceptable_zero_ones[:, 1])]
				
				#Do we need some filter here to check and make sure that we're only constructing max-length patterns?
				
				#Union between all zero - ones, drop ones from one_twos because they're already checked, add two-twos
				#all_acceptable = np.union1d(acceptable_zero_ones.flatten(), acceptable_one_twos[:, 1])
				
				#if all_acceptable.shape[0] > 0:
				#	hits = hits[all_acceptable]
				
				#Remove pattern ID
				hits = hits[:, 1:]
			
				#What we want to do here is to find the furthest away zero for each one value and the furthest away two for each one value,
				#But, if the same zero and two appear on either side of the one, then we want to choose the one that maximizes 0-1 distance?
				
				
				#For each distinct 2, the last time a specific zero value appears must be the max 0 -> 1 distance with a valid two, because ones are ascending
				#For each distinct 1, the first time a specific zero appears 
				
				if acceptable_zero_ones.shape[0] > 0 and acceptable_one_twos.shape[0] > 0:
					#sine_candidates = []
					patterns = {}
					pp = []
					
					#Effective grouping mechanism
					left = np.searchsorted(acceptable_one_twos[:,0], acceptable_zero_ones[:, 1], side = 'left')
					right = np.searchsorted(acceptable_one_twos[:,0], acceptable_zero_ones[:, 1], side = 'right')
					
					
					for zero, one, l, r in zip(acceptable_zero_ones[:, 0], acceptable_zero_ones[:, 1], left, right):
						#if zero not in patterns:
						#	patterns[zero] = {}
						for two in acceptable_one_twos[l:r, 1]:
							#patterns[zero][one] = two
							pp.append((zero, one, two,))
					
					pp = np.array(pp)
					
					if forward:
						#Distance from 0 -> 1
						dist1 = pp[:, 1] - pp[:, 0]
						dist2 = pp[:, 2] - pp[:, 1]
					else:
						#We would be searching in reverse here, so we'd favor a_box -> b_box distance anyway
						dist1 = pp[:, 2] - pp[:, 1]
						dist2 = pp[:, 1] - pp[:, 0]
						
					#Sort to maximize a_box to b_box distance, then by b_box to polyAT distance
					pp = pp[np.lexsort((dist2, dist1,))]
					
					#Find first appearance of each 0, 1, and 2 pattern in order
					unique_elements, first_indices = np.unique(pp[:, 0], return_index=True)
					pp = pp[first_indices]
					unique_elements, first_indices = np.unique(pp[:, 1], return_index=True)
					pp = pp[first_indices]
					unique_elements, first_indices = np.unique(pp[:, 2], return_index=True)
					pp = pp[first_indices]
						
					#Fix ordering for pretty print
					pp = pp[np.argsort(pp[:, 0])]
					
					
					sine_candidates = []	
					for row in pp:
						z, o, t = row[0], row[1], row[2]
						next_candidate = np.concatenate([hits[z], hits[o], hits[t]],)
						sine_candidates.append(next_candidate)
					
					sine_candidates = np.array(sine_candidates)
					#print(sine_candidates[0:75])
					#print(sine_candidates.shape)
					#print('##########')
					
			else:
				hits = None
		else:
			hits = None
		
		return sine_candidates

	#Prepare matches for further processing
	def clean_and_group_matches(self, hits):
		#Convert recovered hits to numpy
		hits_array = np.array(hits)
	
		#sort by regex pattern, start, end location - lexsort reverses the ordering of these components
		hits_array = hits_array[np.lexsort([-hits_array[:,2], hits_array[:,1], hits_array[:,0]])]
				
		#Find (polyAT) items which are not matched at max length and remove them
		hits_array = self.remove_ungreedy(hits_array)
		
		#Sort by start location irrespective of pattern
		hits_array = hits_array[hits_array[:,1].argsort()]
		
		#np.savetxt('raw_hits.txt', hits_array, delimiter='\t', fmt = '%d')
		
		#Break into forward and reverse hits; include polyAT in both
		forward_hits = hits_array[hits_array[:, 0] < 3]
		reverse_hits = hits_array[hits_array[:, 0] > 1]
		
		del hits_array
		
		if len(forward_hits) > 0:
			forward_hits = self.find_plausible_sine_patterns(forward_hits, forward = True)
		else:
			forward_hits = None
			
		
		if len(reverse_hits) > 0:
			reverse_hits = self.find_plausible_sine_patterns(reverse_hits, forward = False)
		else:
			reverse_hits = None
			
		return forward_hits, reverse_hits
		
	def extract_match(self, sequence, recovery, description, forward = True):
		writeout = []
		
		for row in recovery:
			if forward:
				a_box_start  = row[0]
				a_box_end    = row[1]
				b_box_start  = row[2]
				b_box_end    = row[3]
				polyAT_start = row[4]
				polyAT_end   = row[5]
				
				left_tsd_region  = a_box_start - 40
					
				right_tsd_region = polyAT_end + 41
								
				tsd1    = sequence[left_tsd_region:a_box_start]
				a_box   = sequence[a_box_start:a_box_end]
				spacer1 = sequence[a_box_end:b_box_start]
				b_box   = sequence[b_box_start:b_box_end]
				spacer2 = sequence[b_box_end:polyAT_start]
				polyAT  = sequence[polyAT_start:polyAT_end]
				tsd2    = sequence[polyAT_end:right_tsd_region]
				
			else:
				polyAT_start = row[0]
				polyAT_end   = row[1]
				b_box_start  = row[2]
				b_box_end    = row[3]
				a_box_start  = row[4]
				a_box_end    = row[5]
				
				left_tsd_region  = polyAT_start - 40
				right_tsd_region = a_box_end + 41
				
				tsd2    = self.revcomp(sequence[left_tsd_region:polyAT_start])
				polyAT  = self.revcomp(sequence[polyAT_start:polyAT_end])
				spacer2 = self.revcomp(sequence[polyAT_end:b_box_start])
				b_box   = self.revcomp(sequence[b_box_start:b_box_end])
				spacer1 = self.revcomp(sequence[b_box_end:a_box_start])
				a_box   = self.revcomp(sequence[a_box_start:a_box_end])
				tsd1    = self.revcomp(sequence[a_box_end:right_tsd_region])
		
			#Assess whether this SINE candidate actually starts and ends with a TSD
			tsds = self.tsd_checker.operate(tsd1, tsd2)
			#print(tsds)
			#If it does, the best hit is the first hit according to tsd_checker with options selected
			if tsds is not None:
				#tsd = (updates[0], updates[1], left_string_start, left_string_end, right_string_start, right_string_end, updates[2], updates[3]+updates[4], )
				left_tsd, right_tsd, ls, le, rs, rend, tsdl, tsd_mm = tsds[0]
			
				tsd1_spacer = tsd1[le:]
				tsd1 = tsd1[ls:le]
				tsd2_spacer = tsd2[:rs]
				tsd2 = tsd2[rs:rend]
				#SINEfinder output header format
				#>chr1_pat F 876520:876774 TSD-len=10;TSD-score=10;TSD-mism=0
				next_sine_header = f'>{description} {"+" if forward else "-"} {left_tsd_region+ls}:{polyAT_end+rend} TSD-len={tsdl};TSD-score={tsdl-tsd_mm};TSD-mism={tsd_mm}'
				
				#SINEfinder sequence format alternates caps and lowercase with TSDs, a + b boxes, polyAT caps and all else lower
				tsd1 = tsd1.upper()
				tsd1_spacer = tsd1_spacer.lower()
				a_box = a_box.upper()
				spacer1 = spacer1.lower()
				b_box = b_box.upper()
				spacer2 = spacer2.lower()
				polyAT = polyAT.upper()
				tsd2_spacer = tsd2_spacer.lower()
				tsd2 = tsd2.upper()
				next_sine_sequence = f'{tsd1}{tsd1_spacer}{a_box}{spacer1}{b_box}{spacer2}{polyAT}{tsd2_spacer}{tsd2}'
				
				writeout.append(next_sine_header)
				writeout.append(next_sine_sequence)
		
		return writeout
		
	def execute_search(self, sequence_id):
		record = self.fa[sequence_id]
		sequence = record.seq.upper()
		seqlen = len(sequence)
		desc = record.description.split()[0]
		hits = []
		
		#print(f'{sequence_id} loaded')
		
		self.hsdb.scan(sequence.encode(encoding='ascii'), match_event_handler=on_match, context = hits)
		
		#print(f'{sequence_id} scanned')
		
		f = None
		r = None
		
		if len(hits) > 0:
			#Definitely need to remove ungreedy polyAT hits
			f, r = self.clean_and_group_matches(hits)

		#print(f'{sequence_id} cleaned')

		'''
		if f is None:
			fseq = []
		else:
			fseq = []
			for row in f:
				a_box_start  = row[0]
				a_box_end    = row[1]
				b_box_start  = row[2]
				b_box_end    = row[3]
				polyAT_start = row[4]
				polyAT_end   = row[5]
				
				left_tsd_region  = a_box_start - 40
					
				right_tsd_region = polyAT_end + 41
				
				myseq = sequence[left_tsd_region:right_tsd_region]
				fseq.append('>'+sequence_id)
				fseq.append(myseq)
				
				
			
		if r is None:
			rseq = []
			
		else:
			rseq = []
			for row in r:
				polyAT_start = row[0]
				polyAT_end   = row[1]
				b_box_start  = row[2]
				b_box_end    = row[3]
				a_box_start  = row[4]
				a_box_end    = row[5]
				
				left_tsd_region  = polyAT_start - 40
				right_tsd_region = a_box_end + 41
				
				myseq = sequence[left_tsd_region:right_tsd_region]
				
				rseq.append('>'+sequence_id)
				rseq.append(myseq)

		return fseq, rseq
		'''

		del hits
				
		if f is not None:
			#print(f'{len(f)} forward hits in {sequence_id}')
			f_writeout = self.extract_match(sequence, f, desc, forward = True)
			#print(f'{int(len(f_writeout)/2)} fwd remain in {sequence_id}')
		else:
			f_writeout = None
		
		if r is not None:
			#print(f'{len(f)} reverse hits in {sequence_id}')
			r_writeout = self.extract_match(sequence, r, desc, forward = False)
			#print(f'{int(len(r_writeout)/2)} rev remain in {sequence_id}')
		else:
			r_writeout = None
		
		return f_writeout, r_writeout
		
def initialize_sine_worker(genome_file):
	global sinefinder
	sinefinder = hyperSINEfinder(genome_file)
	sinefinder.prep()
		
def worker(seqid):
	f, r = sinefinder.execute_search(seqid)
	return f, r, seqid
	
def jesus_give_me_a_SINE(genome_file, output = None, threads = 1):
	fa = pyfastx.Fasta(genome_file, build_index = True)
	#ids = list(range(len(fa)))
	ids = list(fa.keys())
	#ids = [0]
	
	ok_procs = min([threads, len(ids)])
	
	if output is None:
		import io
		out = io.StringIO('')
	else:
		out = open(output, 'w')
	
	with multiprocessing.Pool(ok_procs, initializer = initialize_sine_worker, initargs = (genome_file, )) as pool:
		fwd = []
		rev = []
		#for f, r, s in pool.map(worker, ids):				
		for f, r, s in pool.imap_unordered(worker, ids):				
			if f is not None:
				for rec in f:
					print(rec, file = out)
			if r is not None:
				for rec in r:
					print(rec, file = out)

	if output is None:
		output = out

	return output


#jesus_give_me_a_SINE(sys.argv[1], sys.argv[2], 10)


#oot = jesus_give_me_a_SINE(f, None, 10)

#print(oot.getvalue())
