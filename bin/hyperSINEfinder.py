import sys
import os
import multiprocessing
import sqlite3
import pyfastx
import hyperscan
import numpy as np

from itertools import product

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
		
		#All possible SINE forward patterns exploded out to explicit regex with no ambiguity
		self.patterns = None
		self.pattern_reverser = None
		#Genome sequence to search
		self.sequence = None
		#Hyperscan database of regex patterns
		self.hsdb = None
	
	
	def revcomp(self, string):
		switcher = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
		revcmp = ''.join([switcher[c] if c in switcher else 'N' for c in string[::-1]])
		return revcmp
		
	def set_patterns(self):
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
	
	#okay, so there's a lot of math here to basically replicate the spacer-matching component of the original regex method
	
	def find_plausible_sine_patterns(self, hits, forward = True):			
		'''
		Where the original regex explicitly searched for spacers, we only have start/stop indices for hyperscan
		
		So, what we do here is we select instances where the correct ordering of SINE elements is found within a max distance of the 
		relevant spacer of the next OK element, maximizing distance between OK elements up to spacer size
		
		Achieved by
		Find instances where hit[ok_order[0]] end within spacers[0] or hit[ok_order[1]] start and
							 hit[ok_order[1]] end within spacers[1] or hit[ok_order[2]] start
							 
		We need to know where the nearest OK hit is for each item
		
		
		'''
		
		#Find rows in hits where regex pattern changes - this is cases where it's plausible a match could occur
		#We do not know which of a run of matches is the best candidate from this, only if it's possible that a run of candidates might be a match
		
		#Run length encoding of regex pattern matches
		#columns: 0 = RLE index, 1 = value, 2 = position, 3 = run_length
		rle_array = self.rle(hits[:, 0])
		
		#Index original hits with an RLE ID matching that in the rle_array, but repeated to the sizes of runs
		#Can we decode from positions instead? Probably
		rle_groups_hits = np.repeat(rle_array[:,0], rle_array[:,3])

		#filterer = np.zeros(rle_array.shape[0], dtype = np.bool_)
		
		#regex pattern order is SINE-like; either a_box -> b_box fwd OR polyAT -> b_box rev
		first_ok = (rle_array[:-2, 1] + 1) == rle_array[1:-1, 1]
		#regex pattern order 2 is SINE-like; either b_box -> polyAT fwd OR b_box -> a_box rev
		second_ok = (rle_array[1:-1, 1] + 1) == rle_array[2:, 1]
		
		#This finds the indices of first patterns where the next two patterns are correct
		correct_pattern_order = np.where(np.logical_and(first_ok, second_ok))[0]
		
		del first_ok
		del second_ok
		
		#Add the ensuing two indices so that we can select the entire runs of correct patterns
		correct_pattern_order = (correct_pattern_order[:, None] + np.arange(3)).flatten()
		
		ok_rle_ids = rle_array[correct_pattern_order, 0]
		
		del correct_pattern_order
		
		hits = hits[np.isin(rle_groups_hits, ok_rle_ids)]
		
		del rle_groups_hits
		del ok_rle_ids
		
		sine_candidates = []
		#It's possible no remaining SINE candidates survive, which would cause RLE to break and also be unnecessary
		if len(hits) > 0:
			rle_array = self.rle(hits[:, 0])
			rle_groups_hits = np.repeat(rle_array[:,0], rle_array[:,3])
			
			#Split by regex pattern index; this will always be 0, 1, 2, 0, 1, 2 .... or 2, 3, 4, 2, 3, 4... depending on forward or reverse
			#Because we are already splitting on the regex index, we do not need the pattern index in the data here
			groups = np.split(hits[:, 1:], np.unique(rle_groups_hits, return_index=True)[1][1:])
			#groups = np.split(hits, np.unique(rle_groups_hits, return_index=True)[1][1:])
			
			if forward:
				spacer_1_low = 25
				spacer_1_high = 50
				spacer_2_low = 20
				spacer_2_high = 500
			else:
				spacer_1_low = 20
				spacer_1_high = 500
				spacer_2_low = 25
				spacer_2_high = 50
			
			
			
			#Triplets of these groups represent a plausible SINE element
			for i in range(0, len(groups), 3):
				first_element = groups[i]
				second_element = groups[i+1]
				third_element = groups[i+2]

				#This one represented regex pattern + start + end indices
				#combinations = np.array(list(product(first_element, second_element, third_element))).reshape(-1, 9)
				
				#this one represents start and end indices only
				#Really, this could just be first_element[end], second element[start, end], third element[start]
				#But that's too much work, not better in any way
				combinations = np.array(list(product(first_element, second_element, third_element))).reshape(-1, 6)
				
				first_to_second_ok = combinations[:,2] - combinations[:, 1]
				second_to_third_ok = combinations[:,4] - combinations[:, 3]
				
				#Instances where the spacer distances are OK between the correctly ordered elements
				valid_sine_trio = np.logical_and(np.logical_and(first_to_second_ok >= spacer_1_low, first_to_second_ok <= spacer_1_high), 
												 np.logical_and(second_to_third_ok >= spacer_2_low, second_to_third_ok <= spacer_2_high))
				
				combinations = combinations[valid_sine_trio]
				
				#Implicitly - zero length hits are intrinsically skipped by not having length >= 1
				
				#Exactly one hit; it must be the longest hit by definition
				if combinations.shape[0] == 1:
					sine_candidates.append(combinations)
				#More than one hit; return all unique; 
				if combinations.shape[0] > 1:
					#We're doing the simple thing here and taking the single longest first_element start : last_element end; 
					#it is most likely to capture a real SINE element in its entirety
					#If tied longest elements, argmax prefers the first
					combinations = combinations[np.argmax(combinations[:,5] - combinations[:, 0])]
					
					sine_candidates.append(combinations)
		
		#Nice numpy arrays are nice
		if len(sine_candidates) > 0:
			sine_candidates = np.vstack(sine_candidates)
		else:
			sine_candidates = None
		
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
		
	def execute_search(self, sequence_id):
		record = self.fa[sequence_id]
		sequence = record.seq
		seqlen = len(sequence)
		desc = record.description
		hits = []
		self.hsdb.scan(sequence.encode(encoding='ascii'), match_event_handler=on_match, context = hits)
		
		f = None
		r = None
		
		if len(hits) > 0:
			#Definitely need to remove ungreedy polyAT hits
			f, r = self.clean_and_group_matches(hits)

		del hits
		
		f_writeout = []
		r_writeout = []

		if f is not None:
			for row in f:
				a_box_start = row[0]
				a_box_end = row[1]
				b_box_start = row[2]
				b_box_end = row[3]
				polyAT_start = row[4]
				polyAT_end = row[5]
				
				left_tsd_region = a_box_start - 40
					
				right_tsd_region = polyAT_end + 41
								
				tsd1 = sequence[left_tsd_region:a_box_start]
				a_box = sequence[a_box_start:a_box_end]
				spacer1 = sequence[a_box_end:b_box_start]
				b_box = sequence[b_box_start:b_box_end]
				spacer2 = sequence[b_box_end:polyAT_start]
				polyAT = sequence[polyAT_start:polyAT_end]
				tsd2 = sequence[polyAT_end:right_tsd_region]
				

				next_sine = (tsd1, a_box, spacer1, b_box, spacer2, polyAT, tsd2,)
				f_writeout.append(next_sine)
				
		
		if r is not None:
			for row in r:
				pass

				
		
		return desc, f_writeout, r_writeout
		
		
	
def initialize_sine_worker(genome_file):
	global sinefinder
	sinefinder = hyperSINEfinder(genome_file)
	sinefinder.prep()
		
def worker(seqid):
	id, f, r = sinefinder.execute_search(seqid)
	return id, f, r
	
def jesus_give_me_a_SINE(genome_file, output = None, threads = 1):
	fa = pyfastx.Fasta(genome_file, build_index = True)
	ids = list(range(len(fa)))
	
	ok_procs = min([threads, len(ids)])
	
	with multiprocessing.Pool(ok_procs, initializer = initialize_sine_worker, initargs = (genome_file, )) as pool:
		with open(output, 'w') as out:
			print('sequence', 'tsd1', 'a_box', 'spacer1', 'b_box', 'spacer2', 'polyAT', 'tsd2', sep = '\t', file = out)
			for id, f, r in pool.imap_unordered(worker, ids):
				for res in f:
					print(id, *res, sep = '\t', file = out)
				
		
f = '../TEtrimmer/try2/bTaeGut7v0.4_MT_rDNA.fa'
jesus_give_me_a_SINE(f, 'example_zebrafinch_sines.txt', 10)