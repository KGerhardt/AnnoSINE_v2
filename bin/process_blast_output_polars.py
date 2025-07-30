import multiprocessing
import sys
import os

import polars as pl
import numpy as np

import pyfastx

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

#the equivalent to the original AnnoSINE_v2's process_blast_output_1 function
#The core behavior of this script
def clean_blast_record(args):
	#Define file format for the reader
	blast_schema = {"qname":pl.String,
				"tname":pl.String,
				"percent_ident":pl.Float32,
				"alen":pl.Int32,
				"nonmatch":pl.Int32,	
				"gap_openings":pl.Int32,
				"qstart":pl.Int32, #These are with reference to short query sequences, 32 bit int is fine
				"qend":pl.Int32,
				"tstart":pl.Int64, #These are with reference to genome position, 64 bit is required
				"tend":pl.Int64,
				"evalue":pl.Float32, #We don't even use most of these, so size isn't an issue
				"bitscore":pl.Float32}

	these_cols = [0, 1, 3, 6, 7]

	file, factor_copy, factor_length, min_copy_number, max_gap, max_shift = args

	df = pl.read_csv(file, has_header = False, separator = '\t', schema = blast_schema, columns = these_cols)
	
	#Instead of trying to clean at the paf2blast outputs, we filter unique over hits here because these files are made to be contiguous
	#Ah, but the paf2blast has overlaps so the end of a previous record is possibly > the start of the next record.
	#So... how do we handle that?
	#df = pl.df.unique(subset['qname', 'tname', 'tstart', 'tend'], maintain_order = True)
	

	#BLAST is 1-indexed, fix this for python by subtracting 1 from all starts
	df = df.with_columns(
		(pl.col('qstart') - 1).alias('qstart')
	)
	
	#Select only the records which are an acceptable percent of their master record's length
	df = df.filter(
		(abs(pl.col("qend") - pl.col('qstart') + 1) >= factor_length * (pl.col('alen').first() - 200))
		.over(pl.col('qname'))
	)
	
	family_counts = df.group_by('qname').agg(pl.max('qend'), pl.len(), pl.col('alen').first())
	#Figure out which queries have enough counts
	family_counts = family_counts.filter(pl.col('len') >= min_copy_number)
	#Get the repeat numbers for every family
	family_counts = family_counts.with_columns((pl.col('len') * factor_copy).ceil().cast(pl.Int32).alias('repeat_number'))

	#Reduce kept data
	df = df.select(['qname', 'tname', 'qstart', 'qend'])
	
	#Filter dataframe to only those records with enough family count to continue
	df = df.join(other = family_counts, on='qname', how = 'inner').select(df.columns)
	
	updates = []
	
	#Surprisingly time intensive data op here.
	for seqid, data in df.select(['qname', 'qstart', 'qend']).group_by('qname'):
		seqid = seqid[0] #comes as a tuple
		
		#Select the metadata for this ID
		getme = family_counts.filter(pl.col('qname') == seqid).select('repeat_number', 'qend', 'alen')
		#Metadata for this record
		repeat_no, maxval, record_length = getme.item(0, 0), getme.item(0, 1), getme.item(0, 2)
		
		hit_count = len(data)
		
		#Get the depth of coverage by position in the query sequence
		arr = np.zeros(maxval, dtype = np.int32)
		for row in data.iter_rows():
			arr[row[1]:row[2]] += 1
				
		#Discount items where the coverage depth is inadequate
		arr[arr < repeat_no] = 0
		
		#Group by runs of nonzero values
		nonzero = arr != 0
		#Get runs of 0 or nonzero
		runs, position_starts, values = rle(nonzero)
		#Subset runs and positions to nonzero only (values is an array of true/false)
		runs = runs[values]
		position_starts = position_starts[values]
		
		#Calculate ends
		position_ends = position_starts + runs - 1
		
		num_gaps = position_starts.shape[0] - 1
		
		#If there is more than one contiguous group of sequences, clean the position start and end records
		if num_gaps > 0:
			'''
			The code here is designed to handle cases where there are positions in a pileup of BLAST alignments with insufficient coverage to be included
			e.g. [1, 2, 3, 3, 3, 2, 4, 4, 4, 2] with acceptable coverage = 3 would result in a cleaned array of
				 [0, 0, 3, 3, 3, 0, 4, 4, 4, 0], 
			
			or, interpreted differently, as runs of acceptable values from indices [2, 5) and [6, 9) with a gap of size 1 at position 5
			
			How this is to be handled:
			
			(1) In cases with multiple gaps, adjacent records with acceptable gap sizes should be considered as one contiguous group
			(2) If there is only one remaining contiguous group after the above logic, then the start and end of that group are used to update records
			(3) If there are still multiple contiguous groups left, the largest contiguous group according to total number of positions covered not including gaps is kept, 
				other groups are discarded
			(4) As a result of (2) or (3), the sole remaining contiguous group's start and end are compared to the whole length of the query sequence as 
				dist_to_start = (group_start - query_start) and dist_to_end = (query_end - group_end)
			(5) Both dist_to_start and dist_to_end must be <= max_shift to update the record with the new start + end values
			
			Note - query_start is always zero because the group start is already the position relative to the start of the query.
			'''
			
			#Find contiguous groupings separated by <= gap_size
			#Gap sizes betweed adjacent records have to be less than max_gap to be included
			gap_sizes = position_starts[1:] - position_ends[0:-1]
			#Initialize a dict of contiguous groups; pre-load group 0 as the first member of the first grouping
			groupings = {0:[0]}
			current_group = 0
			#Group together adjacent records with gaps <= max_gap
			for i in range(1, position_starts.shape[0]):
				if position_starts[i] - position_ends[i-1] <= max_gap:
					#If the current chunk is close enough to the previous one, add it to the previous one.
					groupings[current_group].append(i)
				else:
					#Otherwise, start a new group with the current chunk as the first member.
					current_group += 1
					groupings[current_group] = [i]
			
			#If all groups are contiguous, then there is only one group in groupings and the below purge code is not needed.
			#If there are more than one contiguous groups remaining, identify the longest one
			if len(groupings) > 1:
				#Calculate number of covered bases less gaps for each contiguous group
				contiguous_group_num_bases_covered = {}
				for g in groupings:
					total_covered_bases = 0
					for i in groupings[g]:
						total_covered_bases += (position_ends[i] - position_starts[i] + 1)
					contiguous_group_num_bases_covered[g] = total_covered_bases
					
				#The winning group is whichever group has the most total covered bases in the query
				winning_group = max(contiguous_group_num_bases_covered, key=contiguous_group_num_bases_covered.get)
				
				#Select only the start and end records in the winning group; the others are removed
				position_starts = [position_starts[i] for i in groupings[winning_group]]
				position_ends = [position_ends[i] for i in groupings[winning_group]]

	
		#The first position is always zero, so there is no math needed here
		left_shift = position_starts[0]
		#Length of the query - final position covered adequately
		right_shift = record_length - position_ends[-1]
		
		if left_shift < max_shift and right_shift < max_shift:
			updates.append((seqid, int(position_starts[0]), int(position_ends[-1]), hit_count, record_length))
		else:
			updates.append((seqid, None, None, hit_count, record_length))
	
	return file, updates

#Parallel processing for retrieval of many genome subsequences in little memory.
def ref_sequences_to_partial_file(args):
	index, seqid, ref_genome, my_output, data = args
	
	ref = pyfastx.Fasta(ref_genome, build_index = True)
	
	sequence = ref[seqid]
	seqlen = len(sequence)
	
	with open(my_output, 'w') as outfasta:
		for row in data.iter_rows():
			#Reference name goes unused
			full_description, query_name, classifier, old_start, old_end, tsd_start, tsd_end, start_adjustment, end_adjustment, blast_hit_ct, record_length, reference_name_dupe = row
			if start_adjustment is not None:
				tsd_start = start_adjustment
				tsd_end = end_adjustment
				
				#Sanity fixes
				if tsd_start < 0:
					tsd_start = 0
				if tsd_end > seqlen:
					tsd_end = seqlen
					
			if blast_hit_ct is None:
				blast_hit_ct = 0
			if record_length is None:
				record_length = 0
								
			header = f'>{full_description}|blast_s:{tsd_start}|blast_e:{tsd_end}|blast_count:{blast_hit_ct}|blast_l:{record_length}'
			subseq = sequence[tsd_start:tsd_end]
			
			print(header, file = outfasta)
			print(subseq, file = outfasta)
			
	return index, my_output

#The manager function of this script
def get_updates_from_blasts(output_directory, in_genome_assembly_path, factor_copy = 0.15, factor_length= 0.3, min_copy_number = 1, max_gap = 10, max_shift = 50, threads = 1):
	
	changes_file = os.path.join(output_directory, 'Step_3_tsd_update_ledger.txt')
	completion_marker = os.path.join(output_directory, 'Step_3_blast_processing_log.txt')
	
	input_files = [f for f in os.listdir(output_directory) if 'Step3_blast_output.out_part' in f and not f.endswith('.fa')]
	completed_files = []
	
	#Check previously processed blast inputs
	if os.path.exists(completion_marker):
		with open(completion_marker) as inf:
			for line in inf:
				completed_record = os.path.basename(line.strip())
				completed_files.append(completed_record)
		
		#Purge records that have already been processed			
		completed_files = set(completed_files)
		input_files = [f for f in input_files if f not in completed_files]
		
	#Process remaining blast inputs
	if len(input_files) > 0:
		print('Processing BLAST outputs for changes to TSDs')
		
		commands = [(os.path.join(output_directory, file), factor_copy, factor_length, min_copy_number, max_gap, max_shift) for file in input_files]
		remaining_chunks = len(commands)
		
		print(f'{remaining_chunks} BLAST chunks left to process.')
		
		log = open(completion_marker, 'a')
		with open(changes_file, 'w') as out:
			with multiprocessing.Pool(threads) as pool:
				for r, update_tuples in pool.imap_unordered(clean_blast_record, commands):
					for u in update_tuples:
						print(*u, sep = '\t', file = out)
					remaining_chunks -= 1
					print(f'{r} has been processed! {remaining_chunks} BLAST outputs remain.')
					print(os.path.basename(r), file = log)
	else:
		print('No remaining BLAST outputs to be processed;')
		print('\tYou are probably seeing this because script was previously run and results have simply been loaded.')
		
	#Spacing
	print('')
		
	#Reference genome assemble
	refgen = os.path.normpath(in_genome_assembly_path)
	if not os.path.exists(f'{refgen}.fxi'):
		#This process is pretty quick, a few sec/GB of reference genome
		print('Indexing reference genome...')
	ref = pyfastx.Fasta(refgen, build_index = True)
	
	#The TSD sequences to collect
	sequences_to_update = os.path.join(output_directory, 'Step2_extend_blast_input_rename.fa')
	#Collect update sequences
	if not os.path.exists(f'{sequences_to_update}.fxi'):
		print('Indexing TSD sequences...')
	else:
		print('Collecting TSD sequence information')
	
	fa = pyfastx.Fasta(sequences_to_update, build_index = True)

	#The fasta headers are formatted somewhat strangely - this fixes them
	descriptions = []
	for sequence in fa:
		full_desc = sequence.description
		seqid, classifier, tsd_info = full_desc.split()
		tsd_groups = tsd_info.split('|')
		tsd_groups = [t.split(':') for t in tsd_groups]
		start, end = int(tsd_groups[0][0]), int(tsd_groups[0][1])
		tsd_l = int(tsd_groups[1][1])
		tsd_start = int(tsd_groups[2][1])
		tsd_end = int(tsd_groups[3][1])
		next_item = (full_desc, seqid, classifier, start, end, tsd_start, tsd_end,)
		descriptions.append(next_item)
	
	#Convert list of tuples to polars dataframe
	descriptions_schema = {'full_desc':pl.String,
							'qname':pl.String,
							'classifier':pl.String,
							'start':pl.Int64,
							'end':pl.Int64,
							'tsd_start':pl.Int64,
							'tsd_end':pl.Int64}
							
	descriptions = pl.DataFrame(descriptions, schema = descriptions_schema, orient = 'row')
		
	print('Loading updates to make...')
	#Load the record of changes to make from blast file processing
	changes_schema = {'qname':pl.String,
				'start':pl.Int32,
				'end':pl.Int32,
				'hit_count':pl.Int32,
				'record_length':pl.Int32}
					
	#null_values lets polars interpret Python None writeouts from processing step as nulls, which have meaning to polars
	changes = pl.read_csv(changes_file, has_header = False, separator = '\t', schema = changes_schema, null_values = 'None')
	
	#Join the dataframes to assosciate update data with original data
	descriptions = descriptions.join(other = changes, on = 'qname', how = 'left')
	
	#Free space
	del changes
	
	#Info for printout
	num_tsds = len(descriptions)
	
	#Get the name of the reference sequence for grouping purposes
	descriptions = descriptions.with_columns(
					pl.col('qname').str.replace(r'_\d+$', '').alias('reference_sequence')
					)
	
	#Usually the same as the number of contigs in the reference genome, but could possibly differ if no TSDs found in one contig.
	total_sequences_to_parse = descriptions.select(pl.col('reference_sequence').n_unique()).item()
	
	#Prepare for parallel procesing
	idx = 0
	tsds_by_file = {}
	sequence_retrieval = []
	for refseq, data in descriptions.group_by('reference_sequence'):
		refseq = refseq[0] #Polars encodes group by as tuples in case there are multiple cols used
		output = os.path.join(output_directory, f'Step3_blast_process_output_chunk_{idx}.fa')
		#index, seqid, ref_genome, data
		next_arg = (idx, refseq, refgen, output, data)
		tsds_by_file[idx] = len(data)
		sequence_retrieval.append(next_arg)
		idx += 1

	final_output = os.path.join(output_directory, 'Step3_blast_process_output.fa')
	done_count = 0
	print(f'Making sequence updates. {total_sequences_to_parse} reference genome sequences to process with {num_tsds} total subseqeunces to retrieve.')
	#Process data chunks in parallel; load one contig from the reference genome per parallel process and slice it to each request
	with open(final_output, 'wb') as outfile:
		with multiprocessing.Pool(threads) as pool:
			for index, result_file in pool.imap_unordered(ref_sequences_to_partial_file, sequence_retrieval):
				#Fast file concatenate results to a single final output
				with open(result_file, 'rb') as infile:
					shutil.copyfileobj(infile, outfile)
				
				#Cleanup
				os.remove(result_file)
				done_count += 1
				processed_tsds = tsds_by_file[index]
				num_tsds -= processed_tsds
				#User reports
				print(f'{done_count} genome fragments of {total_sequences_to_parse} complete. {processed_tsds} sequences retrieved in this batch, {num_tsds} remain.')
					
