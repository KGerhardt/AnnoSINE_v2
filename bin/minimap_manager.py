import sys
import os
import multiprocessing
import subprocess
import pyfastx
import polars as pl

import sqlite3

from genomeSplitter import genomeSplitter

from paf2_blast_stream import process_file


import re

'''
#genomesplitter usage
	mn = genomeSplitter(genome_file = genome_file, 
				output_directory = output_dir, 
				chunk_size = c,
				overlap_size = o,
				procs = p,
				smart = s,
				post_index = index_outputs,
				verbose = v,			
				overwrite = ow)
				
	output_files = mn.run()
'''

class minimap_manager:
	def __init__(self, reference_genome, query_sequences, working_directory, k_size = 10, threads = 1):
		self.rg = reference_genome
		self.qs = query_sequences
		
		self.wd = working_directory
		self.thds = threads
		
		self.mm2_k = k_size
		
		self.seqlen_dataframe = None
		
		self.paf_format = {
				"qname":pl.String,
				"qlen":pl.Int32,
				"qstart":pl.Int64,
				"qend":pl.Int64,
				"strand":pl.String,
				"tname":pl.String,
				"tlen":pl.Int64,
				"tstart":pl.Int64,
				"tend":pl.Int64,
				"nmatch":pl.Int32,
				"alen":pl.Int32,
				"mapq":pl.Int32,
				'NM':pl.String, 
				'ms':pl.String, 
				'AS':pl.String, 
				'nn':pl.String, 
				'tp':pl.String, 
				'cm':pl.String,
				's1':pl.String, 
				'de':pl.String, 
				'rl':pl.String, 
				'cg':pl.String, 
				'extra':pl.String
			}
			
		self.blast_format = {
			"qname":pl.String,
			"tname":pl.String,
			"percent_ident":pl.Float32,
			"alen":pl.Int32,
			"nonmatch":pl.Int32,
			"gap_openings":pl.Int32,
			"qstart":pl.Int32,
			"qend":pl.Int32,
			"tstart":pl.Int64,
			"tend":pl.Int64,
			"evalue":pl.Float32,
			"bitscore":pl.Float32
			}
		
		#chr1_179	301	2	114	-	chr17;;0	53461100	13048107	13048234	104	127	5	NM:i:23	ms:i:146	AS:i:130	nn:i:0	tp:A:P	cm:i:8	
		#s1:i:42	s2:i:45	de:f:0.1034	rl:i:10	cg:Z:5M1D37M11D3M2D16M1D51M
		
	def prep_dir(self, directory):
		if not os.path.exists(directory):
			os.makedirs(directory, exist_ok = True)
		
		
	def make_mm_index(self, args):
		file, index = args
		print(f'Indexing {file}...')
		minimap2_index_command = f'minimap2 -t {thread_chunk} -k {self.mm2_k} -d {index} {file}'
		minimap2_index_command = minimap2_index_command.split()
		subprocess.run(minimap2_index_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		os.remove(file)
		return index

		
	def prepare_genome(self):
		rg_splits = os.path.join(self.wd, 'minimap_wrangler', 'minimap_indices')
		
		completion_file = os.path.join(rg_splits, 'indexing_complete.txt')
		
		if not os.path.exists(completion_file):
			self.prep_dir(rg_splits)
			gs = genomeSplitter(genome_file = self.rg,
								output_directory = rg_splits,
								chunk_size = 500_000_000,
								overlap_size = 1_000_000,
								procs = self.thds,
								verbose = False,
								overwrite = False)
			gs.run()
			
			seqlen_dict = gs.seqs_and_lengths
			
			os.remove(os.path.join(rg_splits, 'genomeSplitter.log'))
			
			chunks = []
			
			args = []
			print('Creating minimap2 indices...')
			for f in os.listdir(rg_splits):
				if f.endswith('.fasta'):
					genome_chunk = os.path.join(rg_splits, f)
					idx_file = f.replace('.fasta', '.minimap2.idx')
					index = os.path.join(rg_splits, idx_file)
					args.append((genome_chunk, index,))
				
			
			with multiprocessing.Pool(max([self.thds // thread_chunk, 1])) as pool:
				for mmidx in pool.imap_unordered(self.make_mm_index, args):
					chunks.append(mmidx)
				
			out = open(completion_file, 'w')
			out.close()
			
		else:
			gs = genomeSplitter(genome_file = self.rg,
								output_directory = rg_splits,
								chunk_size = 500_000_000,
								overlap_size = 1_000_000,
								procs = self.thds,
								verbose = False,
								overwrite = False)
								
			gs.get_seqlens()
			seqlen_dict = gs.seqs_and_lengths
			chunks = [os.path.join(rg_splits, f) for f in os.listdir(rg_splits) if f.endswith('minimap2.idx')]
			
		#print(seqlen_dict)
		sld = {'tname':[], 'tlen':[]}
		for seq in seqlen_dict:
			sld['tname'].append(seq)
			sld['tlen'].append(seqlen_dict[seq])
		
		self.seqlen_dataframe = pl.from_dict(sld, schema = {'tname':pl.String, 'tlen':pl.Int64})
		self.seqlen_dataframe = self.seqlen_dataframe.lazy()
						
		return chunks
		
	def prepare_queries(self, query_group_size = 100):
		query_splits = os.path.join(self.wd, 'minimap_wrangler', 'split_queries')
		
		self.prep_dir(query_splits)
		split_check = os.path.join(query_splits, 'pyfastx_split_complete.txt')
		
		if not os.path.exists(split_check):
			query_base = os.path.basename(self.qs)
			split_index = 1
			counter = 0
			with open(self.qs) as inf:
				dat = inf.read()
				
			dat = dat.split('>')
			name = os.path.join(query_splits, f'{query_base}.{split_index}.fa')
			out = open(name, 'w')
			for record in dat:
				out.write(f'>{record}')
				counter += 1
				if counter == query_group_size:
					counter = 0
					out.close()
					split_index += 1
					name = os.path.join(query_splits, f'{query_base}.{split_index}.fa')
					out = open(name, 'w')
					
			out.close()

			out = open(split_check, 'w')
			out.close()

		chunks = [os.path.join(query_splits, f) for f in os.listdir(query_splits) if f != 'pyfastx_split_complete.txt']
		
		return chunks
		
		
	def manage_minimap_alignment(self, pair):
		query, target = queries[pair[0]], targets[pair[1]]
		
		qname = os.path.basename(query)
		tname = os.path.basename(target)
		
		#Run alignment
		#print(f'Searching {qname} vs {tname}')
		my_outfile = os.path.join(working_outputs, f'{qname}_vs_{tname}.paf')
		command = f'minimap2 -c -t {thread_chunk} {target} -p 0.01 -N 10000 {query} -o {my_outfile}'
		command = command.split()
		subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		
		#Convert minimap to BLAST format
		blast_out = f'{my_outfile}.blast'
		#print(f'Converting {my_outfile} to blast')
		success = process_file(my_outfile, blast_out)
		
		#Remove minimap aln
		os.remove(my_outfile)
		
		if not success:
			blast_out = None

		return query, blast_out
			
	#Post-blast conversion
	def merge_query_chunks_blast(self, query_file, query_vs_target_alns):
		qname = os.path.basename(query_file)
		comb_out = os.path.join(working_outputs, f'{qname}_alignments.blast.txt')
		cleaned_q_v_t = [f for f in query_vs_target_alns if f is not None]
		#Concat + sort PAF with polars
		
		if len(cleaned_q_v_t) > 0:
			#Target length will always be wrong for long chromosomes > 500 milion bp, but that's OK we don't use it those
			my_targets = pl.scan_csv(cleaned_q_v_t,
									has_header = False,
									separator = "\t",
									schema = self.blast_format)
			
			#Sort results
			my_targets = my_targets.sort(by = ['qname', 'tname', 'tstart'])

			#Deduplicate rows that appeared in overlap
			my_targets = my_targets.unique(subset=['qname', 'tname', 'tstart'], maintain_order=True)

			try:
				#Write output
				my_targets.sink_csv(path = comb_out,
									include_header = False,
									separator = '\t')
									
				
					
				megabytes = 1024 * 1024
				if os.path.getsize(comb_out) > 500 * megabytes:
					new_output_list = []
					max_lines = 5_000_000 #500 MB/output file
					partial_index = 0
					line_count = 0
					file_handle = os.path.join(working_outputs, f'{qname}_{partial_index}_alignments.blast.txt')
					new_output_list.append(file_handle)
					outfile = open(file_handle, 'w')
					with open(comb_out) as inf:
						for line in inf:
							ele=line.split()
							sid = ele[0]
							
							if line_count >= max_lines:
								# Close the current file and open a new one with an incremented suffix

								#Continue writing lines beyond max until there's a new record so that records are grouped
								if sid==prev:
									outfile.write(line)
									continue

								outfile.close()
								partial_index += 1
								file_handle = os.path.join(working_outputs, f'{qname}_{partial_index}_alignments.blast.txt')
								outfile = open(file_handle, "w")
								new_output_list.append(file_handle)
								line_count = 0  # Reset line count for the new file
							
							outfile.write(line)
							line_count += 1
							prev=sid
							
					outfile.close()
					
					os.remove(comb_out)
					
					comb_out = new_output_list
				else:
					comb_out = [comb_out]
							
			except Exception as e:
				print(f'{comb_out} could not write an output for reason {e}')
				print('This probably is not a problem')
				#Usually the result of an empty file
				comb_out = None
				
		else:
			comb_out = None
		
		#Clean up partial alignments
		for o in cleaned_q_v_t:
			os.remove(o)
			
		return comb_out
			
			
	def run(self):
		thread_group_size = 4
		
		global thread_chunk
		
		if self.thds <= thread_group_size:
			thread_chunk = threads
		else:
			thread_chunk = thread_group_size
	
		global working_outputs
		working_outputs = os.path.join(self.wd, 'minimap_wrangler', 'partial_alignment_results')
		self.prep_dir(working_outputs)
	
		global targets
		#Get minimap indices
		targets = self.prepare_genome()
		num_targets = len(targets)
		
		global queries
		#Get queries split into groups of some size
		#How many files could this be? Potentially thousands, which isn't good
		#Is there a good way too split out a few at a time?
		queries = self.prepare_queries(100)
		query_record = {}
		for q in queries:
			query_record[q] = []
		
		args = []
		for i in range(0, len(queries)):
			for j in range(0, num_targets):
				args.append((i, j,))
		
		
		final_results = []
		#We probably want to change parallelization strategy to search one query chunk against all refgen chunks, combine, then move to the next batch of chunks
		with multiprocessing.Pool(max([self.thds // thread_chunk, 1])) as pool:
			for qname, partial_out in pool.imap_unordered(self.manage_minimap_alignment, args):
				query_record[qname].append(partial_out)
				if len(query_record[qname]) == num_targets:
					#print(f'Merging partial results for {qname}')
					output = self.merge_query_chunks_blast(qname, query_record[qname])
					removed = query_record.pop(qname)
					removed = None
					if output is not None:
						final_results.extend(output)
					
		final_results.sort()
					
		return final_results
		
def run_map(ref_genome, queries, output_dir, mmseqs_k, threads = 1):
	mn = minimap_manager(ref_genome, queries, output_dir, mmseqs_k, threads)
	results = mn.run()
	
	return results

'''
g = sys.argv[1]
q = sys.argv[2]
w = sys.argv[3]
k = 10
t = int(sys.argv[4])

mn = minimap_manager(g, q, w, k, t)
mn.run()
'''