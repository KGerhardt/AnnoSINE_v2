import sys
import os
import multiprocessing
import subprocess
import pyfastx
import polars as pl

import sqlite3

from genomeSplitter import genomeSplitter

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
		pfx = f'pyfastx split {self.qs} -c {query_group_size} --out-dir {query_splits}'
		pfx = pfx.split()
		subprocess.run(pfx)
		chunks = [os.path.join(query_splits, f) for f in os.listdir(query_splits)]
		
		return chunks
		
		
	def manage_minimap_alignment(self, pair):
		query, target = queries[pair[0]], targets[pair[1]]
		
		qname = os.path.basename(query)
		tname = os.path.basename(target)
		
		print(f'Searching {qname} vs {tname}')
		my_outfile = os.path.join(working_outputs, f'{qname}_vs_{tname}.paf')
		command = f'minimap2 -c -t {thread_chunk} {target} -p 0.01 -N 10000 {query} -o {my_outfile}'
		command = command.split()
		subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		
		return query, my_outfile
		
	def merge_query_chunks(self, query_file, query_vs_target_alns):
		qname = os.path.basename(query_file)
		comb_out = os.path.join(working_outputs, f'{qname}_alignments.paf')
		#Concat + sort PAF with polars
		
		#Target length will always be wrong for long chromosomes > 500 milion bp, but that's OK we don't use it those
		my_targets = pl.scan_csv(query_vs_target_alns,
								has_header = False,
								separator = "\t",
								schema = self.paf_format)
		
		#Extract the target sequence offset
		my_targets = my_targets.with_columns(
			pl.col("tname")
			.str.extract(r';;(\d+)', 1)  # Extract the first match of digits
			.cast(pl.Int64)            # Convert to integer
			.alias("offset")
		)

		#Adjust alignment indices accordingly, truncate target name to its non-offset prefix
		my_targets = my_targets.with_columns(
			pl.col("tname").str.replace(r';;\d+', '').alias('tname'),
			(pl.col('offset') + pl.col('tstart')).alias('tstart'),
			(pl.col('offset') + pl.col('tend')).alias('tend')
		)

		#Drop the now superfluous offset indicator
		my_targets = my_targets.drop('offset')

		#Correct sequence lengths to match original genome lengths
		my_targets = my_targets.update(self.seqlen_dataframe, on = 'tname', how = 'left')
		
		#Sort results
		my_targets = my_targets.sort(by = ['qname', 'tname', 'tstart'])

		#Deduplicate rows that appeared in overlap
		my_targets = my_targets.unique(subset=['qname', 'tname', 'tstart'], maintain_order=True)

		#Write output
		my_targets.sink_csv(path = comb_out,
							include_header = False,
							separator = '\t')
		
		#Clean up partial alignments
		for o in query_vs_target_alns:
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
		queries = self.prepare_queries(50)
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
					print(f'Merging partial results for {qname}')
					output = self.merge_query_chunks(qname, query_record[qname])
					removed = query_record.pop(qname)
					removed = None
					final_results.append(output)
					
		final_results.sort()
					
		return final_results
		
		
def run_map(ref_genome, queries, output_dir, mmseqs_k, threads = 1):
	mn = minimap_manager(g, q, w, k, t)
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