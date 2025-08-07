import pyfastx
import mappy
import multiprocessing
import os

#Indices are loaded in advance of mapping each time, so this approach will not work.
#It really is going to have to be either a database approach or let it ride in main
def init_workers(genome_file, mmidx, max_hits, minprob):
	global fa
	global indices
	global maxhit
	global prob
	fa = pyfastx.Fasta(genome_file)
	indices = [mappy.Aligner(fn_idx_in = idx, best_n = max_hits) for idx in mmidx]
	prob = minprob
	
def paf2blast(record):
	pass

#Probably makes more sense to batch sequences to do several seq vs. aln in a row.
#Also makes sense to convert to blast format directly here, skip paf output
def one_sequence_map(seqids):
	'''
	"qname",
	"qlen",
	"qstart",
	"qend",
	"strand",
	"tname",
	"tlen",
	"tstart",
	"tend",
	"nmatch",
	"alen",
	"mapq",
	'NM', 'ms', 'AS', 'nn', 'tp', 'cm', 's1', 'de', 'rl', 'cg', 'extra'
	'''

	#CP002687.1_1	343	0	343	+	CP002687.1;;0	18585056	601935	602278	343	343	60	NM:i:0	ms:i:686	AS:i:686	nn:i:0	tp:A:P	cm:i:57	s1:i:312	s2:i:0	de:f:0	rl:i:0	cg:Z:343M
	for idx in indices:
		for seqid in seqids:
			sequence = fa[seqid].seq
			for i in range(0, 100):
				for hit in idx.map(sequence, MD=True):
					record = [seqid, len(sequence), 
							hit.q_st, hit.q_en, 
							hit.strand,
							hit.ctg, hit.ctg_len,
							hit.r_st, hit.r_en,
							hit.mlen, hit.blen,
							hit.mapq,
							hit.NM,
							hit.cigar_str]
					print(*record, sep = '\t')
			
	
def run(sequences_to_map, minimap_indices, output_file, maxhit = 50000, minprob=0.01, threads = 3):
	fa = pyfastx.Fasta(sequences_to_map, build_index = True)
	seqs = list(fa.keys())
	
	sequence_chunks = [seqs]
	
	with multiprocessing.Pool(threads, initializer = init_workers, initargs = (sequences_to_map, minimap_indices, maxhit, minprob, )) as pool:
		for result in pool.imap_unordered(one_sequence_map, sequence_chunks):
			pass

'''
indir = 'splitme'
for f in os.listdir(indir):
	#command = f'minimap2 -t 1 -k {k_size} -I 1G -d {mmidx} {infile}'
	comm = f'minimap2 -t 10 -k 10 -d {os.path.join(indir, f+'_mmidx')} {os.path.join(indir, f)}'
	os.system(comm)
'''

d = 'splitme'
indices = [os.path.join(d, f) for f in os.listdir(d) if f.endswith('_mmidx')]

#print(indices)

run('output/Step2_extend_blast_input_rename.fa', minimap_indices = indices, output_file = 'blah.txt')

