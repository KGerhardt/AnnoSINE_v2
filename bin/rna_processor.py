import sys
import os
import subprocess

def blast_rna(out_genome_assembly_path, script_dir, threads = 1):
	parent = os.path.dirname(script_dir)
	q = os.path.join(out_genome_assembly_path, 'Step4_rna_input.fasta')
	rna_database = os.path.join(parent, 'Input_Files', 'blast_database', 'rna_database.fa')
	o = os.path.join(out_genome_assembly_path, 'Step4_rna_output.out')
	
	query_is_ok = os.path.getsize(q) > 0
	if query_is_ok:
		command = f'blastn -query {q} -db {rna_database} -out {o} -evalue 1 -num_alignments 50000 -word_size 7 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -outfmt 6 -num_threads {threads}'
		#print(command)
		command = command.split()
		subprocess.run(command, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
	else:
		print('There are no sequences to check for RNA heads')
	
	return query_is_ok


def process_rna_better(out_genome_assembly_path, script_dir):
	parent = os.path.dirname(script_dir)
	rna_label_file = os.path.join(parent, 'Input_Files', 'rna_labels.txt')
	query_input = os.path.join(out_genome_assembly_path, 'Step4_rna_input.fasta')
	
	#Tabular BLAST outputs in -outfmt 6
	input_f = os.path.join(out_genome_assembly_path, 'Step4_rna_output.out')
	
	#Original query sequences
	output_f = os.path.join(out_genome_assembly_path, 'Step4_rna_output.fasta')
	
	#RNA labels is a file that comes with AnnoSINE - it is a tab-sep file which contains the appropriate label by RNA reference sequence
	rna_labels = {}
	with open(rna_label_file) as fh:
		for line in fh:
			segs = line.strip().split('\t')
			target = segs[0]
			label = segs[1]
			rna_labels[target] = label
	
	#e-values: The goal of the program is to find hits of intermediate quality, i.e. those which aren't exact matches but aren't entirely awful either
	#BLAST above handles filtering e > 1; this code also removes e < 1e-15
	#As we read through the alignments which remain, we record the label associated with the best non-removed match by lowest e-value
	keepers = {}
	with open(input_f) as fh:
		for line in fh:
			segs = line.strip().split('\t')
			query = segs[0]
			target = segs[1]
			e = float(segs[10])
			if e >= 1e-15:
				if target in rna_labels:
					label = rna_labels[target]
				else:
					label = 'Unknown'
					
				if query not in keepers:
					keepers[query] = {'e':e, 'label':label}
				else:
					if e < keepers[query]['e']:
						keepers[query] = {'e':e, 'label':label}
	
	#Open the query sequences and output
	do_write = False
	with open(output_f, 'w') as out:
		with open(query_input) as fh:
			for line in fh:
				if line.startswith('>'):
					line = line.strip()
					seqid = line.split()[0][1:]
					#If the sequence had an 1e-15 <= E <= 1 match, write the record with the new label attached and its sequence
					if seqid in keepers:
						do_write = True
						print(f'{line}|{keepers[seqid]['label']}', file = out)
					#Else skip the line and its sequence
					else:
						do_write = False
				else:
					if do_write:
						out.write(line)
					
#script_dir = os.path.dirname(os.path.abspath(__file__))
#outdir = 'output_newopts'
 
#blast_rna(outdir, script_dir, threads = 10)
#process_rna_better(outdir, script_dir)
