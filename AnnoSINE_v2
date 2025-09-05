import time
import os
import math
import re
from collections import Counter
import argparse
from argparse import RawTextHelpFormatter
import matplotlib.pyplot as plt
import operator
import glob
import uuid
import concurrent.futures

from paf2blast6_chunking2_parallel3_ordered_parts import parallel_process_file

import multiprocessing
import subprocess
import pyfastx
import sqlite3
import shutil
from genomeSplitter import genomeSplitter
from process_blast_output_polars import get_updates_from_blasts

print('Example: python3 AnnoSINE.py 2 ../Input_Files/test.fasta ../Output_Files', flush=True)
parser = argparse.ArgumentParser(description="SINE Annotation Tool for Plant Genomes",
								 formatter_class=RawTextHelpFormatter)

# positional arguments
parser.add_argument("mode", type=int,
					help="[1 | 2 | 3]\n"
					"Choose the running mode of the program.\n"
					"\t1--Homology-based method;\n"
					"\t2--Structure-based method;\n"
					"\t3--Hybrid of homology-based and structure-based method.")
parser.add_argument("input_filename", type=str, help="input genome assembly path")
parser.add_argument("output_filename", type=str, help="output files path")

# optional arguments
parser.add_argument("-e", "--hmmer_evalue", metavar='',type=float, default=1e-10,
					help="Expectation value threshold for saving hits of homology search (default: 1e-10)")
parser.add_argument("-v", "--blast_evalue", metavar='',type=float, default=1e-10,
					help="Expectation value threshold for sequences alignment search (default: 1e-10)")
parser.add_argument("-l", "--length_factor", metavar='', type=float, default=0.3,
					help="Threshold of the local alignment length relative to the the BLAST query length (default: 0.3)")
parser.add_argument("-c", "--copy_number_factor", metavar='', type=float, default=0.15,
					help="Threshold of the copy number that determines the SINE boundary (default: 0.15)")
parser.add_argument("-s", "--shift",  metavar='', type=int, default=50,
					help="Maximum threshold of the boundary shift (default: 50)")
parser.add_argument("-g", "--gap",  metavar='', type=int, default=10,
					help="Maximum threshold of the truncated gap (default: 10)")
parser.add_argument("-minc", "--copy_number", metavar='', type=int, default=1,
					help="Minimum threshold of the copy number for each element (default: 20)")
parser.add_argument("-numa", "--num_alignments", metavar='', type=int, default=50000,
					help="--num_alignments value for blast alignments (default: 50000)")

#parser.add_argument("-maxb", "--base_copy_number", type=int, default=1,
					#help="Maximum threshold of copy number for the first and last base (default: 1)")
parser.add_argument("-a", "--animal",metavar='', type=int, default=0,
					help='If set to 1, then Hmmer will search SINE using the animal hmm files from Dfam. If set to 2, then Hmmer will search SINE using both the plant and animal hmm files. (default: 0)')
parser.add_argument("-b", "--boundary", metavar='', type=str, default='msa',
					help="Output SINE seed boundaries based on TSD or MSA (default: msa)")
parser.add_argument("-f", "--figure", metavar='', type=str, default='n',
					help="Output the SINE seed MSA figures and copy number profiles (y/n). Please note that this step may take a long time to process. (default: n)")
parser.add_argument("-temd", "--temp_dir", metavar='',type=str, default='',
					help="The temp dir used by paf2blast6 script. If not set, will use /tmp folder automatically.")
parser.add_argument("-auto", "--automatically_continue", action = 'store_true',
					help="If set, the program will skip finished steps and continue unifinished steps for a previously processed output dir.")
parser.add_argument("-r", "--non_redundant", metavar='', type=str, default='y',
					help="Annotate SINE in the whole genome based on the non-redundant library (y/n) (default: y)")
parser.add_argument("-t", "--threads", metavar='', type=int, default=36,
					help="Threads for each tool in AnnoSINE (default: 36)")
parser.add_argument("-irf", "--irf_path", metavar='', type=str, default='',
					help="Path to the irf program (default: '')")
parser.add_argument("-rpm", "--RepeatMasker_enable", metavar='', type=int, default=1,
					help="If set to 0, then will not run RepearMasker (Step 8 for the code). (default: 1)")
					
parser.add_argument('--stop_at_tsd', action='store_true',
					help='Run HMMsearch, SINEfinder, and TSD searcher, then exit. Useful to developers.')
args = parser.parse_args()

# obtain program paths
script_dir = os.path.dirname(os.path.abspath(__file__)) #shujun
#work_dir = os.getcwd()
work_dir=args.output_filename

def hmm_worker(args):
	splitarg = args.split()
	subprocess.run(splitarg)
	return None

def hmm_predict(genome_assembly_path, cpus, script_dir, work_dir,input_ani,input_hmm_e_value): #shujun
	if input_ani==0:
		db='Family_Seq'
	elif input_ani==1:
		db='Dfam_hmm'
	elif input_ani==2:
		db1='Family_Seq'
		db2='Dfam_hmm'
	if not input_ani==2:
		dir_hmm = os.listdir(script_dir + '/../'+db+'/') #shujun
	else:
		d1=os.listdir(script_dir + '/../'+db1+'/')
		d2=os.listdir(script_dir + '/../'+db2+'/')
		dir_hmm=d1+d2
	#print(genome_assembly_path)
	#exit()

	os.system('mkdir ' + work_dir + '/HMM_out > /dev/null 2>&1') #shujun
	command_args = []
	for num_dir_hmm in range(len(dir_hmm)):
		if dir_hmm[num_dir_hmm] != '.DS_Store':
			# Clear the content of exist output
			if os.path.exists(work_dir + '/HMM_out/' + dir_hmm[num_dir_hmm] + '.out'): #shujun
				clear_filename = work_dir + '/HMM_out/' + dir_hmm[num_dir_hmm] + '.out' #shujun
				with open(clear_filename, "r+") as clear_f:
					clear_f.seek(0)
					clear_f.truncate()
			if not input_ani==2:
				command = f'nhmmer --cpu 4 -E {input_hmm_e_value} -o {work_dir}/HMM_out/{dir_hmm[num_dir_hmm]}.out {script_dir}/../{db}/{dir_hmm[num_dir_hmm]}/{dir_hmm[num_dir_hmm]}.hmm {genome_assembly_path}'
				#os.system(
				#	# output to work directory instead of program directory, shujun
				#	'nhmmer --cpu ' + str(cpus) +' -E '+str(input_hmm_e_value)+ ' -o ' + work_dir + '/HMM_out/' + dir_hmm[num_dir_hmm] + '.out ' 
				#	+ script_dir + '/../'+db+'/' + dir_hmm[num_dir_hmm] + '/' + dir_hmm[num_dir_hmm] + '.hmm '
				#	+ genome_assembly_path)
				command_args.append(command)
			else:
				if os.path.exists(script_dir +'/../'+db1+'/' + dir_hmm[num_dir_hmm]):
					db=db1
				if os.path.exists(script_dir +'/../'+db2+'/' + dir_hmm[num_dir_hmm]):
					db=db2
				
				command = f'nhmmer --cpu 4 -E {input_hmm_e_value} -o {work_dir}/HMM_out/{dir_hmm[num_dir_hmm]}.out {script_dir}/../{db}/{dir_hmm[num_dir_hmm]}/{dir_hmm[num_dir_hmm]}.hmm {genome_assembly_path}'
				command_args.append(command)
				#os.system(
				#	'nhmmer --cpu ' + str(cpus) + ' -E '+str(input_hmm_e_value)+' -o ' + work_dir + '/HMM_out/' + dir_hmm[num_dir_hmm] + '.out '
				#	+ script_dir + '/../'+db+'/' + dir_hmm[num_dir_hmm] + '/' + dir_hmm[num_dir_hmm] + '.hmm '
				#	+ genome_assembly_path)
	
	#This is better - nEfficient is the long term fix.
	with multiprocessing.Pool(cpus // 4) as pool:
		for result in pool.imap_unordered(hmm_worker, command_args):
			pass
	
	

def read_genome_assembly(genome_assembly_path):
	# ============== Read genome sequences according to sequence ID ===============
	def parse_seq_id(seq_header):
		for b in range(len(seq_header)):
			return seq_header.split()[b][1:]

	with open(genome_assembly_path) as genome_f:
		genome_sequences = {}
		lines = genome_f.readlines()
		# buffer list to store SINE segments to be added to dict
		sequence_buffer = []
		for line in lines:
			if line[0] == '>':
				if len(sequence_buffer) > 0:
					# has SINE segments in buffer to be added to dict
					genome_sequences[seq_id] = ''.join(sequence_buffer)  # reset buffer
					sequence_buffer = []  # reset buffer
				# parse new seq_id
				seq_id = parse_seq_id(line)
			else:
				sequence_buffer.append(line.strip().upper())
	if len(sequence_buffer) > 0:
		# add remaining SINE segments in buffer
		genome_sequences[seq_id] = ''.join(sequence_buffer)
		del sequence_buffer
	return genome_sequences

'''
(1) Read the file, check if hits exist
(2) If hits exist, read them and don't read alignment
(3) Filter hits by e value

Solving problems like this is exactly why nEfficient is useful
'''
def process_hmm_output_1(out_file, threshold_hmm_e_value, script_dir):
	# ============================ HMM prediction start and end annotation =======================
	hmm_predict_record_unsort = []
	hmm_predict_record_sort = []
	hmm_predict_family_number = 0
	ani=0
	#Check opening lines, remove '-' characters
	with open(out_file) as predict_f: #shujun
		lines = predict_f.readlines()
		for line in lines[15:]:
			if 'E-value' in line: ani=1
			if 'inclusion threshold' in line or 'No hits detected' in line or line == '\n':
				break
			else:
				check=re.sub('-','',line)
				check=check.strip()
				if not check=='':
					hmm_predict_record_unsort.append(line.split())
	
	#This terminates upon reaching the alignment section, but surely that can be gotten around by just changing the output format to tabular
	if ani==1:
		hmm_predict_record_unsort = []
		#Why reread the file? We already have the lines read into a list
		#with open(out_file) as predict_f:
		#	lines = predict_f.readlines()
		for line in lines[17:]:
			if 'inclusion threshold' in line or 'No hits detected' in line or line == '\n':
				break
			if not line:break
			else:
				check=re.sub('-','',line) 
				check=check.strip()
					
				if not check=='':
					hmm_predict_record_unsort.append(line.split())


	if [] not in hmm_predict_record_unsort:
		#print(hmm_predict_record_unsort)
		out_data = sorted(hmm_predict_record_unsort, key=lambda x: int(x[4]))
		#print(out_data)
		for i in range(len(out_data)):
			if float(out_data[i][0]) < threshold_hmm_e_value:
				if int(out_data[i][4]) < int(out_data[i][5]):
					hmm_predict_record_sort.append({'start': int(out_data[i][4]) - 1,
													'end': int(out_data[i][5]),
													'e_value': float(out_data[i][0]),

													'family': out_file.split('/', 1)[0],
													'id': out_data[i][3],
													'strand': '+'})
					if float(out_data[i][0]) <= 1:
						hmm_predict_family_number += 1
				else:
					hmm_predict_record_sort.append({'start': int(out_data[i][5]) - 1,
													'end': int(out_data[i][4]),
													'e_value': float(out_data[i][0]),
													'family': out_file.split('/', 1)[0],
													'id': out_data[i][3],
													'strand': 'C'})
					if float(out_data[i][0]) <= 1:
						hmm_predict_family_number += 1
	#print(hmm_predict_record_sort)
	#exit()
	return hmm_predict_record_sort, hmm_predict_family_number


'''
(1) Retrieve all outputs from hmmsearch
(2) Call process_hmm_output_1 on each output file to select only the hits - why not output in tabular format and save the alignment space if they aren't used?
(3) Count the number of famies with hits and how many hits per family
(4) Create a list of all the results in total? What is going on here?
'''
def process_hmm_output_2(threshold_hmm_e_value, script_dir):
	print('Processing the hmm prediction ...', flush=True)
	family_count = {}
	family_name = []
	update_hmm_record = []
#	dir_file = os.listdir(work_dir + '/HMM_out/') #shujun
	out_file = glob.glob(work_dir + '/HMM_out/' + '*.out') #shujun
	for a in range(len(out_file)):
		if out_file[a] != '.DS_Store':
			list_pre = process_hmm_output_1(out_file[a], threshold_hmm_e_value, script_dir)[0] #shujun
			for num_pre in range(len(list_pre)):
				if list_pre[num_pre]['e_value'] <= threshold_hmm_e_value:
					family_name.append(os.path.splitext(out_file[a])[0]) #shujun
					if os.path.splitext(out_file[a])[0] not in family_count:
						family_count[os.path.splitext(out_file[a])[0]] = 1 #shujun
					else:
						family_count[os.path.splitext(out_file[a])[0]] += 1 #shujun
			#Iterate through the list again? There's no filtering condition			
			#for num_return_broken in range(len(list_pre)):
			#	update_hmm_record.append(list_pre[num_return_broken])
			update_hmm_record.extend(list_pre)
				
	return update_hmm_record, family_name, family_count


def merge_same_hmm_output(hmm_output_record):
	def is_overlapping(record1, record2):
		if record1['strand'] != record2['strand']:
			return False
		#return (record1['start'] >= record2['start'] and record1['end'] <= record2['end']) or (record1['start'] <= record2['start'] and record1['end'] > record2['end']) or (record1['start'] < record2['start'] and record1['end'] == record2['end']) or ((record1['start'] < record2['start']<record1['end'] < record2['end'] and abs(record1['start']-record2['end'])>=0.5 * abs(record1['start'] - record1['end'])) or (record1['start'] < record2['start']<record1['end'] < record2['end'] and abs(record2['start']-record1['end'])>=0.5 * abs(record1['start'] - record1['end'])))
		return record1['start'] <= record2['end'] and record1['end'] >= record2['start']

	# Sort records by id and start position
	hmm_output_record = sorted(hmm_output_record, key=lambda x: (x['id'], x['start']))

	update_positions = {}
	last_merged_record = None

	for num_record in hmm_output_record:
		record_id = num_record['id']

		# If this is a new record_id or the current record doesn't overlap with the last merged record,
		# then simply append the current record
		if record_id not in update_positions or not is_overlapping(num_record, last_merged_record):
			update_positions.setdefault(record_id, []).append(num_record)
			last_merged_record = num_record
			continue

		# If the current record overlaps with the last merged record, then merge them
		last_merged_record['start'] =min(num_record['start'], last_merged_record['start'])
		last_merged_record['end'] = max(num_record['end'], last_merged_record['end'])
		last_merged_record['family'] = f"{num_record['family']}/{last_merged_record['family']}"
		last_merged_record['e_value'] = num_record['e_value']  # This assumes the current record's e_value should overwrite the last one
	res=[record for records in update_positions.values() for record in records]

	return [record for records in update_positions.values() for record in records]

def merge_same_hmm_output_raw(hmm_output_record):
	print('Merging the same hmm prediction ...', flush=True)
	update_positions = {}
	# check=0
	# check2=0
	# check3=0
	for num_record in hmm_output_record:
		record_id = num_record['id']
		# if record_id=='chr16':
		#	 check+=1

		if record_id in update_positions:
			add_pos = True
			for index in range(len(update_positions[record_id])):
				update_pos = update_positions[record_id][index]
				if (num_record['start'] >= update_pos['start'] and num_record['end'] <= update_pos['end']) or \
						(num_record['start'] <= update_pos['start'] and num_record['end'] > update_pos['end']) or \
						(num_record['start'] < update_pos['start'] and num_record['end'] == update_pos['end']) or \
						(((update_pos['start'] < num_record['start'] < update_pos['end'] < num_record['end'] and
						   abs(num_record['start'] - update_pos['end']) >= 0.5 * abs(
									num_record['start'] - num_record['end'])) or
						  (num_record['start'] < update_pos['start'] < num_record['end'] < update_pos['end']) and
						  abs(update_pos['start'] - num_record['end']) >= 0.5 * abs(
									num_record['start'] - num_record['end']))):
					if num_record['strand'] == update_pos['strand']:
						update_positions[record_id][index] = {'start': min(num_record['start'], update_pos['start']),
													   'end': max(num_record['end'], update_pos['end']),
													   'id': num_record['id'],
													   'strand': num_record['strand'],
													   'family': num_record['family'] + '/' + update_pos['family'],
													   'e_value': num_record['e_value']}
						add_pos = False
						# if record_id == 'chr16':
						#	 print('check2-numrecord-', check2, ':', num_record)
						#	 print('check2-posrecord-', check2, ':', update_pos)
						#	 check2+=1
						break
			if add_pos:
				update_positions[record_id].append(num_record)
				# if record_id == 'chr16':
				#	 check3+=1

		else:
			update_positions[record_id] = [num_record]
	# print('Check:',check,check2,check3)
	# exit()
	# print('Output merged array...')
	# t1=time.time()
	res = [e for r in update_positions for e in update_positions[r]]
	t2=time.time()
	# print('Output uses ',t2-t1)
	return res




'''
(1) Read entire genome into mem. Aggressive.
(2) Call process_hmm_output_3 to get (hmm_records, families_with_hits, counts_of_hits_by_family,)
(3) Merge overlapping domains in a wierd way
(4) Retrieve genome sequences
(5) Output modified seqid fasta
'''
def process_hmm_output_3(threshold_hmm_e_value, in_genome_assembly_path, pattern, out_genome_assembly_path):
	count = 0
	seq_ids = {}
	start = {}
	end = {}
	e_value = {}
	pre_strand = {}
	input_tsd_sequences = {}
	
	t1=time.time()
	output_genome_sequence = read_genome_assembly(in_genome_assembly_path)
	t2=time.time()
	print('process_hmm_output_3:read_genome_assembly uses',t2-t1,flush=True)
	
	t1=time.time()
	update_hmm_record, family_name, family_count = process_hmm_output_2(threshold_hmm_e_value, script_dir)
	t2=time.time()
	print('process_hmm_output_3:process_hmm_output uses',t2-t1,flush=True)

	t1=time.time()
	update_hmm_output = merge_same_hmm_output(update_hmm_record)
	t2=time.time()
	print('process_hmm_output_3:merge_same_hmm_output uses',t2-t1,flush=True)
	
	#Collect sequences from the genome file based on the output.
	#This makes much more sense to do by sequence ID loading either 1 seq/thread or just one at a time instead of keeping the whole genome in mem.
	t1=time.time()
	for num_return_pos in range(len(update_hmm_output)):
		pre_num = count
		start[pre_num] = update_hmm_output[num_return_pos]['start']
		end[pre_num] = update_hmm_output[num_return_pos]['end']
		seq_ids[pre_num] = update_hmm_output[num_return_pos]['id']
		num_id = update_hmm_output[num_return_pos]['id']
		e_value[pre_num] = update_hmm_output[num_return_pos]['e_value']
		input_tsd_sequences[pre_num] = output_genome_sequence[num_id][update_hmm_output[num_return_pos]['start']-30:
																	  update_hmm_output[num_return_pos]['end']+50]
		pre_strand[pre_num] = update_hmm_output[num_return_pos]['strand']
		count += 1
	t2=time.time()
	print('process_hmm_output_3:for_loop uses',t2-t1,flush=True)

	new_start = list(start.values())
	new_end = list(end.values())
	new_sequence_id = list(seq_ids.values())

	#This should only be called for pattern == 1 or 3
	#Looks like this is just writing the sequences to a file using a modified seqid and fasta format.
	if pattern == 1 or pattern == 3:
		#Is this just opening files? No - it's functionally deleting old files but isn't just deleting them. Why?
		if os.path.exists(out_genome_assembly_path+'/Step1_extend_tsd_input_1.fa'):
			modify_text(out_genome_assembly_path+'/Step1_extend_tsd_input_1.fa')
		
		save_to_fna_1(out_genome_assembly_path+'/Step1_extend_tsd_input_1.fa', input_tsd_sequences, pre_strand,
					  new_sequence_id, 0, 0, new_start, new_end)


#There has to be a more elegant way of doing this. Probably pyfastx is better.
def merge_tsd_input(pattern, out_genome_assembly_path):
	if os.path.exists(out_genome_assembly_path+'/Step1_extend_tsd_input.fa'):
		modify_text(out_genome_assembly_path+'/Step1_extend_tsd_input.fa')
	if pattern == 1:
		with open(out_genome_assembly_path+'/Step1_extend_tsd_input_1.fa', 'r') as f1:
			lines1 = f1.readlines()
			lines = lines1
	if pattern == 2:
		with open(out_genome_assembly_path+'/Step1_extend_tsd_input_2.fa', 'r') as f2:
			lines2 = f2.readlines()
			lines = lines2
	elif pattern == 3:
		with open(out_genome_assembly_path+'/Step1_extend_tsd_input_1.fa', 'r') as f1:
			with open(out_genome_assembly_path+'/Step1_extend_tsd_input_2.fa', 'r') as f2:
				lines1 = f1.readlines()
				lines2 = f2.readlines()
				if pattern == 1:
					lines = lines1
				elif pattern == 2:
					lines = lines2
				elif pattern == 3:
					lines = lines1 + lines2
	with open(out_genome_assembly_path+'/Step1_extend_tsd_input.fa', 'w') as f3:
		for line in lines:
			if line[0] == '>':
				head=line
				#f3.write(line)
			else:
				# replace non-ACGT characters to Ns, shujun
				line = line.rstrip('\n')  # Remove newline character
				cleaned_line = ''.join(['N' if c not in 'ACGTacgt' else c for c in line])
				cleaned_line=cleaned_line.strip()
				if not cleaned_line=='':
					f3.write(head+cleaned_line + '\n') # Add newline character back

#Looks like chunk genome to sequences for TSD searcher
def split_file(input_file, output_directory, max_lines_per_chunk=1000):
	# Create the output directory if it doesn't exist
	os.makedirs(output_directory, exist_ok=True)

	# Initialize variables
	chunk_count = 0
	line_count = 0

	# Open the input file in text mode for reading
	with open(input_file, 'r') as file:
		# Create the first output file
		chunk_file = open(os.path.join(output_directory, f'TSD_input_chunk_{chunk_count+1}.fa'), 'w')

		# Iterate over the lines in the input file
		for line in file:
			# Write the line to the current chunk file
			chunk_file.write(line)

			# Increment the line count
			line_count += 1

			# Check if the maximum line count per chunk is reached
			if line_count >= max_lines_per_chunk:
				# Close the current chunk file
				chunk_file.close()

				# Increase the chunk count
				chunk_count += 1

				# Create a new chunk file
				chunk_file = open(os.path.join(output_directory, f'TSD_input_chunk_{chunk_count+1}.fa'), 'w')

				# Reset the line count
				line_count = 0

		# Close the last chunk file
		chunk_file.close()
		
		#Check if the last chunk has no lines and remove it if it is empty
		if line_count == 0:
			os.remove(chunk_file)

	print(f'Successfully split the file into {chunk_count+1} chunks, with a maximum of {max_lines_per_chunk} lines per chunk.')

def check_file_has_rows(input_file):
	# Open the input file in text mode for reading
	with open(input_file, 'r') as file:
		# Read the first line from the input file
		first_line = file.readline()

		# Check if the first line is not empty
		if first_line.strip():
			return True  # File has at least one row
		else:
			return False  # File has no rows

def run_tsd_search(infile,outdir,script_dir,outfile):
	command = f'node {os.path.join(script_dir, 'TSD_Searcher_multi.js')} {outdir} {infile} {outfile}'
	os.system(command)
	retf = os.path.join(outdir, outfile)
	return retf

#It would be better to rewrite this in python, but barring that we can just make checkpointing code
def search_tsd(out_genome_assembly_path, script_dir, cpus=8):
	was_split_checker = os.path.join(out_genome_assembly_path, 'Step1_check_tsd_was_split.txt')
	#Do not re-split files
	if not os.path.exists(was_split_checker):
		split_file(os.path.join(out_genome_assembly_path, 'Step1_extend_tsd_input.fa'), out_genome_assembly_path, 1000)
		#Make a marker that this step was completed.
		with open(was_split_checker, 'w') as temp:
			pass
	
	ifiles=[f for f in os.listdir(out_genome_assembly_path) if 'TSD_input_chunk' in f]
	pres=[f.replace('input', 'output') for f in ofiles]
	
	
	
	ofiles = []
	
	tsd_outf = os.path.join(out_genome_assembly_path, 'Step2_tsd.txt')
	
	#with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as executor:
	with multiprocessing.Pool(cpus) as pool:
		for inf, pre in zip(ifiles, pres):
			#pre = filename.replace('input','output')
			#ofiles.append(executor.submit(run_tsd_search, filename, out_genome_assembly_path, script_dir, pre))
			ofiles.append(executor.submit(run_tsd_search, inf, out_genome_assembly_path, script_dir, pre))
		
		if os.path.exists(tsd_outf):
			needs_header = False
		else:
			needs_header = True
		with open(tsd_outf, 'ab') as out:
			for ef in concurrent.futures.as_completed(ofiles):
				tsd_chunk_output = ef.result()
				with open(tsd_chunk_output, 'rb') as inf:
					header = inf.readline()#Skip first line if it's happened once
					if needs_header:
						out.write(header)
						needs_header = False
					shutil.copyfileobj(inf, out)
					
				os.remove(tsd_chunk_output)
				os.remove(tsd_chunk_output.replace('output', 'input'))

#It would be better to rewrite this in python, but barring that we can just make checkpointing code
def search_tsd(out_genome_assembly_path, script_dir, cpus=8):
	was_split_checker = os.path.join(out_genome_assembly_path, 'Step1_check_tsd_was_split.txt')
	#Do not re-split files
	if not os.path.exists(was_split_checker):
		split_file(os.path.join(out_genome_assembly_path, 'Step1_extend_tsd_input.fa'), out_genome_assembly_path, 1000)
		#Make a marker that this step was completed.
		with open(was_split_checker, 'w') as temp:
			pass
	
	ifiles=[f for f in os.listdir(out_genome_assembly_path) if 'TSD_input_chunk' in f]
	pres=[f.replace('input', 'output') for f in ofiles]
	
	
	
	ofiles = []
	
	tsd_outf = os.path.join(out_genome_assembly_path, 'Step2_tsd.txt')
	
	with concurrent.futures.ProcessPoolExecutor(max_workers=cpus) as executor:
		'''
		for filename in os.listdir(out_genome_assembly_path):
			if not re.search('TSD_input_chunk', filename):
				continue
			if not check_file_has_rows(os.path.join(out_genome_assembly_path, filename)):
				continue
		'''
		for inf, pre in zip(ifiles, pres):
			#pre = filename.replace('input','output')
			#ofiles.append(executor.submit(run_tsd_search, filename, out_genome_assembly_path, script_dir, pre))
			ofiles.append(executor.submit(run_tsd_search, inf, out_genome_assembly_path, script_dir, pre))
		
		if os.path.exists(tsd_outf):
			needs_header = False
		else:
			needs_header = True
		with open(tsd_outf, 'ab') as out:
			for ef in concurrent.futures.as_completed(ofiles):
				tsd_chunk_output = ef.result()
				with open(tsd_chunk_output, 'rb') as inf:
					header = inf.readline()#Skip first line if it's happened once
					if needs_header:
						out.write(header)
						needs_header = False
					shutil.copyfileobj(inf, out)
					
				os.remove(tsd_chunk_output)
				os.remove(tsd_chunk_output.replace('output'))
				
					
				#tofiles.append(ef.result())
	
	#tofiles_sorted = sorted(tofiles, key=lambda x: int(re.search(re.compile(r"TSD_output_chunk_(\d+).fa"), x).group(1)))

	
	#Non os-system remove
	c1 = os.path.join(out_genome_assembly_path, 'TSD_output_chunk_*')
	c2 = os.path.join(out_genome_assembly_path, 'TSD_input_chunk_*')
	for file in glob.glob(c1):
		os.remove(file)
	for file in glob.glob(c2):
		os.remove(file)

def is_at_seq(seq, tolerance):
	"""
	Test if a sequence consists with 'A' ('a') and 'T' ('t').
	Occurrence of other bases cannot be more than `tolerance`.

	:param seq: sequence to test
	:param tolerance: Maximum number of occurrence of other bases.

	:return: Whether the tested sequence matches the condition.
	"""
	base_dict = Counter(seq.lower())
	if 'a' in base_dict:
		base_dict.pop('a')
	if 't' in base_dict:
		base_dict.pop('t')
	return sum(base_dict.values()) <= tolerance

def process_tsd_output(in_genome_assembly_path, out_genome_assembly_path):
	#Allow previous checkpointing code to take over again
	was_split_checker = os.path.join(out_genome_assembly_path, 'Step1_check_tsd_was_split.txt')
	if os.path.exists(was_split_checker):
		os.remove(was_split_checker)
		
	tsd_output_file = os.path.join(out_genome_assembly_path, 'Step1_extend_tsd_input.fa')
	title = []
	sequences = []
	hmm_start = []
	hmm_end = []
	hmm_id = []
	with open(tsd_output_file) as tsd_file:
		lines = tsd_file.readlines()
		hmm_tsd = []
		finder_tsd = []
		filename = os.path.join(out_genome_assembly_path, 'Step2_tsd.txt')
		for line in lines:
			if line[0] == '>':
				title.append(line)
				hmm_start.append(int(line.split()[2].split(':')[0]))
				hmm_end.append(int(line.split()[2].split(':')[1]))
				hmm_id.append(line.split()[0].replace('>', ''))
			else:
				sequences.append(line)
				if not line[0]=='':
				  hmm_tsd.append(len(line.split()[0]))
				else:
				  hmm_tsd.append(0)
			  

	hmm_pos = []
	record_tsd = []
	with open(filename) as f2:
		lines = f2.readlines()
		for line in lines:
			if line[0] != '>':
				hmm_pos.append(line.strip().split())
				if len(line) != 1:
					left_tsd = line.split()[0]
					right_tsd = line.split()[2]
					if is_at_seq(left_tsd.replace('-', ''), 0) or is_at_seq(right_tsd.replace('-', ''), 0):
						record_tsd.append(0)
					else:
						record_tsd.append(1)
				else:
					record_tsd.append(0)

	starts = []
	ends = []
	tsd_info = []
	input_seq = []

	output_genome_sequence = read_genome_assembly(in_genome_assembly_path)
	for t in range(len(hmm_pos)):
		if record_tsd[t] == 0:
			tsd_info.append('tsd not exist')
			start = hmm_start[t] - 100
			end = hmm_end[t] + 100
			seq_id = hmm_id[t]

			input_seq.append(output_genome_sequence[seq_id][start:end])
			starts.append(hmm_start[t] + 30)
			ends.append(hmm_end[t] - 50)
		else:
			tsd_info.append(len(hmm_pos[t][0]))
			start = hmm_start[t] - 30 + int(hmm_pos[t][1].split('-')[1].replace(')', '')) - 100
			end = hmm_start[t] - 30 + int(hmm_pos[t][3].split('-')[0].replace('(', '')) + 100
			seq_id = hmm_id[t]

			input_seq.append(output_genome_sequence[seq_id][start:end])
			starts.append(
				hmm_start[t] - 30 + int(hmm_pos[t][1].split('-')[1].replace(')', '')))  # tsd boundary
			ends.append(hmm_start[t] - 30 + int(hmm_pos[t][3].split('-')[0].replace('(', '')))
	
	tsd_outfa = os.path.join(out_genome_assembly_path, 'Step2_tsd_output.fa')
	if os.path.exists(tsd_outfa):
		modify_text(tsd_outfa)
	save_to_fna_2(tsd_outfa, input_seq, title, tsd_info, starts, ends)

	s2_extend_blast_fa = os.path.join(out_genome_assembly_path, 'Step2_extend_blast_input.fa')
	if os.path.exists(s2_extend_blast_fa):
		modify_text(s2_extend_blast_fa)
	count=0
	
	with open(tsd_outfa, 'r')as tsd_output_file:
		with open(s2_extend_blast_fa, 'w') as blast_input_file:
			flag = False
			for line in tsd_output_file:
				#print(line)
				if line[0] == '>':
					tem=''
					if 'tsd not exist' not in line.split('|')[1]:
						flag = True
						tem=line
					else:
						flag = False
				if flag and not line[0]=='>':
					line=line.strip()
					if line:
						blast_input_file.write(tem)
						blast_input_file.write(line+'\n')
						count+=1
	blast_input_file.close()
	
	rename_fa = os.path.join(out_genome_assembly_path,'Step2_extend_blast_input_rename.fa')
	
	f=open(s2_extend_blast_fa, 'r')
	o=open(rename_fa,'w+')
	c=1
	while True:
		line=f.readline().strip()
		if not line:break
		if re.search('>',line):
			ele=line.split()
			name=ele[0]+'_'+str(c)+' '+' '.join(ele[1:])
			o.write(name+'\n')
			c+=1
		else:
			o.write(line+'\n')
	o.close()
	
	#Create index for later use.
	pyfastx.Fasta(rename_fa, build_index = True)
		
	if count==0:
		print('No sequences in Step2_extend_blast_input.fa! Please check! Exit.')
		exit()

def modify_text(modify_name):
	with open(modify_name, "r+") as f:
		f.truncate()


def save_to_fna_1(filename, sequences, strands, ids, left_offset, right_offset, starts, ends):
	header = '>{} {} {}:{}\n'
	index = 0
	payload = []
	for seq in sequences:
		seq_header = header.format(ids[seq], strands[seq], starts[seq] - left_offset, ends[seq] + right_offset)
		payload.append(seq_header)
		payload.append(sequences[seq] + '\n')
		index += 1
	with open(filename, 'a') as file:
		file.writelines(payload)
	if len(payload)==0:
		print('Warning! No SINE can be detected by HMM!')
		#exit()


def save_to_fna_2(filename, sequences, input_title, input_tsd, input_start, input_end):
	HEADER = '{}'.strip() + '|tsd_l:{}|tsd_s:{}|tsd_e:{}'.strip() + '\n'
	payload = []
	for seq in range(len(sequences)):
		seq_header = HEADER.format(input_title[seq].strip(), input_tsd[seq], input_start[seq], input_end[seq])
		payload.append(seq_header.strip() + '\n')
		payload.append(sequences[seq].strip() + '\n')
	with open(filename, 'a') as file:
		file.writelines(payload)


def is_file_size_exceeded(file_path, max_size):
	file_size = os.path.getsize(file_path)
	file_size_mb = file_size / (1024 * 1024)
	if file_size_mb > max_size:
		return True
	else:
		return False

#Improved readability
'''
'''
def mmindex(infile, k_size):
	mmidx = infile.replace('.fasta', '.mmseqs2db')
		
	command = f'minimap2 -t 1 -k {k_size} -I 2G -d {mmidx} {infile}'
	command = command.split()
	subprocess.run(command, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
	
	return mmidx

def multiple_sequence_alignment(e_value, in_genome_assembly_path, out_genome_assembly_path, cpus, input_num_alignments, input_temd):
	print('Minimap2 againist the genome assembly ...', flush=True)
	#workspace = os.path.join(out_genome_assembly_path, "minimap2")
	#if not os.path.exists(workspace):
	#	os.mkdir(workspace)
	
	out_genome_assembly = os.path.join(out_genome_assembly_path, "Step2_extend_blast_input_rename.fa")
	is_big_file = is_file_size_exceeded(out_genome_assembly, 10)
	if is_big_file:
		ksize = 10
	else:
		ksize = 8
	
	'''
	gs = genomeSplitter(genome_file = in_genome_assembly_path, 
				output_directory = workspace, 
				chunk_size = 2_000_000_000, #Basically allow any size of genome that will fit into mmseqs2 to minimize long sequence breakup. Mostly this just chunks whole genomes.
				overlap_size = 0, #This is a fundamentally inelegant solution for which solutions are also inelegant
				procs = cpus,
				smart = False,
				post_index = False,
				verbose = False,			
				overwrite = False,
				quiet = True)
	
	output_subsets = gs.run()
	
	
	
	ok_num = min([len(output_subsets), cpus])
	args = [(genome_subset, ksize,) for genome_subset in output_subsets]
	
	print("Creating mmseqs index for search...")
	with multiprocessing.Pool(ok_num) as pool:
		for s in pool.starmap(mmindex, args):
			print(f'Created minimap index {s}')
			
	mmseqs_indices = [os.path.join(workspace, f) for f in os.listdir(workspace) if 'mmseqs2db' in f]

	outputs = [idx.replace('mmseqs2db', 'mmseqs.paf') for idx in mmseqs_indices]
	'''
	
	
	mmidx = os.path.join(out_genome_assembly_path, 'Step2_mimimap_index.idx')
	
	comm = f'minimap2  -t {cpus} -k {ksize} -I 1G -d {mmidx} {in_genome_assembly_path}'
	comm = comm.split()
	print('Creating minimap2 index...')
	subprocess.run(comm, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
	
	output_paf = os.path.join(out_genome_assembly_path, "Step3_blast_output.paf")

	comm = f'minimap2 -c -t {cpus} -p 0.01 -K 500M -N {input_num_alignments} {mmidx} {out_genome_assembly} -o {output_paf}'

	comm = comm.split()
	comm = subprocess.run(comm)	
		
	chunk_size=100_000 #Read inputs 500k lines at a time using pandas
	max_lines = 5_000_000 #~500 MB output file size

	#Directly call the python behavior instead of calling by system
	parallel_process_file(output_paf, out_genome_assembly_path, chunk_size, max_lines, cpus, input_temd)
	


	
#Iterate through a blast input one
class blast_input_parser:
	def __init__(self, blast_file):
		self.bf = blast_file
		
	def record_iterator(self):
		with open(self.bf) as inf:
			starts, ends = [], []
			segs = inf.readline().strip().split('\t')
			previous_record_id = segs[0]
			previous_parent_contig = segs[1]
			previous_length = int(segs[3])
			start, end = int(segs[6]), int(segs[7])
			starts.append(start)
			ends.append(end)
			for line in inf:
				segs = line.strip().split('\t')
				group_id = segs[0]
				parent_contig = segs[1]
				start, end = int(segs[6]), int(segs[7])
				if previous_record_id != group_id:
					yield previous_record_id, previous_parent_contig, previous_length, starts, ends
					previous_record_id = group_id
					previous_parent_contig = parent_contig
					previous_length = int(segs[3])
					starts = []
					ends = []
				starts.append(start)
				ends.append(end)
				
		yield previous_record_id, previous_parent_contig, previous_length, starts, ends

'''
What does this do?

(1) Load the entire genome assembly. Ouch.
(2) Concurrent futures run process_file on each output

Process file:
	(1) Read every single line into memory and split lines on tabs; set segment[0] to ''
	(2) 
'''
def process_blast_output_1(in_genome_assembly_path, factor_length, factor_copy, max_shift, max_gap, min_copy_number,
						   pos, out_genome_assembly_path, bound, figure, cpus):
	print('Processing the BLAST output ...', flush=True)
	# ================================ Read BLAST Output ==========================
	step3_files = [f for f in os.listdir(out_genome_assembly_path) if re.search('Step3_blast_output.out', f)]
	pcount = len(step3_files)
	
	#this is definitely something to do with pyfastx instead
	#output_genome_sequence = read_genome_assembly(in_genome_assembly_path)
	output_genome_sequence = pyfastx.Fasta(in_genome_assembly_path)
	
	#This whole function is likely to produce incorrect behavior because of chunking issues. 
	#Blast outputs breaking in the middle of a run of items will result in either double counts or misses for some sequences
	def process_file(filename):
		if not re.search('Step3_blast_output.out', filename):
			return
			
		fa_path = os.path.join(out_genome_assembly_path, f'{filename}.fa')
		
		if os.path.exists(fa_path) and os.path.getsize(fa_path) > 0:
			print(f'Skipping {filename} as {fa_path} already exists', flush=True)
			return
		
		#Do we actually need to do this?
		blast_chunk = os.path.join(out_genome_assembly_path, filename)
		
		#We shouldn't need to do this.
		'''
		dp=set()
		with open(blast_chunk, 'r') as f:
			#Read whole file in a very strange way - there definitely is NOT a need to do this.
			for line in f:
				line=f.readline().strip()
				if not line:break
				ele=line.split('\t')
				dp.add(ele[0])
		'''
		
		previous_output = os.path.join(out_genome_assembly_path, 'Step2_extend_blast_input_rename.fa')
		fa = pyfastx.Fasta(previous_output)
		#dused = set(fa.keys())
		#dused = dused.intersection(dp)
		
		#Get the sequence descriptions
		sine_info = []
		seqids = []
		previous_start = []
		previous_end = []
		previous_id = []
		tsd_length = []
		tsd_bound_start = []
		tsd_bound_end = []
		sinf_dict = {}
		for seq in fa:
			sine_info.append(seq.description)
			seqid, classifier, tsd_info = seq.description.split()
			tsd_groups = tsd_info.split('|')
			tsd_groups = [t.split(':') for t in tsd_groups]
			start, end = int(tsd_groups[0][0]), int(tsd_groups[0][1])
			tsd_l = int(tsd_groups[1][1])
			tsd_start = int(tsd_groups[2][1])
			tsd_end = int(tsd_groups[3][1])
			
			sinf_dict[seqid] = (start, end, tsd_l, tsd_start, tsd_end,)
			
			seqids.append(seqid)
			previous_start.append(start)
			previous_end.append(end)
			tsd_length.append(tsd_l)
			tsd_bound_start.append(tsd_start)
			tsd_bound_end.append(tsd_end)
			
			
		
		'''
		#Create a unique ID for a temp file
		blast_out_filename = uuid.uuid1().hex
		
		dused=set()
		f2=open(previous_output,'r')
		#Open an output file
		ou=open(blast_out_filename,'w+')
		
		for line in f2:
			if not line:break
			if len(dp)==dused:break
			line = line.strip()
			#Check if '>' is in a line and replace if if so
			if line.startswith('>'):
				rep = line.split()
				pr = rep[0]
				pr=pr[1:]
				
				#Split line on '_' and remove the final split from this
				if pr in dp:
					po = line.split('_')
					po = '_'.join(po[:-1])
					
					#line=re.sub(pr,po,line)
					rep[0] = f'>{po}'
					line = ' '.join(rep)
					
					ou.write(line+'\n')
					dused.add(pr)
					
					go=True
				else:
					go=False
			else:
				if go:
					ou.write(line+'\n')
		ou.close()
		'''
		#Up to here - this is filter the extend_blast_input file to only those sequences which appear some time in the blast output.
		#Completely unneeded at this point, pyfastx solves this issue much better anyway
				
		'''
		#Previous seq never actually gets used in this script and so is just wasted time and mem.
		with open(blast_out_filename) as blast_input_file:
			for line in blast_input_file:
				if line[0] != '>':
					previous_seq.append(line)
				else:
					previous_start.append(int(line.split()[2].split('|')[0].split(':')[0]))
					previous_end.append(int(line.split()[2].split('|')[0].split(':')[1]))
					previous_id.append(line.split()[0].replace('>', ''))
					tsd_length.append(int(line.split('|')[1].split(':')[1]))
					tsd_bound_start.append(int(line.split('|')[2].split(':')[1]))
					tsd_bound_end.append(int(line.split('|')[3].split(':')[1]))
					sine_info.append(line)
		
		os.system(f'rm {blast_out_filename}')
		'''
		

		#-------------- Pre-load finished
		#One blast output file.
		fpre=filename
		blast_inter = []
		length = []
		
		#Why save this 
		filename_full_path = os.path.join(out_genome_assembly_path, filename)
				
		tem=[]
		dused=set()
		
		#Okay, so this is all top hits. And therefore, this filtering ought to be done by the chunking/conversion code.
		
		'''
		(1) Read blast input file
		(2) Group by subsequence id
			For each group
			(a) First element determines length
			(b) Additional records are valid IFF their length is > factor_length * (length(first_record) - 200), by default 30%
			(c) Count the number of passing records as family_count
			(d) Get depth of coverage by genomic position for the hit, plot if requested
			(e) Check if there are enough copies of the TE to continue
			(f) If there are enough copies, subset to high coverage positions from the hit
			(g) For those high coverage positions, check if they are contiguous and within acceptable boundaries of the start and end of the alignment window
			(h) For each sequence, if it needs updated, update sequence, collect, and output to fasta
			
			
			
		'''
		#This will be the architecture for the whole per-group logic.
		bp = blast_input_parser(filename_full_path)
		for seqid, parent, group_length, starts, ends in bp.record_iterator():
			threshold = factor_length * (group_length - 200)
			#print(seqid, parent, starts, ends)
		
		#Surely this is better done by database or pandas or polars
		
		#Series of dicts encoding the blast output in memory.
		#Oh, I see - tem is doing a group_by, so it's groups of blast records by, what, gene? Something like that.
		with open(filename_full_path) as blast_output_file_1:
			for line in blast_output_file_1:
				ele=line.split()
				if ele[0] not in dused:
					dused.add(ele[0])
					if len(tem)==0:
						tem.append({'start':int(ele[6])-1,'end':int(ele[7]),'id':ele[1]})
						length.append(int(ele[3]))
					else:
						blast_inter.append(tem)
						length.append(int(ele[3]))
						tem=[]
						tem.append({'start':int(ele[6])-1,'end':int(ele[7]),'id':ele[1]})
				else:
					#I see, this is the add to existing group idea. Length is only the first item, but all hits are recorded
					tem.append({'start':int(ele[6])-1,'end':int(ele[7]),'id':ele[1]})
		
		#Final iteration
		blast_inter.append(tem)
		
		blast_res_0 = []
		blast_res = []
		
		#For each 'gene?' group, do this stuff
		#This can be applied to each tem group as we're loading them instead of after they're all loaded.
		for num_blast in range(len(blast_inter)):
			#Check if other records are an acceptable fraction of the length of the first record
			threshold = factor_length * (length[num_blast] - 200)
			blast_res_0.append([])
			
			#Iterate through records in the list
			for seq in range(len(blast_inter[num_blast])):
				#Always add the first record / besthit
				if seq == 0:
					blast_res_0[num_blast].append(blast_inter[num_blast][seq])
					
				start = blast_inter[num_blast][seq]['start']
				end = blast_inter[num_blast][seq]['end']
				seq_id = blast_inter[num_blast][seq]['id']
				#if subsequent records are long enough compared to the first, add them. 
				#Ah, I see - the first hit is recorded twice
				if abs(end - start + 1) >= threshold:
					blast_res_0[num_blast].append(blast_inter[num_blast][seq])
		
		#Count the number of times each hit appears after length filtering and excluding the first hit.
		#Actually, this includes the first hit because it's added twice. It will, by definition, pass the threshold
		family_count = []
		for num_blast in range(len(blast_res_0)):
			family_count.append(len(blast_res_0[num_blast])-1)
			blast_res.append(blast_res_0[num_blast][1:])

		# ========================= Plot BLAST Sequences Alignment ==============
		def plot_line(y, x1, x2):
			plt.hlines(y, x1, x2, linewidth=0.2)

		if figure == 'y':
			figpath = os.path.join(out_genome_assembly_path, 'Figures')
			if not os.path.exists(figpath):
				os.makedirs(figpath)
			for num in range(len(blast_res)):
				plt.yticks(fontsize=8)
				plt.xticks(fontsize=8)
				y_axis = len(blast_res[num])
				plt.ylabel('No. alignment', fontsize=8)
				plt.xlabel('No. base of each SINE', fontsize=8)
				for j in range(len(blast_res[num])):
					y_axis -= 1
					plot_line(y_axis, blast_res[num][j]['start'], blast_res[num][j]['end'])
				# plt.title('MSA_'+str(num))
				# plt.legend()
				# plt.show()
				plt.savefig(os.path.join(out_genome_assembly_path, 'Figures', 'MSA_{num}.png'))
				plt.close()

		# ================== Statistical Hit Number of Each Position ================
		def get_keys(d):
			sort_key_dic_instance = dict(sorted(d.items(), key=operator.itemgetter(0)))
			return sort_key_dic_instance

		#
		#For each collection of results not including the top hit (why?), get the sum of times each genome position was covered.
		#Would be better as numpy.zeros(shape = (max_end-min_start + 1)), sum indices over start - min_start : end-min_start
		#instead of this dict approach
		hit_list = []
		x_value = []
		y_value = []
		for i in range(len(blast_res)):
			hit_list.append({})
			for j in range(len(blast_res[i])):
				for pos in range(blast_res[i][j]['start'], blast_res[i][j]['end']):
					#print('Check_pos: ',pos)
					if pos in hit_list[i]:
						hit_list[i][pos] += 1
					else:
						hit_list[i][pos] = 1
			
			#These are used only if the
			if figure == 'y':
				sort_hit_list = get_keys(hit_list[i])
				x_value.append(list(get_keys(sort_hit_list).keys()))
				y_value.append(list(get_keys(sort_hit_list).values()))

		if figure == 'y':
			for num_value in range(len(x_value)):
				# ================================== Bar Plot==================================
				plt.plot(x_value[num_value], y_value[num_value], linewidth=1.5, marker='o', markersize=2)
				plt.yticks(fontsize=8)
				plt.xticks(fontsize=8)
				# plt.ylim(0, 150)
				plt.xlabel('No. base of each SINE', fontsize=8)
				plt.ylabel('Copy number', fontsize=8)
				# plt.title(f'Copy number profile_{num_value}')
				plt.tight_layout()
				
				plt.savefig(os.path.join(out_genome_assembly_path+'Figures', 'profile_'+str(num_value)+'.png'))
				# plt.grid(linestyle='dashed')
				plt.close()

		# ================= Choose Repeat Number and Decide Boundary ===============

		#Redefining a function with the same name but different behavior is odd
		def get_keys(d, value):
			return sorted(k for k, v in d.items() if v >= value)

		def split_by_continuity(positions):
			result = []
			pos_start = None
			pos_end = None
			for index in range(len(positions) - 1):
				inter_pos = positions[index]
				next_pos = positions[index + 1]
				if pos_start is None:
					pos_start = inter_pos
				pos_end = inter_pos
				if next_pos != inter_pos + 1:
					# if is not continuous
					result.append({
						'start': pos_start,
						'end': pos_end
					})
					pos_start = None
					pos_end = None
			if pos_start is not None:
				result.append({
					'start': pos_start,
					'end': positions[-1]
				})
			return result

		
		update_pos = []
		blast_count = []
		prob_copy = []
		num = 0
		max_gap_list = []
		print('Copy Number array for your input: ',family_count,'| Required min copy number:',min_copy_number)
		
		#For each set of coverage counts
		for m in range(len(hit_list)):
			#Find the repeat number and record it
			repeat_number = int(math.ceil(family_count[m] * factor_copy))
			prob_copy.append(round(repeat_number, 3))
			blast_count.append(family_count[m])
			#If there are enough copies of the hit
			if family_count[m] >= min_copy_number:
				#I guess make sure that min copy number wasn't 0 or negative? Should probably just happen somewhere in opts
				#Two if checks here for no real reason
				if family_count[m] >= 1:
					#res gets the positions in the genome which are covered at least repeat_number times, equivalent to np.where(counts_arr > repeat_number)
					res = get_keys(hit_list[m], repeat_number)
					#Find starts, ends for indices where multiple contiguous positions are covered at acceptable depth
					blast_pre = split_by_continuity(res)
					# first and last base copy number
					tid=length[m] - 1
					#print(tid)
					#This would be somehow the first or last position isn't covered at least once. That should never happen
					if tid not in hit_list[m] or 0 not in hit_list[m]:continue
					
					#This is *supposed* to be a check that the first and last positions have at least 1 coverage, which they'd have to if the origianl record is kept
					if hit_list[m][0] <= pos and hit_list[m][length[m] - 1] <= pos:
						print('Shift value for your input:', abs(blast_pre[0]['start'] - 100), abs(blast_pre[0]['end'] - (length[m] - 100)),' | Maximum tolerable shift value:',max_shift)
						if len(blast_pre) == 1:
							# shift
							if max(abs(blast_pre[0]['start'] - 100),
								   abs(blast_pre[0]['end'] - (length[m] - 100))) >= max_shift:
								update_pos.append({})
							else:
								update_pos.append(blast_pre)
						#If the coverage isn't fully contiguous...
						elif len(blast_pre) > 1:
							#If there is exactly one gap?
							if len(blast_pre) - 2 == 0:
								# gap
								#If the gap between the first contiguous widow is larger than the max acceptable gap AND
								#The the end of the first contiguous window is within 100 of the end of the whole thing AND
								#The start of the second contiguous window is > 100
								
								#These rules are wild and don't make any sense
								if abs(blast_pre[0]['end'] - blast_pre[1]['start']) >= max_gap and \
									blast_pre[0]['end'] < length[m]-100 and \
										blast_pre[1]['start'] > 100:
									#Then add the item
									max_gap_list.append(m)
									del blast_pre[0]
								#If the first contiguous window's start is within 100+max_shift of 0 and the final contiguous window's is within end > 100+max_shift of 0,
								if max(abs(blast_pre[0]['start'] - 100),
									   abs(blast_pre[-1]['end'] - (length[m] - 100))) >= max_shift:
									update_pos.append({})
								#Keep the position. So basically check if the start start and end end are close to full sized
								else:
									update_pos.append(blast_pre)
							#If there is more than 1 gap:
							#The logic for 1+ gap should work for exactly one gap, too. IDK why these are separate ifs
							else:
								#This doesn't feel like it can be right 
								#It looks like it stole the behavior from above and didn't didn't handle it appropriately
								#Every comparsion between gappy segments is done with the first and second segments, never between, say, segments 2 & 3
								new_blast_pre = []
								#For each gap
								for num_cut in range(0, len(blast_pre) - 2):
									#Check to make sure each gap is smaller than max_gap
									#Check if the end of the first contiguous segnment is within 100 of the end
									#Check if the start of the second segment is < 100
									if abs(blast_pre[num_cut]['end'] - blast_pre[num_cut + 1]['start']) < max_gap and \
										blast_pre[0]['end'] < length[m]-100 and blast_pre[1]['start'] > 100:
										new_blast_pre.append(blast_pre[num_cut])
								
								#This is just checking if the first start and last end of the updated value are within shift, basically like the fiest check
								if len(new_blast_pre) == 0 or \
										max(abs(new_blast_pre[0]['start'] - 100),
											abs(new_blast_pre[-1]['end'] - (length[m] - 100))) >= max_shift:
									update_pos.append({})
								else:
									update_pos.append(new_blast_pre)
						else:
							update_pos.append({})
					else:
						update_pos.append({})
				else:
					update_pos.append({})
			else:
				update_pos.append({})
			num += 1
		p=0
		for e in update_pos:
			if not  e=={}:
				p=1
		if p==0 and pcount==1:
			print('Nothing left after Step3... Program exit.')
			exit()
		
		
		# ================================ Correct Boundary Results ========================
		finder_seq = []
		finder_start = []
		finder_end = []
		#The finder seq logic makes sense to do with pyfastx, but it should be streamed to an output instead of grouped for writes
		for t in range(len(update_pos)):
			#If the position doesn't need an update
			if len(update_pos[t]) == 0:
				# or abs(update_pos[t][0]['start']-update_pos[t][-1]['end']) <= 90:
				seq_id = previous_id[t]
				finder_seq.append(output_genome_sequence[seq_id][previous_start[t]:previous_end[t]])
				finder_start.append(0)
				finder_end.append(0)
			else:
				#Bound is a step parameter
				old_pre_start = tsd_bound_start[t]
				old_pre_end = tsd_bound_end[t]
				
				new_pre_start = old_pre_start + update_pos[t][0]['start'] - 100
				new_pre_end = old_pre_end + update_pos[t][-1]['end'] + 100
				
				finder_start.append(new_pre_start)
				finder_end.append(new_pre_end)
				seq_id = previous_id[t]
				finder_seq.append(output_genome_sequence[seq_id][new_pre_start:new_pre_end])

		fprepath = os.path.join(out_genome_assembly_path, f'{fpre}.fa')
		print("fprepath")
		print(fprepath)
		if os.path.exists(fprepath):
			modify_text(fprepath)
		
		save_to_fna_3(fprepath, finder_seq, sine_info, finder_start, finder_end, blast_count, length)
					  
	with concurrent.futures.ThreadPoolExecutor(max_workers=min(cpus, len(step3_files))) as executor:
		executor.map(process_file, step3_files)
	
	output_wildcard = os.path.join(out_genome_assembly_path, 'Step3_blast_output.out*.fa')
	combined_output = os.path.join(out_genome_assembly_path, 'Step3_blast_process_output.fa')
	chunk_fastas = glob.glob(output_wildcard)
	
	print('chunk fastas')
	print(chunk_fastas)
	
	with open(combined_output, 'wb') as out:
		for f in chunk_fastas:
			with open(f, 'rb') as inf:
				shutil.copyfileobj(inf, out)

	

def save_to_fna_3(filename, sequences, title, bs, be, num, l):
	header = '{}|blast_s:{}|blast_e:{}|blast_count:{}|blast_l:{}'.strip()
	payload = []
	for seq in range(len(sequences)):
		seq_header = header.format(title[seq].strip(), str(bs[seq]), str(be[seq]), num[seq], l[seq])
		payload.append(seq_header.strip() + '\n')
		payload.append(sequences[seq].strip() + '\n')
	with open(filename, 'a') as file:
		file.writelines(payload)


def process_blast_output_2(out_genome_assembly_path):
	input_f1 = out_genome_assembly_path+'/Step3_blast_process_output.fa'
	input_f2 = out_genome_assembly_path+'/Step4_rna_input.fasta'
	with open(input_f1, 'r')as f1:
		with open(input_f2, 'w') as f2:
			flag = False
			for line in f1:
				if line[0] == '>':
					if 'blast_s:0' not in line:
						flag = True
					else:
						flag = False
				if flag:
					f2.write(line)


def blast_rna(out_genome_assembly_path, cpus, script_dir):
	os.system('blastn -query '+out_genome_assembly_path+'/Step4_rna_input.fasta '
			  '-subject ' + script_dir + '/../Input_Files/rna_database.fa ' #shujun
			  '-out '+out_genome_assembly_path+'/Step4_rna_output.out '
			  '-evalue 1 '
			  '-num_alignments 50000 '
			  '-word_size 7 '
			  '-gapopen 5 '
			  '-gapextend 2 '
			  '-penalty -3 '
			  '-reward 2 '
			  '-num_threads '+ str(cpus))


def process_rna(out_genome_assembly_path):
	rna_database = []
	with open(script_dir + '/../Input_Files/rna_database.fa') as f: #shujun
		num = 0
		lines = f.readlines()
		for line in lines:
			if line[0] == '>':
				rna_database.append([])
				if num in range(954):
					rna_database[0].append(line.replace('>', '').strip())
				elif num in range(954, 1666):
					rna_database[1].append(line.replace('>', '').strip())
				elif num in range(1666, 1729):
					rna_database[2].append(line.replace('>', '').strip())
				num += 1

	filter_out = []
	min_e_values = []
	E_VALUE_PATTERN = re.compile(r'Expect = (?P<e_value>\S+)')
	input_f = out_genome_assembly_path+'/Step4_rna_output.out'
	with open(input_f)as f:
		rna_id = []
		title = None
		num = -1
		for line in f:
			if 'Query= ' in line:
				rna_id.append([])
				min_e_values.append(None)
				num += 1
				title = line.split('|')[0]
				filter_out.append(title.strip('\n'))
			if '> ' in line:
				rna_id[num].append(line.split()[1])
			if num > -1 and min_e_values[num] is None:
				# if current e_value is not found, try parsing e_value
				# num > -1 ensures only starting to parse e_value after the query seq is identified
				match = E_VALUE_PATTERN.search(line)
				if match:
					min_e_values[num] = float(match.group('e_value'))

	# remove titles whose e_value is smaller than 1e-20
	filter_out = [title for index, title in enumerate(filter_out)
				  if min_e_values[index] is None or min_e_values[index] >= 1e-15]

	hit_record = []
	for num in range(len(rna_id)):
		if rna_id[num]:
			if rna_id[num][0] in rna_database[0]:
				hit_record.append('tRNA')
			elif rna_id[num][0] in rna_database[1]:
				hit_record.append('5S_rRNA')
			elif rna_id[num][0] in rna_database[2]:
				hit_record.append('7SL_RNA')
			else:
				hit_record.append('Unknown')
		else:
			hit_record.append('Unknown')
	if os.path.exists(out_genome_assembly_path+'/Step4_rna_output.fasta'):
		modify_text(out_genome_assembly_path+'/Step4_rna_output.fasta')
	input_f1 = out_genome_assembly_path+'/Step4_rna_input.fasta'
	input_f2 = out_genome_assembly_path+'/Step4_rna_output.fasta'
	with open(input_f1, 'r')as rna_input_file:
		with open(input_f2, 'w')as rna_output_file:
			flag = False
			num = -1
			for line in rna_input_file:
				if line[0] == '>':
					"""if model == 'rule':"""
					if ('Query= ' + line.split('|')[0].replace('>', '')).split(';')[0].rsplit(' ', 1)[0] in filter_out:
						flag = True
					else:
						flag = False
					"""elif model == 'model':
						if ('Query= ' + line.split('|')[0].replace('>', '')) in filter_out:
							flag = True
						else:
							flag = False"""
				if flag:
					if '>' in line:
						num += 1
						#new_line = line.strip() + '|' + hit_record[num].strip().replace('RNA', '') + '\n'
						#print(line.strip())
						'''
						ele=line.strip().split()
						#print(ele)
						ne=ele[0].split()
						ne[0]=ne[0]+'#SINE/'+hit_record[num].strip()
						ele[0]=' '.join(ne)
						
						line=' '.join(ele)
						#print(line)
						'''
						new_line = line.strip() + '|' + hit_record[num].strip() + '\n'
						rna_output_file.write(new_line)
					else:
						rna_output_file.write(line.strip() + '\n')
			#exit()


def tandem_repeat_finder(uid,out_genome_assembly_path):
	#path = os.path.abspath(os.path.dirname(os.getcwd()))
	#path = os.path.split(os.path.abspath(__file__))[0]
	path=os.getcwd()
	os.system('trf '
			  + uid + '/Step4_rna_output.fasta '
			  '2 5 7 80 10 10 2000 -d -h -l 6')
	if os.path.exists(path+'/Step4_rna_output.fasta.2.5.7.80.10.10.2000.dat'):
		os.system('mv '+path+'/Step4_rna_output.fasta.2.5.7.80.10.10.2000.dat '+out_genome_assembly_path)



def process_trf(input_trf_prob, out_genome_assembly_path, work_dir):
	trf_file = os.path.join(work_dir, 'Step4_rna_output.fasta.2.5.7.80.10.10.2000.dat')
	with open(trf_file, 'r')as trf_f:
		trf_list = []
		num = -1
		trf_lines = trf_f.readlines()
		flag_1 = False
		flag_2 = False
		for line in trf_lines[8:]:
			if 'Parameters: ' in line:
				num += 1
				trf_list.append([])
				flag_1 = True
			if len(line.split()) == 15 and len(line.split()) != 0 and flag_1:
				flag_2 = True
				trf_list[num].append(line.strip())
			if 'Sequence: ' in line and flag_1 and flag_2:
				flag_1 = False
				flag_2 = False
				
	input_f1 = os.path.join(out_genome_assembly_path, 'Step4_rna_output.fasta')
	input_f2 = os.path.join(out_genome_assembly_path, 'Step5_trf_output.fasta')
	
	with open(input_f1)as f:
		seq_length = []
		for line in f:
			if line[0] != '>':
				seq_length.append(len(line.strip()))
				
	with open(input_f1, 'r')as f1:
		with open(input_f2, 'w') as f2:
			counter = -1
			for line in f1:
				if line[0] == '>':
					flag = True
					counter += 1
					trf_length = 0
					if len(trf_list[counter]) != 0:
						for num_irf in range(len(trf_list[counter])):
							trf_length += len((trf_list[counter][num_irf].split()[-1]))
							if trf_length >= input_trf_prob * seq_length[counter]:
								flag = False
								break
					else:
						flag = True
				if flag:
					f2.write(line)


def save_to_fna(filename, input_sequences, input_title):
	header = '{}'.strip()
	index_1 = 0
	payload = []
	for seq in input_sequences:
		seq_header = header.format(input_title[index_1])
		payload.append(seq_header.strip() + '\n')
		payload.append(input_sequences[index_1] + '\n')
		index_1 += 1
	with open(filename, 'a') as file:
		file.writelines(payload)


def extend_seq(in_genome_assembly_path, out_genome_assembly_path):
	output_genome_sequence = read_genome_assembly(in_genome_assembly_path)
	filename_1 = out_genome_assembly_path+'/Step5_trf_output.fasta'

	title = []
	seq = []
	with open(filename_1, 'r')as f1:
		for line in f1:
			if line[0] == '>':
				title.append(line.strip())
				seq_id = line.split()[0].replace('>', '')
				#print(seq_id)
				#seq_id= re.sub('#.*','',seq_id) # Herui New add for header
				sine_start = int(line.split('|')[2].split(':')[1])
				sine_end = int(line.split('|')[3].split(':')[1])
				if abs(sine_end-sine_start) >= 200:
					seq.append(output_genome_sequence[seq_id][sine_start:sine_start+100].strip() +
							   output_genome_sequence[seq_id][sine_end-100:sine_end])
				elif 200 > abs(sine_end-sine_start) >= 100:
					seq.append(output_genome_sequence[seq_id][sine_start:sine_start+50].strip() +
							   output_genome_sequence[seq_id][sine_end-50:sine_end])
				else:
					seq.append(output_genome_sequence[seq_id][sine_start:sine_end])

	filename_2 = out_genome_assembly_path+'/Step6_irf_input.fasta'
	if os.path.exists(filename_2):
		modify_text(filename_2)
	save_to_fna(filename_2, seq, title)


def inverted_repeat_finder(out_genome_assembly_path, irf_path):
	path = os.path.split(os.path.abspath(__file__))[0]
	if irf_path=='': 
		irf_path=path+'/irf308.linux.exe'

	os.system(irf_path +' '+ out_genome_assembly_path +'/Step6_irf_input.fasta ' #shujun
		'2 3 5 80 10 20 500000 10000 -d -h -t4 74 -t5 493 -t7 10000')

	if os.path.exists('Step6_irf_input.fasta.2.3.5.80.10.20.500000.10000.dat'):
		os.system('mv Step6_irf_input.fasta.2.3.5.80.10.20.500000.10000.dat '+out_genome_assembly_path)
	


def process_irf(out_genome_assembly_path):
	irf_file = out_genome_assembly_path+'/Step6_irf_input.fasta.2.3.5.80.10.20.500000.10000.dat'

	with open(irf_file)as irf_f:
		irf_list = []
		num = -1
		irf_lines = irf_f.readlines()
		flag_1 = False
		flag_2 = False
		for line in irf_lines[8:]:
			if 'Parameters: ' in line:
				num += 1
				irf_list.append([])
				flag_1 = True
			if len(line.split()) == 19 and len(line.split()) != 0 and flag_1:
				flag_2 = True
				irf_list[num].append(line.strip())
			if 'Sequence: ' in line and flag_1 and flag_2:
				flag_1 = False
				flag_2 = False

	input_f1 = out_genome_assembly_path+'/Step5_trf_output.fasta'
	input_f2 = out_genome_assembly_path+'/Step6_irf_output.fasta'

	with open(input_f1)as f:
		seq_length = []
		tsd_length = []
		for line in f:
			if line[0] != '>':
				seq_length.append(len(line.strip()))
			else:
				tsd_length.append(int(line.split('|')[1].split(':')[1]))

	with open(input_f1, 'r')as f1:
		with open(input_f2, 'w') as f2:
			counter = -1
			flag = False
			for line in f1:
				if line[0] == '>':
					counter += 1
					if len(irf_list[counter]) != 0:
						for num_irf in range(len(irf_list[counter])):
							irf_left_s = int(irf_list[counter][num_irf].split()[0])
							irf_left_e = int(irf_list[counter][num_irf].split()[1])

							irf_length = int(irf_list[counter][num_irf].split()[2])

							irf_right_s = int(irf_list[counter][num_irf].split()[3])
							irf_right_e = int(irf_list[counter][num_irf].split()[4])
							if irf_length >= 10:
								flag = True
								break
					else:
						flag = False
				if not flag:
					f2.write(line)

#Improved readability, path stability
def cluster_sequences(out_genome_assembly_path,cpus):
	path = os.path.abspath(os.path.dirname(os.getcwd()))
	
	incdh  = os.path.join(out_genome_assembly_path, 'Step6_irf_output.fasta')
	outcdh = os.path.join(out_genome_assembly_path, 'Step7_cluster_output.fasta')
	
	cdhit_comm = f'cd-hit-est -i {incdh} -o {outcdh} -c 0.8 -M 0 -T {cpus} > /dev/null 2>&1'
	os.system(cdhit_comm)
	
	#input_fasta = os.path.join(out_genome_assembly_path, 'Step7_cluster_output.fasta')
	output_fasta = os.path.join(out_genome_assembly_path, 'Seed_SINE.fa')
	
	with open(outcdh, 'r')as f_1:
		with open(output_fasta, 'w')as f_2:
			num = 0
			for line in f_1:
				if line[0] == '>':
					new_line = f'>SINE_{num} ' + line.replace('>', '')
					f_2.write(new_line)
					num += 1
				else:
					f_2.write(line)

def re_process_figure(out_genome_assembly_path):
	comp_1 = []
	comp_2 = []
	
	cf1 = os.path.normpath(os.path.join(out_genome_assembly_path, 'Step2_extend_blast_input.fa'))
	cf2 = os.path.normpath(os.path.join(out_genome_assembly_path, 'Seed_SINE.fa'))
	#We wouldn't need to read through the whole files here if the headers were reformatted to have no spaces; just use pyfastx.Fasta.keys()
	#compf1seqs = pyfastx.Fasta()
	
	with open(cf1) as comp_f1:
		for comp_line1 in comp_f1:
			if comp_line1[0] == '>':
				comp_1.append(comp_line1.replace('>', '').replace('|', ' ').strip('\n'))
				
	with open(cf2) as comp_f2:
		for comp_line2 in comp_f2:
			if comp_line2[0] == '>':
				comp_2.append(comp_line2.replace(comp_line2.split()[0]+' ', '')
							  .split('|blast_s:')[0].replace('|', ' ').strip('\n'))
							  
	for num in range(len(comp_1)):
		if comp_1[num] not in comp_2:
			os.remove(+ '/Figures/' + f'MSA_{num}.png')
			os.remove(out_genome_assembly_path + '/Figures/' + f'profile_{num}.png')

def genome_annotate(in_genome_assembly_path, out_genome_assembly_path, in_nonredundant, rm_cpus): #reduce cpu number for RepeatMasker to avoid overutilization
	#path = os.path.abspath(os.path.dirname(os.getcwd()))
	print('Genome file: ' + in_genome_assembly_path, flush=True)
	print('Annotation results: ' + out_genome_assembly_path + '/RepeatMasker/', flush=True)
	if in_nonredundant == 'y':
		os.system('RepeatMasker -e ncbi -pa ' + str(rm_cpus) + ' -q -no_is -norna -nolow -div 40 '
				  '-lib  ' + out_genome_assembly_path + '/Seed_SINE.fa '
				  '-cutoff 225 ' + in_genome_assembly_path + ' '
				  '-dir ' + out_genome_assembly_path + '/RepeatMasker/ > /dev/null 2>&1')
	elif in_nonredundant == 'n':
		os.system('RepeatMasker -e ncbi -pa ' + str(rm_cpus) + ' -q -no_is -norna -nolow -div 40 '
				  '-lib  ' + out_genome_assembly_path + '/Step7_cluster_output.fasta '
				  '-cutoff 225 ' + in_genome_assembly_path + ' '
				  '-dir ' + out_genome_assembly_path + '/RepeatMasker/ > /dev/null 2>&1')


def sine_finder(genome_assembly_path, script_dir):
	#main_func()
	os.system('python3 ' + script_dir + '/SINEFinder.py ' + genome_assembly_path) #shujun


def save_to_fna_4(filename, input_sequences, input_id, input_direct, input_start, input_end):
	header = '{} {} {}:{}'.strip()
	index_1 = 0
	payload = []
	for seq in input_sequences:
		seq_header = header.format(input_id[index_1], input_direct[index_1], input_start[index_1], input_end[index_1])
		payload.append(seq_header.strip() + '\n')
		payload.append(input_sequences[index_1] + '\n')
		index_1 += 1
	with open(filename, 'a') as file:
		file.writelines(payload)


def process_sine_finder(genome_assembly_path, sine_finder_out, out_genome_assembly_path, pattern):
	output_genome_sequence = read_genome_assembly(genome_assembly_path)
	with open(sine_finder_out, 'r')as f1:
		finder_seq = []
		id = []
		direct = []
		start_position = []
		end_position = []
		lines = f1.readlines()
		flag = False
		for line in lines:
			if line[0] == '>':
				# mis_tsd = int(line.split()[3].split(';')[2].split('=')[1])
				id.append(line.split()[0])
				direct.append(line.split()[1])
				flag = True
				seq_id = line.split()[0].replace('>', '')
				s = int(line.split()[2].split(':')[0])
				e = int(line.split()[2].split(':')[1])
				tsd = int(line.split()[3].split(';')[0].split('=')[1])
				if s <= e:
					start = s + tsd - 30
					end = e - tsd + 50
				else:
					start = e + tsd - 30
					end = s - tsd + 50

			else:
				flag = False
			if flag:
				start_position.append(start+30)
				end_position.append(end-50)
				seq = output_genome_sequence[seq_id][start:end]
				finder_seq.append(seq)
	#print('sine_finder_out: ',sine_finder_out,' out_genome_assembly_path: ',out_genome_assembly_path)
	bname=os.path.basename(sine_finder_out)
	if not os.path.exists(out_genome_assembly_path+'/'+bname):
		#os.system('rm '+out_genome_assembly_path+'/'+bname)
		os.system('mv '+sine_finder_out+' '+out_genome_assembly_path)
	if pattern == 2 or pattern == 3:
		if os.path.exists(out_genome_assembly_path+'/Step1_extend_tsd_input_2.fa'):
			modify_text(out_genome_assembly_path+'/Step1_extend_tsd_input_2.fa')
		save_to_fna_4(out_genome_assembly_path+'/Step1_extend_tsd_input_2.fa', finder_seq, id, direct, start_position, end_position)



def ensure_path(path):
	if not os.path.exists(path):
		os.mkdir(path)


def check_hmm_finished(pre,idir):
	a=True
	if os.path.exists(idir):
		for filename in os.listdir(idir):
			if not os.path.getsize(idir+'/'+a) == 0:
				a=False
				break
	if a==False:
		print(pre+' already finished! Will skip to the next step!')
	return a

def check_finished(pre,arr,at):
	a=True
	for a in arr:
		if os.path.exists(a):
			if not os.path.getsize(a) == 0:
				a=False
			else:
				a=True
				break
		else:
			a=True
			break
	if a==False and at==False:
		print(pre+' already finished! Will skip to the next step!')
	if a==False and at==True:
		print(pre+' already finished! Will regenerate the result because -auto is not set.')
	if a==True and at==False:
		print(pre+' not finished! Will regenerate the result!')
	return a

def convert_ingenome(ingenome,odir):
	uid=uuid.uuid1().hex
	fdir=os.path.dirname(ingenome)
	ig=os.path.basename(ingenome)
	name, ext = os.path.splitext(ig)
	nf=name+'_'+uid+ext
	nd=odir+'/'+nf
	'''
	if fdir=='':
		nd=nf
	else:
		nd=fdir+'/'+nf
	'''
	#print('seqtk seq '+ingenome+' > '+nd)
	#exit()
	os.system('seqtk seq '+ingenome+' > '+nd)
	res=nd
	return res
			

def main_function():
	print('Please input the path of genomic sequence', flush=True) # print out message immediately
	input_pattern = args.mode
	input_genome_assembly_path = args.input_filename
	output_genome_assembly_path = args.output_filename
	if not os.path.exists(output_genome_assembly_path):
	  os.makedirs(output_genome_assembly_path)
	
	input_genome_assembly_path= convert_ingenome(input_genome_assembly_path, output_genome_assembly_path)
	
	ensure_path(output_genome_assembly_path)
	pre=os.path.splitext(input_genome_assembly_path)[0]

	input_sine_finder = pre+'-matches.fasta'

	input_hmm_e_value = args.hmmer_evalue
	input_blast_e_value = args.blast_evalue
	input_factor_length = args.length_factor
	input_factor_copy_number = args.copy_number_factor
	input_max_shift = args.shift
	input_max_gap = args.gap
	input_min_copy_number = args.copy_number
	input_num_alignments=args.num_alignments
	input_ani=args.animal
	input_bound = args.boundary
	input_figure = args.figure
	input_temd=args.temp_dir
	
	input_non_redundant = args.non_redundant
	input_auto=args.automatically_continue
	
	stop_at_tsd_searcher = args.stop_at_tsd
	
	if input_auto:
		at=False
	else:
		at=True

	cpus = args.threads
	rm_cpus = int(cpus // 4) #cpu number for RepeatMasker since each -pa value will invoke 4x rmblast processes

	# obtain program paths
	irf_path = args.irf_path
	rpm=args.RepeatMasker_enable
	start_time = time.time()
	
	print('************************************************************************', flush=True)
	print('*************************** AnnoSINE START! ****************************', flush=True)
	print('************************************************************************', flush=True)
	
	if input_pattern == 1 or input_pattern == 3:
		print('================ HMMER prediction has begun ==================', flush=True)
		step_1_1_file = os.path.join(output_genome_assembly_path, 'Step1_extend_tsd_input_1.fa')
		#if check_hmm_finished('S1_hmm_predict',work_dir+'/HMM_out'):
		if check_finished('Step1_process_hmm',[step_1_1_file], at) or at:
			t1=time.time()
			hmm_predict(input_genome_assembly_path, cpus, script_dir, work_dir,input_ani,input_hmm_e_value)
			t2=time.time()
			print('Step 1 mode-1::hmm_predict uses ',t2-t1,' s',flush=True)
			
			t1=time.time()
			process_hmm_output_3(input_hmm_e_value, input_genome_assembly_path, input_pattern, output_genome_assembly_path)
			t2=time.time()
			print('Step 1 mode-1::process_hmm_output_3 uses ',t2-t1,' s',flush=True)
			
			if os.path.exists(step_1_1_file):
				if os.path.getsize(step_1_1_file) == 0:
					print('HMM can not find any matched SINE! Program exit.')
					exit()
	
	if input_pattern == 2 or input_pattern == 3:
		print('================ Structure search has begun ==================', flush=True)
		step_1_2_file = os.path.join(output_genome_assembly_path, 'Step1_extend_tsd_input_2.fa')
		if check_finished('Step1_sine_part', [step_1_2_file], at) or at:
			t1=time.time()
			sine_finder(input_genome_assembly_path, script_dir)
			t2=time.time()
			print('Step 1 mode-2::sine_finder uses ',t2-t1,' s',flush=True)
			t1=time.time()
			process_sine_finder(input_genome_assembly_path, input_sine_finder, output_genome_assembly_path, input_pattern)
			t2=time.time()
			print('Step 1 mode-2::process_sine_finder uses ',t2-t1,' s',flush=True)
			
			if os.path.exists(step_1_2_file):
				if os.path.getsize(step_1_2_file) == 0:
					print('SINEfinder can not find any SINE! Program exit.')
					exit()
	
	#Data prep ends here
	
	step_1_complete_file = os.path.join(output_genome_assembly_path, 'Step1_extend_tsd_input.fa')
	if check_finished('Step1_merge_part', [step_1_complete_file], at) or at:
		t1=time.time()
		merge_tsd_input(input_pattern, output_genome_assembly_path)
		t2=time.time()
		print('Step 1::merge_tsd_input uses ',t2-t1,' s',flush=True)
	t2=time.time()
	print('Step 1 totally uses ',t2-start_time, ' s',flush=True)
	print('\n======================== Step 1 has been done ========================\n\n', flush=True)

	#here's the part of the program where crashes begin
	print('================ Step 2: TSD identification has begun ================', flush=True)
	t1=time.time()
	if check_finished('Step2_search_tsd',[output_genome_assembly_path+'/Step2_tsd.txt'], at) or at:
		search_tsd(output_genome_assembly_path, script_dir,cpus)
		#search_tsd(output_genome_assembly_path, script_dir,cpus)
		#print('')
	if check_finished('Step2_process_tsd',[output_genome_assembly_path+'/Step2_tsd_output.fa','Step2_extend_blast_input.fa'],at) or at:
		process_tsd_output(input_genome_assembly_path, output_genome_assembly_path)
	t2=time.time()
	print('Step 2 uses ',t2-t1,' s',flush=True)
	print('\n======================== Step 2 has been done ========================\n\n', flush=True)

	if stop_at_tsd_searcher:
		print('AnnoSINE was instructed to stop at TSD searching results.')
		print('This will break EDTA if used in that context; this option was meant as a developer aide and kept for its marginal utility.')
		return None

	'''
	I think the process blast output part is the rate limiting step, but memory use by minimap is also a concern
	'''
	print('================ Step 3: MSA implementation has begun ================', flush=True)
	t1=time.time()
	if check_finished('Step3_MSA', [output_genome_assembly_path+'/Step3_blast_output.out'], at) or at:
		multiple_sequence_alignment(input_blast_e_value, input_genome_assembly_path, output_genome_assembly_path,cpus,input_num_alignments, input_temd)
		
	if check_finished('Step3_process_MSA', [output_genome_assembly_path+'/Step3_blast_process_output.fa'], at) or at:
		#New code for faster, lower memory chunk processing
		get_updates_from_blasts(output_genome_assembly_path, input_genome_assembly_path, factor_copy = input_factor_copy_number, factor_length = input_factor_length, 
								min_copy_number = input_min_copy_number, max_gap = input_max_gap, max_shift = input_max_shift, threads = int(cpus))
		
	
	if check_finished('Step3_process_MSA_p2',[output_genome_assembly_path+'/Step4_rna_input.fasta'],at) or at:
		process_blast_output_2(output_genome_assembly_path)
	
	t2=time.time()
	print('Step 3 uses ',t2-t1,' s',flush=True)

	print('\n======================== Step 3 has been done ========================\n\n', flush=True)
	#exit()

	print('========= Step 4: RNA derived head identification has begun ==========', flush=True)
	t1=time.time()
	if check_finished('Step4_blast_rna',[output_genome_assembly_path+'/Step4_rna_output.out'],at) or at:
		blast_rna(output_genome_assembly_path, cpus, script_dir)
	if check_finished('Step4_process_blast_rna',[output_genome_assembly_path+'/Step4_rna_output.fasta'],at) or at:
		process_rna(output_genome_assembly_path)
	t2=time.time()
	print('Step 4 uses ',t2-t1,' s',flush=True)
	print('\n========================= Step 4 has been done =======================\n\n', flush=True)
	uid=uuid.uuid1().hex 
	if not os.path.exists(uid):
		os.makedirs(uid)
	os.system('cp '+output_genome_assembly_path+'/Step4_rna_output.fasta '+uid)
	print('=============== Step 5: Tandem repeat finder has begun ===============', flush=True)
	t1=time.time()
	if check_finished('Step5_trf',[work_dir+'/Step4_rna_output.fasta.2.5.7.80.10.10.2000.dat'],at) or at:
		#tandem_repeat_finder(output_genome_assembly_path)
		tandem_repeat_finder(uid,output_genome_assembly_path)
	if check_finished('Step5_process_trf',[output_genome_assembly_path+'/Step5_trf_output.fasta'],at) or at:
		process_trf(0.5, output_genome_assembly_path, work_dir)
	t2=time.time()
	print('Step 5 uses ',t2-t1,' s',flush=True)
	os.system('rm -rf '+uid)
	print('\n======================== Step 5 has been done ========================\n\n', flush=True)

	print('=============== Step 6: Inverted repeat finder has begun =============', flush=True)
	t1=time.time()
	if check_finished('Step6_extend_seq',[output_genome_assembly_path+'/Step6_irf_input.fasta'],at) or at:
		extend_seq(input_genome_assembly_path, output_genome_assembly_path)
	if check_finished('Step6_irf',[output_genome_assembly_path+'/Step6_irf_input.fasta.2.3.5.80.10.20.500000.10000.dat'],at) or at:
		inverted_repeat_finder(output_genome_assembly_path, irf_path)
	if check_finished('Step6_process_irf',[output_genome_assembly_path+'/Step6_irf_output.fasta'],at) or at:
		process_irf(output_genome_assembly_path)
	t2=time.time()
	print('Step 6 uses ',t2-t1,' s',flush=True)
	print('\n========================= Step 6 has been done =======================\n\n', flush=True)

	print('=============== Step 7: Sequences clustering has begun ===============', flush=True)
	t1=time.time()
	if check_finished('Step7_cluster_seq',[output_genome_assembly_path+'/Step7_cluster_output.fasta'],at) or at:
		cluster_sequences(output_genome_assembly_path,cpus)
	t2=time.time()
	print('Step 7 uses ',t2-t1,' s',flush=True)
	print('\n======================== Step 7 has been done ========================\n\n', flush=True)

	if rpm==1:
		print('================= Step 8: Genome annotation has begun ================', flush=True)
		t1=time.time()
		if input_figure == 'y':
			dirs = output_genome_assembly_path+'/Figures/'
			if not os.path.exists(dirs):
				os.makedirs(dirs)
			re_process_figure(output_genome_assembly_path)
		genome_annotate(input_genome_assembly_path, output_genome_assembly_path, input_non_redundant, rm_cpus) #shujun
		t2=time.time()
		print('Step 8 uses ',t2-t1,' s',flush=True)
		print('\n========================= Step 8 has been done =======================\n\n', flush=True)

	end_time = time.time()
	print('Total running time: ', end_time - start_time, 's', flush=True)
	print('************************************************************************', flush=True)
	print('************************** AnnoSINE COMPLETE! **************************', flush=True)
	print('************************************************************************', flush=True)
	os.system('rm '+input_genome_assembly_path)
	
	return None


if __name__ == '__main__':
	main_function()
