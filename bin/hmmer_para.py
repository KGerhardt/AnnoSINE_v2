import sys
import os

import multiprocessing
import subprocess
import shutil

import sqlite3


def prep_hmmsearch(genome_file, hmm_directory, working_dir = '.', threads = 1):
	base = os.path.join(working_dir, 'annosine_hmmsearch')
	gendir = os.path.join(base, 'genomes')
	hmmdir = os.path.join(base, 'hmmsearch_output')
	#work_dir = os.path.join(working_dir, 'annosine_hmmsearch', 'genomes')
	if not os.path.exists(base):
		os.mkdir(base)
	
	if not os.path.exists(gendir):
		os.mkdir(gendir)
	else:
		for f in os.listdir(gendir):
			os.remove(os.path.join(gendir, f))
		
	if not os.path.exists(hmmdir):
		os.mkdir(hmmdir)
		
	print(f'Chunking input genome into {threads} partitions...')
	pfx_split = f'pyfastx split {genome_file} -n {threads} -o {gendir}'
	subprocess.run(pfx_split.split())
	print('Chunking complete!')
	
	conn = sqlite3.connect(f'{genome_file}.fxi')
	curs = conn.cursor()
	
	#This is needed to adjust the e-value scaling of the HMM searches later on
	total_megabases = curs.execute('SELECT sum(slen) FROM seq').fetchone()[0] / 1_000_000
	
	curs.close()
	conn.close()
		
	outs = [os.path.join(gendir, f) for f in os.listdir(gendir) if os.path.getsize(os.path.join(gendir, f)) > 0]
		
	hmm_models = [os.path.join(hmm_directory, h) for h in os.listdir(hmm_directory) if h.endswith('.hmm')]
	hmm_models = sorted(hmm_models, key = os.path.getsize, reverse = True)
	
	return outs, hmm_models, total_megabases
	
#Use of the original HMMer method is possibly unneeded with better pyhmmer implementation but
#I want to check pyhmmer nhmmer and comparative speed more carefully first_loc
#There will be some cases where this approach is faster, can't be sure about which always.
def run_one_hmm(args):
	hmm_mod, ingen, mythreads, z_adj = args
	hmm_base = os.path.basename(hmm_mod)
	hmm_base = hmm_base.replace('.hmm', '')
	
	outfile = ingen.replace('genomes', 'hmmsearch_output')
	outfile = f'{outfile}_{hmm_base}.txt'
	
	hmm_one = f'nhmmer --cpu {mythreads} --dna -Z {z_adj} --tblout {outfile} {hmm_mod} {ingen}'
	hmm_one = hmm_one.split()
	
	subprocess.run(hmm_one, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
	#subprocess.run(hmm_one)
	
	return outfile
	
def parse_hmmfile(file, outhandle):
	with open(file) as infh:
		for line in infh:
			if not line.startswith('#'):
				outhandle.write(line)
	
	os.remove(file)

def run_hmmsearch(queries, targets, output, z, threads = 1):
	hmmer_thread_saturation = 4
	
	if threads <= hmmer_thread_saturation:
		thread_chunk = threads
	else:
		thread_chunk = hmmer_thread_saturation
		
	args = []
	for t in targets:
		for q in queries:
			args.append((t, q, thread_chunk, z,))
	
	print(f'Executing HMM search using {threads//thread_chunk} blocks of {thread_chunk} processors.')
	total_chunks = len(args)
	print(f'There are {total_chunks} chunks to process.')
	
	print(f'The search is arranged from largest HMM models to smallest in order to maximize CPU efficiency')
	print(f'Progress reporting will generally speed up as this step proceeds.')
	
	chunk_size = max([int(round(total_chunks/100)), 1])
	processed_chunks = 0
	
	with open(output, 'w') as out:
		#Between 1 and # threads / 6
		with multiprocessing.Pool(threads // thread_chunk) as pool:
			for result in pool.imap_unordered(run_one_hmm, args):
				processed_chunks += 1
				if processed_chunks % chunk_size == 0:
					print(f'HMM search is {round(100 * processed_chunks / total_chunks)}% complete')
				
				parse_hmmfile(result, out)
		
def hmm_output_cleaner(hmm_results_file, threshold_hmm_e_value = 1e-10):
	
	'''
	print('Processing the hmm prediction ...', flush=True)
	family_count = {}
	family_name = []
	update_hmm_record = []
	dir_file = os.listdir(work_dir + '/HMM_out/') #shujun
	out_file = glob.glob(work_dir + '/HMM_out/' + '*.out') #shujun
	for a in range(len(out_file)):
		if out_file[a] != '.DS_Store':
			#list_pre = process_hmm_output_1(out_file[a], threshold_hmm_e_value, script_dir)[0] #shujun
			def process_hmm_output_1(threshold_hmm_e_value, script_dir):
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
	'''
	family_count = {}
	update_hmm_record = []
		
	#This is basically process_hmm_output_1
	with open(hmm_results_file) as fh:
		for line in fh:
			segs = line.strip().split()
			evalue = float(segs[12])
			
			if evalue < threshold_hmm_e_value:
				query_seq = segs[0]
				family = segs[2]
				family = family.replace('-aln-stockholm', '')
				strand = segs[11] == "+"
				
				first_loc = int(segs[6])
				second_loc = int(segs[7])
				
				if strand:
					next_record = {'start': first_loc - 1,
								'end': second_loc,
								'e_value': evalue,

								'family': family,
								'id': query_seq,
								'strand': '+'}
				else:
					next_record = {'start': second_loc - 1,
								'end': first_loc,
								'e_value': evalue,

								'family': family,
								'id': query_seq,
								'strand': '-'}
					
				#And this is basically process_hmm_output_2
				update_hmm_record.append(next_record)
				if family not in family_count:
					family_count[family] = 0
				family_count[family] += 1
		
	#Need to sort output by hit family, sequence id, start index to resemble previous method
	family_name = []
	update_hmm_record = sorted(update_hmm_record, key=lambda x: (x['family'], x['id'],  x['start']))
	for r in update_hmm_record:
		family_name.append(r['family'])

	return update_hmm_record, family_name, family_count
	
					
#f = sys.argv[1]
#hmm_dir = sys.argv[2]

#inputs, hmms, z_adj = prep_hmmsearch(f, hmm_dir, 'test_pfx', 20)
#run_hmmsearch(inputs, hmms, 'test_pfx/hmm_mp.txt', z_adj, 20)

