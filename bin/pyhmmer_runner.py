import sys
import os
import pyfastx
import pyhmmer
import multiprocessing
import io

class pyhmmer_manager:
	def __init__(self, genome_file):
		self.fa = pyfastx.Fasta(genome_file)
		self.hmms = None
		self.sequence = None
		
	def load_model_from_file(self, hmm_file):
		ar = open(hmm_file, mode = 'rb')
		hmm_text = ar.read()
		hmm_io = io.BytesIO(hmm_text)
	
		#Faster to load the file via binary read and then load pyhmmer in memory than to directly use pyhmmer method. Odd, but it works.
		with pyhmmer.plan7.HMMFile(hmm_io) as fh:
			self.hmms = list(fh)
		
	#Use pyfastx to load a sequence and digitize it for search
	def load_sequence(self, sequence_id):
		seqid = sequence_id.encode(encoding = "ascii")
		self.sequence = pyhmmer.easel.TextSequence(name = seqid, sequence = self.fa[sequence_id].seq)
		self.sequence = self.sequence.digitize(pyhmmer.easel.Alphabet.dna())
		self.sequence = [self.sequence]
			
	#Run nhmmer, write results
	def nhmmer(self):
		results = io.BytesIO()
		for hit in pyhmmer.hmmer.nhmmer(self.hmms, self.sequence, cpus = 1):
			hit.write(results, format = 'targets', header = False)
			
		return results.getvalue().decode()
	
def hmm_init(gf, hmms):
	global genome_file
	global hmm_paths
	genome_file = gf
	hmm_paths = hmms	
	
def one_process(args):
	seqid, hmm = args
	print(f'Process {seqid} vs {hmm} begin')
	mn = pyhmmer_manager(genome_file)
	mn.load_model_from_file(hmm_paths[hmm])
	mn.load_sequence(sequence_id = seqid)
	results = mn.nhmmer()
	print(f'Process {seqid} vs {hmm} end')
	return results
		
def process_pyhmmer(genome_file, hmm_model_dir, output_file, threads = 1):
	if not os.path.exists(f'{genome_file}.fxi'):
		print('Indexing genome_file')
		
	fa = pyfastx.Fasta(genome_file, build_index = True)
	sequences = list(fa.keys())	
	
	#Sort hmm models from longest to shortest as an approximation for runtime.
	hmm_models = [os.path.join(hmm_model_dir, h) for h in os.listdir(hmm_model_dir) if h.endswith('.hmm')]
	hmm_models = sorted(hmm_models, key = os.path.getsize, reverse = True)
	hmm_models = {os.path.basename(h):h for h in hmm_models}
	
	#This could end up being a quite large list.
	args = []
	for h in hmm_models:
		for s in sequences:
			args.append((s, h,))
			
	with open(output_file, 'w') as out:
		with multiprocessing.Pool(threads, initializer = hmm_init, initargs = (genome_file, hmm_models,)) as pool:
			for result in pool.imap_unordered(one_process, args):
				if len(result) > 0:
					out.write(result)
					
		
def hmm_output_cleaner(hmm_results_file, threshold_hmm_e_value = 1e-10):
	#fa = pyfastx.Fasta(genome_file)
	
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
								'strand': 'C'}
					
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
	
'''	
md = 'AnnoSINE_v2/hmm_family_seq_easy/'
output = 'test_pyhmmer.txt'
genome_file = sys.argv[1]
threads = 10

#process_pyhmmer(genome_file, md, output, threads)

hmm_output_cleaner('test_pyhmmer.txt')
'''