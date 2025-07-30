import pandas as pd
import os
import sqlite3


class blastdb:
	def __init__(self, path):
		self.dbpath = path
		self.conn = None
		self.curs = None
		
	def open(self):
		if self.conn is None:
			self.conn = sqlite3.connect(self.dbpath)
			self.curs = self.conn.cursor()
		
	def close(self):
		if self.conn is not None:
			self.curs.close()
			self.conn.close()
			
	def init(self):
		'''
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
		'''
	
		blast_schema = ["qname TEXT",
			"tname TEXT",
			"alen INTEGER",
			"qstart INTEGER",
			"qend INTEGER",
			"tstart INTEGER", 
			"tend INTEGER",
			"bitscore REAL"]
		blast_schema = ', '.join(blast_schema)
		table_create = f'CREATE TABLE blast_records ({blast_schema}, PRIMARY KEY (qname, tname, tstart, tend))'
		
		if self.curs is not None:
			self.curs.execute('DROP TABLE IF EXISTS blast_records')
			self.conn.commit()
			self.curs.execute(table_create)
			self.conn.commit()
			
	def insert(self, records):
		bindings = ', '.join(['?'] * 8)
		self.curs.executemany(f'INSERT OR IGNORE INTO blast_records VALUES ({bindings})', records)
		
			

def blast_outputs_to_database(output_directory = '.'):
	dbout = os.path.join(output_directory, 'Step3_blast_sqlitedb.db')
	dbm = blastdb(dbout)
	dbm.open()
	dbm.init()
	
	blast_schema_pandas = {'qname':'str',
							'tname':'str',
							'pctid':'float32',
							'alen':'int32',
							'nonmatch':'int32',
							'gaps':'int32',
							'qstart':'int32',
							'qend':'int32',
							'tstart':'int32',
							'tend':'int32',
							'evalue':'float32',
							'bitscore':'float32'}
							
	#Drop pctid, nonmatch, gaps, evalue
	cols_to_use = [0, 1, 3, 6, 7, 8, 9, 11]
	
	insert_command = ''
	
	files = [os.path.join(output_directory, f) for f in os.listdir(output_directory) if 'olap' in f]
	for f in files:
		for df in pd.read_csv(f, sep = '\t', header = None, dtype = blast_schema_pandas, usecols = cols_to_use, chunksize = 5_000_000):
			df = list(df.itertuples(index=False, name=None))
			dbm.insert(df)
		
	dbm.conn.commit()
	dbm.close()
	
blast_outputs_to_database()