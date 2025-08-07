import sys
import os
import re

nonmatch_pattern = re.compile(r"NM:i:(\d+)")
cigar_pattern = re.compile(r"cg:Z:([A-Za-z0-9]+)")

for line in sys.stdin:
	#Get input data from paf output streaming
	qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, alen, mapq, NM, ms, AS, nn, tp, cm, s1, de, rl, cg, extra = line.strip().split('\t')
	strand = strand=='+'
	qlen = int(qlen)
	qstart = int(qstart)
	qend = int(qend)
	tlen = int(tlen)
	tstart = int(tstart)
	nmatch = int(nmatch)
	alen = int(alen)
	