
#/local/python/3.3.2/bin/python
import sys
import os
import shlex, subprocess
from subprocess import Popen, PIPE, STDOUT,call



def get_read(sequencefile,offset, fasta):
	if fasta:
		offset*=2
	else:
		offset*=4
	sequencefile.seek(offset)
	read=""
	line=sequencefile.readline()
	if not line: 
		print("cannot read read at offset", offset)
		exit(1)
	read+=line#include header
	read+=sequencefile.readline()#include sequence
	if fasta: return read
	read+=sequencefile.readline()#include header2
	read+=sequencefile.readline()#include quality
	return read

def index_bank_offsets(bank_file_name, namesToOffset):
	read_offsets= []

	if "gz" in bank_file_name:
		sequencefile=gzip.open(bank_file_name,"r")
	else: 
		sequencefile=open(bank_file_name,"r")
	line=sequencefile.readline()
	if not line: 
		print("Can't open file", bank_file_name)
		exit(1)
	if line[0]!='@' and line[0]!='>': 
		print("File", bank_file_name, "not correctly formatted")
		exit(1)

	linesperread=2 #fasta by default
	if line[0]=='@': linesperread=4 # fastq
	
	sequencefile.seek(0)
	while True:
		offset=sequencefile.tell()
		line=sequencefile.readline()
		index = "_".join(line.rstrip().split('_')[:2])[1:] 
		if not line: break
		read_offsets.append(offset)
		offsetToName[index] = len(read_offsets) - 1
		for i in range(linesperread-1): line=sequencefile.readline()
	sequencefile.close()
	return read_offsets

if len(sys.argv) > 2:
	clustersFileName = sys.argv[1]
	readsFileName = sys.argv[2]
	path=""
	minSize = 1
	if len(sys.argv) > 3:
		path = sys.argv[3] + "/"
		if not os.path.isdir(path):
			cmdMkdir = 	"mkdir " + path
			subprocess.check_output(['bash','-c', cmdMkdir])
		else:
			cmdRm = "rm " + path + "cluster_*.fasta"
			try:
				DEVNULL = open(os.devnull, 'r+b', 0)
				p = subprocess.call(['bash','-c', cmdRm], stderr=DEVNULL)
			except subprocess.CalledProcessError:
				pass
		if len(sys.argv) > 4:
			minSize = int(sys.argv[4])
	clustersFile = open(clustersFileName, 'r')
	readsFile = open(readsFileName, 'r')
	readsFile.seek(0)
	line=readsFile.readline()
	fasta = True
	if not line: 
		print("cannot read read file")
		exit(1)
	if line[0]!='>':
		if line[0]=='@':
			fasta = False
		else:
			print("reads in wrong format")
			exit(1)
	offsetToName = dict()
	nameToRead = dict()
	# indexing read file
	reads_offsets = index_bank_offsets(readsFileName, offsetToName)
	
	numCluster = 0
	for line in clustersFile:  # 1 line = 1 read and its cluster
		readsIndex = line.rstrip().split(' ')
		if len(readsIndex) >= minSize:
			out = open(path + "cluster_" + str(numCluster) + ".fasta", 'w')
			for r in readsIndex:
				read = get_read(readsFile, int(r), fasta)
				out.write(read)
			out.close()
			numCluster += 1
		    
else:
	print("Usage : python3 CARNAC_to_fasta.py <CARNAC_cluters.txt> <original_read_file.fa> [/path/to/write/files] [cluster_min_size]")
