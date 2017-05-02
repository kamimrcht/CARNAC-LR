
#/local/python/3.3.2/bin/python
import sys
import os
import shlex, subprocess


def get_read(sequencefile,offset):
	sequencefile.seek(offset)
	read=""
	line=sequencefile.readline()
	if not line: 
		print("cannot read read at offset", offset)
		exit(1)
	read+=line#include header
	read+=sequencefile.readline()#include sequence
	if read[0]=='>': return read
	if read[0]!='@': 
		print("read offset", offset, "does not start with @ or >")
		exit(1)
	read+=sequencefile.readline()#include header2
	read+=sequencefile.readline()#include quality
	return read

def index_bank_offsets(bank_file_name, namesToOffset):
	read_offsets= []

	if "gz" in bank_file_name:
		sequencefile=gzip.open(bank_file_name,"r")
	else: 
		sequencefile=open(bank_file_name,"r")
	# i=1
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
	# t=0
	while True:
		offset=sequencefile.tell()
		#~ print(offset)
		line=sequencefile.readline()
		index = "_".join(line.rstrip().split('_')[:2])[1:] + "_"
		#~ print(index, offset)
		
		if not line: break
		# t+=1
		read_offsets.append(offset)
		namesToOffset[index] = len(read_offsets) - 1
		for i in range(linesperread-1): line=sequencefile.readline()
	# print "max=",t
	sequencefile.close()
	return read_offsets

if len(sys.argv) > 2:
	clustersFileName = sys.argv[1]
	readsFileName = sys.argv[2]
	clustersFile = open(clustersFileName, 'r')
	readsFile = open(readsFileName, 'r')
	namesToOffset = dict()
	# indexing read file
	reads_offsets = index_bank_offsets(readsFileName, namesToOffset)
	# getting reads from clusters, 1 file per cluster
	
	i = 0
	out = open("reads_cluster2fasta.fa", 'w')
	out2 = open("clusters.truth", 'w')
	for line in clustersFile:  # 1 line = 1 cluster
		readsList = line.rstrip().split(' ')
		
		if len(readsList) > 2 and len(readsList) < 350:
			toW = ""
			for read in readsList:
				toW += str(i) + " "
				#~ print("*", read, namesToOffset[read])
				#~ print ()
				#~ print(reads_offsets[namesToOffset[read]])
				out.write(get_read(readsFile, reads_offsets[namesToOffset[read]]))
				i += 1
				#~ out.write(get_read(readsFile, reads_offsets[int(read)]))
			out2.write(toW + '\n')
	out.close()
	out2.close()
		    
else:
	print("Usage : ./clusters_to_fasta.py <cluster_file> <read_file>")
