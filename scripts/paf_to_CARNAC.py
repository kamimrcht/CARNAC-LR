#!/usr/bin/env python3
import sys
import os

# example of paf format
#ch114_read985_template_fail_BYK_CB_ONT_1_FAF05104_A     1692    1       1686    +       ch114_read985_template_fail_BYK_CB_ONT_1_FAF05104_A     1692    1       1686    1685    1685    255     cm:i:303
#ch114_read985_template_fail_BYK_CB_ONT_1_FAF05104_A     1692    59      1598    +       ch487_read8370_template_fail_BYK_CB_ONT_1_FAF05104_A    3100    1194    2784    434     1590    255     cm:i:39
# 0:50-25-50.000000 27-32-64.000000, ignores the positions information



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
	while True:
		offset=sequencefile.tell()
		line=sequencefile.readline()
		index = line.rstrip()[1:] 
		if not line: break
		read_offsets.append(offset)
		
		namesToOffset[index] = len(read_offsets) - 1
		for i in range(linesperread-1): line=sequencefile.readline()
	sequencefile.close()
	return read_offsets




if (len(sys.argv) > 3):
	readToRecruited = dict()
	namesToOffset = dict()
	index_bank_offsets(sys.argv[2], namesToOffset)
	with open(sys.argv[1]) as infile:  # paf file
		for rline in infile:
			line = rline.rstrip()
			read = line.split("\t")[0]
			index = namesToOffset[read]
			readRecruited = line.split("\t")[5]
			indexRecruited = namesToOffset[readRecruited]
			if index in readToRecruited:
				if index != indexRecruited:
					readToRecruited[index].add(indexRecruited)
			else:
				if index != indexRecruited:
					readToRecruited.setdefault(index, set())
					readToRecruited[index].add(indexRecruited)
	out = open(sys.argv[3], 'w')
	for read, recruited in readToRecruited.items():
		toWrite = str(read) + ":"
		for r in recruited:
			toWrite += str(r) + "-0-0.0 "
		toWrite += "\n"
		out.write(toWrite)
else:
	print("Usage: python3 paf_to_CARNAC.py file.paf reads.fasta/q(or gz) output.txt")
