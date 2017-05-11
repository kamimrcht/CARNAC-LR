#!/usr/bin/env python3
import sys
import os

#ch114_read985_template_fail_BYK_CB_ONT_1_FAF05104_A     1692    1       1686    +       ch114_read985_template_fail_BYK_CB_ONT_1_FAF05104_A     1692    1       1686    1685    1685    255     cm:i:303
#ch114_read985_template_fail_BYK_CB_ONT_1_FAF05104_A     1692    59      1598    +       ch487_read8370_template_fail_BYK_CB_ONT_1_FAF05104_A    3100    1194    2784    434     1590    255     cm:i:39
# 0:50-25-50.000000 27-32-64.000000, ignores the positions information

if (len(sys.argv) > 2):
	readToIndex = dict()
	index = 0
	readToRecruited = dict()
	with open(sys.argv[1]) as infile:  # reads fasta file
		for rline in infile:
			if '@' in rline or '>' in rline:
				line = rline.rstrip().split("t")
				read = line[0][1:-1]
				#~ print(read)
				readToIndex[read] = index
				index += 1
	with open(sys.argv[2]) as infile:  # minimap output
		for rline in infile:
			line = rline.rstrip()
			read = line.split("\t")[0].split('t')[0][:-1]
			recruited = line.split("\t")[5].split('t')[0][:-1]
			index = readToIndex[read]
			indexRecruited = readToIndex[recruited]
			if index in readToRecruited:
				if index != indexRecruited:
					readToRecruited[index].add(indexRecruited)
			else:
				if index != indexRecruited:
					#~ print(index, indexRecruited)
					readToRecruited.setdefault(index, set())
					readToRecruited[index].add(indexRecruited)
	out = open("minimap_src_format.txt", 'w')
	out.write("#\n#\n#\n")
	for read, recruited in readToRecruited.items():
		toWrite = str(read) + ":"
		for r in recruited:
			toWrite += str(r) + "-0-0.0 "
		toWrite += "\n"
		out.write(toWrite)
