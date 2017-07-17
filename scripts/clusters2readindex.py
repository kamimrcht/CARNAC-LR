#!/usr/bin/env python3
import sys
import os

if (len(sys.argv) > 2):

	with open(sys.argv[1]) as infile:  # read file
		lineToRead = dict()
		lineNumber = 0
		for rline in infile:
			line = rline.rstrip()
			if '>' in line:
				read = "_".join(line.rstrip().split('_')[:2])[1:] + "_"
				lineToRead[lineNumber] = read
				lineNumber += 1
		
	out = open("final_g_clusters_reads_id.txt", 'w')
	with open(sys.argv[2]) as infile:
		readIndexToCluster = dict()
		cluster = 0
		for rline in infile:
				toW = ""
				line = rline.rstrip()
				reads = line.split(' ')
				for r in reads:
					toW += lineToRead[int(r)] + " "
				out.write(toW + '\n')
else:
	print("python3 clusters2readindex.py reads.fa final_g_clusters.txt")
		
		
