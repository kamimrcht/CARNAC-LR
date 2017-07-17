#!/usr/bin/env python3
import sys
import os

if (len(sys.argv) > 2):


	readToTruePairs = dict()
	readToResultPairs = dict()
	seen = set()
	#~ a11 = 0  # when the two clustering put a pair in the same partition
	#~ a01 = 0  # when the true clustering place a pair in the same partition, but the result clustering separates nodes of the pair in different partitions
	#~ a10 = 0  # when the result clustering places a pair in the same partition, but the result clustering separates nodes of the pair in different partitions
	continu = False
	pairsConnectedInResult = set()
	with open(sys.argv[2]) as infile:
		# result clusters : parse and associate to a read its pairs
		for rline in infile:
				line = rline.rstrip()
				reads = line.split(' ')
			#~ if len(reads) > 10 :
				for r in reads:
					seen.add(r)
					for pair in reads:
						seen.add(pair)
						if r != pair:  # r and pair are two nodes
							pairsConnectedInResult.add((min(r,pair),max(r,pair,)))
							if not r in readToResultPairs.keys():
								readToResultPairs[r] = set()
							readToResultPairs[r].add(pair)
							if not pair in readToResultPairs.keys():
								readToResultPairs[pair] = set()
							readToResultPairs[pair].add(r)
	pairsConnectedInTrue = set()
	pairsSeen = set()
	readsInRealClusters = set()  # TO KEEP ONLY READS FROM UNIQUE LOCI
	with open(sys.argv[1]) as infile:
		# true clusters : parse and associate to a read its pairs
		for rline in infile:
				line = rline.rstrip()
				reads = line.split(' ')
			#~ if len(reads) > 10 :
				for r in reads:
					readsInRealClusters.add(r)  # TO CHECK ONLY THIS SUBSET OF READS AFTERWARDS
					for pair in reads:
						if r != pair:
							a = (min(r,pair),max(r,pair,))
							pairsConnectedInTrue.add(a) 
							if a in pairsConnectedInResult:
									#~ a11 += 1
									pairsSeen.add((min(r,pair),max(r,pair,)))
							#~ else:
								#~ a10 +=1
					continu = False
					for pair in reads:
						if r in seen or pair in seen:
							continu = True
							break
					for pair in reads:
						if continu:
							if r != pair:
								if not r in readToTruePairs.keys():
									readToTruePairs[r] = set()
								readToTruePairs[r].add(pair)
								if not pair in readToTruePairs.keys():
									readToTruePairs[pair] = set()
								readToTruePairs[pair].add(r)
	#~ a01 = len(pairsConnectedInTrue) - len(pairsSeen)

	
	# a11
	# a10
	# a01
	
	recall = []
	precision = []
	for read in readToTruePairs.keys():
		#~ if read in readsInRealClusters:   # WE WORK ONLY WITH READS FROM UNIQUE LOCI
			localRecall = 0
			localPrecision = 0
			if read in readToResultPairs:
				if len(readToTruePairs[read]) > 0:
						if len(readToResultPairs[read]) > 0:
							localRecall = len(readToResultPairs[read].intersection(readToTruePairs[read]))*1.0/len(readToTruePairs[read])  #  recall is how many pairs were found among all to find
							localPrecision = len(readToResultPairs[read].intersection(readToTruePairs[read]))*1.0/len(readToResultPairs[read])  # precision if how much predicted pairs are true pairs
						else:
							localRecall = 0
				else:  # when a read is truly alone
					localRecall = 1
					
					if (len(readToResultPairs[read]) != 0):
						localPrecision = 0
				recall.append(localRecall)
				precision.append(localPrecision)
			else:  # when no read was linked to this read in the result
				recall.append(localRecall)

	print("recall ", sum(recall)*1.0/float(len(recall)), " precision ", sum(precision)*1.0/float(len(precision)))
	#~ print("Jaccard index ", a11*1.0/(a11+a01+a10))
else:
	print("validate_with_pairs clusters.truth clusters.result")
		
		
