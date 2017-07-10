#!/usr/bin/env python3
import sys
import os

if (len(sys.argv) > 2):
	threshold = - 1
	if (len(sys.argv) > 3):
		threshold = int(sys.argv[3])
		

	readToTruePairs = dict()
	readToResultPairs = dict()
	seen = set()
	a11 = 0  # when the two clustering put a pair in the same partition
	a01 = 0  # when the true clustering place a pair in the same partition, but the result clustering separates nodes of the pair in different partitions
	a10 = 0  # when the result clustering places a pair in the same partition, but the result clustering separates nodes of the pair in different partitions
	a00 = 0  # pairs of nodes that both clustering do not put in the same partition


	pairsConnectedInTrue = set()
	pairsUnconnectedInTrue = set()
	pairsUnconnectedInResult = set()
	#~ nbNodes = 0
	nodes = set()
	print("Getting true clusters partition")
	with open(sys.argv[1]) as infile:  # first we need to know how many nodes there are
		
		for rline in infile:
			line = rline.rstrip()
			reads = line.split(' ')
			if len(reads) > threshold :
				for r in reads:
					nodes.add(r)
					for pair in reads:
							if r != pair:
								a = (min(r,pair),max(r,pair,))
								pairsConnectedInTrue.add(a)

	pairsUnconnectedInTrue = set()
	for n in nodes:
		for n2 in nodes:
			if n < n2:
				a = (n,n2,)
				if not a in pairsConnectedInTrue:
					pairsUnconnectedInTrue.add((n,n2,))
				pairsUnconnectedInResult.add((n,n2,))
				

	#~ pairsUnconnectedInResult = pairsUnconnectedInTrue
	#~ print("Getting true clusters partition")
	#~ with open(sys.argv[1]) as infile:
		#~ # true clusters : parse and associate to a read its pairs (all nodes are present in this file)
		#~ for rline in infile:
			#~ line = rline.rstrip()
			#~ reads = line.split(' ')
			#~ if len(reads) > threshold :
				#~ for r in reads:
					
					#~ for pair in reads:
						#~ if r != pair:
							#~ a = (min(r,pair),max(r,pair,))
							#~ pairsConnectedInTrue.add(a)
							#~ pairsUnconnectedInTrue.discard(a)


	
	print("Getting result clusters partition and comparison")
	pairsConnectedInResult = set()
	with open(sys.argv[2]) as infile:
		# result clusters : parse and associate to a read its pairs
		for rline in infile:
			line = rline.rstrip()
			reads = line.split(' ')
			if len(reads) > threshold:
				for r in reads:
					seen.add(r)
					for pair in reads:
						seen.add(pair)
						if r != pair:  # r and pair are two nodes
							a = (min(r,pair),max(r,pair,))
							pairsConnectedInResult.add(a)
							pairsUnconnectedInResult.discard(a)
							if a in pairsConnectedInTrue:
									a11 += 1
									#~ pairsSeen.add((min(r,pair),max(r,pair,)))
							else:
								a01 +=1
	for pair in pairsUnconnectedInResult:
		if pair in pairsUnconnectedInTrue:
			a00 += 1
		else:
			a10 += 1
								

	#~ a11 = len(pairsConnectedInTrue.intersection(pairsConnectedInResult))
	#~ a01 = len(pairsUnconnectedInTrue.difference(pairsUnconnectedInResult))
	#~ a00 = len(pairsUnconnectedInTrue.intersection(pairsUnconnectedInResult))
	#~ a10 = len(pairsUnconnectedInResult.difference(pairsUnconnectedInTrue))
			
	
	#~ pairsSeen = set()
	
	#~ a01 = len(pairsConnectedInTrue.difference(pairsSeen))

	
	# a11
	# a10
	# a01
	
	#~ recall = []
	#~ precision = []
	#~ for read in readToTruePairs.keys():
		#~ localRecall = 0
		#~ localPrecision = 0
		#~ if read in readToResultPairs:
			#~ if len(readToTruePairs[read]) > 0:
					#~ if len(readToResultPairs[read]) > 0:
						#~ localRecall = len(readToResultPairs[read].intersection(readToTruePairs[read]))*1.0/len(readToTruePairs[read])  #  recall is how many pairs were found among all to find
						#~ localPrecision = len(readToResultPairs[read].intersection(readToTruePairs[read]))*1.0/len(readToResultPairs[read])  # precision if how much predicted pairs are true pairs
					#~ else:
						#~ localRecall = 0
			#~ else:  # when a read is truly alone
				#~ localRecall = 1
				
				#~ if (len(readToResultPairs[read]) != 0):
					#~ localPrecision = 0
			#~ recall.append(localRecall)
			#~ precision.append(localPrecision)
		#~ else:  # when no read was linked to this read in the result
			#~ recall.append(localRecall)

	#~ print("recall ", sum(recall)*1.0/float(len(recall)), " precision ", sum(precision)*1.0/float(len(precision)))
	print (" a11 " , a11, " a01 " , a01, " a10 ", a10, " a00 ", a00)
	print("Jaccard index ", a11*1.0/(a11+a01+a10))
	print("Rand index ", (a11+a00)*1.0/(a11+a01+a10+a00))
else:
	print("validate_with_pairs clusters.truth clusters.result")
		
		
