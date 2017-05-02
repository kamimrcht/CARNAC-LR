
# -*- coding: utf-8 -*-

import argparse

def collecteLaVerite(clusterDeVerite):
	nbPairesVerite = dict()
	reads=dict()	
	with open(clusterDeVerite, "r") as truth:
		i=0		
		for line in truth:
			i+=1
			if len(line.split())>1:
				nbPairesVerite[i] = len(line.split())*(len(line.split())-1)/2
				for e in line.split():
					reads[e] = str(i)	
	return(reads, nbPairesVerite)


def testTool(clusterVerite, sortieOutil,ext):
	reads = clusterVerite[0]
	nbPairesVerite = clusterVerite[1]
	nbPairesOutil = dict()
	PairesDejaVues = set()
	TP = 0
	FP = 0

	def counting():
		nonlocal TP, FP
		#Tester si reads deja vus, on ordonne les reads pour etre sur de se comparer au meme tuple (A,B)!=(B,A)	
		if (min(readA, readB), max(readA, readB)) not in PairesDejaVues:
			PairesDejaVues.add((min(readA, readB), max(readA, readB)))
			clusA = reads[readA]
			clusB = reads[readB]

			if clusA == clusB:
				TP+=1
				if clusA in nbPairesOutil:
					nbPairesOutil[clusA]+=1
				else:
					#Initialization at 1 if the cluster has never been seen before
					nbPairesOutil[clusA]=1 
			else:	
				FP+=1
					
	with open(sortieOutil, "r") as tool:
		if ext == "minimap":
			for line in tool:
				ls = line.split("\t")
				readA = ls[0]
				readB = ls[5]
				counting()	
	
		elif ext == "graphmap" or ext == "mhap":
			for line in tool:
				readA = line.split()[0]
				readB = line.split()[1]
				counting()

		elif ext == "src":
			tool.readline()
			tool.readline()
			tool.readline()
			for line in tool:
				readA = line.split(":")[0]
				for r in line.split(":")[1].split():
					readB = r.split("-")[0]
					if readA != readB:
						counting()
	
	
	#Precision pour les reads
	if TP + FP != 0:
		precision = TP/(TP+FP)
	else:
		precision = 0
	#Recall par cluster pour les paires de reads
	if sum(nbPairesVerite.values()) != 0:
		recall = sum((nbPairesOutil[clus_id]) for clus_id in nbPairesOutil)/sum(nbPairesVerite.values())
	else:
		recall = 0
	#F-measure
	if precision != 0 or recall != 0:
		f = 2 * precision * recall / (precision + recall)
	else:
		f = 0
	return precision, recall, f



def writing(fichierSortie, metriques, tool):
	o=open(fichierSortie, "w", encoding='ascii')
	o.write("#Tested tool\t%s\n" %tool)	
	o.write('#Precision\t%s\n#Recall\t%s\n#F-measure\t%s\n' %metriques)
	o.close()
			
def main():
	parser=argparse.ArgumentParser(description="Computing recall and precision")
	parser.add_argument('-i', help='input file')
	parser.add_argument('-v', help='truth input file')
	parser.add_argument('-t', help='tested tool', choices=["src", "minimap", "mhap", "graphmap"])
	parser.add_argument('-o', help='output file') 
	args = parser.parse_args()
	
	verite = collecteLaVerite(args.v)
	measures = testTool(verite, args.i, args.t)	
	writing(args.o, measures, args.t)	

main()
