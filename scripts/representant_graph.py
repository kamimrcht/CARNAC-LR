#!/usr/bin/env python
import sys
import pydot


def convert_file_to_clusters(file_name, clusters):  # useful to parse both final_g_clusters file and "true" clusters
		clusterfile = open(file_name, "r")
		i = 0
		for line in clusterfile.readlines():
			line = line.rstrip()
			query_read_ids = line.split(' ')
			for read in query_read_ids:
				if i in clusters:
					clusters[i].add(int(read))
				else:
					clusters[i] = set()
					clusters[i].add(int(read))
			i += 1
		clusterfile.close()

def compute_R_vectors(genes, clusters, allClusterVectors):  # presence/absence of reads from each gene in each cluster
		seen = False
		unclassified = 0
		sumR_ij = 0
		for gene_j in genes.keys():
			#~ if len(genes[gene_j]) > 10 :
				for read in genes[gene_j]:
					for cluster_i in clusters:
						if read in clusters[cluster_i]:
							seen = True
							sumR_ij += 1
							if cluster_i in allClusterVectors:
								allClusterVectors[cluster_i][gene_j] += 1
							else:
								allClusterVectors[cluster_i] = [0] * len(genes)
								allClusterVectors[cluster_i][gene_j] += 1
					if not seen:
						unclassified += 1
					else:
						seen = False
		return (sumR_ij, unclassified)

def compute_max_i_Rij(allClusterVectors, max_i):
		geneToReprCluster = dict()
		for j in range(len(max_i)):  # for each gene, search the best result cluster that represent it
			for result_clust in allClusterVectors.keys():
				if allClusterVectors[result_clust][j] > max_i[j]:
					max_i[j] = allClusterVectors[result_clust][j]
					geneToReprCluster[j] = result_clust
			if max_i[j] == 0:
				geneToReprCluster[j] = -1
		return geneToReprCluster

def compute_max_j_Rij(allClusterVectors, max_j, clusters):
		clusterToReprGene = dict()
		for result_clust in allClusterVectors.keys():
			if len(clusters[result_clust]) > 1 :
				maxj = 0
				for gene in range(len(allClusterVectors[result_clust])):
					if allClusterVectors[result_clust][gene] > maxj:
						maxj = allClusterVectors[result_clust][gene]
						clusterToReprGene[result_clust] = gene
		return clusterToReprGene
		



genes = dict()
clusters = dict()
allClusterVectors = dict()

if (len(sys.argv) > 2):
		genesFile = sys.argv[1]  # true clusters
		clustersFile = sys.argv[2]  # result clusters
		convert_file_to_clusters(genesFile, genes)  # true clusters are in dict genes
		convert_file_to_clusters(clustersFile, clusters)  # result clusters are in dict clusters
		max_i = [0] * len(genes)
		max_j = [0] * len(clusters)
		sumR_ij, unclassified = compute_R_vectors(genes, clusters, allClusterVectors)  #  for each result cluster, allClusterVector stores for each gene the number of reads in the result from this gene
		geneToReprCluster = compute_max_i_Rij(allClusterVectors, max_i)
		clusterToReprGene = compute_max_j_Rij(allClusterVectors, max_j, clusters)


		val = 20

		graph = pydot.Dot(graph_type='digraph', layout='circo')
		for gene in geneToReprCluster.keys():
			if len(genes[gene]) > val:
				node_g = pydot.Node("v"+str(gene), style="filled", fillcolor="red",  xlabel = '"%s"' %(len(genes[gene])))
				#~ node_g = pydot.Node("v"+str(gene), style="filled", fillcolor="red")
				graph.add_node(node_g)
		
		#~ for cluster in clusterToReprGene.keys():
			#~ node_r = pydot.Node("r"+str(cluster), style="filled", fillcolor="green", xlabel = '"%s"' %(len(clusters[cluster])))
			#~ graph.add_node(node_r)

		lookedAt = set()
		resultVisited = set()
		
		for gene in geneToReprCluster.keys():
			if len(genes[gene]) > val:
				if not geneToReprCluster[gene] == - 1:
					if geneToReprCluster[gene] in clusterToReprGene.keys():
						node_r = pydot.Node("r"+str(geneToReprCluster[gene]), style="filled", fillcolor="green", xlabel = '"%s"' %(len(clusters[geneToReprCluster[gene]])))
						resultVisited.add(geneToReprCluster[gene])
						graph.add_node(node_r)
						graph.add_edge(pydot.Edge("v"+str(gene), "r"+str(geneToReprCluster[gene])))
				lookedAt.add(gene)
			
		for cluster in clusterToReprGene.keys():
			if (clusterToReprGene[cluster] in lookedAt):
				if not cluster in resultVisited:
					 node_r = pydot.Node("r"+str(cluster), style="filled", fillcolor="green", xlabel = '"%s"' %(len(clusters[cluster])))
					 graph.add_node(node_r)
				graph.add_edge(pydot.Edge("r"+str(cluster), "v"+str(clusterToReprGene[cluster])))
		graph.write_png('graphe_representants.png')

else:
	print("Outputs metrics for clusters")
	print("Usage: python3 validate_clustering.py true_clusters.txt final_g_clusters.txt")
    
    

