import sys

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
		for j in range(len(max_i)):
			for vector in allClusterVectors.values():
				if vector[j] > max_i[j]:
					max_i[j] = vector[j]  # cluster that have the more nodes of a gene

def compute_max_j_Rij(allClusterVectors, max_j):
		for i in range(len(max_j)):
			max_j[i] = max(allClusterVectors[i])


def compute_precision(max_j, sumR_ij):
		if sumR_ij > 0:
			return sum(max_j)*1.0/sumR_ij
		else:
			return 0

def compute_sensitivity(max_i, sumR_ij, unclassified):
		if sumR_ij > 0:
			return sum(max_i)*1.0/(sumR_ij + unclassified)
		else:
			return 0
    



genes = dict()
clusters = dict()
allClusterVectors = dict()

if (len(sys.argv) > 2):
		genesFile = sys.argv[1]
		clustersFile = sys.argv[2]
		convert_file_to_clusters(genesFile, genes)
		convert_file_to_clusters(clustersFile, clusters)
		#~ print(len(genes), len(clusters))
		max_i = [0] * len(genes)
		max_j = [0] * len(clusters)
		sumR_ij, unclassified = compute_R_vectors(genes, clusters, allClusterVectors)
		#~ print(len(allClusterVectors), sumR_ij, unclassified)
		compute_max_i_Rij(allClusterVectors, max_i)
		compute_max_j_Rij(allClusterVectors, max_j)
		#~ print(max_i, max_j)
		precision = compute_precision(max_j, sumR_ij)
		print("precision =", precision)
		sensitivity = compute_sensitivity(max_i, sumR_ij, unclassified)
		print("sensitivity =", sensitivity)
		print("F-measure =", 2 * precision * sensitivity / (precision + sensitivity))
else:
	print("Outputs metrics for clusters")
	print("Usage: python3 validate_clustering.py true_clusters.txt final_g_clusters.txt")
    
    

