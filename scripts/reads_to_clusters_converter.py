import sys

def convert_file_to_clusters(file_name, indexToClusters):  # clusters file format must be "ch100_read3288", space separators
	clusterfile = open(file_name, "r")
	i = 0
	for line in clusterfile.readlines():
		line = line.rstrip()
		query_read_ids = line.split(' ')
		toWrite = ""
		for read in query_read_ids:
			indexToClusters[read] = i
		i += 1
	clusterfile.close()
	return i
	

def convert_reads_to_clusters(file_name, indexToClusters, clusters):  # >ch100_read3288_template_fail_BYK_CB_ONT_1_FAF05104_A, only FASTA
	readsfile = open(file_name, "r")
	indexRead = 0
	for line in readsfile.readlines():
		if ">" in line:
			line = line.rstrip()
			query_read_ids = '_'.join(line.split('_')[:2])[1:]  # ch100_read3288
			#~ print(query_read_ids)
			cluster = indexToClusters[query_read_ids]
			#~ print(cluster, indexRead)
			clusters[cluster].append(indexRead)
			indexRead += 1


def write_true_clusters(clusters):
	for clust in clusters:
		toWrite = ""
		#~ print("********", clust)
		for r in clusters[clust]:
			#~ print(r)
			toWrite += str(r) + " "
		print(toWrite)

if (len(sys.argv) > 2):
	readsFile = sys.argv[1]
	clustersFile = sys.argv[2]
	indexToClusters = dict()
	nbClusters = convert_file_to_clusters(clustersFile, indexToClusters)
	clusters = {k:[] for k in range(nbClusters)}
	convert_reads_to_clusters(readsFile, indexToClusters, clusters)
	#~ print(len(clusters[0]))
	#~ print(len(clusters[1]))
	#~ print(len(clusters[2]))
	write_true_clusters(clusters)
	
