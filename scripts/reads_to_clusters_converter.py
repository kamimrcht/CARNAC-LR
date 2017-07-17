import sys

def convert_file_to_clusters(file_name, readToCluster):  # clusters file format must be "ch100_read3288 clustxx", space separators
	clusterfile = open(file_name, "r")
	i = 0
	clusters = dict()
	for line in clusterfile.readlines():
		line = line.rstrip()
		if not line : break
		if len(line.split(' ')) > 1:  # some reads have no cluster i don't know why
			readName = line.split(' ')[0]
			cluster = line.split(' ')[1]

			read = '_'.join(line.split('_')[:2])[0:]
			toWrite = ""
			if not cluster in clusters.keys():
				clusters[cluster] = i
				i += 1
				
				
			#~ print(read,cluster)
			readToCluster[read] = clusters[cluster]
	clusterfile.close()
	return i
	

def convert_reads_to_clusters(file_name, readToCluster, clusters):  # >ch100_read3288_template_fail_BYK_CB_ONT_1_FAF05104_A, only FASTA
	readsfile = open(file_name, "r")
	indexRead = 0
	readsList = []
	for line in readsfile.readlines():
		if ">" in line:
			line = line.rstrip()
			read = '_'.join(line.split('_')[:2])[1:]  # ch100_read3288
			if read in readToCluster.keys():
				#~ print("*", readToCluster[read], len(clusters))
				clusters[int(readToCluster[read])].append(indexRead)
			indexRead += 1
	


def write_true_clusters(clusters):
	for clust in clusters:
		toWrite = ""
		if len(clusters[clust]) != 0:
		#~ print("********", clust)
			for r in clusters[clust]:
				#~ print(r)
				toWrite += str(r) + " "
			print(toWrite)

if (len(sys.argv) > 2):
	readsFile = sys.argv[1]
	clustersFile = sys.argv[2]
	readToCluster = dict()
	nbClusters = convert_file_to_clusters(clustersFile,readToCluster)
	clusters = {k:[] for k in range(nbClusters)}
	convert_reads_to_clusters(readsFile, readToCluster, clusters)
	#~ print(len(clusters[0]))
	#~ print(len(clusters[1]))
	#~ print(len(clusters[2]))
	write_true_clusters(clusters)
	
