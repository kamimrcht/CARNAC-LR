import sys

def convert_file_to_clusters(file_name, clusters):  # useful to parse both final_g_clusters file and "true" clusters
	clusterfile = open(file_name, "r")
	i = 0
	for line in clusterfile.readlines():
		if ">" in line:  #>1_referenceNumber:89 alternativeNumber:0 lowExpression
			line = line.rstrip()
			read = int(line.split('_')[0][1:]) - 1
			clust = int(line.split(':')[1].split(' ')[0])
			if clust in clusters:
				clusters[clust].append(str(read))
			else:
				clusters[clust]=[str(read)]
		#~ line = line.rstrip()
		#~ query_read_ids = line.split(' ')
		#~ toWrite = ""
		#~ for read in query_read_ids:
			#~ toWrite += str(i) + " "
			#~ i += 1
		#~ print(toWrite)
	clusterfile.close()

def print_clusters(clusters):
	for clust in clusters:
		tow = ""
		for r in clusters[clust]:
			tow += r + " "
		print(tow)

if len(sys.argv) > 1:
	clustersFile = sys.argv[1]
	clusters = dict()
	convert_file_to_clusters(clustersFile, clusters)
	print_clusters(clusters)
else:
	print("Converts clusters with reads indexes in strings to integer indexes")
	print("Usage: python3 converter_index.py clusters.txt")
    
