#!/usr/bin/env python
import sys
import os
import pydot


if (len(sys.argv) > 5):
	# compute egdes (non redundant)
	with open(sys.argv[1]) as infile:  # src output
		nodes = dict()
		for rline in infile:
			if not '#' in rline:
				read = set()
				line = rline.rstrip().split(":")
				line2 = line[1].split(" ")
				for x in line2:
					neigh = x.split("-")[0]
					if not neigh in nodes:
						if neigh != line[0]:
							read.add(neigh)
							nodes[neigh] = [set(), '']
					else:
						if neigh != line[0]:
							if not line[0] in nodes[neigh]:
								nodes[neigh][0].add(line[0])
				if not line[0] in nodes:
					nodes[line[0]] = [read, '']


	readToIndex = dict()
	index = 0
	with open(sys.argv[2]) as infile:  # fasta file
		for rline in infile:
			if '>' in rline:
				read = rline.rstrip().split("t")[0][1:-1]
				readToIndex[read] = str(index)
				index += 1
				
	colors = {0:'circle', 1:'box'}
	indexColor = 0
	with open(sys.argv[3]) as infile:  # one cluster per line
		for rline in infile:
			cluster = rline.rstrip().split(" ")
			for c in cluster:
				index = readToIndex[c]
				if index in nodes:
					nodes[index][1] = colors[indexColor]  # one color for each cluster
			indexColor += 1


	graph = pydot.Dot(graph_type='graph')
	#~ edge = pydot.Edge('1' , '2')
	#~ graph.add_edge(edge)

	metrics = dict()
	with open(sys.argv[4]) as infile:  # node metrics in format: node cc deg, one node per line
		for rline in infile:
			metricsNode = rline.rstrip().split(" ")
			metrics[metricsNode[0]] = [metricsNode[1], metricsNode[2]]  # metrics[node] = [cc, deg]
			
			
	listEdges = set()
	for n, attr in nodes.items():
		label =  ','.join(metrics[n])
		print(label)
		node = pydot.Node(n, xlabel = '"%s,%s"' %(metrics[n][0][0:4], metrics[n][1]))
		#~ node = pydot.Node(n, xlabel = "0.9,5")
		#~ node = pydot.Node(n + " " + " ".join(metrics[n]))
		node.set("shape", attr[1])
		graph.add_node(node)
		for neighbor in attr[0]:
			node = pydot.Node(neighbor, xlabel = '"%s,%s"' %(metrics[neighbor][0][0:4], metrics[neighbor][1]))
			node.set("shape", nodes[neighbor][1])
			#~ graph.add_node(node)
			#~ listEdges.add((min(n, neighbor), max(n, neighbor)))
			listEdges.add((min(n, neighbor), max(n, neighbor)))
			#~ listEdges.add((n, neighbor))





	#~ graph.write_png('example1_graph.png')
	notRemovedEdges = set()
	with open(sys.argv[5]) as infile:  # cliqueness clustering output
		for rline in infile:
			line = rline.rstrip().split(" ")  # one cluster
			for node1 in line:
				for node2 in line:
					if node1 != node2:
						#~ edge = (node1, node2)
						edge = (min(node1, node2), max(node1, node2))
						if edge in listEdges:
							notRemovedEdges.add(edge)
	for edges in listEdges:
		if not edges in notRemovedEdges:
			edge = pydot.Edge(edges, color = 'red')
			graph.add_edge(edge)
			#~ pass
		else:
			edge = pydot.Edge(edges)
			graph.add_edge(edge)
	print("Printing graph")			
	graph.write_png('example1_graph.png')
else:
	print("Usage: ./print_graphs_cliqueness.py src_outuput_file reads_fasta_file true_clusters_file nodes_metrics_output_by_clustering clusters_ouput_by_clustering")
