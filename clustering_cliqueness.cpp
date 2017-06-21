#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <map>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sstream>
#include <mutex>
#include <functional>
#include <utility>
#include "clustering_cliqueness.hpp"



using namespace std;

struct Node_hash {
    inline std::size_t operator()(const Node& n) const {
        return pow(n.index, 31) + n.degree;
    }
};


inline bool operator == (Node const& n1, Node const& n2)
{
    return (n1.index == n2.index); 
}


vector<uint> removeDuplicates(vector<uint>& vec){
	sort(vec.begin(), vec.end());
	vec.erase(unique(vec.begin(), vec.end()), vec.end());
	return vec;
}

vector<double> removeDuplicatesCC(vector<double>& vec){
	sort(vec.begin(), vec.end());
	vec.erase(unique(vec.begin(), vec.end()), vec.end());
	return vec;
}


// compare nodes first by degree, then by CC
//~ bool CompareNodes(const Node& a, const Node& b) 
//~ { 
   //~ if (a.degree < b.degree) return true;
   //~ if (b.degree < a.degree) return false;

   //~ if (a.CC < b.CC) return true;
   //~ if (b.CC < a.CC) return false;
   //~ return false;
//~ } 


void DFS(uint n, vector<Node>& vecNodes, unordered_set<uint>& visited, set<uint>& nodesInConnexComp, bool& above, double cutoff){
	if (not visited.count(n)){
		if (vecNodes[n].CC >= cutoff){
			above = true;
		}
		visited.insert(n);
		nodesInConnexComp.insert(n);
		for (auto&& neigh : vecNodes[n].neighbors){
			DFS(neigh, vecNodes, visited, nodesInConnexComp, above, cutoff);
		}
	}
}

vector<string> split(const string &s, char delim){
	stringstream ss(s);
	string item;
	vector<string> elems;
	while (getline(ss, item, delim)) {
		elems.push_back(move(item)); 
	}
	return elems;
}


//  parser for infile at SRC format.
// fills a vector of Nodes
// Nodes are reads, they remember their neighbors in the graphs
// Edges exist when a similarity is reported in the infile
// Weights info also come from infile
// in the case infile reports:
// read1: read2 read3
// read2: read3
// edges will be (1,2) (1,3) (2,3) (non oriented) which means a similarity becomes an edge even if it is reported in only one direction (read 1 to 2 for instance)
// weights:

void parsingSRC(ifstream & refFile, vector<Node>& vecNodes, bool weighted){
	string listNodes;
	// header
	getline(refFile, listNodes);
	getline(refFile, listNodes);
	getline(refFile, listNodes);
	vector<uint>  neighbs;
	vector<vector<uint>> clust;
	vector<string> splitted1, splitted2, splitted3;
	uint read, target;
	unordered_map <uint, uint> seenNodes;
	unordered_map <uint, double> neighbToWeight;
	double weight(1);
	while (not refFile.eof()){
		getline(refFile, listNodes);
		splitted1 = split(listNodes, ':');
		if (splitted1.size() > 1){
			splitted2 = split(splitted1[1], ' ');
			target = stoi(splitted1[0]);  // target read's index
			if (not splitted2.empty()){
				for (uint i(0); i < splitted2.size(); ++i){
					splitted3 = split(splitted2[i], '-');
					read = stoi(splitted3[0]);  // recruited read
					if (weighted){
						weight = 1/ (stof(splitted3[1]));
					}
					if (read != target){
						if (not seenNodes.count(target)){ // new node not already in the vector of nodes
							clust = {}; neighbs = {}; neighbToWeight = {};
							Node t({target, 0, 0, clust, neighbs, neighbToWeight});
							vecNodes.push_back({t});  // store in vecNodes
							vecNodes.back().neighbToWeight.insert({read, weight});
							seenNodes.insert({target, vecNodes.size() - 1}); // remember this node has been pushed index -> place in the vector
						}
						if (seenNodes.count(read)){ // this neighbour is already in the vector of nodes
							vecNodes[seenNodes[target]].neighbors.push_back(seenNodes[read]);  // add read as neighbor of target
							if (not vecNodes[seenNodes[target]].neighbToWeight.count(read)){
								vecNodes[seenNodes[target]].neighbToWeight.insert({read, weight});
							}
							vecNodes[seenNodes[target]].neighbors.push_back(seenNodes[read]);  // add read as neighbor of target
							vecNodes[seenNodes[read]].neighbors.push_back(seenNodes[target]);  // add target as neighbor of read
							if (not vecNodes[seenNodes[read]].neighbToWeight.count(target)){
								vecNodes[seenNodes[read]].neighbToWeight.insert({target, weight});
							}
						} else {  // new neighbor not already in the vector of nodes
							clust = {}; neighbs = {}; neighbToWeight = {};
							Node r({read, 0, 0, clust, neighbs, neighbToWeight});
							vecNodes.push_back({r});
							uint position(vecNodes.size() - 1);
							seenNodes.insert({read, position});
							vecNodes[seenNodes[target]].neighbors.push_back(vecNodes.size() - 1);  // store read as neighbor of target
							vecNodes[seenNodes[read]].neighbors.push_back(seenNodes[target]);  // add target as neighbor of read
							vecNodes[seenNodes[read]].neighbToWeight.insert({target, weight}); 
							vecNodes[seenNodes[target]].neighbToWeight.insert({read, weight});  
						}
					}
				}
			}
		}   
	}
}


//compute the CC of a node with its set of neighbors
double getCC(unordered_set<uint>& neighbors, vector<Node>& vecNodes){
	double pairs(0), clusteringCoef(0);
	uint totalPairs;
	for (auto&& neigh : neighbors){  // for each neighbor of the node
		for (auto&& neigh2 : vecNodes[neigh].neighbors){ // for each neighbor of a neighbor
			if (neighbors.count(neigh2)){  // if a neighbor of a neighbor is also a neighbor of the current node = pair of connected neighbors
				++pairs;
			}
		}
	}
	totalPairs = neighbors.size() * (neighbors.size() - 1);
	if (totalPairs > 0){
		clusteringCoef = pairs/totalPairs;
	}
	return clusteringCoef;
}


// the CC values are compared before and after the removal of a sub set of nodes
// if the value is > 0: the removal impacted negatively the cliqueness
// if the value == 0:  no change
// if the value < 0: the cliqueness is better without the nodes

// todo: change = we should take the cluster's CC instead of the cutoff
int getDeltaCC(set<uint>& toRemove, set<uint>& clust1, vector<Node>& vecNodes, double cutoff){
	int deltaCC(0);
	unordered_set<uint> clust1Without;
	for (auto&& i : clust1){
		if (not toRemove.count(i)){
			clust1Without.insert(i);
		}
	}
	double CC1(getCC(clust1Without, vecNodes));
	deltaCC = cutoff - CC1;
	return deltaCC;
}


// compute clustering coefficient and degree for each node
// CC will be stored in ClCo and degrees in vector degrees, and sorted from higher to lower (with no duplicates)
// the CC are stored by decreasing values in the vector as we will start with nodes of highest CC later in the code
// the quantile values are computed before removing the duplicates of CC
// the graph is not oriented so the degree of a node is the nb of edges going through this node
void computeCCandDeg(vector<Node>& vecNodes, vector<double>& ClCo, vector<uint>& degrees, float& lowerCC){
	double clusteringCoef;
	unordered_set<uint> neighbors;
	// start by removing double occurrences in neighbors (a neighbor of a node that would be stored twice)
	for (uint n(0); n < vecNodes.size(); ++n){
		vecNodes[n].neighbors = removeDuplicates(vecNodes[n].neighbors); 
		vecNodes[n].degree = vecNodes[n].neighbors.size();
		degrees.push_back(vecNodes[n].degree);
	}
	
	// fill vector of CC
	for (uint n(0); n < vecNodes.size(); ++n){
		if (vecNodes[n].neighbors.size() > 1){
			neighbors = {};
			copy(vecNodes[n].neighbors.begin(), vecNodes[n].neighbors.end(), inserter(neighbors, neighbors.end()));
			clusteringCoef = getCC(neighbors, vecNodes);
			if (clusteringCoef != 0){
				 vecNodes[n].CC = clusteringCoef;
				ClCo.push_back(clusteringCoef);
			}
		} else {
			ClCo.push_back(0);
		}
	}

	// compute the quantiles of the CC and get the CC value that represents the first quantile
	lowerCC = quantileCC(ClCo, 1, 1000);
	
	// remove duplicated values of CC (so that we don't compute twice a cutoff value later)
	// and sort CC in vector by decreasing values
	ClCo = removeDuplicatesCC(ClCo);
	sort(ClCo.begin(), ClCo.end());
	reverse(ClCo.begin(), ClCo.end());
}

// sort nodes by decreasing degrees
// update vecNodes (neighbours...) as index of nodes may have changed
void sortVecNodes(vector<Node>& vecNodes, vector<uint>& nodesInOrderOfCC){
	vector<Node> vecNodesCpy(vecNodes);
	sort(vecNodesCpy.begin(), vecNodesCpy.end());  // based on degrees then CC
	reverse(vecNodesCpy.begin(), vecNodesCpy.end());
	unordered_map <uint, uint> indexReadsAf;
	for (uint i(0); i < vecNodesCpy.size(); ++i){
		nodesInOrderOfCC.push_back(i);
	}
}



// compute sets around seed nodes
// seed nodes are those whose CC is above the current cutoff
// they gather their direct neighbors in sets that are expected quasi cliques
// nodes which degree is above the last quantile and CC is below the last quantile (even if it is above the current cutoff) are not seeds
// the function is multi threaded with one thread by cutoff
// each Nodes stores a vector of vector (.cluster) )of the size of the length of the list of cutoffs
// in each vector of this vector, if the node is a seed, a number of set is stored (the same number is stored for nodes included in the set)
// this way, at a given cutoff, a node can know if it is in different sets


//todo : change this
// finally we sort the nodes in decreasing order of CC
void computePseudoCliques(vector<double>& cutoffs, vector<Node>& vecNodes, uint nbThreads, vector<uint>& nodesInOrderOfCC, uint higherDegree, float lowerCC){
	vector<uint> v;
	vector<vector<uint>> vec(cutoffs.size());
	// for each node, at each cutoff value we will store a vector that sums up the number of sets the node belongs to
	for (uint i(0); i < vecNodes.size(); ++i){
		vecNodes[i].cluster = vec;
	}
	uint c(0);
	vector<unordered_set<uint>> temp(cutoffs.size());  // sets identifiers for each cutoff value
	#pragma omp parallel num_threads(nbThreads)
	{
		#pragma omp for
		for (c = 0; c < cutoffs.size(); ++c){  // descending cutoffs
			unordered_set<uint> s;
			double cutoff = cutoffs[c];
			for (uint i(0); i < vecNodes.size(); ++i){
				if (vecNodes[i].CC >= cutoff and not (vecNodes[i].degree >= higherDegree and vecNodes[i].CC <= lowerCC)){  // if the node is a seed
					vecNodes[i].cluster[c].push_back(i);  // store a set identifier for this node
					s.insert(i);
					for (auto&& neigh : vecNodes[i].neighbors){
						vecNodes[neigh].cluster[c].push_back(i);  // also include direct neighbors in set
						s.insert(neigh);
					}
				}
			}
			temp[c] = s;
		}
	}
	//~ unordered_set<uint> s;
	//~ for (uint i(0); i < temp.size(); ++i){  // for each cutoff (first higher values)
		//~ for (auto&& n : temp[i]){  // for each set identifier
			//~ if (not s.count(n)){
				//~ s.insert(n);
				//~ nodesInOrderOfCC.push_back(n);  // nodes are ordered from belonging to a cluster of higher CC to lower
			//~ }
		//~ }
	//~ }
}



// generalized local CC computation for a set of nodes (instead of the direct neighboring of a given node)
double computeUnionCC(set<uint>& unionC, vector<Node>& vecNodes){
	double cardUnion(0);
	for (auto&& n : unionC){
		for (auto&& neigh : vecNodes[n].neighbors){
			if (unionC.count(neigh)){
				++cardUnion;
			}
		}
	}
	return cardUnion / ( unionC.size() * (unionC.size() - 1));
}



// transfer nodes from a set to another cancelled set
void transfer(uint tf, uint te, set<uint>& toFill, set<uint>& toEmpty, vector<Node>& vecNodes, vector<set<uint>>& clusters, uint ind){
	vector<uint> vec;
	// remove the  cancelled set identifier from the nodes of the set to cancel: we have to compute again the new list of sets without the cancelled one
	for (auto&& index : toEmpty){
		vec = {};
		for (auto && clust : vecNodes[index].cluster[ind]){
			if (not (clust == te)){
				vec.push_back(clust);
			}
		}
		vec.push_back(tf);
		vecNodes[index].cluster[ind] = removeDuplicates(vec);
	}
	// same operation for the nodes from the nodes in the set to keep
	for (auto&& index : toFill){
		vec = {};
		for (auto && clust : vecNodes[index].cluster[ind]){
			if (not (clust == te)){
				vec.push_back(clust);
			}
		}
		vecNodes[index].cluster[ind] = removeDuplicates(vec);
	}
}


// merge two sets: the smaller is added to the bigger
void merge(uint i1, uint i2, set<uint>& clust1, set<uint>& clust2, vector<set<uint>>& clusters,  vector<Node>& vecNodes, uint ind){
	if (clust1.size() > clust2.size()){  // merge in clust1
		clusters[i1].insert(clusters[i2].begin(), clusters[i2].end());
		transfer(i1, i2, clust1, clust2, vecNodes, clusters, ind);
		clusters[i2] = {};
	} else {  // merge in clust2
		clusters[i2].insert(clusters[i1].begin(), clusters[i1].end());
		transfer(i2, i1, clust2, clust1, vecNodes, clusters, ind);
		clusters[i1] = {};
	}
}



vector<set<uint>> assignNewClusters(set<uint>& clust, vector<Node>& vecNodes, double cutoff){
	bool above(false);
	unordered_set<uint> visited;
	set<uint> nodesInConnexComp;
	vector<set<uint>> newClust;
	for (auto&& node : clust){
		if (not visited.count(node)){
			above = false;
			nodesInConnexComp = {};
			DFS(node, vecNodes, visited, nodesInConnexComp, above, cutoff);
			if (above){
				newClust.push_back(nodesInConnexComp);
			}
		}
	}
	return newClust;
}


// when a set is split, update the nodes' clusters index values (remove the split set index)
void removeSplittedElements(uint index, vector<set<uint>>& clusters, set<uint>& interC){
	set<uint> clust;
	for (auto && elt : clusters[index]){
		if (not interC.count(elt)){
			clust.insert(elt);
		}
	}
	clusters[index] = clust;
}



//compute each cut for each set: the number of edges that connects nodes of the set that are not in the intersection to nodes that are in the intersection
// weights are used to compute the cuts
void getCutsPairSets(vector<Node>& vecNodes, set<uint>& interC, double& cut1, double& cut2, set<uint>& clust1, set<uint>& clust2){
	for (auto&& node : interC){
		for (auto&& neigh : vecNodes[node].neighbors){
			if (clust1.count(neigh) and (not interC.count(neigh))){
				cut1 += vecNodes[node].neighbToWeight[vecNodes[neigh].index];
			}
			if (clust2.count(neigh) and (not interC.count(neigh))){
				cut2 += vecNodes[node].neighbToWeight[vecNodes[neigh].index];
			}
		}
	}
}


// in case of split and ex aequo of the two cuts
// we compute the delta CC (with and without the nodes of the intersection) for each sets, and nodes are let in the set with the smaller deltaCC
void splitExAequo(uint i1, uint i2, set<uint>& clust1, set<uint>& clust2, vector<set<uint>>& clusters,  vector<Node>& vecNodes, set<uint>& interC, uint cutoff, uint ind, double& cut1, double& cut2, double& cut){
	int deltaCC1(getDeltaCC(interC, clust1, vecNodes, cutoff));
	int deltaCC2(getDeltaCC(interC, clust2, vecNodes, cutoff));
	if (deltaCC1 <= deltaCC2){  // keep the intersection in clust1
		transfer(i1, i2, clust1, interC, vecNodes, clusters, ind);
		removeSplittedElements(i2, clusters, interC);
		cut = cut2;
	} else {  // keep the intersection in clust2
		transfer(i2, i1, clust2, interC, vecNodes, clusters, ind);
		removeSplittedElements(i1, clusters, interC);
		cut = cut1;
	}
}



// split a set and keep the node of the intersection in the other set
// when a set is split, it can happen that it becomes disconnected
// we check if that happens, if so we build new clusters with the new  connected subgraphs
void splitProcedure(uint i1, uint i2, set<uint>& clust1, set<uint>& clust2, vector<set<uint>>& clusters,  vector<Node>& vecNodes, set<uint>& interC, uint cutoff, uint ind,  double& cut2, double& cut, vector<set<uint>>& newClust){
	bool more(findArticulPoint(clust2, vecNodes, interC));  // to avoid doing too much computation, if at least one articulation point is found we do a DFS afterwards
	transfer(i1, i2, clust1, interC, vecNodes, clusters, ind);
	removeSplittedElements(i2, clusters, interC);
	cut = cut2;   // todo * 2 ?
	if (more){  // if there was an articulation point in the intersection, it is likely that the split set will be disconnected (if not we are sure it stay connected)
		newClust = assignNewClusters(clust2, vecNodes, cutoff);  // do a DFS to find 1 or more connected components
		if (newClust.size() > 1){  // if there are several connected components do new clusters
			for (uint i(0); i < newClust.size(); ++i){
				clusters.push_back(newClust[i]);
				for (auto&& nodes: newClust[i]){
					vecNodes[nodes].cluster[ind].push_back(clusters.size() - 1);
				}
				transfer(clusters.size() - 1, i2, newClust[i], clust2, vecNodes, clusters, ind);
			}
			clusters[i2].clear();
		}
	}
}



// in case sets are not merged, one of the two sets of the pair is split gives its nodes belonging to the intersection to the other set
// the set that has the most edges connected to nodes of the intersection keeps them, in order to minimize the cut
double splitClust(uint i1, uint i2, set<uint>& clust1, set<uint>& clust2, vector<set<uint>>& clusters,  vector<Node>& vecNodes, set<uint>& interC, uint cutoff, uint ind){
	double cut1(0), cut2(0), cut(0);
	// compute each cut for each set in the pair
	getCutsPairSets(vecNodes, interC, cut1, cut2, clust1, clust2);

	if (clust1.size() == interC.size()){  // if set 1 is exactly the intersection: separate the intersection from set 2
		transfer(i1, i2, clust1, interC, vecNodes, clusters, ind);
		removeSplittedElements(i2, clusters, interC);
		cut = cut2;
	} else if (clust2.size() == interC.size()){ // if set 2 is exactly the intersection: separate the intersection from set 1
		transfer(i2, i1, clust2, interC, vecNodes, clusters, ind);
		removeSplittedElements(i1, clusters, interC);
		cut = cut1;
	} else {
		unordered_set <uint> neighbors;
		vector<set<uint>> newClust;
		if (cut1 == cut2){ 
			// in case of ex aequo choose using the delta CC
			splitExAequo(i1, i2, clust1, clust2, clusters, vecNodes, interC, cutoff, ind, cut1, cut2, cut);

		} else if (cut1 > cut2){  // split clust 2
			splitProcedure(i1, i2, clust1, clust2, clusters, vecNodes, interC, cutoff, ind, cut2, cut, newClust);
		} else {
			// split clust1
			splitProcedure(i2, i1, clust2, clust1, clusters, vecNodes, interC, cutoff, ind, cut1, cut, newClust);
		}
	}
	return cut;
}



// computation of the cut
// any edge that links a node in a cluster to another node which is not in the same cluster increases the cut
// any node (singleton) which is is no cluster increases the cut by its number of edge
double getCut(vector<Node>& vecNodes, vector<set<uint>>& clusters, uint ind){
	double cut = 0;
	vector <uint> seen(vecNodes.size(), 0);
	for (uint i(0); i < vecNodes.size(); ++i){
		if (vecNodes[i].cluster[ind].empty() and (not vecNodes[i].neighbors.empty())){
			for (auto&& ne : vecNodes[i].neighbors){
				cut += vecNodes[i].neighbToWeight[vecNodes[ne].index];
			}
		} else {
			if (not vecNodes[i].cluster[ind].empty()){
				for (auto&& ne : vecNodes[i].neighbors){
					if (not (clusters[vecNodes[i].cluster[ind][0]].count(ne))){
						cut += vecNodes[i].neighbToWeight[vecNodes[ne].index];
					}
				}
			} 
		}
	}
	return cut;
}



// performs merges or splits according to the compared values of the cutoff/generalized CC
// also get a temporary cut value
void mergeOrSplitProcedures(double cutoff, vector<Node>& vecNodes, vector<set<uint>>& clusters, uint ind, double prevCut, vector<uint>& nodesInOrderOfCC, set<uint>& clust1, set<uint>&  clust2, set<uint>&  unionC, set<uint>&  interC, uint& i1, uint& i2, double& cut, uint i){
	// comparison of clusters  by pairs: decisions are taken first for the sets associated with the higher CCs
	double unionCC;
	i1 = vecNodes[i].cluster[ind][0];  // sets are sorted by decreasing seed's CC value, so we compare the two first
	i2 = vecNodes[i].cluster[ind][1];  // that are associated with the two higher CC
	clust1 = clusters[i1];
	clust2 = clusters[i2];
	interC = {};
	set_intersection(clust1.begin(), clust1.end(), clust2.begin(), clust2.end(), inserter(interC, interC.begin()));  // intersection of the two sets
	if (interC.size() == clust1.size() and clust1.size() == clust2.size()){  // clust1 and clust2 are the same
		transfer(i1, i2, clust1, interC,  vecNodes, clusters, ind);  // keep only one and cancel the other
		clusters[i2] = {};
	} else {
		unionC = {};
		set_union(clust1.begin(), clust1.end(), clust2.begin(), clust2.end(), inserter(unionC, unionC.begin()));  // union of the two sets
		unionCC = computeUnionCC(unionC, vecNodes);  // get CC generalized to the union of nodes
		if (unionCC >= cutoff){  // merge operation
			merge(i1, i2, clust1, clust2, clusters, vecNodes, ind);
		} else {  // split operation
			cut += splitClust(i1, i2, clust1, clust2, clusters, vecNodes, interC, cutoff, ind);
		}
	}
}



// for a given cutoff, from the original set of sets of nodes creates from seeds, refine those sets (by merging/splitting strategies) to obtain clusters
// compute the cut associated to those operations
double computeClustersAndCut(double cutoff, vector<Node>& vecNodes, vector<set<uint>>& clusters, uint ind, double prevCut, vector<uint>& nodesInOrderOfCC){
	double cut(0), sCut(0);
	set<uint> clust1, clust2, unionC, interC;
	uint i1, i2;
	// compute list of current sets to be refined in vector clusters
	for (uint n(0); n < vecNodes.size(); ++n){
		if (vecNodes[n].cluster[ind].size() > 0){
			for (auto&& c : vecNodes[n].cluster[ind]){
				clusters[c].insert(n);  // node at index n is in cluster c
			}
		}
	}

	// for each node (by decreasing CC value)
	for (auto&& i : nodesInOrderOfCC){
		cut = 0;
		if (vecNodes[i].cluster[ind].size() > 1){  // if node is in several clusters: choices to make to have only one cluster in the end
			while (vecNodes[i].cluster[ind].size() > 1){
				mergeOrSplitProcedures(cutoff, vecNodes, clusters, ind, prevCut, nodesInOrderOfCC, clust1, clust2, unionC, interC, i1, i2, cut, i);
			}
		} else {
			sCut += cut;
			if (ind > 0 and sCut > prevCut){
				// if the cut is already higher than a precedant cut = we don't reach the minimal cut for this cutoff so we can stop already
				return sCut;
			}
		}
	} // when this procedure stops, each node remains in only one set which is a cluster
	// computation of the cut
	cut = getCut(vecNodes, clusters, ind);
	return cut;
}




// compute a vector of Nodes specific for a given connected component: re-write indexes
void getVecNodes(vector<Node>& vecNodes, vector<Node>& vecNodesGlobal, set<uint>& nodesInConnexComp){
	uint ii(0);
	unordered_map<uint, uint> indexReads;
	for (auto&& val: nodesInConnexComp){
		vecNodes.push_back(vecNodesGlobal[val]);
		indexReads.insert({val, ii});
		++ii;
	}
	vector<uint> vec;
	for (uint i(0); i < vecNodes.size(); ++i){
		vec = {};
		for (auto&& n : vecNodes[i].neighbors){
			vec.push_back(indexReads[n]);
		}
		vecNodes[i].neighbors = vec;
	}
}






bool findArticulPoint(set<uint>& cluster, vector<Node>& vecNodes, set<uint>& interC){
	Graph graph(vecNodes.size());
	unordered_set<uint> visited;
	for (auto&& i : cluster){
		visited.insert(i);
		for (auto&& neigh : vecNodes[i].neighbors){
			if ((not visited.count(neigh)) and cluster.count(neigh)){
				graph.addEdge((int)i, neigh);
			}
		}
	}
	vector<bool> ap; // To store articulation points
    bool b(graph.APBool(ap, interC));
    return b;
}


// disconnects articulations points from the total graph
// each time an articulation point is removed (disconnected), the nb of connected comp. in the graph increases
// it allows to fragment the graph in more connected components and isolate problematic nodes that are articulation points
// articulation points are searched using a DFS then disconnected from the graph
void preProcessGraph(vector<Node>& vecNodes, double cutoff=1.1){
	Graph graph(vecNodes.size());
	unordered_set<uint> visited;
	for (uint i(0); i < vecNodes.size(); ++i){
		visited.insert(i);
		for (auto&& neigh : vecNodes[i].neighbors){
			if (not visited.count(neigh)){
				graph.addEdge((int)i, neigh);
			}
		}
	}
	vector<uint> vec;
	vector<bool> ap; // To store articulation points
    graph.AP(ap); // get articulation points via DFS
    // disconnect nodes:
    for (uint i = 0; i < vecNodes.size(); i++){
        if (ap[i] == true and vecNodes[i].CC < cutoff){
            for (auto&& j : vecNodes[i].neighbors){
				vec = {};
				for (auto&& jj : vecNodes[j].neighbors){
					if (i != jj){
						vec.push_back(jj);
					}
				}
				vecNodes[j].neighbors = vec;
				vecNodes[j].degree = vecNodes[j].neighbors.size();
			}
			vecNodes[i].neighbors = {};
			vecNodes[i].degree = 0;
		}
	}
}



void preProcessGraphQuantiles(vector<Node>& vecNodes, double cutoffCC, uint cutoffEdges){
	vector<uint> vec;
	for (uint i = 0; i < vecNodes.size(); i++){
        if (vecNodes[i].degree >= cutoffEdges and vecNodes[i].CC <= cutoffCC){
            for (auto&& j : vecNodes[i].neighbors){
				vec = {};
				for (auto&& jj : vecNodes[j].neighbors){
					if (i != jj){
						vec.push_back(jj);
					}
				}
				vecNodes[j].neighbors = vec;
				vecNodes[j].degree = vecNodes[j].neighbors.size();
			}
			vecNodes[i].neighbors = {};
			vecNodes[i].degree = 0;
		}
	}
}



uint quantileEdges(vector<uint>&degrees, uint no, uint q){
	double e;
	e = degrees.size()*((double)no/q);
	return (uint)e;
}


double quantileCC(vector<double>&CC, uint no, uint q){
	double cc;
	cc = CC.size()*((float)no/q);
	return cc;
}






void parseArgs(int argc, char** argv, bool& approx, bool& preprocessing, bool& weighted, string& fileName, string& outFileName, uint& nbThreads, uint& granularity){
	approx = false; preprocessing = false; weighted = false;
	outFileName = "final_g_clusters.txt"; fileName = "";
	nbThreads = 2;
	int c;
	granularity = 10;
	while ((c = getopt (argc, argv, "f:o:c:i:pw")) != -1){
		switch(c){
			case 'o':
				outFileName=optarg;
				break;
			case 'f':
				fileName=optarg;
				break;
			case 'c':
				nbThreads=stoi(optarg);
				break;
			case 'i':
				approx = true;
				granularity=stoi(optarg);
				break;
			case 'p':
				preprocessing = true;
				break;
			case 'w':
				weighted = true;
				break;
		}
	}
}


void printHelpCmd(bool help){
	if (help){
		cout << "Usage : ./clustering_cliqueness -f input_file (-o output_file -w -i 10 -p -c nb_cores)" << endl;
		cout << "-f is mandatory"  << endl << "-i performs inexact and speeder research (10 is mandatory value)" << endl << "-p performs pre processing step" << endl << "-c gets the number of threads (default 2)" << endl;
		cout << "-w performs weighted clustering" << endl;
		cout << "Output written in final_g_clusters.txt by default (-o to change output name)" << endl;
	}
}


// find connected components with a DFS, each is stored in vector  nodesInConnexComp
void findConnectedComponents(vector<Node>& vecNodesGlobal, vector<set<uint>>& nodesInConnexComp){
	unordered_set<uint> visited;
	bool b(false);
	for (uint n(0); n < vecNodesGlobal.size(); ++n){
		if (not (visited.count(n))){
			set<uint> s;
			DFS(n, vecNodesGlobal, visited, s, b, 0);
			nodesInConnexComp.push_back(s);
		}
	}
}



// computes a list of cutoffs (stored in vecCC)
// if -i option is not set, the list of cutoff is exactly the list of CC
// else, for connected components with a lot of different CC (more than 100 distinct values), CC values are rounded according to a certain granularity and these rounded values are cutoffs
// then the list of cutoffs is more restrained than the original list of CC and the space to explore is smaller
void computeCutoffs(bool approx, vector<double>& vecCC, vector<double>& ClCo, uint granularity){
	double prev(1.1), cutoffTrunc;
	uint value;
	if (approx){
		prev = 1.1;
		if (ClCo.size() > 100){
			value = granularity;
		} else {
			value = 0;
		}
		for (auto&& cutoff: ClCo){
			if (value != 0){
				cutoffTrunc = trunc(cutoff * value)/value;
			} else {
				cutoffTrunc = cutoff;
			}
			if (cutoffTrunc < prev){
				prev = cutoffTrunc;
				vecCC.push_back(cutoffTrunc);
			}
		}
	} else {
		for (auto&& cutoff: ClCo){
			vecCC.push_back(cutoff);
		}
	}
	
}



bool execute(int argc, char** argv){
	bool printHelp(true);
	bool approx, preprocessing, weighted;
	string outFileName, fileName;
	uint nbThreads, granularity;
	// parsing command line
	parseArgs(argc, argv, approx, preprocessing, weighted, fileName, outFileName, nbThreads, granularity);
	if (not (fileName.empty())){
		printHelp = false;
		cout << "Command line was: " ;
		for (int a(0);  a < argc; ++a){
			cout << argv[a] << " ";
		}
		cout << endl;
		ifstream refFile(fileName);
		vector<Node> vecNodesGlobal;
		cout << "Parsing infile..." << endl;
		// parse similarity information from infile
		parsingSRC(refFile, vecNodesGlobal, weighted);
		if (preprocessing){
			// pre -processing by removing articulation points
			cout << "preprocessing of the graph" << endl;
			preProcessGraph(vecNodesGlobal);
		}
		vector<set<uint>> nodesInConnexComp;
		// decompose graph in connected components
		findConnectedComponents(vecNodesGlobal, nodesInConnexComp);
		cout << "Connected components: " << nodesInConnexComp.size() << endl;

		ofstream out(outFileName);
		ofstream outm("nodes_metrics.txt");
		mutex mm;
		vector<vector<uint>> finalClusters;
		vector<Node> vecNodes;
		vector<double>ClCo, vecCC;
		vector<uint> degrees, nodesInOrderOfCC;
		double minCut(0);
		vector<set<uint>> clustersToKeep;
		uint ccc(0), round(0), higherDegree;
		float lowerCC(0);
		// loop over each connected component
		for (uint c(0); c < nodesInConnexComp.size(); ++c){
			cout << "Connected Component " << c << " size " << nodesInConnexComp[c].size() << endl;
			vecNodes = {};
			
			// compute a vector of nodes for the given connected component
			getVecNodes(vecNodes, vecNodesGlobal, nodesInConnexComp[c]);
			if (preprocessing){
				// pre -processing by removing articulation points
				preProcessGraph(vecNodes);
			}
			ClCo = {};
			// compute CC and degree for each node
			computeCCandDeg(vecNodes, ClCo, degrees, lowerCC); // sorted clustering coefficients
			// compute quantiles for degree distribution
			higherDegree = quantileEdges(degrees, 999, 1000);
			// write nodes metrics
			for (auto&& node : vecNodes){
				outm << node.index << " " << node.CC << " " << node.neighbors.size() << endl;
			}
			vecCC = {}; nodesInOrderOfCC = {};
			// compute a list of cutoffs to loop over
			computeCutoffs(approx, vecCC, ClCo,granularity); // todo: lower bound
			minCut = 0; ccc = 0; round = 0;
			clustersToKeep = {};
			cout << "Computing pseudo cliques" << endl;
			// compute sets originated from seed nodes for each cutoff value
			computePseudoCliques(vecCC, vecNodes, nbThreads, nodesInOrderOfCC, higherDegree, lowerCC);  // todo :  check if we realise a partition
			sortVecNodes(vecNodes, nodesInOrderOfCC);
			cout <<  vecCC.size() << " clustering coefficients to check" << endl;
			bool compute(true);

			// one thread by cutoff
			#pragma omp parallel num_threads(nbThreads)
			{
				#pragma omp for
				for (ccc = 0; ccc < vecCC.size(); ++ccc){
					double cut, prevCut, cutoff(vecCC[ccc]);
					if (ccc != 0){
						if (approx and cutoff == 0 and ccc == vecCC.size() - 1){  // in this case, using the higher cutoff we got cliques, so there is nothing to cut
							mm.lock();
							compute = false;
							mm.unlock();
						}
						if (cut > prevCut){  // we store  non minimal cuts that permit to stop the computation in func computeClustersAndCut anytime an even higher cut is found
							prevCut = cut;
						}
					}
					if (compute){
						vector<Node> vecNodesCpy = vecNodes;
						vector<set<uint>> clusters(vecNodesCpy.size());
						vector<uint> nodesInOrderOfCCcpy = nodesInOrderOfCC;
						cout << "Computing clusters" << endl;
						// refine sets using Clustering coeffs to obtain clusters
						cut = computeClustersAndCut(cutoff, vecNodesCpy, clusters, ccc, prevCut, nodesInOrderOfCCcpy);
						mm.lock();
						cout << round + 1 << "/" << vecCC.size() << " cutoff " << cutoff << " cut " << cut << endl;
						++round;
						mm.unlock();
						// keep the minimal cut and associated clusters:
						if (ccc == 0){
							mm.lock();
							minCut = cut;
							clustersToKeep = clusters;
							if (not weighted){
								if (minCut == 0 and cutoff == 1){
									compute = false;  // clique => stop
								}
							}
							mm.unlock();
						} else {
							if (not weighted){
								if (cut < minCut and cut > 0){
									mm.lock();
									minCut = cut;
									clustersToKeep = clusters;
									mm.unlock();
								}
							} else {
								if (cut < minCut){
									mm.lock();
									minCut = cut;
									clustersToKeep = clusters;
									mm.unlock();
								}
							}
						}
					}
				}
			}
			// print clusters associated to the minimal cut over all cutoff values
			for (uint i(0); i < clustersToKeep.size(); ++i){
				if (not clustersToKeep[i].empty()){
					for (auto&& n : clustersToKeep[i]){
						out << vecNodes[n].index << " " ;
					}
					out << endl;
				}
			}
		}
		cout << "Done." << endl;
	}
	return printHelp;
}
