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

using namespace std;



struct Node{
	uint index;
	uint degree;
	float CC;
	vector<vector<uint>> cluster;
	vector<uint> neighbors;
	bool operator <(const Node&n) const
    {
        return (degree < n.degree);
    }
};

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

vector<float> removeDuplicatesCC(vector<float>& vec){
	sort(vec.begin(), vec.end());
	vec.erase(unique(vec.begin(), vec.end()), vec.end());
	return vec;
}

void DFS(uint n, vector<Node>& vecNodes, unordered_set<uint>& visited, set<uint>& nodesInConnexComp, bool& above, float cutoff){
	unordered_set<uint> neighbors;
	if (not visited.count(n)){
		if (vecNodes[n].CC >= cutoff){
			above = true;
		}
		visited.insert(n);
		nodesInConnexComp.insert(n);
		neighbors = {};
		copy(vecNodes[n].neighbors.begin(), vecNodes[n].neighbors.end(), inserter(neighbors, neighbors.end()));
		for (auto&& neigh : neighbors){
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


void parsingSRC(ifstream & refFile, vector<Node>& vecNodes){
	string listNodes;
	// header
	getline(refFile, listNodes);
	getline(refFile, listNodes);
	getline(refFile, listNodes);
	vector<uint>  neighbs;
	vector<vector<uint>> clust;
	vector<string> splitted1, splitted2, splitted3;
	uint read, target, maxNbNodes(0);
	unordered_map <uint, uint> seenNodes;
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
					if (read != target){
						if (not seenNodes.count(target)){ // new node not already in the vector of nodes
							clust = {}; neighbs = {};
							Node t({target, 0, 0, clust, neighbs});
							vecNodes.push_back({t});  // store in vecNodes
							seenNodes.insert({target, vecNodes.size() - 1}); // remember this node has been pushed index -> place in the vector
						}
						if (seenNodes.count(read)){ // this neighbour is already in the vector of nodes
							vecNodes[seenNodes[target]].neighbors.push_back(seenNodes[read]);  // add read as neighbor of target
							vecNodes[seenNodes[read]].neighbors.push_back(seenNodes[target]);  // add target as neighbor of read
						} else {  // new neighbor not already in the vector of nodes
							clust = {}; neighbs = {};
							Node r({read, 0, 0, clust, neighbs});
							vecNodes.push_back({r});
							uint position(vecNodes.size() - 1);
							seenNodes.insert({read, position});
							vecNodes[seenNodes[target]].neighbors.push_back(vecNodes.size() - 1);  // store read as neighbor of target
							vecNodes[seenNodes[read]].neighbors.push_back(seenNodes[target]);  // add target as neighbor of read
						}
					}
				}
			}
		}   
	}
}


//~ void computeCCandDeg(vector<Node>& vecNodes, set<float, std::greater<float>>& ClCo){
void computeCCandDeg(vector<Node>& vecNodes, vector<float>& ClCo){
	float pairs, clusteringCoef;
	uint totalPairs;
	unordered_set<uint> neighbors;
	// start by removing double occurrences in neighbors
	for (uint n(0); n < vecNodes.size(); ++n){
		vecNodes[n].neighbors = removeDuplicates(vecNodes[n].neighbors);
		vecNodes[n].degree = vecNodes[n].neighbors.size();
	}
	for (uint n(0); n < vecNodes.size(); ++n){
		if (vecNodes[n].degree > 1){
			pairs = 0;
			neighbors = {};
			copy(vecNodes[n].neighbors.begin(), vecNodes[n].neighbors.end(), inserter(neighbors, neighbors.end()));
			for (auto&& neigh : neighbors){  // for each neighbor of the node
				for (auto&& neigh2 : vecNodes[neigh].neighbors){ // for each neighbor of a neighbor
					if (neighbors.count(neigh2)){  // if a neighbor of a neighbor is also a neighbor of the current node = pair of connected neighbors
						++pairs;
					}
				}
			}
			totalPairs = vecNodes[n].neighbors.size() * (vecNodes[n].neighbors.size() - 1);
			if (totalPairs > 0){
				clusteringCoef = pairs/totalPairs;
				vecNodes[n].CC = clusteringCoef;
				ClCo.push_back(clusteringCoef);
				//~ ClCo.insert(clusteringCoef);
			} else {
				//~ ClCo.insert(0);
				ClCo.push_back(0);
			}
		} else {
			//~ ClCo.insert(0);
			ClCo.push_back(0);
		}
	}
	ClCo = removeDuplicatesCC(ClCo);
	sort(ClCo.begin(), ClCo.end());
	reverse(ClCo.begin(), ClCo.end());
}


void sortVecNodes(vector<Node>& vecNodes){
	unordered_map <uint, uint> indexReadsBef;
	for (uint i(0); i < vecNodes.size(); ++i){
		indexReadsBef.insert({vecNodes[i].index, i});
	}
	sort(vecNodes.begin(), vecNodes.end());
	reverse(vecNodes.begin(), vecNodes.end());
	unordered_map <uint, uint> indexReadsAf;
	for (uint i(0); i < vecNodes.size(); ++i){
		indexReadsAf.insert({indexReadsBef[vecNodes[i].index], i});  // former and new index in vecNodes
	}
	
	vector<uint> vec;
	
	for (uint i(0); i < vecNodes.size(); ++i){
		vec = {};
		for (auto&& n : vecNodes[i].neighbors){
			vec.push_back(indexReadsAf[n]);
		}
		vecNodes[i].neighbors = vec;
	}
}


void computePseudoCliques(vector<float>& cutoffs, vector<Node>& vecNodes){
	vector<uint> v;
	float cutoff;
	vector<vector<uint>> vec(cutoffs.size());
	for (uint i(0); i < vecNodes.size(); ++i){
		vecNodes[i].cluster = vec;
	}
	uint c(0), nv(0);
	unordered_set<uint> done;
	
	#pragma omp parallel num_threads(20)
	{
		#pragma omp for
		for (c = 0; c < cutoffs.size(); ++c){
			cutoff = cutoffs[c];
			for (uint i(0); i < vecNodes.size(); ++i){
				if (vecNodes[i].CC >= cutoff){
							vecNodes[i].cluster[c].push_back(i);
							for (auto&& neigh : vecNodes[i].neighbors){
								vecNodes[neigh].cluster[c].push_back(i);
							}
				}
			}
		}
		
		#pragma omp for
		for (nv = 0; nv < vecNodes.size(); ++nv){
			for (uint ii(0); ii < vecNodes[nv].cluster.size(); ++ii){
				vecNodes[nv].cluster[ii] = removeDuplicates(vecNodes[nv].cluster[ii]);
			}
		}
	}
}


float computeUnionCC(set<uint>& unionC, vector<Node>& vecNodes){
	float cardUnion(0);
	unordered_set<uint> neighbors;
	for (auto&& n : unionC){
		neighbors = {};
		copy(vecNodes[n].neighbors.begin(), vecNodes[n].neighbors.end(), inserter(neighbors, neighbors.end()));
		for (auto&& neigh : vecNodes[n].neighbors){
			if (unionC.count(neigh)){
				++cardUnion;
			}
		}
	}
	return cardUnion / ( unionC.size() * (unionC.size() - 1));
}


void transfer(uint tf, uint te, set<uint>& toFill, set<uint>& toEmpty, vector<Node>& vecNodes, vector<set<uint>>& clusters, uint ind){
	vector<uint> vec;
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



vector<set<uint>> assignNewClusters(set<uint>& clust, vector<Node>& vecNodes, float cutoff){
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


void removeSplittedElements(uint index, vector<set<uint>>& clusters, set<uint>& interC){
	set<uint> clust;
	for (auto && elt : clusters[index]){
		if (not interC.count(elt)){
			clust.insert(elt);
		}
	}
	clusters[index] = clust;
}


uint splitClust(uint i1, uint i2, set<uint>& clust1, set<uint>& clust2, vector<set<uint>>& clusters,  vector<Node>& vecNodes, set<uint>& interC, uint cutoff, uint ind){
	uint cut1(0), cut2(0), cut(0);
	for (auto&& node : interC){
		for (auto&& neigh : vecNodes[node].neighbors){
			if (clust1.count(neigh)){
				++cut1;
			}
			if (clust2.count(neigh)){
				++cut2;
			}
		}
	}
	if (clust1.size() == interC.size()){
		transfer(i1, i2, clust1, interC, vecNodes, clusters, ind);
		removeSplittedElements(i2, clusters, interC);
		cut = cut2;
	} else if (clust2.size() == interC.size()){
		transfer(i2, i1, clust2, interC, vecNodes, clusters, ind);
		removeSplittedElements(i1, clusters, interC);
		cut = cut1;
	} else {
		unordered_set <uint> neighbors;
		vector<set<uint>> newClust;
		if (cut1 >= cut2){  // todo = treat equality case
			transfer(i1, i2, clust1, interC, vecNodes, clusters, ind);
			removeSplittedElements(i2, clusters, interC);
			cut = cut2;   // todo * 2 ?
			newClust = assignNewClusters(clust2, vecNodes, cutoff);
			if (newClust.size() > 1){
				for (uint i(0); i < newClust.size(); ++i){
					clusters.push_back(newClust[i]);
					for (auto&& nodes: newClust[i]){
						vecNodes[nodes].cluster[ind].push_back(clusters.size() - 1);
					}
					transfer(clusters.size() - 1, i2, newClust[i], clust2, vecNodes, clusters, ind);
				}
				clusters[i2] = {};
			}
		} else {
			// split clust1
			transfer(i2, i1, clust2, interC, vecNodes, clusters, ind);
			removeSplittedElements(i1, clusters, interC);
			cut = cut1;
			newClust = assignNewClusters(clust1, vecNodes, cutoff);
			
			if (newClust.size() > 1){
				for (uint i(0); i < newClust.size(); ++i){
					clusters.push_back(newClust[i]);
					for (auto&& nodes: newClust[i]){
						vecNodes[nodes].cluster[ind].push_back(clusters.size() - 1);
					}
					transfer(clusters.size() - 1, i1, newClust[i], clust1, vecNodes, clusters, ind);
				}
				clusters[i1] = {};
			}
		}
	}
	return cut;
}


uint computeClustersAndCut(float cutoff, vector<Node>& vecNodes, vector<set<uint>>& clusters, uint ind){
	uint cut(0), cuthalf(0);
	set<uint> clust1, clust2, unionC, interC;
	uint i1, i2;
	float unionCC;
	for (uint n(0); n < vecNodes.size(); ++n){
		//~ vecNodes[n].cluster[ind] = removeDuplicates(vecNodes[n].cluster[ind]);
		if (vecNodes[n].cluster[ind].size() > 0){
			for (auto&& c : vecNodes[n].cluster[ind]){
				clusters[c].insert(n);  // node at index n is in cluster c
			}
		}
	}
	//~ cout << "clusters  at ind" << ind << ":" << endl;
	//~ for (uint i(0); i < clusters.size(); ++i){
		//~ if (not clusters[i].empty()){
			//~ for (auto&& c: clusters[i]){
				//~ cout << vecNodes[c].index << " ";
			//~ }
			//~ cout << endl;
		//~ }
		
	//~ }
	//~ cin.get();
	for (uint i(0); i < vecNodes.size(); ++i){
		if (vecNodes[i].cluster[ind].size() > 1){  // node is in several clusters
			while (vecNodes[i].cluster[ind].size() > 1){
				//~ cout << "node " << vecNodes[i].index << " in " << vecNodes[i].cluster[ind].size() << " clust at ind " << ind << endl;
				
				i1 = vecNodes[i].cluster[ind][0];
				i2 = vecNodes[i].cluster[ind][1];
				clust1 = clusters[i1];
				clust2 = clusters[i2];
				interC = {};
				set_intersection(clust1.begin(), clust1.end(), clust2.begin(), clust2.end(), inserter(interC, interC.begin()));
				if (interC.size() == clust1.size() and clust1.size() == clust2.size()){  // clust1 and clust2 are the same
					transfer(i1, i2, clust1, interC,  vecNodes, clusters, ind);
					//~ cout << "transfer" << endl;
					clusters[i2] = {};
					//~ cout << "node " << vecNodes[i].index << " in " << vecNodes[i].cluster[ind].size() << " clust at ind " << ind << endl;
				} else {
					unionC = {};
					set_union(clust1.begin(), clust1.end(), clust2.begin(), clust2.end(), inserter(unionC, unionC.begin()));
					unionCC = computeUnionCC(unionC, vecNodes);
					if (unionCC >= cutoff){  // merge
						//~ cout << "merge" << endl;
						merge(i1, i2, clust1, clust2, clusters, vecNodes, ind);
					} else {  // split
						//~ cout << "split" << endl;
						cut += splitClust(i1, i2, clust1, clust2, clusters, vecNodes, interC, cutoff, ind);
					}
				}
				//~ cin.get();
			}
		}
	}

	cut = 0;

	

	uint halfcut(0);
	for (uint clust1(0); clust1 < clusters.size(); ++clust1){
		if (not clusters[clust1].empty()){
			for (auto && i : clusters[clust1]){
				for (auto&& n : vecNodes[i].neighbors){
					if (not (clusters[clust1].count(n))){
						if (vecNodes[n].cluster[ind].empty()){
							++cut;
							//~ cout << "appen" << endl;
						} else {
							++halfcut;
							//~ cout << "happen" << endl;
						}
					}
				}
			}
		}
	}
	
	//~ for (uint i(0); i < clusters.size(); ++i){
		//~ if (not clusters[i].empty()){
			//~ for (auto&& node: clusters[i]){
				
			//~ }
		//~ }
		
	//~ }
	//~ cout << "clusters  at end of ind" << ind << ":" << endl;
	//~ for (uint i(0); i < clusters.size(); ++i){
		//~ if (not clusters[i].empty()){
			//~ for (auto&& c: clusters[i]){
				//~ cout << vecNodes[c].index << " ";
			//~ }
			//~ cout << endl;
		//~ }
		
	//~ }
	//~ cin.get();
	return cut + halfcut/2;
}

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


int main(int argc, char** argv){

	if (argc > 1){
		bool approx(false);
		string outFileName("final_g_clusters.txt");
		if (argc > 2){
			if ((string)argv[2] == "--approx"){
				approx = true;
			} else {
				outFileName = argv[2];
			}
			if (argc > 3){
				if ((string)argv[3] == "--approx"){
					approx = true;
				}
			}
		}
		
		string fileName(argv[1]);
        ifstream refFile(fileName);
		vector<Node> vecNodesGlobal;
		cout << "Parsing..." << endl;
		parsingSRC(refFile, vecNodesGlobal);
		unordered_set<uint> visited;
		vector<set<uint>> nodesInConnexComp;
		bool b(false);
		for (uint n(0); n < vecNodesGlobal.size(); ++n){
			if (not (visited.count(n))){
				set<uint> s;
				DFS(n, vecNodesGlobal, visited, s, b, 0);
				nodesInConnexComp.push_back(s);
			}
		}
		cout << "Connected components: " << nodesInConnexComp.size() << endl;

		ofstream out(outFileName);
		mutex mm;

		vector<vector<uint>> finalClusters;
		for (uint c(0); c < nodesInConnexComp.size(); ++c){
			mm.lock();
			cout << "Connected Component " << c << " size " << nodesInConnexComp[c].size() << endl;
			mm.unlock();
			vector<Node> vecNodes;
			getVecNodes(vecNodes, vecNodesGlobal, nodesInConnexComp[c]);
			vector<float> ClCo;
			computeCCandDeg(vecNodes, ClCo);
			vector<float>vecCC;
			if (approx){
				float prev(1.1), cutoffTrunc;
				uint value;
				if (ClCo.size() > 10000){
					value = 1000;
				} else if (ClCo.size() > 1000){
					value = 100;
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
			uint minCut(0);
			vector<set<uint>> clustersToKeep;
			uint ccc(0);
			sortVecNodes(vecNodes);  // sort by decreasing degree
			cout << "Computing pseudo cliques" << endl;
			computePseudoCliques(vecCC, vecNodes);
			#pragma omp parallel num_threads(20)
			{
				#pragma omp for
				for (ccc = 0; ccc < vecCC.size(); ++ccc){
					uint cut;
					float precCutoff = -1;
					float cutoff(vecCC[ccc]);
					if (ccc != 0){
						precCutoff = vecCC[ccc - 1];
					}
					vector<Node> vecNodesCpy = vecNodes;
					vector<set<uint>> clusters(vecNodesCpy.size());
					cut = computeClustersAndCut(cutoff, vecNodesCpy, clusters, ccc);
					mm.lock();
					cout << "cutoff " << cutoff << " cut " << cut << endl;
					mm.unlock();
					if (ccc == 0){
						mm.lock();
						minCut = cut;
						clustersToKeep = clusters;
						mm.unlock();
					} else {
						if (cut < minCut and cut > 0){
							mm.lock();
							minCut = cut;
							clustersToKeep = clusters;
							mm.unlock();
						}
					}
				}
			}
			
			vector <uint> v;
			for (uint i(0); i < clustersToKeep.size(); ++i){
				v = {};
				if (not clustersToKeep[i].empty()){
					for (auto&& n : clustersToKeep[i]){
						v.push_back(vecNodes[n].index);
					}
					mm.lock();
					finalClusters.push_back(v);
					mm.unlock();
				}
			}
		}
		cout << "Writing clusters in output" << endl;
		for (uint i(0); i < finalClusters.size(); ++i){
			for (auto&& n : finalClusters[i]){
				out << n << " ";
			}
			out << endl;
		}
	} else {
		cout << "Usage : ./clustering_cliqueness (src_output) (--approx)" << endl;
		cout << "Output written in final_g_clusters.txt" << endl;
	}
}
