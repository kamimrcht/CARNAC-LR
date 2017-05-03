
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
#include <list>

using namespace std;


void DFS(uint n, unordered_map <uint, unordered_set<uint>>& nodeToNeighbors, unordered_set<uint>& visited, unordered_set<uint>& nodesInConnexComp){
	unordered_set<uint> neighbours;
	if (not visited.count(n)){
		visited.insert(n);
		nodesInConnexComp.insert(n);
		neighbours = nodeToNeighbors[n];
		for (auto&& neigh : neighbours){
			DFS(neigh, nodeToNeighbors, visited, nodesInConnexComp);
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

// awaits  input lines like this 0:50-25-50.000000 27-32-64.000000, ignores the positions information
uint parsingSRC(ifstream & refFile, unordered_map <uint, unordered_set<uint>>& nodeToNeighbors){
	string listNodes;
	// header
	getline(refFile, listNodes);
	getline(refFile, listNodes);
	getline(refFile, listNodes);
	vector<string> splitted1, splitted2, splitted3;
	uint read, target, maxnbNode(0);
	while (not refFile.eof()){
		getline(refFile, listNodes);
		splitted1 = split(listNodes, ':');
		if (splitted1.size() > 1){
			splitted2 = split(splitted1[1], ' ');
			unordered_set<uint> reads;
			target = stoi(splitted1[0]);  // target read
			if (target > maxnbNode){
				maxnbNode = target;
			}
			if (not splitted2.empty()){
				for (uint i(0); i < splitted2.size(); ++i){
					splitted3 = split(splitted2[i], '-');
					read = stoi(splitted3[0]);  // recruited read
					
					if (read != target){
						if (read > maxnbNode){
							maxnbNode = read;
						}
						reads.insert(read);
						if (nodeToNeighbors.count(read)){
							nodeToNeighbors[read].insert(target);
						} else {
							nodeToNeighbors.insert({read, {target}});
						}
					}
				}
			}
			if (nodeToNeighbors.count(target)){
				for (auto&& n : reads){
					nodeToNeighbors[target].insert(n);
				}
			} else {
				nodeToNeighbors.insert({target, reads});
			}
		}   
	}
	return maxnbNode;
}

void computeCCandDeg(unordered_map <uint, unordered_set<uint>>& nodeToNeighbors, unordered_map <uint, pair<float, uint>>& nodeToMetrics, map <uint, vector<uint>, std::greater<uint>>& degToNode, set<float, std::greater<float>>& CC){
	uint deg;
	float clusteringCoef(0), pairs(0);
	for (auto node(nodeToNeighbors.begin()); node != nodeToNeighbors.end(); ++node){
		unordered_set <uint> neighbors(node->second.begin(), node->second.end());
		float clusteringCoef(0);
		float pairs(0);
		deg = neighbors.size();
		if (deg > 1){
			for (auto&& neigh : node->second){  // for each neighbor of the node
				for (auto&& neigh2 : nodeToNeighbors.at(neigh)){ // for each neighbor of a neighbor
					if (neighbors.count(neigh2)) {  // if a neighbor of a neighbor is also a neighbor of the current node = pair of connected neighbors
						++pairs;
					}
				}
			}
		}
		uint totalPairs(deg * (deg - 1));
		if (totalPairs > 0){
			clusteringCoef = pairs/totalPairs;
		}
		nodeToMetrics.insert({node->first, {clusteringCoef, deg}});
		if (degToNode.count(deg)){
			degToNode[deg].push_back(node->first);
		} else {
			degToNode.insert({deg, {node->first}});
		}
		CC.insert(clusteringCoef);
	}
	
}


void computePseudoCliques(float cutoff, map <uint, vector<uint>, std::greater<uint>>& degToNode, vector<set<uint>>& pcliqueToNodes, unordered_map <uint, unordered_set<uint>>& nodeToPCliques, unordered_map <uint, pair<float, uint>>& nodeToMetrics, uint& nbPCliques, unordered_map <uint, unordered_set<uint>>& nodeToNeighbors, unordered_set<uint>& pCliquesAboveThresh, vector<list<uint>>& clustersOfNodes){
	unordered_set<uint> neighbors;
	for (auto deg(degToNode.begin()); deg != degToNode.end(); ++deg){
		for (auto&& node : deg->second){  // for all nodes by decreasing degrees
			if (nodeToMetrics[node].first >= cutoff){  // do a new pseudo clique from this node
				if (nodeToPCliques.count(node)){ 
					nodeToPCliques[node].insert(nbPCliques);
				} else {
					nodeToPCliques.insert({node, {nbPCliques}});
				}
				pcliqueToNodes.push_back({node});
				neighbors = nodeToNeighbors[node];
				//~ cout << node << " " << clustersOfNodes.size() << endl;
				clustersOfNodes[node].push_back(nbPCliques);
				for (auto&& neigh: neighbors){  // include neighbours in the pseudo clique
					if (nodeToPCliques.count(neigh)){
						nodeToPCliques[neigh].insert(nbPCliques);
					} else {
						nodeToPCliques.insert({neigh, {nbPCliques}});
					}
					pcliqueToNodes.back().insert(neigh);
					clustersOfNodes[neigh].push_back(nbPCliques);
				}
				 pCliquesAboveThresh.insert(nbPCliques);
				++ nbPCliques;
			} else {  // nodes under the threshold
				if (not nodeToPCliques.count(node)){  // if the node is not already in a cluster create a singleton
					nodeToPCliques.insert({node, {nbPCliques}});
					pcliqueToNodes.push_back({node});
					//~ cout << node << " " << clustersOfNodes.size() << endl;
					clustersOfNodes[node].push_back(nbPCliques);
					++ nbPCliques;
				}
			}
		}
	}
}



void mergeProcedure(vector<set<uint>>& pcliqueToNodes, uint indexC1, uint indexC2, unordered_map <uint, uint>& nodeToCluster, list<uint>::iterator& it1, list<uint>::iterator& it2, vector<list<uint>>& clustersOfNodes, bool& modif1, uint n, bool merginFirst=false){
	if (pcliqueToNodes[indexC1].size() >= pcliqueToNodes[indexC2].size() or merginFirst){  // merge c2 in c1
		pcliqueToNodes[indexC1].insert(pcliqueToNodes[indexC2].begin(), pcliqueToNodes[indexC2].end());
		pcliqueToNodes[indexC2].clear();
		for (auto&& node : pcliqueToNodes[indexC1]){
			if (nodeToCluster.count(node)){
				nodeToCluster[node] = indexC1;
			} else {
				nodeToCluster.insert({node, indexC1});
			}
			//~ cout << "erase cluster " << *it2  << " from node " << node << endl;
			//~ cout << "add cluster " << indexC1  << " to node " << node << endl;
			if (node != n){
				if (pcliqueToNodes[indexC2].count(node)){
					clustersOfNodes[node].remove(*it2);
					clustersOfNodes[node].push_back(indexC1);
				}
			}
			//~ cout << "done" << endl;
		}
	} else { // merge c1 in c2
		pcliqueToNodes[indexC2].insert(pcliqueToNodes[indexC1].begin(), pcliqueToNodes[indexC1].end());
		pcliqueToNodes[indexC1].clear();
		for (auto&& node : pcliqueToNodes[indexC2]){
			if (nodeToCluster.count(node)){
				nodeToCluster[node] = indexC2;
			} else {
				nodeToCluster.insert({node, indexC2});
			}
			//~ cout << "erase cluster " << *it1 << " from node " << node <<  endl;
			//~ cout << "add cluster " << indexC2  << " to node " << node << endl;
			if (node != n){
				if (pcliqueToNodes[indexC1].count(node)){
					clustersOfNodes[node].remove(*it1);
					clustersOfNodes[node].push_back(indexC2);
				}
			}
			//~ cout << "done" << endl;
			modif1 = true;
		}
	}
	
}


void assignClusterDFS(uint node, unordered_map <uint, set<uint>>& temporaryClusters, unordered_map <uint, uint>& temporaryNodesToClusters, uint& indexCluster, unordered_map <uint, unordered_set<uint>>& nodeToNeighbors, set<uint>& pcliqueOut, set<uint>& pcliqueIn, unordered_set<uint>& pCliquesAboveThresh, unordered_map <uint, pair<float, uint>>& nodeToMetrics, float cutoff){
		bool above(false);
		if (temporaryClusters.count(indexCluster)){
			temporaryClusters[indexCluster].insert(node);
			temporaryNodesToClusters.insert({node, indexCluster});
		} else {
			temporaryClusters.insert({indexCluster, {node}});
			temporaryNodesToClusters.insert({node, indexCluster});
			if (nodeToMetrics[node].first >= cutoff){
				pCliquesAboveThresh.insert(indexCluster);
			}
		}
		for (auto&& neigh : nodeToNeighbors[node]){
			
			if (not temporaryNodesToClusters.count(neigh)){
				if (not (pcliqueOut.count(neigh)) and pcliqueIn.count(neigh)){
					assignClusterDFS(neigh, temporaryClusters, temporaryNodesToClusters, indexCluster, nodeToNeighbors, pcliqueOut, pcliqueIn, pCliquesAboveThresh, nodeToMetrics, cutoff);
				}
			}
		}

}


uint cutProcedure(set<uint>& interC, unordered_map <uint, unordered_set<uint>>& nodeToNeighbors, vector<set<uint>>& pcliqueToNodes, uint indexC1, uint indexC2, unordered_map <uint, uint>& nodeToCluster, unordered_set<uint>& pCliquesAboveThresh, unordered_map <uint, pair<float, uint>>& nodeToMetrics, float cutoff, list<uint>::iterator& it1, list<uint>::iterator& it2, vector<list<uint>>& clustersOfNodes, bool& modif1, uint n){
	uint cut1(0), cut2(0), cut(0);
	for (auto&& node : interC){
		for (auto&& neigh : nodeToNeighbors[node]){
			if (pcliqueToNodes[indexC1].count(neigh)){
				++cut1;
			} else if (pcliqueToNodes[indexC2].count(neigh)){
				++cut2;
			}
		}
	}
	if (cut1 <= cut2){  // transfer nodes in cluster 2
		modif1 = true;
		for (auto&& node : interC){
			pcliqueToNodes[indexC1].erase(node);
			if (nodeToCluster.count(node)){
				nodeToCluster[node] = indexC2;
			} else {
				nodeToCluster.insert({node, indexC2});
			}
			if (n != node){
				//~ cout << "Remove cluster " << *it1 << " from node " << node << endl;
				
				clustersOfNodes[node].remove(*it1);
				//~ cout << "done" << endl;
			}
		}
		// nodes from the cluster that transferred can be disjoint
		unordered_map <uint, set<uint>> temporaryClusters;
		unordered_map <uint,uint> temporaryNodesToClusters;
		uint indexCluster(pcliqueToNodes.size());
		for (auto&& node : pcliqueToNodes[indexC1]){
			if (not temporaryNodesToClusters.count(node)){
				assignClusterDFS(node, temporaryClusters, temporaryNodesToClusters, indexCluster, nodeToNeighbors, interC, pcliqueToNodes[indexC1], pCliquesAboveThresh,nodeToMetrics, cutoff);
				++indexCluster;
			}
			
		}
		if (temporaryClusters.size() > 1){  // several connected component in the remaining c1
			//~ modif = true;
			pcliqueToNodes[indexC1] = {};
			for (auto newClust(temporaryClusters.begin()); newClust != temporaryClusters.end(); ++newClust){
				pcliqueToNodes.push_back(newClust->second);
				for (auto&& nn : newClust->second){
					if (nn != n and not (interC.count(nn))){
						//~ cout << "remove cluster " << *it1 << " from node " << nn << endl;
						clustersOfNodes[nn].remove(*it1);
						//~ cout << "done" << endl;
					}
					//~ cout << "new cluster " << pcliqueToNodes.size() - 1 << " for node " << nn << endl;
					clustersOfNodes[nn].push_back(pcliqueToNodes.size() - 1);
					//~ cout << "done" << endl;
				}
			}
		} else {
			pCliquesAboveThresh.erase(pcliqueToNodes.size());
		}
		cut += cut2;
	} else {  // transfer nodes in cluster 1
		for (auto&& node : interC){
			pcliqueToNodes[indexC2].erase(node);
			if (nodeToCluster.count(node)){
				nodeToCluster[node] = indexC1;
			} else {
				nodeToCluster.insert({node, indexC1});
			}
			if (n != node){
				//~ if (node == 2524){
					//~ for(auto&& c: clustersOfNodes[node]){
						//~ cout << c << " " ;
					//~ }
					//~ cout << endl;
				//~ }
				//~ cout << "Remove cluster " << *it2 << " from node " << node << endl; 
				clustersOfNodes[node].remove(*it2);
				//~ if (node == 2524){
					//~ for(auto&& c: clustersOfNodes[node]){
						//~ cout << c << " " ;
					//~ }
					//~ cout << endl;
				//~ }
				//~ cout << "done" << endl;
			}
		}

		unordered_map <uint, set<uint>> temporaryClusters;
		unordered_map <uint, uint> temporaryNodesToClusters;
		uint indexCluster(pcliqueToNodes.size());
		for (auto&& node : pcliqueToNodes[indexC2]){
			if (not temporaryNodesToClusters.count(node)){
				assignClusterDFS(node, temporaryClusters, temporaryNodesToClusters, indexCluster, nodeToNeighbors, interC, pcliqueToNodes[indexC2], pCliquesAboveThresh, nodeToMetrics, cutoff);
				++indexCluster;
			}
		}
		if (temporaryClusters.size() > 1){  // several connected component in the remaining c2
			//~ modif = true;
			pcliqueToNodes[indexC2] = {};
			for (auto newClust(temporaryClusters.begin()); newClust != temporaryClusters.end(); ++newClust){
				pcliqueToNodes.push_back(newClust->second);
				for (auto&& nn : newClust->second){
					if (nn != n and not (interC.count(nn))){
						//~ cout << "remove cluster " << *it2 << " from node " << nn << endl;
						clustersOfNodes[nn].remove(*it2);
						//~ cout << "done" << endl;
					}
					//~ cout << "new cluster " << pcliqueToNodes.size() - 1 << " for node " << nn << endl;
					clustersOfNodes[nn].push_back(pcliqueToNodes.size() - 1);
					//~ cout << "done" << endl;
				}
			}
		} else {
			pCliquesAboveThresh.erase(pcliqueToNodes.size());
		}
		cut += cut1;
	}
	return cut;
}


uint  computeClustersAndCut(vector<set<uint>>& pcliqueToNodes, unordered_map <uint, unordered_set<uint>>& nodeToNeighbors, float cutoff, unordered_set<uint>& nodeSingletons, unordered_set<uint>& pCliquesAboveThresh,  unordered_map <uint, pair<float, uint>>& nodeToMetrics, vector<list<uint>>& clustersOfNodes){
	float unionCC(0);
	uint cut(0);
	bool modif(false), modif1(false);
	unordered_map <uint, uint> nodeToCluster;
	cout << "Comparing " << pcliqueToNodes.size() << " clusters (cutoff " << cutoff<< ")" <<  endl;
	for (uint i(0); i < clustersOfNodes.size(); ++i){
		if (clustersOfNodes[i].size() > 1){
			//~ cout << "node " << i << " size cl " << clustersOfNodes[i].size()<< endl;
			// todo add while clustersOfNodes[i].size() > 1
			//~ while (clustersOfNodes[i].size() > 1){
				auto it(clustersOfNodes[i].begin());
				uint clust1(*it);
				for (auto it2(clustersOfNodes[i].begin()), end = clustersOfNodes[i].end(); it2 != end;) {
					//~ cout << "clusters " << *it << " " << *it2 << endl;
					modif = false;
					uint clust2(*it2);
					if ((not pcliqueToNodes[clust1].empty()) and (not pcliqueToNodes[clust2].empty()) and pCliquesAboveThresh.count(clust1) and pCliquesAboveThresh.count(clust2) and clust1 != clust2){
						set<uint> unionC;
						set<uint> interC;
						unionCC = 0;
						set_intersection(pcliqueToNodes[clust1].begin(), pcliqueToNodes[clust1].end(), pcliqueToNodes[clust2].begin(), pcliqueToNodes[clust2].end(), inserter(interC, interC.begin()));
						if (interC.size() == pcliqueToNodes[clust1].size() and interC.size() == pcliqueToNodes[clust2].size()){  // the two clusters are exactly the same
							pcliqueToNodes[clust2] = {};
							it2 = clustersOfNodes[i].erase(it2);
							modif = true;
						} else {
							set_union(pcliqueToNodes[clust1].begin(), pcliqueToNodes[clust1].end(), pcliqueToNodes[clust2].begin(), pcliqueToNodes[clust2].end(), inserter(unionC, unionC.begin()));
							float cardUnion(0);
							for (auto&& node : unionC){
								for (auto&& neigh : nodeToNeighbors[node]){
									if (unionC.count(neigh)){
										++cardUnion;
									}
								}
							}
							unionCC = cardUnion / (unionC.size() * (unionC.size() - 1));
							if (unionCC >= 100 * cutoff/100){
								//~ cout << "merge 1" << endl;
								mergeProcedure(pcliqueToNodes, clust1, clust2, nodeToCluster, it, it2 , clustersOfNodes, modif1, i);
								if (modif1){
									 it = clustersOfNodes[i].erase(it);
								} else {
									 it2 = clustersOfNodes[i].erase(it2);
								}
								modif = true;
							} else {
								//~ cout << "cut" << endl;
								cutProcedure(interC, nodeToNeighbors, pcliqueToNodes, clust1, clust2, nodeToCluster, pCliquesAboveThresh, nodeToMetrics, cutoff, it, it2 , clustersOfNodes, modif1, i);
								//~ cutProcedure(interC, nodeToNeighbors, pcliqueToNodes, clust1, clust2, nodeToCluster, pCliquesAboveThresh, nodeToMetrics, cutoff,  nodeToCluster, it,  it2, clustersOfNodes, modif1, i);
								if (modif1){
									it = clustersOfNodes[i].erase(it);
								} else {
									it2 = clustersOfNodes[i].erase(it2);
								}
								modif = true;
							}
						}
					}  else if ((not pcliqueToNodes[clust2].empty()) and pCliquesAboveThresh.count(clust1) and clust1 != clust2){
						set<uint> interC;
						set_intersection(pcliqueToNodes[clust1].begin(), pcliqueToNodes[clust1].end(), pcliqueToNodes[clust2].begin(), pcliqueToNodes[clust2].end(), inserter(interC, interC.begin()));
						//~ if (not interC.empty()){
							modif = true;
							//~ cout << "merge 2" << endl;
							mergeProcedure(pcliqueToNodes, clust1, clust2, nodeToCluster, it, it2 , clustersOfNodes, modif1, i, true);
							if (modif1){
								 it = clustersOfNodes[i].erase(it);
							} else {
								 it2 = clustersOfNodes[i].erase(it2);
							}
						//~ }
						
					} else if ((not pcliqueToNodes[clust2].empty()) and pCliquesAboveThresh.count(clust2) and clust1 != clust2){
						set<uint> interC;
						set_intersection(pcliqueToNodes[clust1].begin(), pcliqueToNodes[clust1].end(), pcliqueToNodes[clust2].begin(), pcliqueToNodes[clust2].end(), inserter(interC, interC.begin()));
						//~ if (not interC.empty()){
							modif = true;
							//~ cout << "merge 3" << endl;
							mergeProcedure(pcliqueToNodes, clust2, clust1, nodeToCluster, it, it2 , clustersOfNodes, modif1, i, true);
							if (modif1){
								 it = clustersOfNodes[i].erase(it);
							} else {
								 it2 = clustersOfNodes[i].erase(it2);
							}
						//~ }
					}
					//~ cout << modif << " " << modif1 << endl;
					if (modif1){
						//~ cout << "modif1" << endl;
						break;
					}
					if (not modif){
						 ++it2;
					}
				}
			}
		//~ }
	}
	//~ for (uint clust1(0); clust1 < pcliqueToNodes.size(); ++clust1){
		//~ for (uint clust2(0); clust2 < pcliqueToNodes.size(); ++clust2){
			//~ if ((not pcliqueToNodes[clust1].empty()) and (not pcliqueToNodes[clust2].empty()) and pCliquesAboveThresh.count(clust1) and pCliquesAboveThresh.count(clust2) and clust1 != clust2){
				//~ set<uint> unionC;
				//~ set<uint> interC;
				//~ unionCC = 0;
				//~ set_intersection(pcliqueToNodes[clust1].begin(), pcliqueToNodes[clust1].end(), pcliqueToNodes[clust2].begin(), pcliqueToNodes[clust2].end(), inserter(interC, interC.begin()));
				//~ if (not interC.empty()){  // some nodes belong to both clusters
					//~ modif = true;
					//~ if (interC.size() ==pcliqueToNodes[clust1].size() and interC.size() == pcliqueToNodes[clust2].size()){
						//~ pcliqueToNodes[clust2] = {};
					//~ } else {
						//~ set_union(pcliqueToNodes[clust1].begin(), pcliqueToNodes[clust1].end(), pcliqueToNodes[clust2].begin(), pcliqueToNodes[clust2].end(), inserter(unionC, unionC.begin()));
						//~ float cardUnion(0);
						//~ for (auto&& node : unionC){
							//~ for (auto&& neigh : nodeToNeighbors[node]){
								//~ if (unionC.count(neigh)){
									//~ ++cardUnion;
								//~ }
							//~ }
						//~ }
						//~ unionCC = cardUnion / (unionC.size() * (unionC.size() - 1));
						//~ if (unionCC >= 100 * cutoff/100){
							//~ mergeProcedure(pcliqueToNodes, clust1, clust2, nodeToCluster);
							//~ modif1 = true;
						//~ } else {
							//~ cutProcedure(interC, nodeToNeighbors,  pcliqueToNodes, clust1, clust2, nodeToCluster, pCliquesAboveThresh, nodeToMetrics, cutoff, modif1);
						//~ }
					//~ }
				//~ }
			//~ } else if ((not pcliqueToNodes[clust2].empty()) and pCliquesAboveThresh.count(clust1) and clust1 != clust2){
				//~ set<uint> interC;
				//~ set_intersection(pcliqueToNodes[clust1].begin(), pcliqueToNodes[clust1].end(), pcliqueToNodes[clust2].begin(), pcliqueToNodes[clust2].end(), inserter(interC, interC.begin()));
				//~ if (not interC.empty()){
					//~ modif = true; 
					//~ mergeProcedure(pcliqueToNodes, clust1, clust2, nodeToCluster, true);
					//~ modif1 = true;
				//~ }
				
			//~ } else if ((not pcliqueToNodes[clust2].empty()) and pCliquesAboveThresh.count(clust2) and clust1 != clust2){
				//~ set<uint> interC;
				//~ set_intersection(pcliqueToNodes[clust1].begin(), pcliqueToNodes[clust1].end(), pcliqueToNodes[clust2].begin(), pcliqueToNodes[clust2].end(), inserter(interC, interC.begin()));
				//~ if (not interC.empty()){
					//~ modif = true;
					//~ mergeProcedure(pcliqueToNodes, clust2, clust1, nodeToCluster, true);
					//~ modif1 = true;
				//~ }
			//~ }
		//~ }
	//~ }

	cut = 0;
	for (uint clust1(0); clust1 < pcliqueToNodes.size(); ++clust1){
		for (auto && node : pcliqueToNodes[clust1]){
			for (auto&& n : nodeToNeighbors[node]){
				if (not (pcliqueToNodes[clust1].count(n))){
					++cut;
				}
			}
		}
	}
	return cut/2;
}


int main(int argc, char** argv){

	if (argc > 1){
		string outFileName("final_g_clusters.txt");
		if (argc > 2){
			outFileName = argv[2];
		}
		
		string fileName(argv[1]);
        ifstream refFile(fileName);
        unordered_map <uint, unordered_set<uint>> nodeToNeighborsGlobal;
        // parse SRC's output, for each line we get a node and its neighbors and fill the map nodeToNeighbors
        cout << "Parsing..." << endl;
		uint maxNbNode(parsingSRC(refFile, nodeToNeighborsGlobal));

		uint nbConnexComp(0);
		unordered_set<uint> visited;
		vector<unordered_set<uint>> nodesInConnexComp;
		for (auto node(nodeToNeighborsGlobal.begin()); node != nodeToNeighborsGlobal.end(); ++node){
			if (not (visited.count(node->first))){
				unordered_set<uint> s;
				nodesInConnexComp.push_back(s);
				DFS(node->first, nodeToNeighborsGlobal, visited, nodesInConnexComp.back());
				++ nbConnexComp;
			}
		}
		cout << "Connected components: " << nbConnexComp << endl;

		ofstream out(outFileName);
		mutex mm;
		uint c(0);
		//~ #pragma omp parallel for
		for (c=0; c < nodesInConnexComp.size(); ++c){
			unordered_map <uint, unordered_set<uint>> nodeToNeighbors;
			for (auto node(nodeToNeighborsGlobal.begin()); node != nodeToNeighborsGlobal.end(); ++node){
				if (nodesInConnexComp[c].count(node->first)){
					nodeToNeighbors.insert({node->first, node->second});
				}
			}
			unordered_map <uint, pair<float, uint>> nodeToMetrics;
			map <uint, vector<uint>, std::greater<uint>> degToNode;
			set <float, std::greater<float>> CC;
			// get CC and degree for each node
			cout << "Computing graph..." << endl;
			computeCCandDeg(nodeToNeighbors, nodeToMetrics, degToNode, CC);
			ofstream outm("nodes_metrics.txt");
			mm.lock();
			for (auto node(nodeToMetrics.begin()); node != nodeToMetrics.end(); ++node){
				outm << node->first << " " << node->second.first << " " << node->second.second << endl; 
			}
			mm.unlock();
			//~ unordered_map <uint, set<uint>> pcliqueToNodes;
			vector< set<uint>> pcliqueToNodes;
			unordered_map <uint, unordered_set<uint>> nodeToPCliques;
			uint nbPCliques(0), cut(0);
			vector<uint> globalCut;
			// compute min cut for each clustering coef value
			unordered_set<uint> nodeSingletons;
			unordered_set<uint> pCliquesAboveThresh;
			
			cout << "Testing cutoff values..." << endl;
			for (auto&& cutoff: CC){
				nbPCliques = 0;
				pcliqueToNodes = {};
				nodeToPCliques = {};
				nodeSingletons = {};
				pCliquesAboveThresh = {};
				vector<list<uint>> clustersOfNodes(maxNbNode + 1);
				computePseudoCliques(cutoff, degToNode, pcliqueToNodes, nodeToPCliques, nodeToMetrics, nbPCliques, nodeToNeighbors, pCliquesAboveThresh, clustersOfNodes);
				cout << "compute cut" << endl;
				cut = computeClustersAndCut(pcliqueToNodes, nodeToNeighbors, cutoff, nodeSingletons, pCliquesAboveThresh, nodeToMetrics, clustersOfNodes);
				globalCut.push_back(cut);
			}
			uint i(0), minCut(0), val(0);
			// get the min cut over all cc values
			for (auto&& cut : globalCut){
				if (i == 0){  // if the graph is composed only of clique the mincut will be zeros since the first cc
					minCut = cut;
					val = i;
				} else {
					if (cut < minCut and cut > 0){
						val = i;
						minCut = cut;
					}
				}
				++i;
			}
			
			uint index(0);
			for (auto&& cutoff: CC){
				if (index == val){
					cout <<"computing final cluster with cutoff:" <<  cutoff << " and min cut:" << minCut << endl;
					pcliqueToNodes = {};
					nodeToPCliques = {};
					nbPCliques = 0;
					nodeSingletons = {};
					pCliquesAboveThresh = {};
					vector<list<uint>> clustersOfNodes(maxNbNode + 1);
					computePseudoCliques(cutoff, degToNode, pcliqueToNodes, nodeToPCliques, nodeToMetrics, nbPCliques, nodeToNeighbors, pCliquesAboveThresh, clustersOfNodes);
					computeClustersAndCut(pcliqueToNodes, nodeToNeighbors, cutoff, nodeSingletons, pCliquesAboveThresh, nodeToMetrics, clustersOfNodes);
					break;
				}
				++index;
			}
			for (uint p(0); p < pcliqueToNodes.size(); ++p){
				if (not pcliqueToNodes[p].empty()){
					mm.lock();
					for (auto&& node : pcliqueToNodes[p]){
						out << node << " " ;
					}
					out << endl;
					mm.unlock();
				}
			}
		}
	} else {
		cout << "Usage : ./clustering_cliqueness src_output" << endl;
		cout << "Output written in final_g_clusters.txt" << endl;
	}
}
