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
#include "findArticulationPoints.hpp"

#ifndef CLUST
#define CLUST

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


vector<uint> removeDuplicates(vector<uint>& vec);

vector<float> removeDuplicatesCC(vector<float>& vec);

bool findBridge(vector<Node>& vecNodes, set<uint>& cluster, set<uint>& toRemove);

void DFS(uint n, vector<Node>& vecNodes, unordered_set<uint>& visited, set<uint>& nodesInConnexComp, bool& above, float cutoff);

vector<string> split(const string &s, char delim);


void parsingSRC(ifstream & refFile, vector<Node>& vecNodes);


float getCC(unordered_set<uint>& neighbors, vector<Node>& vecNodes);


int getDeltaCC(set<uint>& toRemove, set<uint>& clust1, vector<Node>& vecNodes, float cutoff);

void computeCCandDeg(vector<Node>& vecNodes, vector<float>& ClCo);

void sortVecNodes(vector<Node>& vecNodes);

void computePseudoCliques(vector<float>& cutoffs, vector<Node>& vecNodes, uint nbThreads);

float computeUnionCC(set<uint>& unionC, vector<Node>& vecNodes);

void transfer(uint tf, uint te, set<uint>& toFill, set<uint>& toEmpty, vector<Node>& vecNodes, vector<set<uint>>& clusters, uint ind);

void merge(uint i1, uint i2, set<uint>& clust1, set<uint>& clust2, vector<set<uint>>& clusters,  vector<Node>& vecNodes, uint ind);

vector<set<uint>> assignNewClusters(set<uint>& clust, vector<Node>& vecNodes, float cutoff);

void removeSplittedElements(uint index, vector<set<uint>>& clusters, set<uint>& interC);

uint splitClust(uint i1, uint i2, set<uint>& clust1, set<uint>& clust2, vector<set<uint>>& clusters,  vector<Node>& vecNodes, set<uint>& interC, uint cutoff, uint ind);

uint computeClustersAndCut(float cutoff, vector<Node>& vecNodes, vector<set<uint>>& clusters, uint ind);

void getVecNodes(vector<Node>& vecNodes, vector<Node>& vecNodesGlobal, set<uint>& nodesInConnexComp);

void cutBrigdesInConnectedComp(vector<Node>& vecNodes, uint val);


#endif
