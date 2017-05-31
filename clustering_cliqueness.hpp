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
	double CC;
	vector<vector<uint>> cluster;
	vector<uint> neighbors;
	unordered_map<uint, double> neighbToWeight;
	bool operator <(const Node&n) const
    {
        return (degree < n.degree);
    }
};


vector<uint> removeDuplicates(vector<uint>& vec);

vector<double> removeDuplicatesCC(vector<double>& vec);

bool findBridge(vector<Node>& vecNodes, set<uint>& cluster, set<uint>& toRemove);

void DFS(uint n, vector<Node>& vecNodes, unordered_set<uint>& visited, set<uint>& nodesInConnexComp, bool& above, double cutoff);

vector<string> split(const string &s, char delim);


void parsingSRC(ifstream & refFile, vector<Node>& vecNodes);


double getCC(unordered_set<uint>& neighbors, vector<Node>& vecNodes);


int getDeltaCC(set<uint>& toRemove, set<uint>& clust1, vector<Node>& vecNodes, double cutoff);

void computeCCandDeg(vector<Node>& vecNodes, vector<double>& ClCo, vector<uint>& degrees, float& lowerCC);

void sortVecNodes(vector<Node>& vecNodes);

//~ void computePseudoCliques(vector<double>& cutoffs, vector<Node>& vecNodes, uint nbThreads, vector<uint>& nodesInOrderOfCC);
void computePseudoCliques(vector<double>& cutoffs, vector<Node>& vecNodes, uint nbThreads, vector<uint>& nodesInOrderOfCC, uint higherDegree, float lowerCC);

double computeUnionCC(set<uint>& unionC, vector<Node>& vecNodes);

void transfer(uint tf, uint te, set<uint>& toFill, set<uint>& toEmpty, vector<Node>& vecNodes, vector<set<uint>>& clusters, uint ind);

void merge(uint i1, uint i2, set<uint>& clust1, set<uint>& clust2, vector<set<uint>>& clusters,  vector<Node>& vecNodes, uint ind);

vector<set<uint>> assignNewClusters(set<uint>& clust, vector<Node>& vecNodes, double cutoff);

void removeSplittedElements(uint index, vector<set<uint>>& clusters, set<uint>& interC);

//~ uint splitClust(uint i1, uint i2, set<uint>& clust1, set<uint>& clust2, vector<set<uint>>& clusters,  vector<Node>& vecNodes, set<uint>& interC, uint cutoff, uint ind);
double splitClust(uint i1, uint i2, set<uint>& clust1, set<uint>& clust2, vector<set<uint>>& clusters,  vector<Node>& vecNodes, set<uint>& interC, uint cutoff, uint ind);

//~ uint computeClustersAndCut(float cutoff, vector<Node>& vecNodes, vector<set<uint>>& clusters, uint ind, uint prevCut, vector<uint>& nodesInOrderOfCC);
double computeClustersAndCut(double cutoff, vector<Node>& vecNodes, vector<set<uint>>& clusters, uint ind, double prevCut, vector<uint>& nodesInOrderOfCC);

void getVecNodes(vector<Node>& vecNodes, vector<Node>& vecNodesGlobal, set<uint>& nodesInConnexComp);

void cutBrigdesInConnectedComp(vector<Node>& vecNodes, uint val);

bool findArticulPoint(set<uint>& cluster, vector<Node>& vecNodes, set<uint>& interC);

void preProcessGraph(vector<Node>& vecNodes, double cutoff);

void mergeOrSplitProcedures(double cutoff, vector<Node>& vecNodes, vector<set<uint>>& clusters, uint ind, double prevCut, vector<uint>& nodesInOrderOfCC, set<uint>& clust1, set<uint>&  clust2, set<uint>&  unionC, set<uint>&  interC, uint& i1, uint& i2, double& cut, uint i);

uint quantileEdges(vector<uint>&degrees, uint no, uint q);

double quantileCC(vector<double>&CC, uint no, uint q);

#endif
