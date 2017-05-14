// A C++ program to find articulation points in an undirected graph
#include<iostream>
#include <list>

#ifndef PREPROC
#define PREPROC


using namespace std;
 
// A class that represents an undirected graph of reads
class Graph
{

public:
    int nbNodes;  
    list<uint> *edges;    
    //~ void APUtil(uint v, bool visited[], uint disc[], uint low[], uint parent[], bool ap[]);
    void APUtil(uint v, bool visited[], uint disc[], uint low[], uint parent[], vector<bool>& ap);
    bool APUtilBool(uint v, bool visited[], uint disc[], uint low[], uint parent[], vector<bool>& ap, set<uint>& interC);

    Graph(uint nbNodes);   // Constructor
    ~Graph();
    void addEdge(uint v, uint w);   // function to add an edge to graph
    //~ void AP(bool* ap);    // get articulation points
    void AP(vector<bool>& ap);    // get articulation points
    bool APBool(vector<bool>& ap, set<uint>& interC);    // get articulation points
};
 


#endif
