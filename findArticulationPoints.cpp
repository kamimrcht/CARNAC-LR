// A C++ program to find articulation points in an undirected graph
#include <iostream>
#include <list>
#include <vector>
#include <set>
#include "findArticulationPoints.hpp"


using namespace std;
 
// A class that represents an undirected graph of reads
//~ class Graph
//~ {

//~ public:
    //~ int nbNodes;  
    //~ list<uint> *edges;    
    //~ void APUtil(uint v, bool visited[], uint disc[], uint low[], 
                //~ uint parent[], bool ap[]);

    //~ Graph(uint nbNodes);   // Constructor
    //~ void addEdge(uint v, uint w);   // function to add an edge to graph
    //~ void AP(bool* ap);    // get articulation points
//~ };
 
Graph::Graph(uint nbNodes){
    this->nbNodes = nbNodes;
    edges = new list<int>[nbNodes];
}
 
void Graph::addEdge(int v, int w)
{
    edges[v].push_back(w);
    edges[w].push_back(v);  
}

Graph::~Graph(){
	delete [] edges;
}
 
// A recursive function that find articulation points using DFS traversal
// u --> The vertex to be visited next
// visited[] --> keeps tract of visited vertices
// disc[] --> Stores discovery times of visited vertices
// parent[] --> Stores parent vertices in DFS tree
// ap[] --> Store articulation points
//~ void Graph::APUtil(uint u, bool visited[], uint disc[], uint low[], uint parent[], bool ap[]){
void Graph::APUtil(int u, bool visited[], int disc[], int low[], int parent[], vector<bool>& ap){
	
    // A static variable is used for simplicity, we can avoid use of static
    // variable by passing a pointer.
    static int time = 0;
 
    // Count of children in DFS Tree
    int children = 0;
 
    // Mark the current node as visited
    visited[u] = true;
 
    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;
 
    // Go through all vertices adjacent to this
    list<int>::iterator i;
    for (i = edges[u].begin(); i != edges[u].end(); ++i)
    {
        int v = *i;  // v is current adjacent of u
 
        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it
        if (!visited[v])
        {
            children++;
            parent[v] = u;
            APUtil(v, visited, disc, low, parent, ap);
 
            // Check if the subtree rooted with v has a connection to
            // one of the ancestors of u
            low[u]  = min(low[u], low[v]);
 
            // u is an articulation point in following cases
 
            // (1) u is root of DFS tree and has two or more chilren.
            if (parent[u] == -1 && children > 1)
               ap[u] = true;
 
            // (2) If u is not root and low value of one of its child is more
            // than discovery value of u.
            if (parent[u] != -1 && low[v] >= disc[u])
               ap[u] = true;
        }
 
        // Update low value of u for parent function calls.
        else if (v != parent[u])
            low[u]  = min(low[u], disc[v]);
    }
}


bool Graph::APUtilBool(int u, bool visited[], int disc[], int low[], int parent[], vector<bool>& ap, set<uint>& interC){
	bool found(false);
    // A static variable is used for simplicity, we can avoid use of static
    // variable by passing a pointer.
    static int time = 0;
 
    // Count of children in DFS Tree
    int children = 0;
 
    // Mark the current node as visited
    visited[u] = true;
 
    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;
 
    // Go through all vertices adjacent to this
    list<int>::iterator i;
    for (i = edges[u].begin(); i != edges[u].end(); ++i)
    {
        int v = *i;  // v is current adjacent of u
 
        // If v is not visited yet, then make it a child of u
        // in DFS tree and recur for it
        if (!visited[v])
        {
            children++;
            parent[v] = u;
            APUtilBool(v, visited, disc, low, parent, ap, interC);
 
            // Check if the subtree rooted with v has a connection to
            // one of the ancestors of u
            low[u]  = min(low[u], low[v]);
 
            // u is an articulation point in following cases
 
            // (1) u is root of DFS tree and has two or more chilren.
            if (parent[u] == -1 && children > 1)
               ap[u] = true;
 
            // (2) If u is not root and low value of one of its child is more
            // than discovery value of u.
            if (parent[u] != -1 && low[v] >= disc[u])
               ap[u] = true;
        }
		if (ap[u] and interC.count((uint)u)){
			return true;
		}
        // Update low value of u for parent function calls.
        else if (v != parent[u])
            low[u]  = min(low[u], disc[v]);
    }
    return found;
}
 
// The function to do DFS traversal. It uses recursive function APUtil()
//~ void Graph::AP(bool* ap)
void Graph::AP(vector<bool>& ap){
    // Mark all the vertices as not visited
    bool *visited = new bool[nbNodes];
    int *disc = new int[nbNodes];
    int *low = new int[nbNodes];
    int *parent = new int[nbNodes];

    // Initialize parent and visited, and ap(articulation point) arrays
    for (uint i = 0; i < (uint)nbNodes; i++)
    {
        parent[i] = -1;
        visited[i] = false;
        ap.push_back(false);
        //~ ap[i] = false;
    }
 
    // Call the recursive helper function to find articulation points
    // in DFS tree rooted with vertex 'i'
    for (uint i = 0; i < (uint)nbNodes; i++)
        if (visited[i] == false)
            APUtil(i, visited, disc, low, parent, ap);
    delete [] visited;
    delete [] disc;
    delete [] low;
    delete [] parent;
}
 
bool Graph::APBool(vector<bool>& ap,  set<uint>& interC){
	bool found(false);
    // Mark all the vertices as not visited
    bool *visited = new bool[nbNodes];
    int *disc = new int[nbNodes];
    int *low = new int[nbNodes];
    int *parent = new int[nbNodes];

    // Initialize parent and visited, and ap(articulation point) arrays
    for (uint i = 0; i < (uint)nbNodes; i++)
    {
        parent[i] = -1;
        visited[i] = false;
        ap.push_back(false);
        //~ ap[i] = false;
    }
 
    // Call the recursive helper function to find articulation points
    // in DFS tree rooted with vertex 'i'
    for (uint i = 0; i < (uint)nbNodes; i++){
        if (visited[i] == false){
            found = APUtilBool(i, visited, disc, low, parent, ap, interC);
            if (found){
				delete [] visited;
				delete [] disc;
				delete [] low;
				delete [] parent;
				return false;
			}
		}
     }
    delete [] visited;
    delete [] disc;
    delete [] low;
    delete [] parent;
    return found;
}
 


