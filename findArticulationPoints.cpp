// A C++ program to find articulation points in an undirected graph
#include<iostream>
#include <list>
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
    edges = new list<uint>[nbNodes];
}
 
void Graph::addEdge(uint v, uint w)
{
    edges[v].push_back(w);
    edges[w].push_back(v);  
}
 
// A recursive function that find articulation points using DFS traversal
// u --> The vertex to be visited next
// visited[] --> keeps tract of visited vertices
// disc[] --> Stores discovery times of visited vertices
// parent[] --> Stores parent vertices in DFS tree
// ap[] --> Store articulation points
void Graph::APUtil(uint u, bool visited[], uint disc[], 
                                      uint low[], uint parent[], bool ap[])
{
    // A static variable is used for simplicity, we can avoid use of static
    // variable by passing a pointer.
    static int time = 0;
 
    // Count of children in DFS Tree
    uint children = 0;
 
    // Mark the current node as visited
    visited[u] = true;
 
    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;
 
    // Go through all vertices adjacent to this
    list<uint>::iterator i;
    for (i = edges[u].begin(); i != edges[u].end(); ++i)
    {
        uint v = *i;  // v is current adjacent of u
 
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
 
// The function to do DFS traversal. It uses recursive function APUtil()
void Graph::AP(bool* ap)
{
    // Mark all the vertices as not visited
    bool *visited = new bool[nbNodes];
    uint *disc = new uint[nbNodes];
    uint *low = new uint[nbNodes];
    uint *parent = new uint[nbNodes];
    //~ bool *ap = new bool[nbNodes]; // To store articulation points
 
    // Initialize parent and visited, and ap(articulation point) arrays
    for (uint i = 0; i < nbNodes; i++)
    {
        parent[i] = -1;
        visited[i] = false;
        ap[i] = false;
    }
 
    // Call the recursive helper function to find articulation points
    // in DFS tree rooted with vertex 'i'
    for (uint i = 0; i < nbNodes; i++)
        if (visited[i] == false)
            APUtil(i, visited, disc, low, parent, ap);
 
    // Now ap[] contains articulation points, print them
    //~ for (uint i = 0; i < nbNodes; i++)
        //~ if (ap[i] == true)
            //~ cout << i << " ";
}
 
// Driver program to test above function
//~ int main()
//~ {
    //~ // Create graphs given in above diagrams
    //~ cout << "\nArticulation points in first graph \n";
    

    //~ Graph g1(5);
    //~ g1.addEdge(1, 0);
    //~ g1.addEdge(0, 2);
    //~ g1.addEdge(2, 1);
    //~ g1.addEdge(0, 3);
    //~ g1.addEdge(3, 4);
    //~ bool *ap1 = new bool[g1.nbNodes]; // To store articulation points
    //~ g1.AP(ap1);
 
    //~ cout << "\nArticulation points in second graph \n";
    //~ Graph g2(4);
    //~ g2.addEdge(0, 1);
    //~ g2.addEdge(1, 2);
    //~ g2.addEdge(2, 3);
    //~ bool *ap2 = new bool[g2.nbNodes]; // To store articulation points
    //~ g2.AP(ap2);
 
    //~ cout << "\nArticulation points in third graph \n";
    //~ Graph g3(7);
    //~ g3.addEdge(0, 1);
    //~ g3.addEdge(1, 2);
    //~ g3.addEdge(2, 0);
    //~ g3.addEdge(1, 3);
    //~ g3.addEdge(1, 4);
    //~ g3.addEdge(1, 6);
    //~ g3.addEdge(3, 5);
    //~ g3.addEdge(4, 5);
    //~ bool *ap3 = new bool[g3.nbNodes]; // To store articulation points
    //~ g3.AP(ap3);
 
    //~ return 0;
//~ }


