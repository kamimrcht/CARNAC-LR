/*****************************************************************************
 * * * * *   CARNAC: Clustering coefficient-based Acquisition of RNA Communities
 * * * * *
 * * * * *   Authors: Camille Marchet
 * * * * *   Contact: camille.marchet@irisaa.fr, INRIA/IRISA/GenScale, Campus de Beaulieu, 35042 Rennes Cedex, France
 * * * * *   Source: https://github.com/Kamimrcht/CARNAC
 * * * * *
 * * * * *
 * * * * *  This program is free software: you can redistribute it and/or modify
 * * * * *  it under the terms of the GNU Affero General Public License as
 * * * * *  published by the Free Software Foundation, either version 3 of the
 * * * * *  License, or (at your option) any later version.
 * * * * *
 * * * * *  This program is distributed in the hope that it will be useful,
 * * * * *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 * * * * *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * * * * *  GNU Affero General Public License for more details.
 * * * * *
 * * * * *  You should have received a copy of the GNU Affero General Public License
 * * * * *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * * * * * *****************************************************************************/




#include <iostream>
#include <list>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_set>
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
    edges = new vector<int>[nbNodes];
}

void Graph::addEdge(int v, int w)
{
    edges[v].push_back(w);
    edges[w].push_back(v);
}

Graph::~Graph(){
	delete [] edges;
}

//~ // A recursive function that find articulation points using DFS traversal
//~ // u --> The vertex to be visited next
//~ // visited[] --> keeps tract of visited vertices
//~ // disc[] --> Stores discovery times of visited vertices
//~ // parent[] --> Stores parent vertices in DFS tree
//~ // ap[] --> Store articulation points

void Graph::APUtil(int U, bool visited[], int disc[], int low[], int parent[], vector<bool>& ap){
	vector<int> children(ap.size(),0);
	vector<pair<int,int>> mem;
	mem.push_back({U,0});
    static int time = 0;
	while(not mem.empty()){
		pair<int,int> balls(mem[mem.size()-1]);
		int u=balls.first;
		mem.pop_back();
		if(balls.second>=0){
			if(not visited[u]){
				visited[u] = true;
				disc[u] = low[u] = ++time;
			}
			int i=balls.second;
			if( edges[u].empty()){
				continue;
			}
			int v =  edges[u][i];  // v is current adjacent of u
			if((uint)i+1<edges[u].size()){
				mem.push_back({u,i+1});
			}
			if (not visited[v]){
				children[u]++;
				parent[v] = u;

				mem.push_back({u,-v-1});
				mem.push_back({v,0});

			}else if (v != parent[u]){// Update low value of u for parent function calls.
				low[u]  = min(low[u], disc[v]);
			}
		}else{
			int v=-(balls.second+1);
			low[u]  = min(low[u], low[v]);
            if (parent[u] == -1 && children[u]>1)
               ap[u] = true;
            if (parent[u] != -1 && low[v] >= disc[u])
               ap[u] = true;
		}
	}
}



//~ void Graph::APUtil(int u, bool visited[], int disc[], int low[], int parent[], vector<bool>& ap){
    //~ static int time = 0;
    //~ int children = 0;
    //~ visited[u] = true;
    //~ disc[u] = low[u] = ++time;
    //~ for (int i = 0; i != edges[u].size(); ++i){
		//~ int v =  edges[u][i]; // v is current adjacent of u
        //~ if (!visited[v])
        //~ {
            //~ children++;
            //~ parent[v] = u;
            //~ APUtil(v, visited, disc, low, parent, ap);
            //~ low[u]  = min(low[u], low[v]);
            //~ if (parent[u] == -1 && children > 1)
               //~ ap[u] = true;
            //~ if (parent[u] != -1 && low[v] >= disc[u])
               //~ ap[u] = true;
        //~ }
        //~ else if (v != parent[u])
            //~ low[u]  = min(low[u], disc[v]);
    //~ }
//~ }


//~ bool Graph::APUtilBool(int u, bool visited[], int disc[], int low[], int parent[], vector<bool>& ap, set<uint>& interC){
	//~ bool found(false);
    //~ static int time = 0;
    //~ int children = 0;
    //~ visited[u] = true;
    //~ disc[u] = low[u] = ++time;
    //~ for (uint i = 0; i < edges[u].size(); ++i){
		//~ int v =  edges[u][i];
        //~ if (!visited[v])
        //~ {
            //~ children++;
            //~ parent[v] = u;
            //~ APUtilBool(v, visited, disc, low, parent, ap, interC);
            //~ low[u]  = min(low[u], low[v]);
            //~ if (parent[u] == -1 && children > 1)
               //~ ap[u] = true;
            //~ if (parent[u] != -1 && low[v] >= disc[u])
               //~ ap[u] = true;
        //~ }else if (v != parent[u]){
            //~ low[u]  = min(low[u], disc[v]);
		//~ }
		//~ if (ap[u] and interC.count((uint)u)){
			//~ return true;
		//~ }

    //~ }
    //~ return found;
//~ }


bool Graph::APUtilBool(int U, bool visited[], int disc[], int low[], int parent[], vector<bool>& ap, set<uint>& interC){
	vector<int> children(ap.size(),0);
	vector<pair<int,int>> mem;
	mem.push_back({U,0});
    static int time = 0;
	while(not mem.empty()){
		pair<int,int> balls(mem[mem.size()-1]);
		int u=balls.first;
		mem.pop_back();
		if(balls.second>=0){
			if(not visited[u]){
				visited[u] = true;
				disc[u] = low[u] = ++time;
			}
			int i=balls.second;
			if( edges[u].empty()){
				continue;
			}
			int v =  edges[u][i];  // v is current adjacent of u
			if((uint)i+1<edges[u].size()){
				mem.push_back({u,i+1});
			}
			if (not visited[v]){
				children[u]++;
				parent[v] = u;

				mem.push_back({u,-v-1});
				mem.push_back({v,0});

			}else{
				if (v != parent[u]){// Update low value of u for parent function calls.
					low[u]  = min(low[u], disc[v]);
				}
				if (ap[u] and interC.count((uint)u)){
					return true;
				}
			}


		}else{
			int v=-(balls.second+1);
			low[u]  = min(low[u], low[v]);
            if (parent[u] == -1 && children[u]>1)
               ap[u] = true;
            if (parent[u] != -1 && low[v] >= disc[u])
               ap[u] = true;

            if (ap[u] and interC.count((uint)u)){
				return true;
			}
			if (v != parent[u]){// Update low value of u for parent function calls.
				low[u]  = min(low[u], disc[v]);
			}

		}
	}
	return false;
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
        if (visited[i] == false){
			//~ cout<<"go"<<endl;
            APUtil(i, visited, disc, low, parent, ap);
		}
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
