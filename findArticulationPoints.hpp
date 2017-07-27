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



#include<iostream>
#include <list>

#ifndef PREPROC
#define PREPROC


using namespace std;
 
// A class that represents an undirected graph of reads
class Graph
{

public:
    uint nbNodes;  
    list<int> *edges;    
    //~ void APUtil(uint v, bool visited[], uint disc[], uint low[], uint parent[], bool ap[]);
    void APUtil(int v, bool visited[], int disc[], int low[], int parent[], vector<bool>& ap);
    bool APUtilBool(int v, bool visited[], int disc[], int low[], int parent[], vector<bool>& ap, set<uint>& interC);

    Graph(uint nbNodes);   // Constructor
    ~Graph();
    void addEdge(int v, int w);   // function to add an edge to graph
    //~ void AP(bool* ap);    // get articulation points
    void AP(vector<bool>& ap);    // get articulation points
    bool APBool(vector<bool>& ap, set<uint>& interC);    // get articulation points
};
 


#endif
