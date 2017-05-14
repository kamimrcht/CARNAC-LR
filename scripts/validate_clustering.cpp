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




vector<string> split(const string &s, char delim){
	stringstream ss(s);
	string item;
	vector<string> elems;
	while (getline(ss, item, delim)) {
		elems.push_back(move(item)); 
	}
	return elems;
}

void parsing(ifstream & file, unordered_map<uint, uint>& readToCluster, uint& nbClusters, uint& nbReads){
	string cluster;
	vector<string> vecClust;
	nbClusters = 0;
	nbReads = 0;
	uint i(0), nr(0);
	while (not file.eof()){
		getline(file, cluster);
		vecClust = split(cluster, ' ');
		if (not vecClust.empty()){
			for(auto&& r : vecClust){
				readToCluster.insert({stoi(r), i});
				++nr;
			}
			++i;
		}
	}
	nbClusters = i;
	nbReads = nr;
}


void representantReads(unordered_map<uint, uint>& readToTrueCluster, vector<vector<uint>>& representantsOfTrueClusters, unordered_map<uint, uint>& readToResultCluster, vector<vector<uint>>& representantsOfResultClusters, uint& unassignedReads, uint nbTrueClusters, uint nbResultClusters){
	unassignedReads = 0;
	for (auto read(readToTrueCluster.begin()); read != readToTrueCluster.end(); ++read){
		if (not readToResultCluster.count(read->first)){
			++ unassignedReads;
		} else {
			if (representantsOfTrueClusters[read->second].empty()){
				vector<uint> v(nbResultClusters, 0);
				representantsOfTrueClusters[read->second] = v;
			}
			++representantsOfTrueClusters[read->second][readToResultCluster[read->first]];
			if (representantsOfResultClusters[readToResultCluster[read->first]].empty()){
				vector<uint> v(nbTrueClusters, 0);
				representantsOfResultClusters[readToResultCluster[read->first]] = v;
			}
			++representantsOfResultClusters[readToResultCluster[read->first]][read->second];
		}
	}
}


void computePrecisionAndRecall(float& recall, float& precision, vector<vector<uint>>& representantsOfTrueClusters, vector<vector<uint>>& representantsOfResultClusters, uint unassignedReads, uint nbReadsTrue, uint nbReadsResult){
	precision = 0; recall = 0;
	uint max(0);
	uint sum(0);
	for (uint c(0); c < representantsOfTrueClusters.size(); ++c){
		max = 0;
		for(auto&& nbReads : representantsOfTrueClusters[c]){
			if (nbReads > max){
				//~ cout << nbReads << endl;
				max = nbReads;
			}
		}
		//~ cout << sum << endl;
		sum += max;
	}
	//~ cout << sum << " " << nbReadsTrue << " " << unassignedReads << endl;
	recall = (float)sum/(float)(nbReadsResult + unassignedReads);
	sum = 0; max = 0;
	for (uint c(0); c < representantsOfResultClusters.size(); ++c){
		max = 0;
		for(auto&& nbReads : representantsOfResultClusters[c]){
			if (nbReads > max){
				max = nbReads;
			}
		}
		sum += max;
	}
	//~ cout << sum << " " << nbReadsResult << endl;
	precision = (float)sum/(float)nbReadsResult;
}

int main(int argc, char** argv){
	if (argc > 2){
		ifstream trueClustersFile(argv[1]);
		ifstream resultClustersFile(argv[2]);
		unordered_map<uint, uint> resultClusters, trueClusters;
		uint nbTrueClusters, nbResultClusters, nbReadsTrue, nbReadsResult, unassignedReads;
		parsing(trueClustersFile, trueClusters, nbTrueClusters, nbReadsTrue);
		//~ cout << nbTrueClusters << " " << nbReadsTrue << endl;
		parsing(resultClustersFile, resultClusters, nbResultClusters, nbReadsResult);
		//~ cout << nbResultClusters << " " << nbReadsResult << endl;
		vector<vector<uint>> representantsOfTrueClusters(nbTrueClusters);
		vector<vector<uint>> representantsOfResultClusters(nbResultClusters);
		//~ cout << "* " << representantsOfTrueClusters.size() << " " << representantsOfResultClusters.size()  << endl;
		representantReads(trueClusters, representantsOfTrueClusters, resultClusters, representantsOfResultClusters, unassignedReads, nbTrueClusters, nbResultClusters);
		float sensitivity, precision;
		computePrecisionAndRecall(sensitivity, precision, representantsOfTrueClusters, representantsOfResultClusters, unassignedReads, nbReadsTrue, nbReadsResult);
		cout << "Sensitivity " << sensitivity << endl;
		cout << "Precision " << precision << endl;
		cout << "F-measure " <<  2 * precision * sensitivity / (precision + sensitivity) << endl;
	}
}
