#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <unistd.h>
#include <vector>
#include <stdlib.h>
#include <unordered_map>



using namespace std;



char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}



string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}



string getCanonical(const string& str){
	return (min(str,revComp(str)));
}

uint count_str(const string& str,char c){
	uint res(0);
	for(uint i(0);i<str.size();++i){
		if(str[i]==c){
			++res;
		}
	}
	return res;
}



void help(){
	cout<<"./CARNAC_2_fasta reads.fa|q CARNAC_output min_size_cluster"<<endl;
}


int main(int argc, char *argv[]) {
	if(argc<2){
		help();
		exit(0);
	}
	vector<string> reads;
	string readFile(argv[1]);
	string line,number;
	string clusterFile(argv[2]);
	uint minSizeCluster(1);
	if (argc > 2){
		minSizeCluster = stoi(argv[3]);
	}
	string unitig,useless,useless2,useless3,msp;
	ifstream readStream(readFile);
	ifstream  clusterStream(clusterFile);
	if(not clusterStream.good()){
		cerr<<"Problem with cluster file"<<endl;
		return 0;
	}
	if(not readStream.good()){
		cerr<<"Problem with read file"<<endl;
		return 0;
	}
	//~ ofstream out("numbers.txt");
	char c (readStream.peek());
	bool fastaMode(false);
	if(c=='>'){
		fastaMode=true;
	}
	if(fastaMode){
		while(not readStream.eof()){
			getline(readStream,useless);
			getline(readStream,line);
			if(line.size()>2){
				reads.push_back(useless+"\n"+line);
			}
		}
	}else{
		while(not readStream.eof()){
			getline(readStream,useless);
			getline(readStream,useless2);
			getline(readStream,useless3);
			getline(readStream,line);
			if(line.size()>2){
				reads.push_back(useless+"\n"+useless2+"\n"+useless3+"\n"+line);
			}
		}
	}
	uint clusterNumber(0);
	while(not clusterStream.eof()){
		getline(clusterStream,line);
		if(count_str(line,' ')>=minSizeCluster){
			ofstream out(("cluster"+to_string(clusterNumber)+".fasta").c_str(),ofstream::out);
			uint64_t i(1),lasti(0);
			while(i<line.size()){
				if(line[i]==' '){
					number=line.substr(lasti,i-lasti);
					lasti=i+1;
					uint uNumber=stoi(number);
					out<< reads[uNumber]<<"\n";
				}
				++i;
			}
			number=line.substr(lasti,i-lasti);
			if(number.size()>0){
				uint uNumber=stoi(number);
				out << reads[uNumber]<<"\n";
			}
		}
		clusterNumber++;
	}

    return 0;
}
