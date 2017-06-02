#include <iostream>
#include <fstream>
#include "clustering_cliqueness.hpp"



int main(int argc, char** argv){
	bool cmd(true);
	if (argc > 1){
		cmd = execute(argc, argv);
	}
	printHelpCmd(cmd);
	return 0;
}
