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
#include "clustering_cliqueness.hpp"

int main(int argc, char** argv){
	bool cmd(true);
	if (argc > 1){
		cmd = execute(argc, argv);
	}
	printHelpCmd(cmd);
	return 0;
}
