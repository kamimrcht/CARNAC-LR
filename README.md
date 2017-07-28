CARNAC (Clustering coefficient-based Acquisition of RNA Communities)
====================================================================

# Get CARNAC:

	git clone https://github.com/kamimrcht/CARNAC.git

# Compilation and Usage:

	cd CARNAC;
	make;

## Warning:

For MacOS users (clang version < ), the flag -fopenmp must be removed from Makefile before the compilation. In this configuration CARNAC cannot be run with several threads.



## Usage:


	./CARNAC -f input_file (-o output_file -t nb_cores)


## Test:

	./CARNAC -f test/test.input
The expected result is the same than in test/final_g_clusters_to_obtain.txt


# Options:

	./CARNAC
to output the options
* -f is mandatory
* -t gets the number of threads (default 2)
* Output is written in final_g_clusters.txt by default (-o to change output name)
	
	
# I/O:

## Input format:

## Output format:

Reads are attributed indices from 0 to #reads-1 in order of appearance in input file.
The output is a .txt file with one line per cluster, with the indices of its read members separated by spaces.
For instance:
	0 1 2 3
	4 5
The first four reads are in a cluster, the two last reads are in a second cluster.


# Contact:

camille.marchet@irisa.fr
