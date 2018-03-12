CARNAC-LR (Clustering coefficient-based Acquisition of RNA Communities in Long Reads)
====================================================================

# Get CARNAC:

	git clone https://github.com/kamimrcht/CARNAC.git

# Requirements
* C++11 and GCC version from 4.9 / CLANG from 3.9
* Python3 for scripts

# Compilation and Usage:

	cd CARNAC;
	./install;

## Warning:

For MacOS users (clang version < ), the flag -fopenmp must be removed from Makefile before the compilation. In this configuration CARNAC cannot be run with several threads.



## Usage:


	./CARNAC-LR -f input_file (-o output_file -t nb_cores)


## Quick test:

	./CARNAC-LR -f test/test.input
The result clusters should be the same than in test/final_g_clusters_to_obtain.txt


# Options:

	./CARNAC-LR
to output the options
* -f is mandatory
* -t gets the number of threads (default 2)
* Output is written in final_g_clusters.txt by default (-o to change output name)
	
	
# I/O:

## Input format:
* Pairwise mApping Format (PAF) (see for instance https://github.com/lh3/minimap/blob/master/minimap.1) can be converted in CARNAC-LR input format using ./scripts/paf_to_CARNAC.py. The fastq/a file (also .gz) used for the run is also mandatory.
* CARNAC-LR can directly read Short Read Connector Linker (see https://github.com/GATB/short_read_connector) output format

Example with example_file.paf and reads.fa:

	python3 ./scripts/paf_to_CARNAC.py example_file.paf reads.fa input_CARNAC.txt
	./CARNAC-lR -f input_CARNAC.txt -o output_CARNAC.txt

## Output format:

Reads are attributed indices from 0 to #reads-1 in order of appearance in input file.
The output is a .txt file with one line per cluster, with the indices of its read members separated by spaces.
For instance with 6 reads in the input in the previous example:

	less output_CARNAC.txt
	0 1 2 3
	4 5

The first four reads are in one cluster, the two last reads are in a second cluster.
Transform clusters file to separated Fasta files:


	./scripts/CARNAC_to_fasta <original_read_file.fa> <CARNAC_clusters_files> [cluster_min_size]

Mandatory arguments are the output of CARNAC followed by the read file. Clusters are output in fasta format, with a file name that correspond to their order of appearance in CARNAC's output. A A minimum size of the clusters to be written can be set.


# Contact:

camille.marchet@irisa.fr
