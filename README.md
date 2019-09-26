CARNAC-LR (Clustering coefficient-based Acquisition of RNA Communities in Long Reads)
====================================================================

# Get CARNAC:
* Using bioconda:
Thanks to @npavlovikj, CARNAC-LR is now available on bioconda: https://anaconda.org/bioconda/carnac-lr, for both Linux and OSX. 

* Using git clone:
	git clone https://github.com/kamimrcht/CARNAC.git
then compile and install:
	cd CARNAC;
	./install;

# Requirements
* C++11 and GCC version from 4.9 / CLANG from 3.9
* Python3 for scripts

## Warning:

For MacOS users (clang version < ), the flag -fopenmp must be removed from Makefile before the compilation. In this configuration CARNAC cannot be run with several threads.



## Usage:
First compute overlaps between reads using [minimap2](https://github.com/lh3/minimap2). 
It is very likely that minimap2 will return as primary alignments the reads mapped on themselves. To prevent this, and obtain more read connections in minimap2 output, I recommend to use -X option:

	minimap2 reads.fq reads.fq -X > minimap_output.paf

Then convert to CARNAC-LR format:

	python CARNAC-LR/scripts/paf_to_CARNAC.py minimap_output.paf reads.fq input_carnac.txt

Before running CARNAC-LR I recommend to increase the stack size

	ulimit -s unlimited

And then, launch CARNAC-LR:

	./CARNAC-LR -f input_carnac.txt (-o output_file -t nb_cores)

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

camille.marchet@univ-lille.fr
