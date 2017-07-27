CARNAC (Clustering coefficient-based Acquisition of RNA Communities)
====================================================================


Compilation and Usage:
----------------------
	make;

* Warning:
For MacOS users (clang version < ), the flag -fopenmp must be removed from Makefile before the compilation. In this configuration CARNAC cannot be run with several threads.



* Usage :

	./CARNAC -f input_file (-o output_file -t nb_cores)

* Options:
	* -f is mandatory
	* -t gets the number of threads (default 2)
	* Output is written in final_g_clusters.txt by default (-o to change output name)
	
	./CARNAC
to output the options


Test:
-----
./CARNAC -f test/test.input
The expected result is the same than in test/final_g_clusters_to_obtain.txt


Contact:
--------
camille.marchet@irisa.fr
