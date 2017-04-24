Compilation and Usage:
 g++ -std=c++11 clustering_cliqueness2.cpp -o clustering_cliqueness;  ./clustering_cliqueness short_read_connector_res.txt

Results are written in:
 final_g_clusters.txt

Test:
./clustering_cliqueness test/short_read_connector_res.txt
The expected result is the same than in test/final_g_clusters_to_obtain.txt

The executable also provides a file for graph visualization (nodes_metrics.txt).
Small graphs clustering can be seen (with removed edges during the clustering in red); by doing:
cd scripts
python print_graphs_cliqueness.py ../test/short_read_connector_res.txt ../test/two_clusters.fasta ../test/clusters.list ../nodes_metrics.txt ../final_g_clusters.txt
A .png file will be generated in the scripts directory.

