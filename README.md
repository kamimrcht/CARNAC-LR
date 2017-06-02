Compilation and Usage:
make;
./clustering_cliqueness -f infile [-i 10, -c val, -p, -o path]

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

