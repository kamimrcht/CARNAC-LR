Compilation and Usage:
 g++ -std=c++11 clustering_cliqueness2.cpp -o clustering_cliqueness;  ./clustering_cliqueness short_read_connector_res.txt

Results are written in:
 final_g_clusters.txt

Test:
./clustering_cliqueness test/short_read_connector_res.txt
The expected result is the same than in test/final_g_clusters_to_obtain.txt
