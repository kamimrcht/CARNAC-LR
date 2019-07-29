#CC=/usr/bin/g++
CC=g++
#CC=clang++
CFLAGS=  -Wall  -O3 -std=c++11 -march=native -pthread -fopenmp
CFLAGS_SIMPLE=-std=c++11

LDFLAGS= -pthread -fopenmp


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O4   -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif



EXEC=CARNAC-LR scripts/CARNAC_to_fasta

all: $(EXEC)

scripts/CARNAC_to_fasta: scripts/CARNAC_to_fasta.cpp
	$(CC) -o $@ -c $^ $(CFLAGS_SIMPLE)

CARNAC-LR: main.o clustering_cliqueness.o preprocessing.o
	$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.cpp clustering_cliqueness.hpp
	$(CC) -o $@ -c $< $(CFLAGS)

clustering_cliqueness.o: clustering_cliqueness.cpp clustering_cliqueness.hpp findArticulationPoints.hpp
	$(CC) -o $@ -c $< $(CFLAGS)

preprocessing.o: findArticulationPoints.cpp findArticulationPoints.hpp
	$(CC) -o $@ -c $< $(CFLAGS)




clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
