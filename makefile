#CC=/usr/bin/g++
CC=g++
#CC=clang++
CFLAGS=  -Wall  -O3 -std=c++11 -march=native -pthread -fopenmp
LDFLAGS= -pthread -fopenmp


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O4   -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O4 -g
LDFLAGS=-g
endif



EXEC=CARNAC

all: $(EXEC)


CARNAC: main.o clustering_cliqueness.o preprocessing.o
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
