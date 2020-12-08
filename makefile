IDIR=./header
LDIR=/usr/local/lib
CC=clang++
# -lboost_filesystem -lboost_system
#-fopenmp
CFLAGS=-I$(IDIR) -L$(LDIR) -Xpreprocessor -fopenmp -lm -std=c++17 -Wall -O3
DBFLAGS=-I$(IDIR) -L$(LDIR) -Xpreprocessor -fopenmp -lm -std=c++17 -Wall -Wpedantic -g -DDEBUG

SRC = ./src/tracer.cpp ./src/lattice.cpp ./src/site.cpp ./src/wrapper.cpp

$(info $$SRC is [${SRC}])

debug : ./src/lattice_gas.cpp $(SRC)
	$(CC) $(DBFLAGS) $^ -o lattice_gas -lomp

build : ./src/lattice_gas.cpp $(SRC)
	$(CC) $(CFLAGS) $^ -o lattice_gas -lomp
