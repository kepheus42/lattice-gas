IDIR = ./header
LDIR = /usr/local/lib
CC   = clang++

SRC = ./src/tracer.cpp ./src/lattice.cpp ./src/site.cpp ./src/wrapper.cpp

CFLAGS  = -I$(IDIR) -L$(LDIR) -Xpreprocessor -fopenmp -lomp -std=c++17 -Wall -O3
DBFLAGS = -I$(IDIR) -L$(LDIR) -Xpreprocessor -fopenmp -lomp -std=c++17 -Wall -Wpedantic -g -DDEBUG

debug : ./src/lattice_gas.cpp $(SRC)
	$(CC) $(DBFLAGS) $^ -o lattice_gas

build : ./src/lattice_gas.cpp $(SRC)
	$(CC) $(CFLAGS) $^ -o lattice_gas
