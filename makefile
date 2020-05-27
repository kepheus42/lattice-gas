IDIR=./header
LDIR=/usr/local/lib
CC=g++
CFLAGS=-I$(IDIR) -L$(LDIR) -lm -lboost_filesystem -lboost_system -std=c++17 -Wall -Ofast
DBFLAGS=-I$(IDIR) -L$(LDIR) -lm -lboost_filesystem -lboost_system -std=c++17 -Wall -Wpedantic -g -DDEBUG
LLDBFLAGS=-I$(IDIR) -L$(LDIR) -lm -lboost_filesystem -lboost_system -std=c++17 -Wall -Wpedantic -g

SRC = ./src/global.cpp ./src/tracer.cpp ./src/lattice.cpp ./src/output.cpp ./src/wrapper.cpp

$(info $$SRC is [${SRC}])

debug : ./src/lattice_gas.cpp $(SRC)
	$(CC) $(DBFLAGS) $^ -o lattice_gas

lldebug : ./src/lattice_gas.cpp $(SRC)
	$(CC) $(LLDBFLAGS) $^ -o lattice_gas

build : ./src/lattice_gas.cpp $(SRC)
	$(CC) $(CFLAGS) $^ -o lattice_gas
