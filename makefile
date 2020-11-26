IDIR=./header
LDIR=/usr/local/lib
CC=g++
CFLAGS=-I$(IDIR) -L$(LDIR) -lm -lboost_filesystem -lboost_system -std=c++17 -Wall -Ofast
PFLAGS=-I$(IDIR) -L$(LDIR) -lm -lboost_filesystem -lboost_system -std=c++17 -Wall -Wpedantic -g
DBFLAGS=-I$(IDIR) -L$(LDIR) -lm -lboost_filesystem -lboost_system -std=c++17 -Wall -Wpedantic -g -DDEBUG

SRC = ./src/tracer.cpp ./src/lattice.cpp ./src/site.cpp ./src/wrapper.cpp ./src/global.cpp

$(info $$SRC is [${SRC}])

debug : ./src/lattice_gas.cpp $(SRC)
	$(CC) $(DBFLAGS) $^ -o lattice_gas

build : ./src/lattice_gas.cpp $(SRC)
	$(CC) $(CFLAGS) $^ -o lattice_gas

profile : ./src/lattice_gas.cpp $(SRC)
	$(CC) $(PFLAGS) $^ -o lattice_gas
