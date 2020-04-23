IDIR=./header
LDIR=/usr/local/lib
CC=g++
CFLAGS=-I$(IDIR) -L$(LDIR) -lm -lboost_filesystem -lboost_system -std=c++17 -Wall -Ofast

SRC = ./src/global.cpp ./src/tracer.cpp ./src/lattice.cpp ./src/output.cpp ./src/wrapper.cpp

$(info $$SRC is [${SRC}])

lattice_gas : ./src/lattice_gas.cpp $(SRC)
	$(CC) $(CFLAGS) $^ -o $@
