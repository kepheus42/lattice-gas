IDIR=./header
LDIR=/usr/local/lib
CC=g++
CFLAGS=-I$(IDIR) -L$(LDIR) -lm -lboost_filesystem -lboost_system -std=c++11 -Wall -Ofast

SRC = ./src/global.cpp ./src/tracer.cpp ./src/lattice.cpp ./src/output.cpp ./src/wrapper.cpp


$(info $$SRC is [${SRC}])

lattice_gas : ./src/lattice_gas.cpp $(SRC)
	$(CC) $(CFLAGS) $^ -o $@

tracer_test : ./src/tracer_test.cpp $(SRC)
	$(CC) $(CFLAGS) $^ -o $@

rw_stats : ./src/random_walk_stats.cpp $(SRC)
	$(CC) $(CFLAGS) $^ -o $@

lattice_gas_trajectories : ./src/lattice_gas_trajectories.cpp $(SRC)
	$(CC) $(CFLAGS) $^ -o $@

lattice_gas_v2 : ./src/lattice_gas_wtds.cpp $(SRC)
	$(CC) $(CFLAGS) $^ -o $@

lattice_gas_multisim : ./src/lattice_gas_multisim.cpp $(SRC)
	$(CC) $(CFLAGS) $^ -o $@

lattice_gas_graphical_output : ./src/lattice_gas_graphical.cpp $(SRC)
	$(CC) $(CFLAGS) $^ -o $@
#
#lg_traps : ./src/lattice_gas_with_traps.cpp $(SRC)
#	$(CC) $(CFLAGS) $^ -o $@
#
#cr_rw: ./src/crowded_random_walk.cpp $(SRC)
#	$(CC) $(CFLAGS) $^ -o $@
