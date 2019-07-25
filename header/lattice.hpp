#ifndef LATTICE_H
#define LATTICE_H

#include "tracer.hpp"
#include <vector>
#include <random>

class Lattice {

public:
// constructor
Lattice(int,int,int,int);
// logic
// computes the 1d linear coordinate from x,y tuple
int coord(int,int);
// performs a single timestep
void timestep();
void timestep_no_interaction();
// getters
int get_t();
int get_grid_size_x();
int get_grid_size_y();
int get_number_of_tracers_1x1();
int get_number_of_tracers_2x2();
std::vector<int> get_tracer_positions();
std::vector<int> get_occupation_map();

double get_avg_lsquared_1x1();
double get_avg_lsquared_2x2();
// Debugging
void print_occupation_map();
void print_tracer_positions();

private:
int m_t;
int m_grid_size_x;
int m_grid_size_y;
int m_number_of_tracers_1x1;
int m_number_of_tracers_2x2;
std::vector<int> m_movement_order;
std::vector<int> m_occupation_map;
std::vector<Tracer*> m_tracers;

};

#endif
