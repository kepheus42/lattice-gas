#ifndef LATTICE_H
#define LATTICE_H

#include "tracer.hpp"
#include "site.hpp"
#include "global.hpp"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>

class Tracer;
class Tracer_1x1;
class Tracer_2x2;

class Site;
class Site_1x1;
class Site_2x2;

class Lattice {
public:
Lattice(int,int,int,double,double);

void setup_sites();
void setup_tracers();
void setup_movement_selection_list();

void timestep();
void timestep_warmup();
void timestep_no_interaction();

inline int coord(int,int);
inline int coord_1x1(int,int);
inline int coord_2x2(int,int);

inline int coord_to_x(int);
inline int coord_to_y(int);
// - - - - - - - - - - - - - - -
// getters
// - - - - - - - - - - - - - - -
int get_grid_size();
int get_number_of_tracers_1x1();
int get_number_of_tracers_2x2();
// - - - - - - - - - - - - - - -
int get_tracer_size(int);
int get_tracer_last_move(int);
int get_tracer_position_x(int);
int get_tracer_position_y(int);
// - - - - - - - - - - - - - - -
std::vector<Tracer *> get_tracers();
std::vector<Tracer *> get_tracers_1x1();
std::vector<Tracer *> get_tracers_2x2();
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> get_tracer_positions();

void set_neighbor_sites(Site_1x1 *);
void set_neighbor_sites(Site_2x2 *);

// - - - - - - - - - - - - - - - - - - - - - - - - -
// Debugging
// - - - - - - - - - - - - - - - - - - - - - - - - -
void print_tracer_positions();
void print_neighbors();
void print_sites();
// - - - - - - - - - - - - - - - - - - - - - - - - -
private:
// parameters of the simulation
int m_grid_size;
int m_number_of_sites;
int m_number_of_tracers;
int m_number_of_tracers_1x1;
int m_number_of_tracers_2x2;
// variable base step rates for all tracer types (between 0.0 and 1.0)
double m_step_rate_1x1;
double m_step_rate_2x2;
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> m_movement_selector;
int m_movement_selector_length;
int m_step_attempts_per_timestep;
// - - - - - - - - - - - - - - - - - - - - - - - - -
// to keep track of the lattice geometry
std::vector<Site*> m_sites_1x1;
std::vector<Site*> m_sites_2x2;
// to store the tracer objects
std::vector<Tracer*> m_tracers;
std::vector<Tracer*> m_tracers_1x1;
std::vector<Tracer*> m_tracers_2x2;
};

#endif
