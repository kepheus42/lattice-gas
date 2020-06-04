#ifndef LATTICE_H
#define LATTICE_H

#include "tracer.hpp"
#include "sites.hpp"
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>

class Lattice {
public:
// constructor
// three species (1x1,2x2), numbers of tracers, step rates
Lattice(int,int,int,int,int,int,double,double,double,int,int);
// function to create the specified numbers of tracers
// m_number_of_tracers_1x1, m_number_of_tracers_2x2
// with the corresponding step rates
void setup_tracers();
// function to initialize the tracer ID lists for random step order selection
void setup_movement_selection_list();
//
void setup_sites();
//
void timestep();
void timestep_warmup();
void timestep_no_interaction();
//
// logic
// computes the 1d linear coordinate from x,y tuple
inline int coord(int,int);
// performs a single timestep
// warm up lattice before actual simulation
void warm_up();
// - - - - - - - - - - - - - - -
// getters
// - - - - - - - - - - - - - - -
int get_t();
int get_current_step_attempt();
int get_grid_size_x();
int get_grid_size_y();
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
std::vector<Site_1x1 *> get_sites_1x1();
std::vector<Site_2x2 *> get_sites_2x2();

std::vector<Site_1x1 *> get_neighbor_sites_1x1(Site_1x1 *);
std::vector<Site_2x2 *> get_neighbor_sites_2x2(Site_2x2 *);

std::vector<Site_1x1 *> get_blocking_sites_1x1_1x1(Site_1x1 *);
std::vector<Site_2x2 *> get_blocking_sites_2x2_1x1(Site_1x1 *);
std::vector<Site_1x1 *> get_blocking_sites_1x1_2x2(Site_2x2 *);
std::vector<Site_2x2 *> get_blocking_sites_2x2_2x2(Site_2x2 *);
// - - - - - - - - - - - - - - - - - - - - - - - - -
unsigned long get_sum_steps_taken_1x1();
unsigned long get_sum_steps_taken_2x2();
// - - - - - - - - - - - - - - - - - - - - - - - - -
double get_sum_lsqaured_1x1();
double get_sum_lsquared_2x2();
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> get_tracer_placement_1x1();
std::vector<int> get_tracer_placement_2x2();
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> get_tracer_positions();
std::vector<int> get_occupation_map();
// - - - - - - - - - - - - - - - - - - - - - - - - -
// Debugging
// - - - - - - - - - - - - - - - - - - - - - - - - -
void print_occupation_map();
void print_tracer_positions();
// - - - - - - - - - - - - - - - - - - - - - - - - -
private:
// parameters of the simulation
int m_t;
int m_grid_size_x;
int m_grid_size_y;
int m_number_of_sites;
int m_number_of_timesteps;
int m_number_of_tracers_1x1;
int m_number_of_tracers_2x2;
int m_number_of_tracers_total;
// variable base step rates for all tracer types (between 0.0 and 1.0)
double m_step_rate_1x1;
double m_step_rate_2x2;
int m_wtd_max;
int m_wtd_res;
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> m_movement_selector;
int m_movement_selector_length;
int m_step_attempts_per_timestep;
int m_current_step_attempt;
// - - - - - - - - - - - - - - - - - - - - - - - - -
double m_t_increment;
// - - - - - - - - - - - - - - - - - - - - - - - - -
// to keep track of the relations between sites on the lattice
std::vector<Site*> m_sites_1x1;
std::vector<Site*> m_sites_2x2;
// to keep track of the positions of all tracers
std::vector<int> m_occupation_map;
// to store the tracer objects
std::vector<Tracer*> m_tracers;
std::vector<Tracer*> m_tracers_1x1;
std::vector<Tracer*> m_tracers_2x2;
};

#endif
