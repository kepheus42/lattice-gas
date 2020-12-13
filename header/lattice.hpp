#ifndef LATTICE_H
#define LATTICE_H

#include "tracer.hpp"
#include "site.hpp"
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

class Lattice {
public:
Lattice(int,int,int,int,int,double,double,unsigned int);
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void setup_sites();
void setup_tracers();
void setup_movement_selection_list();
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void warmup();
void evolve();
void evolve_no_interaction();
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void timestep();
void timestep_warmup();
void timestep_no_interaction();
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void update_data();
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void set_rng_seed(int);
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
inline int coord(int,int);
inline int coord_to_x(int);
inline int coord_to_y(int);
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
// getters
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
int get_grid_size();
int get_number_of_tracers_1x1();
int get_number_of_tracers_2x2();
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<Tracer *> get_tracers();
std::vector<Tracer *> get_tracers_1x1();
std::vector<Tracer *> get_tracers_2x2();
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<double> get_avg_rate_1x1();
std::vector<double> get_avg_rate_2x2();

std::vector<double> get_avg_diff_1x1();
std::vector<double> get_avg_diff_2x2();

std::vector<double> get_avg_lsq_1x1();
std::vector<double> get_avg_lsq_2x2();
std::vector<double> get_avg_corr_1x1();
std::vector<double> get_avg_corr_2x2();
std::vector<double> get_site_corr_1x1();
std::vector<double> get_site_corr_2x2();
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<int> get_pos_1x1();
std::vector<int> get_pos_2x2();
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
// Debugging
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void print_tracers();
void print_tracer_positions();
void print_site_1x1_states();
void print_site_2x2_states();
void print_sites();
template <typename T>
void db_print_vector(std::vector<T>);
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
private:
// parameters of the simulation
int m_id;
int m_t;
int m_w;
int m_timesteps;
int m_timesteps_w;
int m_grid_size;
int m_number_of_sites;
int m_number_of_tracers;
int m_number_of_tracers_1x1;
int m_number_of_tracers_2x2;
// variable base step rates for all tracer types (between 0.0 and 1.0)
double m_step_rate_1x1;
double m_step_rate_2x2;
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
int m_step_attempts_per_timestep;
int m_movement_selector_length;
std::vector<int> m_movement_selector;
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
// handles data storage
int m_dpoints;
int m_dinterval;
// randomness
std::mt19937 m_rng; // generates the random numbers for stepping
std::uniform_int_distribution<int> m_random_par; // for picking a random particle
std::uniform_int_distribution<int> m_random_dir; // for picking a random direction
// pre computed values
double m_one_over_n_1x1;
double m_one_over_four_n_1x1;
double m_one_over_n_2x2;
double m_one_over_four_n_2x2;
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
// to keep track of the lattice geometry
// std::vector<Site*> m_sites;
std::vector<Site*> m_sites_1x1;
std::vector<Site*> m_sites_2x2;
// to store the tracer objects
std::vector<Tracer*> m_tracers;
std::vector<Tracer*> m_tracers_1x1;
std::vector<Tracer*> m_tracers_2x2;
// to store data
std::vector<double> m_avg_rate_1x1;
std::vector<double> m_avg_rate_2x2;
std::vector<double> m_avg_diff_1x1;
std::vector<double> m_avg_diff_2x2;
std::vector<double> m_avg_lsq_1x1;
std::vector<double> m_avg_lsq_2x2;
std::vector<double> m_site_corr_1x1;
std::vector<double> m_site_corr_2x2;
std::vector<int> m_pos_1x1;
std::vector<int> m_pos_2x2;

std::vector<int> m_site_corr_1x1_counter;
std::vector<int> m_site_corr_2x2_counter;
};

#endif
