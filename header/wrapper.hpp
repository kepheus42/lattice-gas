#ifndef WRAPPER_H
#define WRAPPER_H

#include "lattice.hpp"
#include "global.hpp"
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>

class Wrapper {
public:

Wrapper(int,int,int,int,int,int,double,double);

inline int coord(int,int);

void warmup();
void evolve();
void evolve_no_interaction();

void timestep();
void timestep_warmup();
void timestep_no_interaction();

// wrapper function to call all the updating routines
void update_data();

int get_t();

// compute current ensemble average step rates
double get_avg_rate_1x1();
double get_avg_rate_2x2();

// compute current ensemble average lsq
double get_avg_lsq_1x1();
double get_avg_lsq_2x2();

// compute updated correlations
void update_correlations_1x1();
void update_correlations_2x2();

void update_positions_1x1();
void update_positions_2x2();

// New layout:
// Getters for sim results:
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> get_result_rate_1x1();
std::vector<double> get_result_rate_2x2();

std::vector<double> get_result_lsq_1x1();
std::vector<double> get_result_lsq_2x2();

std::vector<double> get_result_corr_1x1();
std::vector<double> get_result_corr_2x2();

std::vector<double> get_result_pos_1x1();
std::vector<double> get_result_pos_2x2();

// for debugging
std::vector<int> get_raw_pos_1x1();
std::vector<int> get_raw_pos_2x2();
std::vector<unsigned long long> get_raw_corr_1x1();
std::vector<unsigned long long> get_raw_corr_2x2();


private:
// - - - - - - - - - - - - - - - - - - - - - - - -
// stores pointers to all lattices
std::vector<Lattice*> m_lattices;

// stores pointers to all tracers (all & by type)
std::vector<Tracer*> m_tracers;
std::vector<Tracer*> m_tracers_1x1;
std::vector<Tracer*> m_tracers_2x2;

int m_number_of_lattices;
int m_grid_size;
int m_timesteps;
int m_timesteps_w;
int m_t;
int m_w;
int m_number_of_tracers_1x1;
int m_number_of_tracers_2x2;
double m_step_rate_1x1;
double m_step_rate_2x2;
int m_dpoints;
/*
   // handles data storage
   int m_data_points;
   std::vector<int> m_data_storage_intervals;
   int m_next_data_storage_interval;
   std::vector<int> m_steps_taken_divisor;

   int m_number_of_tracers_total_times_number_of_lattices;
   int m_number_of_tracers_1x1_times_number_of_lattices;
   int m_number_of_tracers_2x2_times_number_of_lattices;

   // measured average stepping rates of the tracers
   std::vector<double> m_avg_rate_1x1;
   std::vector<double> m_avg_rate_2x2;

   // measured average mean squared displacement of the tracers
   std::vector<double> m_avg_lsq_1x1;
   std::vector<double> m_avg_lsq_2x2;

   // for counting series of steps with length 2/3/4 for 1x1/2x2/3x3 tracers, sorted by directions
   std::vector<unsigned long> m_correlations_1x1;
   std::vector<unsigned long> m_correlations_2x2;

   //
   std::vector<int> m_positions_1x1;
   std::vector<int> m_positions_2x2;
 */
};

#endif
