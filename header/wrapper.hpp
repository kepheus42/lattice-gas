#ifndef WRAPPER_H
#define WRAPPER_H

#include "lattice.hpp"
#include "global.hpp"
#include <vector>
#include <numeric>

class Wrapper {
public:

Wrapper(int,int,int,int,int,int,double,double,int);

inline int coord(int,int);

void timestep();
void timestep_warmup();

// wrapper function to call all the updating routines
inline void update_data();

int get_t();

// compute current ensemble average step rates
double get_avg_rate_1x1();
double get_avg_rate_2x2();

// compute current ensemble average lsquared
double get_avg_lsquared_1x1();
double get_avg_lsquared_2x2();

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

std::vector<double> get_result_lsquared_1x1();
std::vector<double> get_result_lsquared_2x2();

std::vector<unsigned long> get_result_correlations_1x1();
std::vector<unsigned long> get_result_correlations_2x2();

std::vector<double> get_result_norm_correlations_1x1();
std::vector<double> get_result_norm_correlations_2x2();

std::vector<int> get_result_positions_1x1();
std::vector<int> get_result_positions_2x2();

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
int m_number_of_timesteps;
int m_t;
int m_number_of_tracers_1x1;
int m_number_of_tracers_2x2;
double m_step_rate_1x1;
double m_step_rate_2x2;

// handles data storage
int m_data_points;
std::vector<int> m_data_storage_intervals;
int m_next_data_storage_interval;

int m_number_of_pos_to_save;
int m_pos_saving_interval;

int m_number_of_tracers_total_times_number_of_lattices;
int m_number_of_tracers_1x1_times_number_of_lattices;
int m_number_of_tracers_2x2_times_number_of_lattices;

// measured average stepping rates of the tracers
std::vector<double> m_avg_rate_1x1;
std::vector<double> m_avg_rate_2x2;

// measured average mean squared displacement of the tracers
std::vector<double> m_avg_lsquared_1x1;
std::vector<double> m_avg_lsquared_2x2;

// for counting series of steps with length 2/3/4 for 1x1/2x2/3x3 tracers, sorted by directions
std::vector<unsigned long> m_correlations_1x1;
std::vector<unsigned long> m_correlations_2x2;

//
std::vector<int> m_positions_1x1;
std::vector<int> m_positions_2x2;

};

#endif
