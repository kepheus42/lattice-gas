#ifndef WRAPPER_H
#define WRAPPER_H

#include "lattice.hpp"
#include <vector>
#include <numeric>

class Wrapper {
public:
// - - - - - - - - - - - - - - - - - - - - - - - -
// C O N S T R U C T O R
// - - - - - - - - - - - - - - - - - - - - - - - -
Wrapper(int,int,int,int,int,int,int,double,double,double,int,int);
// - - - - - - - - - - - - - - - - - - - - - - - -
// void init_wtd(int);
void timestep();
//
// wrapper function to call all the updating routines
inline void update_data();
// - - - - - - - - - - - - - - - - - - - - - - - -
unsigned int get_t();
// = = = = = = = = = = = = = = = = = = = = = = = =
// compute current ensemble average step rates
double get_avg_rate_1x1();
double get_avg_rate_2x2();
double get_avg_rate_3x3();
// = = = = = = = = = = = = = = = = = = = = = = = =
// compute current ensemble average lsquared
double get_avg_lsquared_1x1();
double get_avg_lsquared_2x2();
double get_avg_lsquared_3x3();
// = = = = = = = = = = = = = = = = = = = = = = = =
// compute updated correlations and 2-step wtds, by summing up correlations contributions from all lattices
void update_correlations_1x1();
void update_correlations_2x2();
void update_correlations_3x3();
// = = = = = = = = = = = = = = = = = = = = = = = =
// New layout:
// Getters for sim results:
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> get_result_rate_1x1();
std::vector<double> get_result_rate_2x2();
std::vector<double> get_result_rate_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> get_result_lsquared_1x1();
std::vector<double> get_result_lsquared_2x2();
std::vector<double> get_result_lsquared_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<unsigned long> get_result_wtd_1x1();
std::vector<unsigned long> get_result_wtd_2x2();
std::vector<unsigned long> get_result_wtd_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<unsigned long> get_result_correlations_1x1();
std::vector<unsigned long> get_result_correlations_2x2();
std::vector<unsigned long> get_result_correlations_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
private:
std::vector<Lattice*> m_lattices;
unsigned int m_number_of_lattices;
unsigned int m_grid_size_x;
unsigned int m_grid_size_y;
unsigned int m_number_of_timesteps;
unsigned int m_number_of_timesteps_warmup;
unsigned int m_t;
unsigned int m_number_of_tracers_1x1;
unsigned int m_number_of_tracers_2x2;
unsigned int m_number_of_tracers_3x3;
double m_step_rate_1x1;
double m_step_rate_2x2;
double m_step_rate_3x3;
unsigned int m_wtd_max;
unsigned int m_wtd_res;
// - - - - - - - - - - - - - - - - - - - - - - - -
double m_this_number_of_tracers_1x1_times_number_of_lattices;
double m_this_number_of_tracers_2x2_times_number_of_lattices;
double m_this_number_of_tracers_3x3_times_number_of_lattices;
// - - - - - - - - - - - - - - - - - - - - - - - -
// measured average stepping rates of the tracers
std::vector<double> m_avg_rate_1x1;
std::vector<double> m_avg_rate_2x2;
std::vector<double> m_avg_rate_3x3;
// - - - - - - - - - - - - - - - - - - - - - - - -
// measured average mean squared displacement of the tracers
std::vector<double> m_avg_lsquared_1x1;
std::vector<double> m_avg_lsquared_2x2;
std::vector<double> m_avg_lsquared_3x3;
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<unsigned long> m_wtd_1x1;
std::vector<unsigned long> m_wtd_2x2;
std::vector<unsigned long> m_wtd_3x3;
// - - - - - - - - - - - - - - - - - - - - - - - -
// for counting series of steps with length 2/3/4 for 1x1/2x2/3x3 tracers, sorted by directions
std::vector<unsigned long> m_correlations_1x1;
std::vector<unsigned long> m_correlations_2x2;
std::vector<unsigned long> m_correlations_3x3;
// - - - - - - - - - - - - - - - - - - - - - - - -
};

#endif
