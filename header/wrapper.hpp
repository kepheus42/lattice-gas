#ifndef WRAPPER_H
#define WRAPPER_H

#include "lattice.hpp"
#include <vector>

class Wrapper {
public:
// - - - - - - - - - - - - - - - - - - - - - - - -
// C O N S T R U C T O R
// - - - - - - - - - - - - - - - - - - - - - - - -
Wrapper(int,int,int,int,int,int,int,double,double,double,int);
// - - - - - - - - - - - - - - - - - - - - - - - -
// void init_wtd(int);
void timestep();
//
void update_avg_rate();
void update_avg_lsquared();
// - - - - - - - - - - - - - - - - - - - - - - - -
int get_t();
// = = = = = = = = = = = = = = = = = = = = = = = =
// compute current (at T = m_t) ensemble average step rates
double get_avg_rate_1x1();
double get_avg_rate_2x2();
double get_avg_rate_3x3();
// = = = = = = = = = = = = = = = = = = = = = = = =
// compute current (at T = m_t) ensemble average lsquared
double get_avg_lsquared_1x1();
double get_avg_lsquared_2x2();
double get_avg_lsquared_3x3();
// = = = = = = = = = = = = = = = = = = = = = = = =
// get stepping rate distribution for 1x1/2x2/3x3 tracers
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> get_stepping_rate_distribution_1x1();
std::vector<double> get_stepping_rate_distribution_2x2();
std::vector<double> get_stepping_rate_distribution_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
// get avg rate correlations for 1x1/2x2/3x3 tracers
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> get_avg_rate_correlations_1x1();
std::vector<double> get_avg_rate_correlations_2x2();
std::vector<double> get_avg_rate_correlations_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
// get wtds for 1x1 tracers
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> get_wtd_1x1();
std::vector<int> get_wtd_2x2();
std::vector<int> get_wtd_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
// get normalized WTDs
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> get_norm_wtd_1x1();
std::vector<double> get_norm_wtd_2x2();
std::vector<double> get_norm_wtd_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
// T W O  A N D
// T H R E E  S T E P
// C O R R E L A T I O N S
// - - - - - - - - - - - - - - - - - - - - - - - - -
// average two step correlations
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> get_avg_two_step_correlation_1x1();
std::vector<double> get_avg_two_step_correlation_2x2();
std::vector<double> get_avg_two_step_correlation_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
// returns relative two step correlations for all particles
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> get_relative_two_step_correlations_1x1();
std::vector<double> get_relative_two_step_correlations_2x2();
std::vector<double> get_relative_two_step_correlations_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
// average three step correlations
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> get_avg_three_step_correlation_1x1();
std::vector<double> get_avg_three_step_correlation_2x2();
std::vector<double> get_avg_three_step_correlation_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
// returns relative three step correlations for all particles
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> get_relative_three_step_correlations_1x1();
std::vector<double> get_relative_three_step_correlations_2x2();
std::vector<double> get_relative_three_step_correlations_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
// average four step correlations
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> get_avg_four_step_correlation_1x1();
std::vector<double> get_avg_four_step_correlation_2x2();
std::vector<double> get_avg_four_step_correlation_3x3();
// - - - - - - - - - - - - - - - - - - - - - - - -
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

private:
std::vector<Lattice*> m_lattices;
int m_number_of_lattices;
int m_grid_size_x;
int m_grid_size_y;
int m_number_of_timesteps;
int m_number_of_timesteps_warmup;
int m_t;
int m_number_of_tracers_1x1;
int m_number_of_tracers_2x2;
int m_number_of_tracers_3x3;
double m_step_rate_1x1;
double m_step_rate_2x2;
double m_step_rate_3x3;
int m_wtd_max;
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
// time- and ensemble averaged conditional probabilities for series of steps with length 2/3/4 for 1x1/2x2/3x3 tracers
// 1x1
std::vector<double> m_avg_corr_2s_1x1;
std::vector<double> m_avg_corr_3s_1x1;
std::vector<double> m_avg_corr_4s_1x1;
// 2x2
std::vector<double> m_avg_corr_2s_2x2;
std::vector<double> m_avg_corr_3s_2x2;
std::vector<double> m_avg_corr_4s_2x2;
// 3x3
std::vector<double> m_avg_corr_2s_3x3;
std::vector<double> m_avg_corr_3s_3x3;
std::vector<double> m_avg_corr_4s_3x3;
// - - - - - - - - - - - - - - - - - - - - - - - -
};

#endif
