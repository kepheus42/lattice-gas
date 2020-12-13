#ifndef TRACER_H
#define TRACER_H

#include <vector>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "site.hpp"

class Lattice;

class Tracer {
public:
Tracer(int,int,Site*);
//
void step(int);
void step_warmup(int);
void step_unhindered(int);
//
void update_last_step(int);
//
// change mobility state of tracer
void stuck();
void unstuck();
// property getters
// identifier
int get_id();
// position
int get_pos();
int get_x();
int get_y();
// displacement
long get_dx();
long get_dy();
// type
int get_type();
// lsq
double get_lsq();
// not/stuck
bool get_isstuck();
//
//unsigned long long
double get_steps_taken();
//
std::vector<int> get_site_correlation();
//
int get_pos_corr(int);
//
std::vector<long> get_correlations();
std::vector<std::vector<Site*> > get_blocking_sites();
//

protected:
// = = = = = = = = = = = = =
int m_id;
int m_type;
Site* m_site;
// unsigned long m_steps_taken;
// int m_dx;
// int m_dy;
// = = = = = = = = = = = = =
std::vector<int> m_powers_of_four;
std::vector<int> m_last_step_dir;
// Stores the last N moves of the tracer, with m_last_step_dir[0] being the most recent one.
// This is used to compute the conditional likelyhoods for all possible step sequences with length up to N.
// std::vector<int> m_last_step_idx;
// This vector of length m_last_step_dir.size() contains the indices of the last 1-N-step sequences
std::vector<long> m_correlations;
// = = = = = = = = = = = = =
//int m_steps_taken;
// int m_steps_taken;
// Total number of steps taken, gives actual step rate when divided by number of timesteps.
// = = = = = = = = = = = = =
bool m_isstuck;
// true if tracer is caught in a trap (TODO: add trapping process to Lattice class), false if not
// = = = = = = = = = = = = =
};

#endif
