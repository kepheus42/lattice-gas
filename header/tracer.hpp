#ifndef TRACER_H
#define TRACER_H

#include <vector>
#include <chrono>
#include <iostream>
#include <algorithm>
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
inline void update_last_step(int);
//
// change mobility state of tracer
void stuck();
void unstuck();
// property getters
// --- id --
int get_id();
//
int get_position();
// --- x ---
int get_x();
int get_dx();
// --- y ---
int get_y();
int get_dy();
// --- size ---
int get_size();
// --- lsquare ---
double get_lsquared();
// --- stuck or not stuck ---
bool get_isstuck();
//
int get_steps_taken();
//
bool get_last_step();
//
std::vector<long> get_correlations();
//

protected:
// = = = = = = = = = = = = =
int m_id;
int m_type;
Site* m_site;
int m_dx;
int m_dy;
double m_lsquared;
// = = = = = = = = = = = = =
std::vector<int> m_last_step_dir;
// Stores the last N moves of the tracer, with m_last_step_dir[0] being the most recent one.
// This is used to compute the conditional likelyhoods for all possible step sequences with length up to N.
// std::vector<int> m_last_step_idx;
// This vector of length m_last_step_dir.size() contains the indices of the last 1-N-step sequences
std::vector<long> m_correlations;
// = = = = = = = = = = = = =
//int m_steps_taken;
int m_steps_taken;
// Total number of steps taken, gives actual step rate when divided by number of timesteps.
// = = = = = = = = = = = = =
bool m_isstuck;
// true if tracer is caught in a trap (TODO: add trapping process to Lattice class), false if not
// = = = = = = = = = = = = =
};

class Tracer_2x2 : public Tracer {
public:
Tracer_2x2(int,int);
};


#endif
