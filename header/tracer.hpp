#ifndef TRACER_H
#define TRACER_H

#include <vector>
#include <chrono>
#include <iostream>
#include <algorithm>
#include "global.hpp"

class Tracer {
public:
// constructor (x,y,grid_size_x,grid_size_y)
Tracer(int,int,int,int,int,double);
// logic to move the tracer on the 2d lattice
// different versions of these functions are required for the larger square tracers (2x2 and 3x3), therefore virtual
virtual void step(std::vector<int> &, int, double);
virtual void step_warmup(std::vector<int> &, int);
// this one is the same for all types, therefore not virtual
void step_unhindered(int, double);
//
inline void update_last_step(int);
//
inline int coord(int,int);

void reset();
// change mobility state of tracer
void stuck();
void unstuck();
// property getters
// --- id --
int get_id();
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
// --- last move (0,1,2,3,4) of the tracer ---
int get_last_move();
// --- time since last move
double get_time_since_last_move();
// --- total number of steps taken
int get_steps_taken();
// --- 1 if step taken during current timestep, 0 else
bool get_last_step();
//
protected:
// = = = = = = = = = = = = =
int m_id;
// Integer assigned at creation
int m_size;
// size of the tracer in lattice sites
// = = = = = = = = = = = = =
int m_x;
int m_y;
int m_dx;
int m_dy;
// location and distance covered in lattice units along x- and y-direction
// = = = = = = = = = = = = =
double m_step_rate;
// Base step rate (0.0 to 1.0): how many steps does the tracer take per unit time, on average.
// = = = = = = = = = = = = =
bool m_last_step;
// 1 if the tracers most recent attempt to take a step was valid, 0 otherwise.
// = = = = = = = = = = = = =
std::vector<int> m_last_step_dir;
// Stores the last N moves of the tracer, with m_last_step_dir[0] being the most recent one.
// This is used to compute the conditional likelyhoods for all possible step sequences with length up to N.
std::vector<int> m_last_step_idx;
// This vector of length m_last_step_dir.size() contains the indices of the last 1-N-step sequences
// = = = = = = = = = = = = =
int m_last_step_wtd_idx;
int m_wtd_max;
int m_wtd_res;
// = = = = = = = = = = = = =
double m_time_of_last_step;
double m_time_since_last_step;
// time when the last move took place, and time elapsed since then
// = = = = = = = = = = = = =
int m_steps_taken;
// Total number of steps taken, gives actual step rate when divided by number of timesteps.
// = = = = = = = = = = = = =
int m_grid_size_x;
int m_grid_size_y;
// dimensions of the parent grid
// = = = = = = = = = = = = =
double m_lsquared;
// distance moved squared since birth in lattice units
// = = = = = = = = = = = = =
bool m_isstuck;
// true if tracer is caught in a trap (TODO: add trapping process to Lattice class), false if not
// = = = = = = = = = = = = =
};

class Tracer_2x2 : public Tracer {
public:
Tracer_2x2(int,int,int,int,int,double);
void step(std::vector<int> &, int, double);
void step_warmup(std::vector<int> &, int);
};

class Tracer_3x3 : public Tracer {
public:
Tracer_3x3(int,int,int,int,int,double);
void step(std::vector<int> &, int, double);
void step_warmup(std::vector<int> &, int);
};

#endif
