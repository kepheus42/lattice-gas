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
virtual void step(std::vector<int> &, int, double);
virtual void unhindered_step();
//
void update_last_moves(int);
//
void update_two_step_correlation();
void update_three_step_correlation();
void update_four_step_correlation();
//
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
int get_time_since_last_move();
// --- total number of steps taken
int get_steps_taken();
// --- 1 if step taken during current timestep, 0 else
int get_last_step();
// --- stats for the trajectory of the random walker (waiting time distributions for 4 possible movements w.r.t. last move)
std::vector<int> get_wtd();
//
std::vector<double> get_relative_two_step_correlation();
std::vector<double> get_relative_three_step_correlation();
std::vector<double> get_relative_four_step_correlation();
//
protected:
int m_id;
// location and distance covered in lattice units along x- and y-direction
int m_x;
int m_y;
int m_dx;
int m_dy;
// base step rate (0.0 to 1.0)
double m_step_rate;
// last move of the tracers
int m_last_step;
std::vector<int> m_last_step_dir;
int m_current_move;
int m_wtd_index;
// time when the last move took place, and time elapsed since then
double m_time_of_last_move;
double m_time_since_last_move;

int m_steps_taken;
// dimensions of the parent grid
int m_grid_size_x;
int m_grid_size_y;
// distance moved squared since birth in lattice units
double m_lsquared;
// can tracer move?
bool m_isstuck;
// size of the tracer in lattice sites
int m_size;
// waiting times for different move directions w.r.t. the last move
// stored in one vector, as:
// 0 : +x+x
// 1 : +x-x
// 2 : +x+y
// 3 : +x-y
// next 4 with one delta_t waiting time
// 4 : +x+0+x
// 5 : +x+0-x
// 6 : +x+0+y
// 7 : +x+0-y
// etc.
// TODO: move storage to <lattice> or <wrapper>
std::vector<int> m_wtd;
std::vector<int> m_two_step_correlation;
std::vector<int> m_three_step_correlation;
std::vector<int> m_four_step_correlation;
};

class Tracer_2x2 : public Tracer {
public:
Tracer_2x2(int,int,int,int,int,double);
// direction selection in tracer->step function
void step(std::vector<int> &);
// direction selection in tracer->step function, m_t as input
void step(std::vector<int> &, int);
// direction selection in lattice->timestep function, m_t as input
void step(std::vector<int> &, int, int);
void unhindered_step();
};

class Tracer_3x3 : public Tracer {
public:
Tracer_3x3(int,int,int,int,int,double);
void step(std::vector<int> &);
void step(std::vector<int> &, int);
void step(std::vector<int> &, int, int);
void unhindered_step();
};

#endif
