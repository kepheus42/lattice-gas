#ifndef TRACER_H
#define TRACER_H

#include <vector>
#include <iostream>
#include "global.hpp"

class Tracer {
public:
// constructor (x,y,grid_size_x,grid_size_y)
Tracer(int,int,int,int,int);
// logic to move the tracer on the 2d lattice
virtual void step(std::vector<int> &);
virtual void unhindered_step();
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
//
protected:
int m_id;
// location and distance covered in lattice units along x- and y-direction
int m_x;
int m_y;
int m_dx;
int m_dy;
// dimensions of the parent grid
int m_grid_size_x;
int m_grid_size_y;
// distance moved squared since birth in lattice units
double m_lsquared;
// can tracer move?
bool m_isstuck;
// size of the tracer in lattice sites
int m_size;
};

class Tracer_2x2 : public Tracer {
public:
Tracer_2x2(int,int,int,int,int);
void step(std::vector<int> &);
void unhindered_step();
};

class Tracer_3x3 : public Tracer {
public:
Tracer_3x3(int,int,int,int,int);
void step(std::vector<int> &);
void unhindered_step();
};

#endif
