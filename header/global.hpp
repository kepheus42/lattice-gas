#ifndef GLOBAL_H
#define GLOBAL_H

#include <random>
#include <cmath>
#include <vector>

extern std::mt19937 generator;

void set_rng_seed(int);
int random_int(int,int);

int coord(int,int,int,int);

std::vector<int> shuffle_vector(std::vector<int>);
#endif
