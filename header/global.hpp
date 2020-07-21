#ifndef GLOBAL_H
#define GLOBAL_H

#include <random>
#include <cmath>
#include <vector>
#include <algorithm>

extern std::mt19937 generator;

void set_rng_seed(int);
int random_int(int,int);
int random_int_binomial(int,double);

double random_0_to_1();

int coord(int,int,int,int);

// std::vector<int> shuffle_vector(std::vector<int>);
void shuffle_vector(std::vector<int>);

#endif
