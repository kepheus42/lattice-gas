#ifndef GLOBAL_H
#define GLOBAL_H

#include <random>
#include <cmath>
#include <vector>
#include <algorithm>
#include <ctime>
#include <iostream>


extern std::mt19937 generator;

inline void set_rng_seed(int seed){
        if(seed == 0) {
                int time = std::time(NULL);
                generator.seed(time);
        } else {
                generator.seed(seed);
        }
}

inline double random_0_to_1(){
        std::uniform_real_distribution<> rnd_0_to_1(0,1);
        return rnd_0_to_1(generator);
}

inline int random_int(int min, int max){
        std::uniform_int_distribution<> rnd_int(min, max);
        return rnd_int(generator);
}

inline int random_int_binomial(int max, double p){
        std::binomial_distribution<> rnd_int(max,p);
        return rnd_int(generator);
}

inline void shuffle_vector(std::vector<int> vec){
        std::shuffle(vec.begin(),vec.end(),generator);
}

int coord(int,int,int,int);

// std::vector<int> shuffle_vector(std::vector<int>);
void shuffle_vector(std::vector<int>);

// for debugging
template <typename T>
static void print_vector(std::vector<T> vec)
{
        for(T el : vec)
        {
                std::cout << el << std::endl;
        }
}

#endif
