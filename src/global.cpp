#include "global.hpp"
#include <ctime>

std::mt19937 generator;

void set_rng_seed(int seed){
        if(seed == 0) {
                int time = std::time(NULL);
                generator.seed(time);
        } else {
                generator.seed(seed);
        }
}

double random_0_to_1(){
        std::uniform_real_distribution<> rnd_0_to_1(0,1);
        return rnd_0_to_1(generator);
}

int random_int(int min, int max){
        std::uniform_int_distribution<> rnd_int(min, max);
        return rnd_int(generator);
}

int random_int_binomial(int max, double p){
        std::binomial_distribution<> rnd_int(max,p);
        return rnd_int(generator);
}

void shuffle_vector(std::vector<int> vec){
        std::shuffle(vec.begin(),vec.end(),generator);
}

int coord(int x, int y, int grid_size_x, int grid_size_y){
        return ((x+grid_size_x)%grid_size_x)*grid_size_y+((y+grid_size_y)%grid_size_y);
}
