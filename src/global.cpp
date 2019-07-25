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

int random_int(int min, int max){
        std::uniform_int_distribution<> rnd_int(min, max);
        return rnd_int(generator);
}

std::vector<int> shuffle_vector(std::vector<int> vec){
        int j = 0;
        int tmp = 0;
        for(int i = vec.size()-1; i > 0; i--) {
                j = random_int(0,i);
                tmp = vec[i];
                vec[i] = vec[j];
                vec[j] = tmp;
        }
        return vec;
}

int coord(int x, int y, int grid_size_x, int grid_size_y){
        return (x%grid_size_x) * grid_size_y + (y%grid_size_y);
}
