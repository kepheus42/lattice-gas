#include "global.hpp"
#include "output.hpp"
#include <iostream>
#include <ctime>
#include <vector>
#include <boost/format.hpp>
#include <string>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>


int coord(int x, int y, int max_x){
        return x+y*(255);
}

int select_direction(){
        return directions(generator);
}

int main(int ac, char** av){
        int grid_size = 255;
        std::vector<int> grid_vector (grid_size*grid_size,0);
        int time = std::time(NULL);
        int c = 1;
        int x = 128;
        int y = 128;
        int dir_x = 0;
        int dir_y = 0;
        int max_t = atoi(av[1]);
        int seed = atoi(av[2]);
        generator.seed(seed);
        std::vector<int> pos_vector = {x,y};

        // 127 spots to the left, 127 spots to the right
        grid_vector[coord(128,128,grid_size)] = 1;

        for (int t = 1; t < max_t; t++) {
                bool walker_is_on_graph = true;
                while(walker_is_on_graph) {
                        int dir = select_direction();
                        switch(dir) {
                        case 1: dir_x = 1; dir_y = 0; break;
                        case 2: dir_x = 0; dir_y = 1; break;
                        case 3: dir_x =-1; dir_y = 0; break;
                        case 4: dir_x = 0; dir_y =-1; break;
                        }
                        if(grid_vector[coord(x+dir_x,y+dir_y,grid_size)] == 0) {
                                c++;
                                x+=dir_x;
                                y+=dir_y;
                                grid_vector[coord(x,y,grid_size)] = c;
                                walker_is_on_graph = false;
                        }
                        // check if the walker is trapped
                        if(!(grid_vector[coord(x+1,y,grid_size)]==0)&&!(grid_vector[coord(x,y+1,grid_size)]==0)&&!(grid_vector[coord(x-1,y,grid_size)]==0)&&!(grid_vector[coord(x,y-1,grid_size)]==0)) {
                                t = max_t;
                                break;
                        }
                }
                pos_vector.push_back(x);
                pos_vector.push_back(y);
        }
        vector_to_file(pos_vector, boost::str(boost::format("../data/self_avoiding_walk_t%i_%i.bin") % max_t % seed));
}
