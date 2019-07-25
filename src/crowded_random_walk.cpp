#include "global.hpp"
#include "lattice.hpp"
#include "tracer.hpp"
#include "output.hpp"
#include <iostream>
#include <vector>
#include <boost/format.hpp>
#include <string>

#define STORE_POSITION false

std::uniform_int_distribution<> step_directions(1,4);

int coord(int x, int y){
        return x+y*grid_size;
}

int select_direction(){
        return step_directions(generator);
}

int main(int ac, char** av){
        // sim parameters
        int seed = atoi(av[1]);
        int max_t = atoi(av[2]);
        int size_x = atoi(av[3]);
        int size_y = atoi(av[4]);
        double density_1 = atof(av[5]);
        double density_2 = atof(av[6]);
        // fix grid size
        // number of tracers from density:
        int number_of_tracers_1x1 = size_x * size_y * density_1;
        int number_of_tracers_2x2 = size_x * size_y * density_2;
        // initialize lattice:
        Lattice lattice = new Lattice(size_x,size_y,number_of_tracers_1x1,number_of_tracers_2x2);
        set_rng_seed(seed);
        // run for max_t timesteps
        for (int t = 0; t < max_t; t++) {
                lattice.timestep();
        }
        // output

        std::uniform_int_distribution<> lattice_position(0,grid_size-1);
        // set the rng seed for reproducible results
        // 0 sets rng.seed to current time
        //
        std::vector<int> update_order (n_walkers,0);
        std::vector<double> lsq_vector = {0};
        std::cout << grid_vector.size() << std::endl;
        std::vector<int> pos_vector;
        // create a vector containing the desired number of tracers
        for (int n = 0; n < n_walkers; n++) {
                update_order[n] = n;
                int x = lattice_position(generator);
                int y = lattice_position(generator);
                tracers.push_back(new Tracer(x,y));
                // place tracers on the lattice
                grid_vector[coord(x,y)] = 1;
                if(STORE_POSITION) {
                        pos_vector.push_back(x);
                        pos_vector.push_back(y);
                }
        }
        // loop over time range
        for (int t = 1; t < max_t; t++) {
                // temporary
                double lsq_average = 0.0;
                std::vector<double> lsq_vector_temp;
                // loop over tracers
                update_order = shuffle_vector(update_order);
                for (int i : update_order) {
                        Tracer * tr = tracers[i];
                        // time-step logic:
                        // randomly select direction
                        int dir_x, dir_y;
                        int dir = select_direction();
                        switch(dir) {
                        case 1: dir_x = 1; dir_y = 0; break;
                        case 2: dir_x = 0; dir_y = 1; break;
                        case 3: dir_x =-1; dir_y = 0; break;
                        case 4: dir_x = 0; dir_y =-1; break;
                        }
                        // update
                        int new_x = tr->x+dir_x;
                        int new_x_mod = (grid_size+new_x%grid_size)%grid_size;
                        int new_y = tr->y+dir_y;
                        int new_y_mod = (grid_size+new_y%grid_size)%grid_size;
                        // if new location is empty
                        if (grid_vector[coord(new_x_mod,new_y_mod)] == 0) {
                                // set old location to empty, update location at tr*, and set new location to filled
                                grid_vector[coord(tr->x,tr->y)] = 0;
                                tr->x=new_x_mod;
                                tr->y=new_y_mod;
                                tr->dx += dir_x;
                                tr->dy += dir_y;
                                grid_vector[coord(tr->x,tr->y)] = 1;
                        }
                        // statistics computation and storage
                        // calculate lsquared
                        tr->calc_lsquared();
                        // store lsq for averaging
                        lsq_vector_temp.push_back(tr->l_squared);
                        // store tracer position
                        if(STORE_POSITION) {
                                pos_vector.push_back(tr->x);
                                pos_vector.push_back(tr->y);
                        }
                }
                // compute average lsqaured of all tracers
                for (double lsq : lsq_vector_temp) {
                        lsq_average += lsq;
                }
                lsq_average /= lsq_vector_temp.size();
                // push back average lsquared to vector for output
                lsq_vector.push_back(lsq_average);
        }
        density = (double)n_walkers/(double)grid_size/(double)grid_size;
        if(STORE_POSITION) {
                vector_to_file(pos_vector, boost::str(boost::format("../data/crowded_random_walk_positions_%f_n%i_t%i_s%i.bin") % density % n_walkers % max_t % seed));
        }
        vector_to_file(lsq_vector, boost::str(boost::format("../data/crowded_random_walk_statistics_c%f_n%i_t%i_s%i.bin") % density % n_walkers % max_t % seed));
}
