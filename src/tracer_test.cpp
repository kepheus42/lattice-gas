#include "global.hpp"
#include "tracer.hpp"
#include "lattice.hpp"
#include "output.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <boost/format.hpp>
#include <string>

int main(int ac, char** av){
        set_rng_seed(0);
        // - - - - - - - - - - - - -
        // A L L  V A R I A B L E S
        // - - - - - - - - - - - - -
        int grid_size_x = atoi(av[1]);
        int grid_size_y = atoi(av[2]);
        int number_of_tracers_1x1 = atoi(av[3]);
        int number_of_tracers_2x2 = atoi(av[4]);
        int max_number_of_timesteps = atoi(av[5]);
        // - - - - - - - - - - - - -
        Lattice * lattice = new Lattice(grid_size_x,grid_size_y,number_of_tracers_1x1,number_of_tracers_2x2);
        // - - - - - - - - - - - - -
        std::cout << "Timestep + avg lsquared of tracer populations" << std::endl;
        int exponent_counter = 0;
        while(lattice->get_t() < max_number_of_timesteps)
        {
                if(lattice->get_t()/(int)pow(10,exponent_counter+1)) {
                        exponent_counter++;
                }
                if(lattice->get_t()%(int)pow(10,exponent_counter) == 0) {
                        std::cout << lattice->get_t() << "\t" << lattice->get_avg_lsquared_1x1() << "\t" << lattice->get_avg_lsquared_2x2() << "\t" << std::endl;
                }
                lattice->timestep_no_interaction();
        }
        std::cout << lattice->get_t() << "\t" << lattice->get_avg_lsquared_1x1() << "\t" << lattice->get_avg_lsquared_2x2() << std::endl;
}
