#include "global.hpp"
#include "tracer.hpp"
#include "output.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <string>

int main(int ac, char** av){
        set_rng_seed(0);
        // - - - - - - - - - - - - -
        Wrapper * wrapper = new Wrapper(10,10,1,0,0,1.0,1.0,1.0,25,50);
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
