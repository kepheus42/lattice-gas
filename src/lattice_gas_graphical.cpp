#include "global.hpp"
#include "tracer.hpp"
#include "lattice.hpp"
#include "output.hpp"
#include <iostream>
#include <chrono>
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
        int number_of_tracers_3x3 = atoi(av[5]);
        int max_number_of_timesteps = atoi(av[6]);
        int timesteps_warmup = atoi(av[7]);
        int term_output = atoi(av[8]);
        //
        // exit program if too many tracers are called
        if (grid_size_x*grid_size_y < number_of_tracers_1x1 + 4*number_of_tracers_2x2) {
                throw std::invalid_argument("Too many tracers for grid of the given size!");
        }
        std::cout << "X = " << grid_size_x << " Y = " << grid_size_y << " N1 = " << number_of_tracers_1x1 << " N2 = " << number_of_tracers_2x2 << " N3 = " << number_of_tracers_3x3 << " T = " << max_number_of_timesteps << " O = " << term_output << '\n';
        // - - - - - - - - - - - - -
        std::string output_file_name;
        output_file_name = "snapshot_X_"+std::to_string(grid_size_x)+"_Y_"+std::to_string(grid_size_y)+"_N1_"+std::to_string(number_of_tracers_1x1)+"_N2_"+std::to_string(number_of_tracers_2x2)+"_N3_"+std::to_string(number_of_tracers_3x3)+"_T_"+std::to_string(max_number_of_timesteps);
        // - - - - - - - - - - - - -
        std::vector<int> output_snapshot;
        // - - - - - - - - - - - - -
        Lattice * lattice = new Lattice(grid_size_x,grid_size_y,number_of_tracers_1x1,number_of_tracers_2x2,number_of_tracers_3x3);
        // - - - - - - - - - - - - -
        int progress_percent = 0;
        // run warm up to disperse 2x2 and 3x3 tracers
        lattice->warm_up(timesteps_warmup);
        //
        while(lattice->get_t() < max_number_of_timesteps)
        {
                if (term_output)
                {
                        if((int)(100*lattice->get_t()/max_number_of_timesteps) > progress_percent) {
                                progress_percent = (int)(100*lattice->get_t()/max_number_of_timesteps);
                                fflush(stdout);
                                std::cout << progress_percent << "%\r";
                        }
                }
                // timestep
                for(int entry : lattice->get_tracer_positions()) {
                        output_snapshot.push_back(entry);
                }
                lattice->timestep();
        }
        if (term_output)
        {
                std::cout << "100%" << std::endl;
        }
        // final output after last timestep
        for(int entry : lattice->get_tracer_positions()) {
                output_snapshot.push_back(entry);
        }

        vector_to_file(output_snapshot,output_file_name);
}
