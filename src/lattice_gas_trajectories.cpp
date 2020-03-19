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
        // exit program if too many tracers are called
        if (grid_size_x*grid_size_y < number_of_tracers_1x1 + 4*number_of_tracers_2x2) {
                throw std::invalid_argument("Too many tracers for grid of the given size!");
        }
        std::cout << "X = " << grid_size_x << " Y = " << grid_size_y << " N1 = " << number_of_tracers_1x1 << " N2 = " << number_of_tracers_2x2 << " N3 = " << number_of_tracers_3x3 << " T = " << max_number_of_timesteps << " O = " << term_output << '\n';
        // - - - - - - - - - - - - -
        std::vector<double> output_trajectories;
        std::vector<double> output_msd;
        std::string output_file_name_msd;
        std::string output_file_name_trajectories;
        std::string header_string;
        output_file_name_msd = "msd_X_"+std::to_string(grid_size_x)+"_Y_"+std::to_string(grid_size_y)+"_N1_"+std::to_string(number_of_tracers_1x1)+"_N2_"+std::to_string(number_of_tracers_2x2)+"_N3_"+std::to_string(number_of_tracers_3x3)+".bin";
        output_file_name_trajectories = "trajectories_X_"+std::to_string(grid_size_x)+"_Y_"+std::to_string(grid_size_y)+"_N1_"+std::to_string(number_of_tracers_1x1)+"_N2_"+std::to_string(number_of_tracers_2x2)+"_N3_"+std::to_string(number_of_tracers_3x3)+".bin";
        header_string = "Hello!";
        // TODO: search for files with output_file_name in their name
        // set a counter to n+1 where n is in output_file_name_n
        // if none exist, start with n=0
        // if any exist, continue with n++
        // - - - - - - - - - - - - -
        Lattice * lattice = new Lattice(grid_size_x,grid_size_y,number_of_tracers_1x1,number_of_tracers_2x2,number_of_tracers_3x3);
        // - - - - - - - - - - - - -
        int progress_percent = 0;
        // for trajectories output
        int max_number_of_trajectories = 500;
        int max_number_of_trajectories_1x1 = number_of_tracers_1x1;
        int max_number_of_trajectories_2x2 = number_of_tracers_2x2;
        int max_number_of_trajectories_3x3 = number_of_tracers_3x3;
        if(number_of_tracers_1x1 > max_number_of_trajectories)
        {
                max_number_of_trajectories_1x1 = max_number_of_trajectories;
        }
        if(number_of_tracers_2x2 > max_number_of_trajectories)
        {
                max_number_of_trajectories_2x2 = max_number_of_trajectories;
        }
        if(number_of_tracers_3x3 > max_number_of_trajectories)
        {
                max_number_of_trajectories_3x3 = max_number_of_trajectories;
        }
        int total_number_of_trajectories = max_number_of_trajectories_1x1 + max_number_of_trajectories_2x2 + max_number_of_trajectories_3x3;
        // run warm up to disperse 2x2 and 3x3 tracers
        lattice->warm_up(timesteps_warmup);
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
                // save msd vs t
                output_msd.push_back((double)lattice->get_t());
                output_msd.push_back(lattice->get_avg_lsquared_1x1());
                output_msd.push_back(lattice->get_avg_lsquared_2x2());
                output_msd.push_back(lattice->get_avg_lsquared_3x3());
                // timestep
                lattice->timestep();
                for(int n = 0; n < total_number_of_trajectories; n++)
                {
                        int m = n;
                        if((max_number_of_trajectories_2x2 == 100)&&(n > (max_number_of_trajectories_2x2-1)))
                        {
                                m += number_of_tracers_2x2-1-max_number_of_trajectories_2x2;
                        }
                        output_trajectories.push_back(lattice->get_tracer_last_move(m));
                }
        }
        if (term_output)
        {
                std::cout << "100%" << std::endl;
        }
        // save msd vs t
        output_msd.push_back((double)lattice->get_t());
        output_msd.push_back(lattice->get_avg_lsquared_1x1());
        output_msd.push_back(lattice->get_avg_lsquared_2x2());
        output_msd.push_back(lattice->get_avg_lsquared_3x3());
        // write msd and trajectories to disk
        vector_to_file(output_msd,output_file_name_msd);
        vector_to_file(output_trajectories,output_file_name_trajectories);
}
