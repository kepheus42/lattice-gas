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
        int number_of_tracers_2x2 = 0;
        int number_of_tracers_3x3 = 0;
        int max_number_of_timesteps = atoi(av[4]);
        int term_output = atoi(av[5]);
        // exit program if too many tracers are called
        if (grid_size_x*grid_size_y < number_of_tracers_1x1 + 4*number_of_tracers_2x2) {
                throw std::invalid_argument("Too many tracers for grid of the given size!");
        }
        std::cout << grid_size_x << ' ' << grid_size_y << ' ' << number_of_tracers_1x1 << ' ' << number_of_tracers_2x2 << ' ' << max_number_of_timesteps << ' ' << term_output << '\n';
        // - - - - - - - - - - - - -
        std::vector<double> output;
        std::string output_file_name;
        std::string header_string;
        output_file_name = "X_"+std::to_string(grid_size_x)+"_Y_"+std::to_string(grid_size_y)+"_N1_"+std::to_string(number_of_tracers_1x1)+"_N2_"+std::to_string(number_of_tracers_2x2)+"_N3_"+std::to_string(number_of_tracers_3x3)+".bin";
        header_string = "Hello!";
        // TODO: search for files with output_file_name in their name
        // set a counter to n+1 where n is in output_file_name_n
        // if none exist, start with n=0
        // if any exist, continue with n++
        // - - - - - - - - - - - - -
        Lattice * lattice = new Lattice(grid_size_x,grid_size_y,number_of_tracers_1x1,number_of_tracers_2x2,number_of_tracers_3x3);
        // - - - - - - - - - - - - -
        if (term_output)
        {
                std::cout << "Timestep + avg lsquared of tracer populations" << std::endl;
        }
        int exponent_counter = 0;
        while(lattice->get_t() < max_number_of_timesteps)
        {
                if(lattice->get_t()/(int)pow(10,exponent_counter+1)) {
                        exponent_counter++;
                }
                if(lattice->get_t()%(int)pow(10,exponent_counter) == 0) {
                        if (term_output)
                        {
                                std::cout << lattice->get_t() << "\t" << lattice->get_avg_lsquared_1x1() << "\t" << lattice->get_avg_lsquared_2x2() << "\t" << std::endl;
                        }
                        // output.push_back((double)lattice->get_t());
                        // output.push_back(lattice->get_avg_lsquared_1x1());
                        // output.push_back(lattice->get_avg_lsquared_2x2());
                }
                output.push_back((double)lattice->get_t());
                output.push_back(lattice->get_avg_lsquared_1x1());
                output.push_back(lattice->get_avg_lsquared_2x2());
                lattice->timestep_no_interaction();
        }
        if (term_output)
        {
                std::cout << lattice->get_t() << "\t" << lattice->get_avg_lsquared_1x1() << "\t" << lattice->get_avg_lsquared_2x2() << std::endl;
        }
        output.push_back((double)lattice->get_t());
        output.push_back(lattice->get_avg_lsquared_1x1());
        output.push_back(lattice->get_avg_lsquared_2x2());
        vector_to_file(output,output_file_name);
}
