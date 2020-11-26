#include "global.hpp"
#include "lattice.hpp"
#include "tracer.hpp"
#include "output.hpp"
#include "wrapper.hpp"
#include "site.hpp"
#include <iostream>
#include <chrono>
#include <vector>
#include <random>
#include <boost/format.hpp>
#include <string>

// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif
// To add debugging messages, use D(std::cerr << "Debugging message 1 2 3!" << std::endl; )

int main(int ac, char** av){
        using clock = std::chrono::steady_clock;
        clock::time_point start = clock::now();
        printf("START!\n");
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // A L L  I N P U T  V A R I A B L E S
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int rng_seed = atoi(av[1]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int grid_size = atoi(av[2]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int number_of_timesteps = atoi(av[3]);
        int number_of_timesteps_warmup = atoi(av[4]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int number_of_tracers_1x1 = atoi(av[5]);
        int number_of_tracers_2x2 = atoi(av[6]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        double step_rate_1x1 = atof(av[7]);
        double step_rate_2x2 = atof(av[8]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int number_of_lattices = atoi(av[9]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int number_of_pos_to_save = atoi(av[10]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        set_rng_seed(rng_seed);
        printf("X:%4d 1:%6d 2:%6d T:%8d W:%6d L:%10d\n",
               grid_size,
               number_of_tracers_1x1,
               number_of_tracers_2x2,
               number_of_timesteps,
               number_of_timesteps_warmup,
               number_of_lattices);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_msd_1x1;
        std::string output_file_name_msd_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_rate_1x1;
        std::string output_file_name_rate_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_corr_1x1;
        std::string output_file_name_corr_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_pos_1x1;
        std::string output_file_name_pos_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_base_string =
                "XY_"+std::to_string(grid_size)+
                "_N1_"+std::to_string(number_of_tracers_1x1)+
                "_N2_"+std::to_string(number_of_tracers_2x2)+
                "_T_"+std::to_string(number_of_timesteps)+
                "_W_"+std::to_string(number_of_timesteps_warmup)+".bin";
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_msd_1x1 = "msd_1x1_"+output_file_name_base_string;
        output_file_name_msd_2x2 = "msd_2x2_"+output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_rate_1x1 = "rate_1x1_"+output_file_name_base_string;
        output_file_name_rate_2x2 = "rate_2x2_"+output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_corr_1x1 = "corr_1x1_"+output_file_name_base_string;
        output_file_name_corr_2x2 = "corr_2x2_"+output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_pos_1x1 = "pos_1x1_"+output_file_name_base_string;
        output_file_name_pos_2x2 = "pos_2x2_"+output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string header_string;
        header_string = "Hello!";
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // create a wrapper
        D(std::cerr << "Creating Wrapper ..." << std::endl);
        Wrapper * wrapper = new Wrapper(number_of_lattices,
                                        grid_size,
                                        number_of_timesteps,
                                        number_of_timesteps_warmup,
                                        number_of_tracers_1x1,
                                        number_of_tracers_2x2,
                                        step_rate_1x1,
                                        step_rate_2x2,
                                        number_of_pos_to_save);
        D(std::cerr << "Done!" << std::endl);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // perform warm up if number_of_timesteps_warmup
        D(std::cerr << "Starting Warmup ..." << std::endl);
        for(int t = 0; t < number_of_timesteps_warmup; t++) { wrapper->timestep_warmup(); }
        D(std::cerr << "Done!" << std::endl);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // let the simulation run for number_of_timesteps iterations
        int progress = 0;
        D(std::cerr << "Starting Simulation Run ..." << std::endl);
        for(int t = 0; t < number_of_timesteps; t++)
        {
                if(1000*t/number_of_timesteps > progress) {
                        progress = 1000*t/number_of_timesteps;
                        printf("\rProgress:%6.1f",(double)progress/10.0);
                        fflush(stdout);
                }
                wrapper->timestep();
        }
        D(std::cerr << "Done!" << std::endl);
        //printf("RUN:%f4.1\n", (double)progress_percent/10.0);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        if(number_of_tracers_1x1)
        {
                vector_to_file(wrapper->get_result_lsquared_1x1(),
                               output_file_name_msd_1x1);
                vector_to_file(wrapper->get_result_rate_1x1(),
                               output_file_name_rate_1x1);
                vector_to_file(wrapper->get_result_correlations_1x1(),
                               output_file_name_corr_1x1);
                vector_to_file(wrapper->get_result_positions_1x1(),
                               output_file_name_pos_1x1);
        }
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        if(number_of_tracers_2x2)
        {
                vector_to_file(wrapper->get_result_lsquared_2x2(),
                               output_file_name_msd_2x2);
                vector_to_file(wrapper->get_result_rate_2x2(),
                               output_file_name_rate_2x2);
                vector_to_file(wrapper->get_result_correlations_2x2(),
                               output_file_name_corr_2x2);
                vector_to_file(wrapper->get_result_positions_2x2(),
                               output_file_name_pos_2x2);
        }
        // for(double d : wrapper->get_result_sublattice_conc_2x2()) { std::cout << d << std::endl;}
        clock::time_point end = clock::now();
        std::chrono::duration<double> execution_time = end - start;
        printf("\rProgress:  DONE\n");
        long n_step_attempts = (long)(number_of_timesteps+number_of_timesteps_warmup)*(number_of_tracers_1x1*step_rate_1x1+number_of_tracers_2x2*step_rate_2x2)*number_of_lattices;
        printf("Total number of step attempts: %li Time elapsed: %f\n",n_step_attempts,execution_time);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
}
