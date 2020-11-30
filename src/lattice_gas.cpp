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
        printf("X:%4d 1:%6d 2:%6d T:%8d W:%6d L:%10d\n",
               grid_size,
               number_of_tracers_1x1,
               number_of_tracers_2x2,
               number_of_timesteps,
               number_of_timesteps_warmup,
               number_of_lattices);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_base_string =
                "X_"+std::to_string(grid_size)+
                "_N1_"+std::to_string(number_of_tracers_1x1)+
                "_N2_"+std::to_string(number_of_tracers_2x2)+
                "_T_"+std::to_string(number_of_timesteps)+
                "_W_"+std::to_string(number_of_timesteps_warmup)+".bin";
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_msd_1x1 = "msd_1x1_"+output_file_name_base_string;
        std::string output_file_name_msd_2x2 = "msd_2x2_"+output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_rate_1x1 = "rate_1x1_"+output_file_name_base_string;
        std::string output_file_name_rate_2x2 = "rate_2x2_"+output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_corr_1x1 = "corr_1x1_"+output_file_name_base_string;
        std::string output_file_name_corr_2x2 = "corr_2x2_"+output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_pos_1x1 = "pos_1x1_"+output_file_name_base_string;
        std::string output_file_name_pos_2x2 = "pos_2x2_"+output_file_name_base_string;
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
                                        step_rate_2x2);
        D(std::cerr << "Done!" << std::endl);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // perform warm up if number_of_timesteps_warmup
        clock::time_point start_dynamics = clock::now();
        wrapper->warmup();
        wrapper->evolve();
        clock::time_point end_dynamics = clock::now();
        std::chrono::duration<double> duration_dynamics = end_dynamics - start_dynamics;
        long n_step_attempts = ((long)number_of_timesteps+(long)number_of_timesteps_warmup)*((long)number_of_tracers_1x1*(long)step_rate_1x1+(long)number_of_tracers_2x2*(long)step_rate_2x2)*(long)number_of_lattices;
        // for(int t = 0; t < number_of_timesteps_warmup; t++) { wrapper->timestep_warmup(); }
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // let the simulation run for number_of_timesteps iterations
        // for(int t = 0; t < number_of_timesteps; t++) { wrapper->timestep(); }
        /*
           int progress = 0;
           if(1000*t/number_of_timesteps > progress) {
                progress = 1000*t/number_of_timesteps;
                printf("\rProgress:%6.1f",(double)progress/10.0);
                fflush(stdout);
           }
         */
        //printf("RUN:%f4.1\n", (double)progress_percent/10.0);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        D( std::cout << "Data output..." << std::endl );
        if(number_of_tracers_1x1)
        {
                vector_to_file(wrapper->get_result_lsq_1x1(), output_file_name_msd_1x1);
                vector_to_file(wrapper->get_result_rate_1x1(), output_file_name_rate_1x1);
                // vector_to_file(wrapper->get_result_correlations_1x1(), output_file_name_corr_1x1);
                // vector_to_file(wrapper->get_result_positions_1x1(), output_file_name_pos_1x1);
        }
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        if(number_of_tracers_2x2)
        {
                vector_to_file(wrapper->get_result_lsq_2x2(), output_file_name_msd_2x2);
                vector_to_file(wrapper->get_result_rate_2x2(), output_file_name_rate_2x2);
                // vector_to_file(wrapper->get_result_correlations_2x2(),output_file_name_corr_2x2);
                // vector_to_file(wrapper->get_result_positions_2x2(), output_file_name_pos_2x2);
        }
        D( std::cout << "Done!" << std::endl );
        // for(double d : wrapper->get_result_sublattice_conc_2x2()) { std::cout << d << std::endl;}
        // printf("\rProgress:  DONE\n");
        printf("MCS: %li %.2fs\n",n_step_attempts,duration_dynamics);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
}
