
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
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        /* I N P U T   P A R A M E T E R S                         */
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        int rng_seed = atoi(av[1]);
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        int grid_size = atoi(av[2]);
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        int number_of_timesteps = atoi(av[3]);
        int number_of_timesteps_warmup = atoi(av[4]);
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        int number_of_tracers_1x1 = atoi(av[5]);
        int number_of_tracers_2x2 = atoi(av[6]);
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        double step_rate_1x1 = atof(av[7]);
        double step_rate_2x2 = atof(av[8]);
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        int number_of_lattices = atoi(av[9]);
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        long n_step_attempts = ((long)number_of_timesteps+(long)number_of_timesteps_warmup)*((long)number_of_tracers_1x1/(long)(1/step_rate_1x1)+(long)number_of_tracers_2x2/(long)(1/step_rate_2x2))*(long)number_of_lattices;
        printf("START! MCS: %12li\n",n_step_attempts);
        printf("X:%4d 1:%6d 2:%6d T:%8d W:%6d L:%10d\n",
               grid_size,
               number_of_tracers_1x1,
               number_of_tracers_2x2,
               number_of_timesteps,
               number_of_timesteps_warmup,
               number_of_lattices);
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        std::string output_fname_base =
                "X_"   + std::to_string(grid_size) +
                "_N1_" + std::to_string(number_of_tracers_1x1) +
                "_N2_" + std::to_string(number_of_tracers_2x2) +
                "_T_"  + std::to_string(number_of_timesteps) +
                "_W_"  + std::to_string(number_of_timesteps_warmup) + ".bin";
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        // create a wrapper
        D(std::cerr << "Creating Wrapper ..." << std::endl);
        Wrapper * wrapper = new Wrapper(rng_seed,
                                        number_of_lattices,
                                        grid_size,
                                        number_of_timesteps,
                                        number_of_timesteps_warmup,
                                        number_of_tracers_1x1,
                                        number_of_tracers_2x2,
                                        step_rate_1x1,
                                        step_rate_2x2);
        D(std::cerr << "Done!" << std::endl);
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        clock::time_point start_dynamics = clock::now();
        wrapper->warmup();
        wrapper->evolve();
        clock::time_point end_dynamics = clock::now();
        std::chrono::duration<double> duration_dynamics = end_dynamics - start_dynamics;
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        D( std::cout << "Data output..." << std::endl );
        if(number_of_tracers_1x1)
        {
                vector_to_file( wrapper->get_result_diff_1x1(),      "1x1_diff_"  + output_fname_base );
                vector_to_file( wrapper->get_result_rate_1x1(),      "1x1_rate_"  + output_fname_base );
                vector_to_file( wrapper->get_result_step_corr_1x1(), "1x1_steps_" + output_fname_base );
                vector_to_file( wrapper->get_result_pos_1x1(),       "1x1_pos_"   + output_fname_base );
                vector_to_file( wrapper->get_result_site_corr_1x1(), "1x1_site_corr_"   + output_fname_base );

        }
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
        if(number_of_tracers_2x2)
        {
                vector_to_file( wrapper->get_result_diff_2x2(),      "2x2_diff_"  + output_fname_base );
                vector_to_file( wrapper->get_result_rate_2x2(),      "2x2_rate_"  + output_fname_base );
                vector_to_file( wrapper->get_result_step_corr_2x2(), "2x2_steps_" + output_fname_base );
                vector_to_file( wrapper->get_result_pos_2x2(),       "2x2_pos_"   + output_fname_base );
                vector_to_file( wrapper->get_result_site_corr_2x2(), "2x2_site_corr_"   + output_fname_base );
        }
        D( std::cout << "Done!" << std::endl );
        printf("Done! %.2fs\n",duration_dynamics.count());
        /* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
}
