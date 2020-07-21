#include "global.hpp"
#include "tracer.hpp"
#include "lattice.hpp"
#include "output.hpp"
#include "wrapper.hpp"
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

<<<<<<< HEAD
=======

>>>>>>> 0bff648d5adb391c294373620f1565d115c1c73c
int main(int ac, char** av){
        set_rng_seed(0);
        printf("START!\n");
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // A L L  I N P U T  V A R I A B L E S
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int grid_size_x = atoi(av[1]);
        int grid_size_y = atoi(av[2]);
        int number_of_timesteps = atoi(av[3]);
        int number_of_timesteps_warmup = atoi(av[4]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int number_of_tracers_1x1 = atoi(av[5]);
        int number_of_tracers_2x2 = atoi(av[6]);
<<<<<<< HEAD
=======
        int number_of_tracers_3x3 = 0; //atoi(av[5]);
>>>>>>> 0bff648d5adb391c294373620f1565d115c1c73c
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int wtd_max = atoi(av[7]);
        int wtd_res = atoi(av[8]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
<<<<<<< HEAD
        double step_rate_1x1 = atoi(av[9]);
        double step_rate_2x2 = atoi(av[10]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int number_of_lattices = atoi(av[11]);
=======
        double step_rate_1x1 = 1;
        double step_rate_2x2 = 1/2;
        double step_rate_3x3 = 1/3; // atof(av[11]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int number_of_lattices = atoi(av[9]);
>>>>>>> 0bff648d5adb391c294373620f1565d115c1c73c
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        printf("X:%4d Y:%4d 1:%6d 2:%6d T:%8d W:%6d L:%10d\n",grid_size_x,grid_size_y,number_of_tracers_1x1,number_of_tracers_2x2,number_of_timesteps,number_of_timesteps_warmup,number_of_lattices);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_msd_1x1;
        std::string output_file_name_msd_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_rate_1x1;
        std::string output_file_name_rate_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
<<<<<<< HEAD
        std::string output_file_name_sublattice_conc_1x1;
        std::string output_file_name_sublattice_conc_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
=======
>>>>>>> 0bff648d5adb391c294373620f1565d115c1c73c
        std::string output_file_name_wtd_1x1;
        std::string output_file_name_wtd_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_corr_1x1;
        std::string output_file_name_corr_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_base_string = "X_"+std::to_string(grid_size_x)+"_Y_"+std::to_string(grid_size_y)+"_N1_"+std::to_string(number_of_tracers_1x1)+"_N2_"+std::to_string(number_of_tracers_2x2)+"_N3_"+std::to_string(number_of_tracers_3x3)+"_T_"+std::to_string(number_of_timesteps)+"_W_"+std::to_string(number_of_timesteps_warmup)+".bin";
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_msd_1x1 = "msd_1x1_"+output_file_name_base_string;
        output_file_name_msd_2x2 = "msd_2x2_"+output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_rate_1x1 = "rate_1x1_"+output_file_name_base_string;
        output_file_name_rate_2x2 = "rate_2x2_"+output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
<<<<<<< HEAD
        output_file_name_sublattice_conc_1x1 = "sub_c_1x1_"+output_file_name_base_string;
        output_file_name_sublattice_conc_2x2 = "sub_c_2x2_"+output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_neighbor_count_1x1 = "";
        output_file_name_neighbor_count_2x2 = "";
        // - - - - - - - - - - - - - - - - - - - - - - - - -
=======
>>>>>>> 0bff648d5adb391c294373620f1565d115c1c73c
        output_file_name_wtd_1x1 = "wtd_1x1_"+output_file_name_base_string;
        output_file_name_wtd_2x2 = "wtd_2x2_"+output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_corr_1x1 = "corr_1x1_"+output_file_name_base_string;
        output_file_name_corr_2x2 = "corr_2x2_"+output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string header_string;
        header_string = "Hello!";
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // create a wrapper
        D(std::cerr << "Creating Wrapper ..." << std::endl);
        Wrapper * wrapper = new Wrapper(number_of_lattices,
                                        grid_size_x,
                                        grid_size_y,
                                        number_of_timesteps,
                                        number_of_timesteps_warmup,
                                        number_of_tracers_1x1,
                                        number_of_tracers_2x2,
                                        number_of_tracers_3x3,
                                        step_rate_1x1,
                                        step_rate_2x2,
                                        step_rate_3x3,
                                        wtd_max,
                                        wtd_res);
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
<<<<<<< HEAD
                vector_to_file(wrapper->get_result_sublattice_conc_1x1(),
                               output_file_name_sublattice_conc_1x1);
=======
>>>>>>> 0bff648d5adb391c294373620f1565d115c1c73c
                vector_to_file(wrapper->get_result_wtd_1x1(),
                               output_file_name_wtd_1x1);
                vector_to_file(wrapper->get_result_correlations_1x1(),
                               output_file_name_corr_1x1);
        }
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        if(number_of_tracers_2x2)
        {
                vector_to_file(wrapper->get_result_lsquared_2x2(),
                               output_file_name_msd_2x2);
                vector_to_file(wrapper->get_result_rate_2x2(),
                               output_file_name_rate_2x2);
<<<<<<< HEAD
                vector_to_file(wrapper->get_result_sublattice_conc_2x2(),
                               output_file_name_sublattice_conc_2x2);
=======
>>>>>>> 0bff648d5adb391c294373620f1565d115c1c73c
                vector_to_file(wrapper->get_result_wtd_2x2(),
                               output_file_name_wtd_2x2);
                vector_to_file(wrapper->get_result_correlations_2x2(),
                               output_file_name_corr_2x2);
        }
<<<<<<< HEAD
        // for(double d : wrapper->get_result_sublattice_conc_2x2()) { std::cout << d << std::endl;}
=======
>>>>>>> 0bff648d5adb391c294373620f1565d115c1c73c
        printf("\rProgress:  DONE\n");
        // - - - - - - - - - - - - - - - - - - - - - - - - -
}
