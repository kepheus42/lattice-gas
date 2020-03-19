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

int main(int ac, char** av){
        set_rng_seed(0);
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
        int number_of_tracers_3x3 = 0; //atoi(av[5]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int max_wtd = 100;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        double step_rate_1x1 = 1;
        double step_rate_2x2 = 1/2;
        double step_rate_3x3 = 1/3; // atof(av[11]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        int number_of_lattices = atoi(av[7]);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        printf("X:%4d Y:%4d 1:%6d 2:%6d T:%8d W:%6d L:%10d",grid_size_x,grid_size_y,number_of_tracers_1x1,number_of_tracers_2x2,number_of_timesteps,number_of_timesteps_warmup,number_of_lattices);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_avg_msd_1x1;
        std::string output_file_name_avg_msd_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_avg_rate_1x1;
        std::string output_file_name_avg_rate_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_wtd_1x1;
        std::string output_file_name_wtd_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_rate_dist_1x1;
        std::string output_file_name_rate_dist_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_rates_1x1;
        std::string output_file_name_rates_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_two_step_correlation_1x1;
        std::string output_file_name_two_step_correlation_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_three_step_correlation_1x1;
        std::string output_file_name_three_step_correlation_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string output_file_name_four_step_correlation_1x1;
        std::string output_file_name_four_step_correlation_2x2;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_base_string = "X_"+std::to_string(grid_size_x)+"_Y_"+std::to_string(grid_size_y)+"_N1_"+std::to_string(number_of_tracers_1x1)+"_N2_"+std::to_string(number_of_tracers_2x2)+"_N3_"+std::to_string(number_of_tracers_3x3)+"_T_"+std::to_string(number_of_timesteps)+"_W_"+std::to_string(number_of_timesteps_warmup)+".bin";
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_avg_msd_1x1 = "avg_msd_1x1_" + output_file_name_base_string;
        output_file_name_avg_msd_2x2 = "avg_msd_2x2_" + output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_avg_rate = "avg_rate_1x1_" + output_file_name_base_string;
        output_file_name_avg_rate = "avg_rate_2x2_" + output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_wtd_1x1 = "wtd_1x1_" + output_file_name_base_string;
        output_file_name_wtd_2x2 = "wtd_2x2_" + output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_rate_dist_1x1 = "rate_dist_1x1_" + output_file_name_base_string;
        output_file_name_rate_dist_2x2 = "rate_dist_2x2_" + output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_rates_1x1 =  "rates_1x1_" + output_file_name_base_string;
        output_file_name_rates_2x2 =  "rates_2x2_" + output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_two_step_correlation_1x1 = "two_step_corr_1x1_" + output_file_name_base_string;
        output_file_name_two_step_correlation_2x2 = "two_step_corr_2x2_" + output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_three_step_correlation_1x1 = "three_step_corr_1x1_" + output_file_name_base_string;
        output_file_name_three_step_correlation_2x2 = "three_step_corr_2x2_" + output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        output_file_name_four_step_correlation_1x1 = "four_step_corr_1x1_" + output_file_name_base_string;
        output_file_name_four_step_correlation_2x2 = "four_step_corr_2x2_" + output_file_name_base_string;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        std::string header_string;
        header_string = "Hello!";
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // create a wrapper
        Wrapper * wrapper = new Wrapper(number_of_lattices,
                                        grid_size_x,
                                        grid_size_y,
                                        number_of_timesteps,
                                        number_of_tracers_1x1,
                                        number_of_tracers_2x2,
                                        number_of_tracers_3x3,
                                        step_rate_1x1,
                                        step_rate_2x2,
                                        step_rate_3x3,
                                        max_wtd);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // perform warm up if number_of_timesteps_warmup
        for(int t = 0; t < number_of_timesteps_warmup; t++) { wrapper->timestep_warmup(); }
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        // let the simulation run for number_of_timesteps iterations
        int progress_percent = 0;
        for(int t = 0; t < number_of_timesteps; t++)
        {
                if((int)(1000*t/number_of_timesteps) > progress_percent) {
                        progress_percent = (double)(100*t)/(double)number_of_timesteps;
                        fflush(stdout);
                        printf("RUN:%f4.1", progress_percent);
                }
                wrapper->timestep();
        }
        printf("RUN:%f4.1", progress_percent);
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        if(number_of_tracers_1x1)
        {
                vector_to_file(wrapper->get_avg_msd_1x1(),
                               output_file_name_avg_msd_1x1);
                vector_to_file(wrapper->get_avg_stepping_rate_1x1(),
                               output_file_name_avg_rate_1x1);
                vector_to_file(wrapper->get_avg_two_step_correlation_1x1(),
                               output_file_name_two_step_correlation_1x1);
                vector_to_file(wrapper->get_avg_three_step_correlation_1x1(),
                               output_file_name_three_step_correlation_1x1);
                vector_to_file(wrapper->get_avg_four_step_correlation_1x1(),
                               output_file_name_four_step_correlation_1x1);
        }
        // - - - - - - - - - - - - - - - - - - - - - - - - -
        if(number_of_tracers_2x2)
        {
                vector_to_file(wrapper->get_avg_msd_2x2(),
                               output_file_name_avg_msd_2x2);
                vector_to_file(wrapper->get_avg_stepping_rate_2x2(),
                               output_file_name_avg_rate_2x2);
                vector_to_file(wrapper->get_avg_two_step_correlation_2x2(),
                               output_file_name_two_step_correlation_2x2);
                vector_to_file(wrapper->get_avg_three_step_correlation_2x2(),
                               output_file_name_three_step_correlation_2x2);
                vector_to_file(wrapper->get_avg_four_step_correlation_2x2(),
                               output_file_name_four_step_correlation_2x2);
        }
        std::cout << "DONE!" << std::endl;
        // - - - - - - - - - - - - - - - - - - - - - - - - -
}
