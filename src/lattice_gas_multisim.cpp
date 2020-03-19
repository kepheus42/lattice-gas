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
        int number_of_lattices = atoi(av[8]);
        // TODO: find some reasonable maximum number depth for the WTDs
        // * note: at 95% concentration 1x1 tracers, after 1 million timesteps, the depth is only t~800
        // * note: 2x2 tracers might require longer wtd storage
        // * truncate at some high number (e.g. 1k timesteps), should easily be enough to infer shape
        int max_wtd = atoi(av[9]);
        // - - - - - - - - - - - - -
        std::cout << "X = " << grid_size_x << " Y = " << grid_size_y << " N1 = " << number_of_tracers_1x1 << " N2 = " << number_of_tracers_2x2 << " N3 = " << number_of_tracers_3x3 << " T = " << max_number_of_timesteps << '\n';
        // - - - - - - - - - - - - -
        std::vector<int> output_wtd_1_tmp;
        std::vector<int> output_wtd_2_tmp;
        std::vector<int> output_wtd_3_tmp;
        std::vector<double> output_msd;
        std::vector<double> output_rates;
        output_rates.reserve(3);
        // - - - - - - - - - - - - -
        std::string output_file_name_base_string;
        std::string output_file_name_msd;
        std::string output_file_name_rates;
        std::string output_file_name_rate_distribution_1x1;
        std::string output_file_name_rate_distribution_2x2;
        std::string output_file_name_rate_distribution_3x3;
        std::string output_file_name_wtd_1;
        std::string output_file_name_wtd_2;
        std::string output_file_name_wtd_3;
        std::string header_string;
        // - - - - - - - - - - - - -
        output_file_name_base_string =  "_X_"+std::to_string(grid_size_x)
                                       +"_Y_"+std::to_string(grid_size_y)
                                       +"_N1_"+std::to_string(number_of_tracers_1x1)
                                       +"_N2_"+std::to_string(number_of_tracers_2x2)
                                       +"_N3_"+std::to_string(number_of_tracers_3x3)
                                       +"_T_"+std::to_string(max_number_of_timesteps)
                                       +".bin";
        //
        output_file_name_msd   = "multisim_msd"+output_file_name_base_string;
        output_file_name_rates = "multisim_rates"+output_file_name_base_string;
        output_file_name_rate_distribution_1x1 = "multisim_rate_dist_1x1"+output_file_name_base_string;
        output_file_name_rate_distribution_2x2 = "multisim_rate_dist_2x2"+output_file_name_base_string;
        output_file_name_rate_distribution_3x3 = "multisim_rate_dist_3x3"+output_file_name_base_string;
        output_file_name_wtd_1 = "multisim_wtd_1x1"+output_file_name_base_string;
        output_file_name_wtd_2 = "multisim_wtd_2x2"+output_file_name_base_string;
        output_file_name_wtd_3 = "multisim_wtd_3x3"+output_file_name_base_string;
        header_string          = "parameters:";
        // - - - - - - - - - - - - -
        Wrapper * wrapper = new Wrapper(number_of_lattices, grid_size_x, grid_size_y, number_of_tracers_1x1, number_of_tracers_2x2, number_of_tracers_3x3);
        wrapper->init_wtd(max_wtd);
        // - - - - - - - - - - - - -
        int progress_percent = 0;
        // run warm up to disperse 2x2 and 3x3 tracers
        wrapper->warm_up(timesteps_warmup);
        //
        while(wrapper->get_t() < max_number_of_timesteps)
        {
                if((int)(100*wrapper->get_t()/max_number_of_timesteps) > progress_percent) {
                        progress_percent = (int)(100*wrapper->get_t()/max_number_of_timesteps);
                        fflush(stdout);
                        std::cout << " " << progress_percent << "%\r";
                }
                // save msd vs t
                output_msd.push_back((double)wrapper->get_t());
                output_msd.push_back(wrapper->get_avg_lsquared_1x1());
                output_msd.push_back(wrapper->get_avg_lsquared_2x2());
                output_msd.push_back(wrapper->get_avg_lsquared_3x3());
                // timestep
                wrapper->timestep();
                //
        }
        std::cout << "100%" << std::endl;
        // save msd vs t
        output_msd.push_back((double)wrapper->get_t());
        output_msd.push_back(wrapper->get_avg_lsquared_1x1());
        output_msd.push_back(wrapper->get_avg_lsquared_2x2());
        output_msd.push_back(wrapper->get_avg_lsquared_3x3());
        //
        output_rates.push_back(wrapper->get_avg_stepping_rate_1x1());
        output_rates.push_back(wrapper->get_avg_stepping_rate_2x2());
        output_rates.push_back(wrapper->get_avg_stepping_rate_3x3());
        // write msd, rates, and wtds to disk
        vector_to_file(output_msd,output_file_name_msd);
        vector_to_file(output_rates,output_file_name_rates);
        // 3 wtds
        // for non-normalized waiting time distributions, use:
        // - - - - - - - - -
        // vector_to_file(wrapper->get_wtd_1x1(),output_file_name_wtd_1);
        // vector_to_file(wrapper->get_wtd_2x2(),output_file_name_wtd_2);
        // vector_to_file(wrapper->get_wtd_3x3(),output_file_name_wtd_3);
        // - - - - - - - - -
        // save WTDs only if they were actually recorded!
        if(number_of_tracers_1x1 > 0)
        {
                vector_to_file(wrapper->get_norm_wtd_1x1(),output_file_name_wtd_1);
                vector_to_file(wrapper->get_stepping_rate_distribution_1x1(),output_file_name_rate_distribution_1x1);
        }
        if(number_of_tracers_2x2 > 0)
        {
                vector_to_file(wrapper->get_norm_wtd_2x2(),output_file_name_wtd_2);
                vector_to_file(wrapper->get_stepping_rate_distribution_2x2(),output_file_name_rate_distribution_2x2);
        }
        if(number_of_tracers_3x3 > 0)
        {
                vector_to_file(wrapper->get_norm_wtd_3x3(),output_file_name_wtd_3);
                vector_to_file(wrapper->get_stepping_rate_distribution_3x3(),output_file_name_rate_distribution_3x3);
        }
}
