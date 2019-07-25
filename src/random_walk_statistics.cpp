#include "global.hpp"
#include "tracer.hpp"
#include "output.hpp"
#include <iostream>
#include <vector>
#include <boost/format.hpp>
#include <string>

int main(int ac, char** av){
        int max_t = atoi(av[1]);
        int n_walkers = atoi(av[2]);
        int seed = atoi(av[3]);
        set_rng_seed(seed);
        std::vector<double> lsq_vector = {0};
        std::vector<Tracer*> tracers;
        // create a vector containing the desired number of tracers
        for (int n = 0; n < n_walkers; n++) {
                tracers.push_back(new Tracer(0,0));
        }
        // iterate through all tracers, do for t_max timesteps
        for (int t = 1; t < max_t; t++) {
                double lsq_average;
                std::vector<double> lsq_vector_temp;
                for(Tracer * tr : tracers) {
                        // evolve tracer and push back lsquared for later use
                        tr->step();
                        lsq_vector_temp.push_back(tr->l_squared);
                }
                // compute average lsqaured of all tracers
                for (double lsq : lsq_vector_temp) {
                        lsq_average += lsq;
                }
                lsq_average /= lsq_vector_temp.size();
                // push back average lsquared to vector for output
                lsq_vector.push_back(lsq_average);
        }
        vector_to_file(lsq_vector, boost::str(boost::format("../data/random_walk_statistics_n%i_t%i_%i.bin") %n_walkers % max_t % seed));
}
