#include "global.hpp"
#include "tracer.hpp"
#include "output.hpp"
#include <iostream>
#include <vector>
#include <boost/format.hpp>
#include <string>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>


int main(int ac, char** av){
        int max_t = atoi(av[1]);
        int seed = atoi(av[2]);
        generator.seed(seed);
        std::vector<int> pos_vector = {0,0};
        Tracer tr = Tracer (0,0);
        for (int t = 1; t < max_t; t++) {
                tr.step();
                pos_vector.push_back(tr.x);
                pos_vector.push_back(tr.y);
        }
        vector_to_file(pos_vector, boost::str(boost::format("../data/random_walk_t%i_%i.bin") % max_t % seed));
}
