#include "global.hpp"
#include "tracer.hpp"
#include "output.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <boost/format.hpp>
#include <string>

#define STORE_POSITION false

std::uniform_int_distribution<> step_directions(1,4);

int coord(int x, int y){
        return x+y*grid_size;
}

int select_direction(){
        return step_directions(generator);
}

int main(int ac, char** av){
        // sim parameters
        // duration of the simulation in time_steps
        int max_t = atoi(av[1]);
        // number of tracers on the lattice (lattice size scales with n), number of traps
        int n_tracers = atoi(av[2]);
        int n_traps;
        // desired concentrations of tracers + traps
        double c_traps = atof(av[3]);
        // trapping and release probabilities
        double p_on = atof(av[4]);
        dobule one_minus_p_on = 1.0 - p_on;
        // trapping probabilites for tracers touching multiple traps
        double p_on_2 = 2*p_on*one_minus_p_on + p_on*p_on;
        double p_on_3 = 3*p_on*one_minus_p_on*one_minus_p_on + 3*p_on*p_on*one_minus_p_on + p_on*p_on*p_on;
        double p_on_4 = 4*p_on*one_minus_p_on*one_minus_p_on*one_minus_p_on + 6*p_on*p_on*one_minus_p_on*one_minus_p_on + 4*p_on*p_on*p_on*one_minus_p_on + p_on*p_on*p_on*p_on;
        double p_off = atof(av[5]);
        // seed for rng
        int seed = atoi(av[6]);
        // int grid_size is declared as extern in GLOBAL_H and global.cpp
        // grid size scales with particle number
        // higher particle number -> more exact match to desired c_tracers
        grid_size = round(sqrt((double) n_tracers / c_tracers));
        // we want
        n_traps = round(c_traps * grid_size * grid_size);
        std::uniform_int_distribution<> lattice_position(0,grid_size-1);
        std::uniform_real_distribution<> dist(0.0,1.0);
        // set the rng seed for reproducible results
        // warm up rng by running it 10^4 times
        set_rng_seed(seed);
        generator.discard(9999);
        //
        std::vector<double> lsq_vector = {0};
        std::vector<Tracer*> tracers;
        // storage of grid state and trap positions
        // in this model, tracers don't exclude each other from lattice sites, so we dont need the grid_vector
        // * * * *
        // std::vector<int> grid_vector (grid_size*grid_size,0);
        // * * * *
        std::vector<int> trap_vector (grid_size*grid_size,0);
        // storage of tracer positions for output, if deisred
        std::vector<int> pos_vector;
        // populate the lattice with the desired number of traps
        for (int n = 0; n < n_traps; n++) {
                int x = lattice_position(generator);
                int y = lattice_position(generator);
                trap_vector[coord(x,y)] = 1;
        }
        // populate the lattice with the desired number of tracers
        for (int n = 0; n < n_tracers; n++) {
                int x = lattice_position(generator);
                int y = lattice_position(generator);
                tracers.push_back(new Tracer(x,y));
                // grid_vector[coord(x,y)] = 1;
                if(STORE_POSITION) {
                        pos_vector.push_back(x);
                        pos_vector.push_back(y);
                }
        }
        for (Tracer* tr : tracers) {
                // 2x2 tracers have surface_area of 2
                // they map to the lattice sites at:
                // x+0, x+1, x+0, x+1
                // y+0, y+0, y+1, y+1
                // for now, this way of modelling is the best idea i have
                // tr->surface_area = 4;
                // tr->footprint = {0,0,1,0,0,1,1,1};
                // some tracers might start out on a trap site - or rather touching a trap site
                // check for those and set their state to stuck with p_on probability
                int hits = trap_vector[coord(tr->x,tr->y)] + trap_vector[coord(tr->x+1,tr->y)] + trap_vector[coord(tr->x,tr->y+1)] + trap_vector[coord(tr->x+1,tr->y+1)];
                switch(hits) {
                case 1: double d = dist(generator); if (d < p_on) { tr->stuck(); } break;
                case 2: double d = dist(generator); if (d < p_on_2) { tr->stuck(); } break;
                case 3: double d = dist(generator); if (d < p_on_3) { tr->stuck(); } break;
                case 4: double d = dist(generator); if (d < p_on_4) { tr->stuck(); } break;
                default: break;
                }
                // maybe we should allow for immediate unstucking, if the tracer is lucky?
                // if (d < p_on - p_off) {....}
                // * * * *



        }
        // loop over time range
        for (int t = 1; t < max_t; t++) {
                // temporary
                double lsq_average = 0.0;
                std::vector<double> lsq_vector_temp;
                // loop over tracers
                for (Tracer * tr : tracers) {
                        // time-step logic:
                        if(tr->is_stuck) {
                                // if tracer is stuck:
                                // unstuck with probability p_off
                                double d = dist(generator);
                                if(d < p_off) {
                                        tr->unstuck();
                                }
                        }
                        // do newly unstuck tracers get to walk in the same timestep?
                        // in this version of the code, they do
                        // TODO: check back with Ana, if this behaviour is the desired one.
                        if (!tr->is_stuck) {
                                // if tracer is not stuck:
                                // randomly select direction
                                int dir = select_direction();
                                int dir_x, dir_y;
                                switch(dir) {
                                case 1: dir_x = 1; dir_y = 0; break;
                                case 2: dir_x = 0; dir_y = 1; break;
                                case 3: dir_x =-1; dir_y = 0; break;
                                case 4: dir_x = 0; dir_y =-1; break;
                                }
                                // update
                                int new_x = tr->x+dir_x;
                                int new_x_mod = (grid_size+new_x%grid_size)%grid_size;
                                int new_y = tr->y+dir_y;
                                int new_y_mod = (grid_size+new_y%grid_size)%grid_size;
                                tr->x=new_x_mod;
                                tr->y=new_y_mod;
                                // if the periodic boundary in x or y was crossed, increment passage counters
                                if(new_x == grid_size) {
                                        tr->increment_passages_x();
                                } else if (new_x == -1) {
                                        tr->decrement_passages_x();
                                }
                                if(new_y == grid_size) {
                                        tr->increment_passages_y();
                                }else if (new_y == -1) {
                                        tr->decrement_passages_y();
                                }
                                // TODO: rework if to swith/case statement
                                // switch(new_x) {
                                // case grid_size: tr->increment_passages_x(); break;
                                // case -1: tr->decrement_passages_x(); break;
                                // default: break;
                                // }
                                // switch(new_y) {
                                // case grid_size: tr->increment_passages_y(); break;
                                // case -1: tr->decrement_passages_y(); break;
                                // default: break;
                                // }
                                //

                                // if new location is trap
                                if(trap_vector[coord(tr->x, tr->y)] == 1) {
                                        double d = dist(generator);
                                        if (d < p_on) {
                                                tr->stuck();
                                        }
                                }
                        }
                        // statistics computation and storage
                        // calculate lsquared
                        tr->calc_lsquared();
                        // store lsq for averaging
                        lsq_vector_temp.push_back(tr->l_squared);
                        // store tracer position
                        if(STORE_POSITION) {
                                pos_vector.push_back(tr->x);
                                pos_vector.push_back(tr->y);
                        }
                }
                // compute average lsqaured of all tracers
                for (double lsq : lsq_vector_temp) {
                        lsq_average += lsq;
                }
                lsq_average /= lsq_vector_temp.size();
                // push back average lsquared to vector for output
                lsq_vector.push_back(lsq_average);
        }
        c_tracers = (double)n_tracers/(double)grid_size/(double)grid_size;
        c_traps = (double)n_traps/(double)grid_size/(double)grid_size;
        if(STORE_POSITION) {
                vector_to_file(pos_vector, boost::str(boost::format("./data/lgtr_positions_ctraps=%f_ctracers=%f_pon=%f_poff=%f_ntracers=%i_tmax=%i_seed=%i.bin") % c_traps % c_tracers %p_on % p_off % n_tracers % max_t % seed));
        }

        vector_to_file(lsq_vector, boost::str(boost::format("./data/lgtr_statistics_ctraps=%f_ctracers=%f_pon=%f_poff=%f_ntracers=%i_tmax=%i_seed=%i.bin") % c_traps % c_tracers %p_on % p_off % n_tracers % max_t % seed));
}
