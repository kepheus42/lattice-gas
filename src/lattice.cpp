#include <stdexcept>
#include <iostream>
#include <cmath>
#include "lattice.hpp"
#include "global.hpp"
// - - - - - - - - - - - - - - - - - - - - - - - -
Lattice::Lattice(int grid_size_x,
                 int grid_size_y,
                 int number_of_timesteps,
                 int number_of_tracers_1x1,
                 int number_of_tracers_2x2,
                 int number_of_tracers_3x3,
                 double step_rate_1x1,
                 double step_rate_2x2,
                 double step_rate_3x3,
                 int wtd_max,
                 int wtd_res) :
        m_t(0),
        m_grid_size_x(grid_size_x),
        m_grid_size_y(grid_size_y),
        m_number_of_timesteps(number_of_timesteps),
        m_number_of_tracers_1x1(number_of_tracers_1x1),
        m_number_of_tracers_2x2(number_of_tracers_2x2),
        m_number_of_tracers_3x3(number_of_tracers_3x3),
        m_number_of_tracers_total(number_of_tracers_1x1+number_of_tracers_2x2+number_of_tracers_3x3),
        m_step_rate_1x1(number_of_tracers_1x1 > 0 ? step_rate_1x1 : 1.0),
        m_step_rate_2x2(number_of_tracers_2x2 > 0 ? step_rate_2x2 : 1.0),
        m_step_rate_3x3(number_of_tracers_3x3 > 0 ? step_rate_3x3 : 1.0),
        m_movement_selector_list((int)(this->m_number_of_tracers_1x1/this->m_step_rate_2x2/this->m_step_rate_3x3)+
                                 (int)(this->m_number_of_tracers_2x2/this->m_step_rate_1x1/this->m_step_rate_3x3)+
                                 (int)(this->m_number_of_tracers_3x3/this->m_step_rate_1x1/this->m_step_rate_2x2),0),
        m_movement_selector_length(this->m_movement_selector_list.size()),
        m_step_attempts_per_timestep((int)(this->m_number_of_tracers_1x1*this->m_step_rate_1x1)+
                                     (int)(this->m_number_of_tracers_2x2*this->m_step_rate_2x2)+
                                     (int)(this->m_number_of_tracers_3x3*this->m_step_rate_3x3)),
        m_t_increment(1.0/(double)this->m_step_attempts_per_timestep),
        m_occupation_map(grid_size_x*grid_size_y,0),
        m_wtd_max(wtd_max),
        m_wtd_res(wtd_res)
{
        // check if there's enough space on the grid to place the tracers
        if (grid_size_x*grid_size_y < number_of_tracers_1x1 + 4*number_of_tracers_2x2 + 9*number_of_tracers_3x3)
        {
                throw std::invalid_argument("Too many tracers for grid of the given size!");
        }
        else if (((number_of_tracers_2x2 > 0)&&((grid_size_x < 2)||(grid_size_y < 2)))||((number_of_tracers_3x3 > 0)&&((grid_size_x < 3)||(grid_size_y < 3))))
        {
                throw std::invalid_argument("Grid is too narrow for 2x2 or 3x3 tracers!");
        }
        this->setup_movement_selection_list();
        this->setup_tracers();
}
// - - - - - - - - - - - - - - - - - - - - - - - -
void setup_movement_selection_list()
{
        //
        // structure of the selection list vector:
        //  T Y P E 1               T Y P E 2
        // 0 0 0 1 1 1 2 2 2 ... 10 10 11 11 12 12
        // in this case, rate_1/rate_2 = 3/2

        // use the movement selector list to randomly pick particles to move, with probabilities proportional to the ratios of the stepping rates and numbers of particles
        int tmp_idx = 0;
        // The first 1x1 tracer is this->Tracers[n_3x3 + n_2x2], the first 2x2 tracer is this->Tracers[n_3x3]
        int tmp_idx_1x1 = this->m_number_of_tracers_3x3 + this->m_number_of_tracers_2x2;
        int tmp_idx_2x2 = this->m_number_of_tracers_3x3;
        // tmp_m1 - tmp_m3 : how many times each tracer idx for tracers of types 1-3 appears in the list
        int tmp_m1 = (int)(1/this->m_step_rate_2x2/this->m_step_rate_3x3);
        int tmp_m2 = (int)(1/this->m_step_rate_1x1/this->m_step_rate_3x3);
        int tmp_m3 = (int)(1/this->m_step_rate_1x1/this->m_step_rate_2x2);
        // tmp_n1 - tmp_n3 : how many list entires do types 1 - 3 get in the final list
        int tmp_n1 = (int)(this->m_number_of_tracers_1x1/this->m_step_rate_2x2/this->m_step_rate_3x3);
        int tmp_n2 = (int)(this->m_number_of_tracers_2x2/this->m_step_rate_1x1/this->m_step_rate_3x3);
        int tmp_n3 = (int)(this->m_number_of_tracers_3x3/this->m_step_rate_1x1/this->m_step_rate_2x2);
        for(int n = 0; n < tmp_n3; n++)
        {
                this->m_movement_selector_list[idx] = n/tmp_m3;
                tmp_idx++;
        }
        for(int n = 0; n < tmp_n2; n++)
        {
                this->m_movement_selector_list[idx] = tmp_idx_2x2 + n/tmp_m2;
                tmp_idx++;
        }
        for(int n = 0; n < tmp_n1; n++)
        {
                this->m_movement_selector_list[idx] =  tmp_idx_1x1 + n/tmp_m1;
                tmp_idx++;
        }
}
// - - - - - - - - - - - - - - - - - - - - - - - -
void Lattice::setup_tracers()
{
        // create a vector to store pointers to all the tracers
        this->m_tracers.reserve(this->m_number_of_tracers_1x1 + this->m_number_of_tracers_2x2 + this->m_number_of_tracers_3x3);
        // tmp vector to handle start position allocation for 3x3 tracers
        int tmp_id = 1;
        int tmp_x;
        int tmp_y;
        if(this->m_number_of_tracers_3x3)
        {
                this->m_tracers_3x3.reserve(this->m_number_of_tracers_3x3);
                std::vector<int> tmp_start_positions_3x3;
                tmp_start_positions_3x3.reserve(grid_size_x/3 * grid_size_y/3);
                for (int y=0; y < this->m_grid_size_y; y+=3)
                {
                        for (int x=0; x < this->m_grid_size_x; x+=3)
                        {
                                tmp_start_positions_3x3.push_back(x * grid_size_y + y);
                        }
                }
                // create the 3x3 tracers:
                shuffle_vector(tmp_start_positions_3x3);
                for (int n = 0; n < this->m_number_of_tracers_3x3; n++)
                {
                        tmp_y = tmp_start_positions_3x3[n] % this->m_grid_size_y;
                        tmp_x = (tmp_start_positions_3x3[n] - tmp_y)/this->m_grid_size_y;
                        this->m_tracers.push_back(new Tracer_3x3(tmp_id,tmp_x,tmp_y,this->m_grid_size_x,this->m_grid_size_y,this->m_step_rate_3x3));
                        this->m_tracers_3x3.push_back(this->m_tracers.back());
                        this->m_occupation_map[this->coord(tmp_x+0,tmp_y+0)] = tmp_id;
                        this->m_occupation_map[this->coord(tmp_x+1,tmp_y+0)] = tmp_id;
                        this->m_occupation_map[this->coord(tmp_x+2,tmp_y+0)] = tmp_id;
                        this->m_occupation_map[this->coord(tmp_x+0,tmp_y+1)] = tmp_id;
                        this->m_occupation_map[this->coord(tmp_x+1,tmp_y+1)] = tmp_id;
                        this->m_occupation_map[this->coord(tmp_x+2,tmp_y+1)] = tmp_id;
                        this->m_occupation_map[this->coord(tmp_x+0,tmp_y+2)] = tmp_id;
                        this->m_occupation_map[this->coord(tmp_x+1,tmp_y+2)] = tmp_id;
                        this->m_occupation_map[this->coord(tmp_x+2,tmp_y+2)] = tmp_id;
                        tmp_id++;
                }
        }
        if(this->m_number_of_tracers_2x2)
        {
                this->m_tracers_2x2.reserve(this->m_number_of_tracers_2x2);
                // tmp vector to handle start position allocation for 2x2 tracers
                std::vector<int> tmp_start_positions_2x2;
                tmp_start_positions_2x2.reserve(this->m_grid_size_x/2 * this->m_grid_size_y/2);
                for (int y=0; y < this->m_grid_size_y; y+=2)
                {
                        for (int x=0; x < this->m_grid_size_x; x+=2)
                        {
                                // TODO: add check for occupancy
                                tmp_start_positions_2x2.push_back(this->coord(x,y));
                        }
                }
                // create the 2x2 tracers:
                shuffle_vector(tmp_start_positions_2x2);
                for (int n = 0; n < number_of_tracers_2x2; n++)
                {
                        tmp_y = tmp_start_positions_2x2[n] % this->m_grid_size_y;
                        tmp_x = (tmp_start_positions_2x2[n] - tmp_y)/this->m_grid_size_y;
                        this->m_tracers.push_back(new Tracer_2x2(tmp_id,tmp_x,tmp_y,this->m_grid_size_x,this->m_grid_size_y,this->m_step_rate_2x2));
                        this->m_tracers_2x2.push_back(this->m_tracers.back());
                        this->m_occupation_map[this->coord(tmp_x+0,tmp_y+0)] = tmp_id;
                        this->m_occupation_map[this->coord(tmp_x+1,tmp_y+0)] = tmp_id;
                        this->m_occupation_map[this->coord(tmp_x+0,tmp_y+1)] = tmp_id;
                        this->m_occupation_map[this->coord(tmp_x+1,tmp_y+1)] = tmp_id;
                        tmp_id++;
                }
        }
        if(this->m_number_of_tracers_1x1)
        {
                this->m_tracers_1x1.reserve(this->m_number_of_tracers_1x1);
                // tmp vector to handle creation of the 1x1 tracers
                std::vector<int> tmp_start_positions_1x1;
                tmp_start_positions_1x1.reserve(this->m_grid_size_x * this->m_grid_size_y - 4 * this->m_number_of_tracers_2x2 - 9 * this->m_number_of_tracers_3x3);
                for (int i = 0; i < this->m_grid_size_x*this->m_grid_size_y; i++)
                {
                        if(!this->m_occupation_map[i])
                        {
                                tmp_start_positions_1x1.push_back(i);
                        }
                }
                shuffle_vector(tmp_start_positions_1x1);
                for (int n = 0; n < this->m_number_of_tracers_1x1; n++)
                {
                        tmp_y = tmp_start_positions_1x1[n] % this->m_grid_size_y;
                        tmp_x = (tmp_start_positions_1x1[n] - tmp_y)/this->m_grid_size_y;
                        this->m_tracers.push_back(new Tracer(tmp_id,tmp_x,tmp_y,this->m_grid_size_x,this->m_grid_size_y,this->m_step_rate_1x1));
                        this->m_tracers_1x1.push_back(this->m_tracers.back());
                        this->m_occupation_map[this->coord(tmp_x,tmp_y)] = tmp_id;
                        tmp_id++;
                }
        }
}
// - - - - - - - - - - - - - - - - - - - - - - - -
void Lattice::update_wtd()
{
        for(Tracer * tr : this->m_tracers_1x1)
        {
                int idx = tr->get_wtd_index();
                if(idx > 4*this->m_wtd_max) {
                        idx = (idx-1)%4+4*this->m_wtd_max+1;
                }
                if(idx != 0) {
                        this->m_wtd_1x1[idx-1]++;
                }
        }
        for(Tracer * tr : this->m_tracers_2x2)
        {
                int idx = tr->get_wtd_index();
                if(idx > 4*this->m_wtd_max) {
                        idx = (idx-1)%4+4*this->m_wtd_max+1;
                }
                if(idx != 0) {
                        this->m_wtd_2x2[idx-1]++;
                }
        }
        for(Tracer * tr : this->m_tracers_3x3)
        {
                int idx = tr->get_wtd_index();
                if(idx > 4*this->m_wtd_max) {
                        idx = (idx-1)%4+4*this->m_wtd_max+1;
                }
                if(idx != 0) {
                        this->m_wtd_3x3[idx-1]++;
                }
        }
}
// = = = = = = = = = = = = = = = = = = = = = =
// = = = M E M B E R   F U N C T I O N S = = =
// = = = = = = = = = = = = = = = = = = = = = =
inline int Lattice::coord(int x, int y){
        //
        return ((x+this->m_grid_size_x)%this->m_grid_size_x)*this->m_grid_size_y+((y+this->m_grid_size_y)%this->m_grid_size_y);
}
// - - - - - - - - - - - - - - - - - - - - - - - -
void Lattice::warm_up(){
        bool tmp_do_nothing = true;
}
//
void Lattice::timestep(){
        int tmp_rnd; // random number between 0 and 4*length of selector list
        int tmp_par; // which particle moves
        int tmp_dir; // which direction it moves in
        double tmp_time = (double)this->m_t;
        for (int n=0; n < this->m_step_attempts_per_timestep; n++)
        {
                tmp_rnd = random_int(0,4*this->m_movement_selector_length);
                tmp_par = this->m_movement_selector[std::floor(tmp_rnd/4)];
                tmp_dir = 1+tmp_rnd%4;
                // eacth timestep consists of m_step_attempts_per_timestep equal intervals, thus the time at which a given step attempt takes place is t = m_t + n / m_step_attempts_per_timestep
                this->m_tracers[tmp_par]->step(this->m_occupation_map,tmp_dir,tmp_time);
                tmp_time += this->m_t_increment;
        }
        this->m_t++;
}
//
void Lattice::timestep_warmup(){
        int tmp_rnd; // random number between 0 and 4*length of selector list
        int tmp_par; // which particle moves
        int tmp_dir; // which direction it moves in
        for (int n=0; n < this->m_step_attempts_per_timestep; n++)
        {
                tmp_rnd = random_int(0,4*this->m_movement_selector_length);
                tmp_par = this->m_movement_selector[std::floor(tmp_rnd/4)];
                tmp_dir = 1+tmp_rnd%4;
                this->m_tracers[tmp_par]->step_warmup(this->m_occupation_map,tmp_dir);
        }
        this->m_t++;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
void Lattice::timestep_no_interaction(){
        int tmp_rnd; // random number between 0 and 4*length of selector list
        int tmp_par; // which particle moves
        int tmp_dir; // which direction it moves in
        double tmp_time = (double)this->m_t;
        for (int n=0; n < this->m_step_attempts_per_timestep; n++)
        {
                tmp_rnd = random_int(0,4*this->m_movement_selector_length);
                tmp_par = this->m_movement_selector[std::floor(tmp_rnd/4)];
                tmp_dir = 1+tmp_rnd%4;
                this->m_tracers[tmp_par]->unhindered_step(tmp_dir,tmp_time);
                tmp_time += this->m_t_increment;
        }
        this->m_t++;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
// G E T T E R S
// to obtain the values of member properties of the lattice class
// - - - - - - - - - - - - - - - - - - - - - - - -
int Lattice::get_t(){
        return this->m_t;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
int Lattice::get_number_of_tracers_1x1()
{
        return this->m_number_of_tracers_1x1;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
int Lattice::get_number_of_tracers_2x2()
{
        return this->m_number_of_tracers_2x2;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
int Lattice::get_number_of_tracers_3x3()
{
        return this->m_number_of_tracers_3x3;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
int Lattice::get_grid_size_x()
{
        return this->m_grid_size_x;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
int Lattice::get_grid_size_y()
{
        return this->m_grid_size_y;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> Lattice::get_occupation_map()
{
        return this->m_occupation_map;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
int Lattice::get_tracer_size(int id)
{
        return this->m_tracers[id]->get_size();
}
// get the last move of a tracer (0,1,2,3,4)
int Lattice::get_tracer_last_move(int id)
{
        return this->m_tracers[id]->get_last_move();
}
// get dx or dy of a tracer
int Lattice::get_tracer_position_x(int id)
{
        return this->m_tracers[id]->get_dx();
}
// - - - - - - - - - - - - - - - - - - - - - - - -
int Lattice::get_tracer_position_y(int id)
{
        return this->m_tracers[id]->get_dy();
}
// - - - - - - - - - - - - - - - - - - - - - - - -
// get vectors containing pointers to all tracers for wrapper class to operate on them
std::vector<Tracer *> Lattice::get_tracers()
{
        return this->m_tracers;
}
std::vector<Tracer *> Lattice::get_tracers_1x1()
{
        return this->m_tracers_1x1;
}
std::vector<Tracer *> Lattice::get_tracers_2x2()
{
        return this->m_tracers_2x2;
}
std::vector<Tracer *> Lattice::get_tracers_3x3()
{
        return this->m_tracers_3x3;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
// instantaneous average stepping rates for the different tracer types
// - - - - - - - - - - - - - - - - - - - - - - - -
double Lattice::get_avg_rate_1x1()
{
        // returns the current average step rate
        double tmp_sum_steps_taken = 0;
        if(!this->m_number_of_tracers_1x1) { return tmp_sum_steps_taken; }
        for(Tracer * tr : this->m_tracers_1x1)
        {
                tmp_sum_steps_taken += (double)tr->get_last_step();
        }
        // if number of tracers, return average, else return 0
        return(tmp_sum_steps_taken/(double)this->m_number_of_tracers_1x1);
}
// - - - - - - - - - - - - - - - - - - - - - - - -
double Lattice::get_avg_rate_2x2()
{
        // returns the current average step rate
        double tmp_sum_steps_taken = 0;
        if(!this->m_number_of_tracers_2x2) { return tmp_sum_steps_taken; }
        for(Tracer * tr : this->m_tracers_2x2)
        {
                tmp_sum_steps_taken += (double)tr->get_last_step();
        }
        // if number of tracers, return average, else return 0
        return(tmp_sum_steps_taken/(double)this->m_number_of_tracers_2x2);
}
// - - - - - - - - - - - - - - - - - - - - - - - -
double Lattice::get_avg_rate_3x3()
{
        // returns the current average step rate
        double tmp_sum_steps_taken = 0;
        if(!this->m_number_of_tracers_3x3) { return tmp_sum_steps_taken; }
        for(Tracer * tr : this->m_tracers_3x3)
        {
                tmp_sum_steps_taken += (double)tr->get_last_step();
        }
        // if number of tracers, return average, else return 0
        return(tmp_sum_steps_taken/(double)this->m_number_of_tracers_3x3);
}
// - - - - - - - - - - - - - - - - - - - - - - - -
double Lattice::get_avg_lsquared_1x1()
{
        // returns the current average step rate
        double tmp_sum_lsquared = 0;
        if(!this->m_number_of_tracers_1x1) { return tmp_sum_lsquared; }
        for(Tracer * tr : this->m_tracers_1x1)
        {
                tmp_sum_lsquared += tr->get_lsquared();
        }
        // if number of tracers, return average, else return 0
        return(tmp_sum_lsquared/(double)this->m_number_of_tracers_1x1);
}
// - - - - - - - - - - - - - - - - - - - - - - - -
double Lattice::get_avg_lsquared_2x2()
{
        // returns the current average step rate
        double tmp_sum_lsquared = 0;
        if(!this->m_number_of_tracers_2x2) { return tmp_sum_lsquared; }
        for(Tracer * tr : this->m_tracers_2x2)
        {
                tmp_sum_lsquared += tr->get_lsquared();
        }
        // if number of tracers, return average, else return 0
        return(tmp_sum_lsquared/(double)this->m_number_of_tracers_2x2);
}
// - - - - - - - - - - - - - - - - - - - - - - - -
double Lattice::get_avg_lsquared_3x3()
{
        // returns the current average step rate
        double tmp_sum_lsquared = 0;
        if(!this->m_number_of_tracers_3x3) { return tmp_sum_lsquared; }
        for(Tracer * tr : this->m_tracers_3x3)
        {
                tmp_sum_lsquared += tr->get_lsquared();
        }
        // if number of tracers, return average, else return 0
        return(tmp_sum_lsquared/(double)this->m_number_of_tracers_3x3);
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<unsigned int> Lattice::get_wtd_1x1()
{
        std::vector<unsigned int> tmp_wtd(4+16*(this->m_wtd_max*this->m_wtd_res+1),0);

        int tmp_wtd_idx = 0;
        for(Tracer * tr : this->m_tracers_1x1)
        {
                if(tr->get_last_step_wtd_idx())
                {
                        tmp_wtd_idx = tr->get_last_step_wtd_idx()+(int)(this->m_wtd_res*tr->get_time_since_last_step());
                        tmp_wtd[std::min(tmp_wtd_idx,tmp_wtd.size()-1)]++;
                }
        }
        return tmp_wtd;
}
std::vector<unsigned int> Lattice::get_wtd_2x2()
{
        std::vector<unsigned int> tmp_wtd(4+16*(this->m_wtd_max*this->m_wtd_res+1),0);
        for(Tracer * tr : this->m_tracers_2x2)
        {
                if(tr->get_last_step_wtd_idx())
                {
                        tmp_wtd[tr->get_last_step_wtd_idx()]++;
                }
        }
        return tmp_wtd;
}
std::vector<unsigned int> Lattice::get_wtd_3x3()
{
        std::vector<unsigned int> tmp_wtd(4+16*(this->m_wtd_max*this->m_wtd_res+1),0);
        for(Tracer * tr : this->m_tracers_3x3)
        {
                if(tr->get_last_step_wtd_idx())
                {
                        tmp_wtd[tr->get_last_step_wtd_idx()]++;
                }
        }
        return tmp_wtd;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<unsigned int> Lattice::get_correlations_1x1()
{
        std::vector<unsigned int> tmp_correlations(341,0);
        for(Tracer * tr : this->m_tracers_1x1)
        {
                if(tr->get_last_step())
                {
                        for(int idx : tr->get_last_step_idx()) { tmp_correlations[idx]++; }
                }
        }
        return tmp_correlations;
}
std::vector<unsigned int> Lattice::get_correlations_2x2()
{
        std::vector<unsigned int> tmp_correlations(341,0);
        for(Tracer * tr : this->m_tracers_2x2)
        {
                if(tr->get_last_step())
                {
                        for(int idx : tr->get_last_step_idx()) { tmp_correlations[idx]++; }
                }
        }
        return tmp_correlations;
}
std::vector<unsigned int> Lattice::get_correlations_3x3()
{
        std::vector<unsigned int> tmp_correlations(341,0);
        for(Tracer * tr : this->m_tracers_3x3)
        {
                if(tr->get_last_step())
                {
                        for(int idx : tr->get_last_step_idx()) { tmp_correlations[idx]++; }
                }
        }
        return tmp_correlations;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_avg_rate_correlations_1x1()
{
        //
        // in total, we have to keep track of 6 different possible site blocking cofigurations
        // obtain ids for the particles sitting on the 4 nearest neighbor sites
        //   2
        // 3 X 1
        //   4
        std::vector<double> tmp_avg_rate_correlations = {0.,0.,0.};
        for(Tracer * tr : this->m_tracers_1x1)
        {
                int x_0 = tr->get_x();
                int y_0 = tr->get_y();
                int x_1 = x_0+1;
                int x_2 = x_0-1;
                int y_1 = y_0+1;
                int y_2 = y_0-1;
                int neighboring_particle_id_1 = this->m_occupation_map[this->coord(x_1,y_0)];
                int neighboring_particle_id_2 = this->m_occupation_map[this->coord(x_0,y_1)];
                int neighboring_particle_id_3 = this->m_occupation_map[this->coord(x_2,y_0)];
                int neighboring_particle_id_4 = this->m_occupation_map[this->coord(x_0,y_2)];
                if(neighboring_particle_id_1)
                {
                        switch(this->m_tracers[neighboring_particle_id_1 - 1]->get_size())
                        {
                        case 1: tmp_avg_rate_correlations[0]++; break;
                        case 4: tmp_avg_rate_correlations[1]++; break;
                        case 9: tmp_avg_rate_correlations[2]++; break;
                        default: break;
                        }
                }
                if(neighboring_particle_id_2)
                {
                        switch(this->m_tracers[neighboring_particle_id_2 - 1]->get_size())
                        {
                        case 1: tmp_avg_rate_correlations[0]++; break;
                        case 4: tmp_avg_rate_correlations[1]++; break;
                        case 9: tmp_avg_rate_correlations[2]++; break;
                        default: break;
                        }
                }
                if(neighboring_particle_id_3)
                {
                        switch(this->m_tracers[neighboring_particle_id_3 - 1]->get_size())
                        {
                        case 1: tmp_avg_rate_correlations[0]++; break;
                        case 4: tmp_avg_rate_correlations[1]++; break;
                        case 9: tmp_avg_rate_correlations[2]++; break;
                        default: break;
                        }
                }
                if(neighboring_particle_id_4)
                {
                        switch(this->m_tracers[neighboring_particle_id_4 - 1]->get_size())
                        {
                        case 1: tmp_avg_rate_correlations[0]++; break;
                        case 4: tmp_avg_rate_correlations[1]++; break;
                        case 9: tmp_avg_rate_correlations[2]++; break;
                        default: break;
                        }
                }
        }
        return tmp_avg_rate_correlations;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_avg_rate_correlations_2x2()
{
        // - - - - - - - - - - - - - - -
        // how many possible blocking configurations?
        // obtain ids for the particles sitting on the 4 nearest neighbor sites
        //   4 3
        // 5 X X 2
        // 6 X X 1
        //   7 8
        std::vector<double> tmp_avg_rate_correlations = {0.,0.,0.,0.,0.,0.};
        for(Tracer * tr : this->m_tracers_2x2)
        {
                int x_0 = tr->get_x();
                int y_0 = tr->get_y();
                int x_1 = x_0+2;
                int x_2 = x_0+1;
                int x_3 = x_0-1;
                int y_1 = y_0+2;
                int y_2 = y_0+1;
                int y_3 = y_0-1;
                int neighboring_particle_id_1 = this->m_occupation_map[this->coord(x_1,y_0)];
                int neighboring_particle_id_2 = this->m_occupation_map[this->coord(x_1,y_2)];
                int neighboring_particle_id_3 = this->m_occupation_map[this->coord(x_2,y_1)];
                int neighboring_particle_id_4 = this->m_occupation_map[this->coord(x_0,y_1)];
                int neighboring_particle_id_5 = this->m_occupation_map[this->coord(x_3,y_2)];
                int neighboring_particle_id_6 = this->m_occupation_map[this->coord(x_3,y_0)];
                int neighboring_particle_id_7 = this->m_occupation_map[this->coord(x_0,y_3)];
                int neighboring_particle_id_8 = this->m_occupation_map[this->coord(x_2,y_3)];
        }
        return tmp_avg_rate_correlations;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_avg_rate_correlations_3x3()
{
        std::vector<double> tmp_avg_rate_correlations = {0.0};
        return tmp_avg_rate_correlations;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
// Waiting time distributions for the different tracer types
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> Lattice::get_wtd_1x1()
{
        std::vector<unsigned int> tmp_wtd_1x1 {,0};
        for(Tracer * tr : this->m_tracers_1x1)
        {

        }
        return this->m_wtd_1x1;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> Lattice::get_wtd_2x2()
{
        return this->m_wtd_2x2;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> Lattice::get_wtd_3x3()
{
        return this->m_wtd_3x3;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_norm_wtd_1x1()
{
        // If there are no tracers of the requestd type on the lattice, return 0
        if(this->m_number_of_tracers_1x1 == 0)
        {
                std::vector<double> tmp_wtd = {0};
                return tmp_wtd;
        }
        double tmp_wtd_norm = 0;
        std::vector<double> tmp_wtd;
        tmp_wtd.reserve(4*(this->m_wtd_max+1));
        // normalize sum of the resulting wtd to one
        for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
        {
                tmp_wtd_norm += (double)this->m_wtd_1x1[n];
        }
        for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
        {
                tmp_wtd.push_back((double)this->m_wtd_1x1[n]/tmp_wtd_norm);
        }
        return tmp_wtd;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_norm_wtd_2x2()
{
        // If there are no tracers of the requestd type on the lattice, return 0
        if(this->m_number_of_tracers_2x2 == 0)
        {
                std::vector<double> tmp_wtd = {0};
                return tmp_wtd;
        }
        double tmp_wtd_norm = 0;
        std::vector<double> tmp_wtd;
        tmp_wtd.reserve(4*(this->m_wtd_max+1));
        // normalize sum of the resulting wtd to one
        for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
        {
                tmp_wtd_norm += (double)this->m_wtd_2x2[n];
        }
        for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
        {
                tmp_wtd.push_back((double)this->m_wtd_2x2[n]/tmp_wtd_norm);
        }
        return tmp_wtd;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_norm_wtd_3x3()
{
        // If there are no tracers of the requestd type on the lattice, return 0
        if(this->m_number_of_tracers_3x3 == 0)
        {
                std::vector<double> tmp_wtd = {0};
                return tmp_wtd;
        }
        double tmp_wtd_norm = 0;
        std::vector<double> tmp_wtd;
        tmp_wtd.reserve(4*(this->m_wtd_max+1));
        // normalize sum of the resulting wtd to one
        for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
        {
                tmp_wtd_norm += (double)this->m_wtd_3x3[n];
        }
        for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
        {
                tmp_wtd.push_back((double)this->m_wtd_3x3[n]/tmp_wtd_norm);
        }
        return tmp_wtd;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_avg_two_step_correlation_1x1()
{
        std::vector<double> tmp_vec(4,0);
        int tmp_number_of_tracers = 0;
        for(Tracer * tr : this->m_tracers_1x1)
        {
                int idx = 0;
                if(!(tr->get_steps_taken() < 2))
                {
                        tmp_number_of_tracers++;
                        for(double tmp_dbl : tr->get_relative_two_step_correlation())
                        {
                                tmp_vec[idx] += tmp_dbl;
                                idx++;
                        }
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= tmp_number_of_tracers;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_avg_two_step_correlation_2x2()
{
        std::vector<double> tmp_vec(4,0);
        int tmp_number_of_tracers = 0;
        for(Tracer * tr : this->m_tracers_2x2)
        {
                int idx = 0;
                if(!(tr->get_steps_taken() < 2))
                {
                        tmp_number_of_tracers++;
                        for(double tmp_dbl : tr->get_relative_two_step_correlation())
                        {
                                tmp_vec[idx] += tmp_dbl;
                                idx++;
                        }
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= tmp_number_of_tracers;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_avg_two_step_correlation_3x3()
{
        std::vector<double> tmp_vec(4,0);
        int tmp_number_of_tracers = 0;
        for(Tracer * tr : this->m_tracers_3x3)
        {
                int idx = 0;
                if(!(tr->get_steps_taken() < 2))
                {
                        tmp_number_of_tracers++;
                        for(double tmp_dbl : tr->get_relative_two_step_correlation())
                        {
                                tmp_vec[idx] += tmp_dbl;
                                idx++;
                        }
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= tmp_number_of_tracers;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_relative_two_step_correlations_1x1()
{
        std::vector<double> tmp_vec;
        tmp_vec.reserve(4*this->m_number_of_tracers_1x1);
        for(Tracer * tr : this->m_tracers_1x1)
        {
                for(double dbl : tr->get_relative_two_step_correlation())
                {
                        tmp_vec.push_back(dbl);
                }
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_relative_two_step_correlations_2x2()
{
        std::vector<double> tmp_vec;
        tmp_vec.reserve(4*this->m_number_of_tracers_2x2);
        for(Tracer * tr : this->m_tracers_2x2)
        {
                for(double dbl : tr->get_relative_two_step_correlation())
                {
                        tmp_vec.push_back(dbl);
                }
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_relative_two_step_correlations_3x3()
{
        std::vector<double> tmp_vec;
        tmp_vec.reserve(4*this->m_number_of_tracers_3x3);
        for(Tracer * tr : this->m_tracers_3x3)
        {
                for(double dbl : tr->get_relative_two_step_correlation())
                {
                        tmp_vec.push_back(dbl);
                }
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_avg_three_step_correlation_1x1()
{
        std::vector<double> tmp_vec(16,0);
        int tmp_number_of_tracers = 0;
        for(Tracer * tr : this->m_tracers_1x1)
        {
                int idx = 0;
                if(!(tr->get_steps_taken() < 3))
                {
                        tmp_number_of_tracers++;
                        for(double tmp_dbl : tr->get_relative_three_step_correlation())
                        {
                                tmp_vec[idx] += tmp_dbl;
                                idx++;
                        }
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= tmp_number_of_tracers;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_avg_three_step_correlation_2x2()
{
        std::vector<double> tmp_vec(16,0);
        int tmp_number_of_tracers = 0;
        for(Tracer * tr : this->m_tracers_2x2)
        {
                int idx = 0;
                if(!(tr->get_steps_taken() < 3))
                {
                        tmp_number_of_tracers++;
                        for(double tmp_dbl : tr->get_relative_three_step_correlation())
                        {
                                tmp_vec[idx] += tmp_dbl;
                                idx++;
                        }
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= tmp_number_of_tracers;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_avg_three_step_correlation_3x3()
{
        std::vector<double> tmp_vec(16,0);
        int tmp_number_of_tracers = 0;
        for(Tracer * tr : this->m_tracers_3x3)
        {
                int idx = 0;
                if(!(tr->get_steps_taken() < 3))
                {
                        tmp_number_of_tracers++;
                        for(double tmp_dbl : tr->get_relative_three_step_correlation())
                        {
                                tmp_vec[idx] += tmp_dbl;
                                idx++;
                        }
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= tmp_number_of_tracers;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_relative_three_step_correlations_1x1()
{
        std::vector<double> tmp_vec;
        tmp_vec.reserve(16*this->m_number_of_tracers_1x1);
        for(Tracer * tr : this->m_tracers_1x1)
        {
                for(double dbl : tr->get_relative_three_step_correlation())
                {
                        tmp_vec.push_back(dbl);
                }
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_relative_three_step_correlations_2x2()
{
        std::vector<double> tmp_vec;
        tmp_vec.reserve(16*this->m_number_of_tracers_2x2);
        for(Tracer * tr : this->m_tracers_2x2)
        {
                for(double dbl : tr->get_relative_three_step_correlation())
                {
                        tmp_vec.push_back(dbl);
                }
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_relative_three_step_correlations_3x3()
{
        std::vector<double> tmp_vec;
        tmp_vec.reserve(16*this->m_number_of_tracers_3x3);
        for(Tracer * tr : this->m_tracers_3x3)
        {
                for(double dbl : tr->get_relative_three_step_correlation())
                {
                        tmp_vec.push_back(dbl);
                }
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_avg_four_step_correlation_1x1()
{
        std::vector<double> tmp_vec(64,0);
        int tmp_number_of_tracers = 0;
        for(Tracer * tr : this->m_tracers_1x1)
        {
                int idx = 0;
                if(!(tr->get_steps_taken() < 4))
                {
                        tmp_number_of_tracers++;
                        for(double tmp_dbl : tr->get_relative_four_step_correlation())
                        {
                                tmp_vec[idx] += tmp_dbl;
                                idx++;
                        }
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= tmp_number_of_tracers;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_avg_four_step_correlation_2x2()
{
        std::vector<double> tmp_vec(64,0);
        int tmp_number_of_tracers = 0;
        for(Tracer * tr : this->m_tracers_2x2)
        {
                int idx = 0;
                if(!(tr->get_steps_taken() < 4))
                {
                        tmp_number_of_tracers++;
                        for(double tmp_dbl : tr->get_relative_four_step_correlation())
                        {
                                tmp_vec[idx] += tmp_dbl;
                                idx++;
                        }
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= tmp_number_of_tracers;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Lattice::get_avg_four_step_correlation_3x3()
{
        std::vector<double> tmp_vec(64,0);
        int tmp_number_of_tracers = 0;
        for(Tracer * tr : this->m_tracers_3x3)
        {
                int idx = 0;
                if(!(tr->get_steps_taken() < 4))
                {
                        tmp_number_of_tracers++;
                        for(double tmp_dbl : tr->get_relative_four_step_correlation())
                        {
                                tmp_vec[idx] += tmp_dbl;
                                idx++;
                        }
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= tmp_number_of_tracers;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> Lattice::get_tracer_positions()
{
        std::vector<int> tmp_tracer_positions;
        tmp_tracer_positions.reserve(4*this->m_number_of_tracers_total);
        // Prints Tracer information to std output
        for(Tracer * tr : this->m_tracers)
        {
                tmp_tracer_positions.push_back(tr->get_id());
                tmp_tracer_positions.push_back(tr->get_size());
                tmp_tracer_positions.push_back(tr->get_x());
                tmp_tracer_positions.push_back(tr->get_y());
        }
        return tmp_tracer_positions;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> get_relative_particle_placement_1x1()
{
        std::vector<int> tmp_mod_tracer_positions;
        tmp_mod_tracer_positions.reserve(4);
        return tmp_mod_tracer_positions;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> get_relative_particle_placement_2x2()
{
        std::vector<int> tmp_mod_tracer_positions;
        tmp_mod_tracer_positions.reserve(4);
        return tmp_mod_tracer_positions;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
// D E B U G G I N G  O U T P U T
// - - - - - - - - - - - - - - - - - - - - - - - -
void Lattice::print_occupation_map()
{
        // Prints the occupation map vector to std output
        std::cout << std::endl;
        for(int y = this->m_grid_size_y - 1; y >= 0; y--)
        {
                for(int x = 0; x < this->m_grid_size_x; x++)
                {
                        std::cout << this->m_occupation_map[this->coord(x,y)] << " ";
                }
                std::cout << std::endl;
        }
        std::cout << std::endl;
}
// - - - - - - - - - - - - - - - - - - - - - - - -
void Lattice::print_tracer_positions()
{
        // Prints Tracer information to std output
        for(Tracer * tr : this->m_tracers)
        {
                std::cout << "T: " << this->m_t << " ID: " << tr->get_id() << " X: " << tr->get_x() << " Y: " << tr->get_y() << " dx: " << tr->get_dx() << " dy: " << tr->get_dy() << " l^2 " << tr->get_lsquared() << " Size: " << tr->get_size() << std::endl;
        }
}
