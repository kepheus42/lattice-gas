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
        m_wtd_max(wtd_max),
        m_wtd_res(wtd_res),
        m_movement_selector((int)(this->m_number_of_tracers_1x1/this->m_step_rate_2x2/this->m_step_rate_3x3)+
                            (int)(this->m_number_of_tracers_2x2/this->m_step_rate_1x1/this->m_step_rate_3x3)+
                            (int)(this->m_number_of_tracers_3x3/this->m_step_rate_1x1/this->m_step_rate_2x2),0),
        m_movement_selector_length(this->m_movement_selector.size()),
        m_step_attempts_per_timestep((int)(this->m_number_of_tracers_1x1*this->m_step_rate_1x1)+
                                     (int)(this->m_number_of_tracers_2x2*this->m_step_rate_2x2)+
                                     (int)(this->m_number_of_tracers_3x3*this->m_step_rate_3x3)),
        m_t_increment(1.0/(double)this->m_step_attempts_per_timestep),
        m_occupation_map(grid_size_x*grid_size_y,0)

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
void Lattice::setup_movement_selection_list()
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
                this->m_movement_selector[tmp_idx] = n/tmp_m3;
                tmp_idx++;
        }
        for(int n = 0; n < tmp_n2; n++)
        {
                this->m_movement_selector[tmp_idx] = tmp_idx_2x2 + n/tmp_m2;
                tmp_idx++;
        }
        for(int n = 0; n < tmp_n1; n++)
        {
                this->m_movement_selector[tmp_idx] =  tmp_idx_1x1 + n/tmp_m1;
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
                tmp_start_positions_3x3.reserve(this->m_grid_size_x/3 * this->m_grid_size_y/3);
                for (int y=0; y < this->m_grid_size_y; y+=3)
                {
                        for (int x=0; x < this->m_grid_size_x; x+=3)
                        {
                                tmp_start_positions_3x3.push_back(x * this->m_grid_size_y + y);
                        }
                }
                // create the 3x3 tracers:
                shuffle_vector(tmp_start_positions_3x3);
                for (int n = 0; n < this->m_number_of_tracers_3x3; n++)
                {
                        tmp_y = tmp_start_positions_3x3[n] % this->m_grid_size_y;
                        tmp_x = (tmp_start_positions_3x3[n] - tmp_y)/this->m_grid_size_y;
                        this->m_tracers.push_back(new Tracer_3x3(tmp_id,tmp_x,tmp_y,this->m_grid_size_x,this->m_grid_size_y,this->m_step_rate_3x3,this->m_wtd_max,this->m_wtd_res,this->m_step_attempts_per_timestep));
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
                for (int n = 0; n < this->m_number_of_tracers_2x2; n++)
                {
                        tmp_y = tmp_start_positions_2x2[n] % this->m_grid_size_y;
                        tmp_x = (tmp_start_positions_2x2[n] - tmp_y)/this->m_grid_size_y;
                        this->m_tracers.push_back(new Tracer_2x2(tmp_id,tmp_x,tmp_y,this->m_grid_size_x,this->m_grid_size_y,this->m_step_rate_2x2,this->m_wtd_max,this->m_wtd_res,this->m_step_attempts_per_timestep));
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
                        this->m_tracers.push_back(new Tracer(tmp_id,tmp_x,tmp_y,this->m_grid_size_x,this->m_grid_size_y,this->m_step_rate_1x1,this->m_wtd_max,this->m_wtd_res,this->m_step_attempts_per_timestep));
                        this->m_tracers_1x1.push_back(this->m_tracers.back());
                        this->m_occupation_map[this->coord(tmp_x,tmp_y)] = tmp_id;
                        tmp_id++;
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
        return;
}
//
void Lattice::timestep(){
        int tmp_rnd; // random number between 0 and 4*length of selector list
        int tmp_par; // which particle moves
        int tmp_dir; // which direction it moves in
        long tmp_time = this->m_t * this->m_step_attempts_per_timestep;
        for (long n=0; n < this->m_step_attempts_per_timestep; n++)
        {
                tmp_rnd = random_int(0,4*this->m_movement_selector_length-1);
                tmp_par = this->m_movement_selector[tmp_rnd/4];
                tmp_dir = 1+tmp_rnd%4;
                // eacth timestep consists of m_step_attempts_per_timestep equal intervals, thus the time at which a given step attempt takes place is t = m_t + n / m_step_attempts_per_timestep
                this->m_tracers[tmp_par]->step(this->m_occupation_map,tmp_dir,tmp_time+n);
        }
        this->m_t++;
}
//
void Lattice::timestep_warmup(){
        int tmp_rnd; // random number between 0 and 4*length of selector list
        int tmp_par; // which particle moves
        int tmp_dir; // which direction it moves in
        for (long n=0; n < this->m_step_attempts_per_timestep; n++)
        {
                tmp_rnd = random_int(0,4*this->m_movement_selector_length-1);
                tmp_par = this->m_movement_selector[tmp_rnd/4];
                tmp_dir = 1+tmp_rnd%4;
                this->m_tracers[tmp_par]->step_warmup(this->m_occupation_map,tmp_dir);
        }
}
// - - - - - - - - - - - - - - - - - - - - - - - -
void Lattice::timestep_no_interaction(){
        int tmp_rnd; // random number between 0 and 4*length of selector list
        int tmp_par; // which particle moves
        int tmp_dir; // which direction it moves in
        long tmp_time = this->m_t * this->m_step_attempts_per_timestep;
        for (long n=0; n < this->m_step_attempts_per_timestep; n++)
        {
                tmp_rnd = random_int(0,4*this->m_movement_selector_length);
                tmp_par = this->m_movement_selector[tmp_rnd/4];
                tmp_dir = 1+tmp_rnd%4;
                this->m_tracers[tmp_par]->step_unhindered(tmp_dir,tmp_time+n);
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
int Lattice::get_tracer_size(int id)
{
        return this->m_tracers[id]->get_size();
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
// - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> Lattice::get_occupation_map()
{
        return this->m_occupation_map;
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
