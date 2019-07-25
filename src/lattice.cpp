#include <stdexcept>
#include <iostream>

#include "lattice.hpp"
#include "global.hpp"

Lattice::Lattice(int grid_size_x, int grid_size_y, int number_of_tracers_1x1, int number_of_tracers_2x2) : m_t(0),  m_grid_size_x(grid_size_x), m_grid_size_y(grid_size_y), m_number_of_tracers_1x1(number_of_tracers_1x1), m_number_of_tracers_2x2(number_of_tracers_2x2){
        // check if there's enough space on the grid to place the tracers
        if (grid_size_x*grid_size_y < number_of_tracers_1x1 + 4*number_of_tracers_2x2) {
                throw std::invalid_argument("Too many tracers for grid of the given size!");
        }
        // create a vector to handle the randomization of the order of movement every timestep
        this->m_movement_order.reserve(number_of_tracers_1x1 + number_of_tracers_2x2);
        for (int n = 0; n < number_of_tracers_1x1 + number_of_tracers_2x2; n++) {
                m_movement_order.push_back(n);
        }
        // create a vector to store pointers to all the tracers
        this->m_tracers.reserve(number_of_tracers_1x1 + number_of_tracers_2x2);
        // create a vector to store the state of the lattice
        this->m_occupation_map.reserve(grid_size_x * grid_size_y);
        for(int i=0; i < grid_size_x * grid_size_y; i++)
        {
                this->m_occupation_map[i] = 0;
        }
        // tmp vector to handle site allocation for 2x2 tracers
        std::vector<int> tmp_start_positions_2x2;
        tmp_start_positions_2x2.reserve(grid_size_x/2 * grid_size_y/2);
        for (int y=0; y < grid_size_y; y+=2)
        {
                for (int x=0; x < grid_size_x; x+=2)
                {
                        tmp_start_positions_2x2.push_back(x * grid_size_y + y);
                }
        }
        // create the 2x2 tracers:
        tmp_start_positions_2x2 = shuffle_vector(tmp_start_positions_2x2);
        int tmp_x;
        int tmp_y;
        int tmp_id = 1;
        for (int n = 0; n < number_of_tracers_2x2; n++)
        {
                tmp_y = tmp_start_positions_2x2[n] % grid_size_y;
                tmp_x = (tmp_start_positions_2x2[n] - tmp_y)/grid_size_y;
                this->m_tracers.push_back(new Tracer_2x2(tmp_id,tmp_x,tmp_y,grid_size_x,grid_size_y));
                this->m_occupation_map[tmp_x * grid_size_y + tmp_y] = tmp_id;
                this->m_occupation_map[(tmp_x+1)%grid_size_x * grid_size_y + tmp_y] = tmp_id;
                this->m_occupation_map[tmp_x * grid_size_y + (tmp_y+1)%grid_size_y] = tmp_id;
                this->m_occupation_map[(tmp_x+1)%grid_size_x * grid_size_y + (tmp_y+1)%grid_size_y] = tmp_id;
                tmp_id++;
        }
        // tmp vector to handle creation of the 1x1 tracers
        std::vector<int> tmp_start_positions_1x1;
        tmp_start_positions_1x1.reserve(grid_size_x * grid_size_y - 4 * number_of_tracers_2x2);
        for (int i = 0; i < grid_size_x*grid_size_y; i++)
        {
                if(!this->m_occupation_map[i])
                {
                        tmp_start_positions_1x1.push_back(i);
                }
        }
        tmp_start_positions_1x1 = shuffle_vector(tmp_start_positions_1x1);
        for (int n = 0; n < number_of_tracers_1x1; n++)
        {
                tmp_y = tmp_start_positions_1x1[n] % grid_size_y;
                tmp_x = (tmp_start_positions_1x1[n] - tmp_y)/grid_size_y;
                this->m_tracers.push_back(new Tracer(tmp_id,tmp_x,tmp_y,grid_size_x,grid_size_y));
                this->m_occupation_map[tmp_x * grid_size_y + tmp_y] = tmp_id;
                tmp_id++;
        }
}
int Lattice::coord(int x, int y){
        return x*this->m_grid_size_y + y;
}
void Lattice::timestep(){
        this->m_movement_order = shuffle_vector(this->m_movement_order);
        for (int n=0; n < this->m_number_of_tracers_1x1 + this->m_number_of_tracers_2x2; n++)
        {
                this->m_tracers[this->m_movement_order[n]]->step(this->m_occupation_map);
        }
        this->m_t++;
}
void Lattice::timestep_no_interaction(){
        // use this to let all tracers step without interacting
        for (int n=0; n < this->m_number_of_tracers_1x1 + this->m_number_of_tracers_2x2; n++)
        {
                this->m_tracers[n]->unhindered_step();
        }
        this->m_t++;
}
// - - - - - - - - - -
// G E T T E R S
// - - - - - - - - - -
int Lattice::get_t(){
        return this->m_t;
}
int Lattice::get_number_of_tracers_1x1()
{
        return this->m_number_of_tracers_1x1;
}
int Lattice::get_number_of_tracers_2x2()
{
        return this->m_number_of_tracers_2x2;
}
int Lattice::get_grid_size_x()
{
        return this->m_grid_size_x;
}
int Lattice::get_grid_size_y()
{
        return this->m_grid_size_y;
}
std::vector<int> Lattice::get_occupation_map()
{
        return this->m_occupation_map;
}
double Lattice::get_avg_lsquared_2x2()
{
        // Returns the avg lsquared of all 2x2 tracers
        if(this->m_number_of_tracers_2x2 == 0)
        {
                return 0.0;
        }
        double tmp_sum_lsquared = 0;
        // double tmp_avg_lsquared = 0.0;
        for(Tracer * tr : this->m_tracers)
        {
                if(tr->get_size() == 4)
                {
                        tmp_sum_lsquared += tr->get_lsquared();
                }
        }
        return (tmp_sum_lsquared/(double)this->m_number_of_tracers_2x2);
}
double Lattice::get_avg_lsquared_1x1()
{
        // Returns the avg lsquared of all 1x1 tracers
        if(this->m_number_of_tracers_1x1 == 0)
        {
                return 0.0;
        }
        double tmp_sum_lsquared = 0;
        // double tmp_avg_lsquared = 0.0;
        for(Tracer * tr : this->m_tracers)
        {
                if(tr->get_size() == 1)
                {
                        tmp_sum_lsquared += tr->get_lsquared();
                }
        }
        return (tmp_sum_lsquared/(double)this->m_number_of_tracers_1x1);
}
// - - - - - - - - - - - - - - - -
// D E B U G G I N G  O U T P U T
// - - - - - - - - - - - - - - - -
void Lattice::print_occupation_map()
{
        // Prints the occupation map vector to std output
        std::cout << std::endl;
        for(int y = this->m_grid_size_y - 1; y >= 0; y--)
        {
                for(int x = 0; x < this->m_grid_size_x; x++)
                {
                        std::cout << this->m_occupation_map[x*this->m_grid_size_y + y] << " ";
                }
                std::cout << std::endl;
        }
        std::cout << std::endl;
}
void Lattice::print_tracer_positions()
{
        // Prints Tracer information to std output
        for(Tracer * tr : this->m_tracers)
        {
                std::cout << "T: " << this->m_t << " ID: " << tr->get_id() << " X: " << tr->get_x() << " Y: " << tr->get_y() << " dx: " << tr->get_dx() << " dy: " << tr->get_dy() << " l^2 " << tr->get_lsquared() << " Size: " << tr->get_size() << std::endl;
        }
}
