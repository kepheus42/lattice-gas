#include "lattice.hpp"

// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif
// To add debugging messages, use D(std::cerr << "Debugging message 1 2 3!" << std::endl; )

// - - - - - - - - - - - - - - - - - - - - - - - -
Lattice::Lattice(int grid_size,
                 int number_of_tracers_1x1,
                 int number_of_tracers_2x2,
                 double step_rate_1x1,
                 double step_rate_2x2) :
        m_grid_size(grid_size),
        m_number_of_sites(grid_size*grid_size),
        m_number_of_tracers(number_of_tracers_1x1+number_of_tracers_2x2),
        m_number_of_tracers_1x1(number_of_tracers_1x1),
        m_number_of_tracers_2x2(number_of_tracers_2x2),
        m_step_rate_1x1(this->m_number_of_tracers_1x1 > 0 ? step_rate_1x1 : 1.0),
        m_step_rate_2x2(this->m_number_of_tracers_2x2 > 0 ? step_rate_2x2 : 1.0),
        m_movement_selector(this->m_number_of_tracers_1x1*(int)(1/this->m_step_rate_2x2)+this->m_number_of_tracers_2x2*(int)(1/this->m_step_rate_1x1),0),
        m_movement_selector_length(this->m_movement_selector.size()),
        m_step_attempts_per_timestep((int)(this->m_number_of_tracers_1x1*this->m_step_rate_1x1)+(int)(this->m_number_of_tracers_2x2*this->m_step_rate_2x2)),
        m_tracer_locations(this->m_number_of_tracers,0),
        m_occupation_map(2*this->m_number_of_sites,0)
{
        // check if there's enough space on the grid to place the tracers
        if (grid_size*grid_size < number_of_tracers_1x1 + 4*number_of_tracers_2x2)
        {
                throw std::invalid_argument("Too many tracers for grid of the given size!");
        }
        else if ((number_of_tracers_2x2 > 0)&&((grid_size < 2)||(grid_size < 2)))
        {
                throw std::invalid_argument("Grid is too narrow for 2x2 or 3x3 tracers!");
        }

        this->m_sites_1x1.reserve(this->m_number_of_sites);
        this->m_sites_2x2.reserve(this->m_number_of_sites);

        this->m_tracers.reserve(this->m_number_of_tracers);
        this->m_tracers_1x1.reserve(this->m_number_of_tracers_1x1);
        this->m_tracers_2x2.reserve(this->m_number_of_tracers_2x2);

        //D( std::cout << "Capacity of m_sites_1x1: " << this->m_sites_1x1.capacity() << std::endl << "Capacity of m_sites_2x2: " << this->m_sites_2x2.capacity() << std::endl );
        //D( std::cout << "Parameters of lattice: R1 " << this->m_step_rate_1x1 << " R2 " << this->m_step_rate_2x2 << std::endl );
        this->setup_sites();
        this->setup_movement_selection_list();
        this->setup_tracers();

}

void Lattice::setup_sites(){
        int tmp_x;
        int tmp_y;
        int* tmp_occupation_map = this->m_occupation_map.data();
        D( std::cout << "Starting neighbors setup" << std::endl );
        for(int n = 0; n < this->m_number_of_sites; n++) {
                tmp_x = n/this->m_grid_size;
                tmp_y = n%this->m_grid_size;
                std::vector< std::vector< int *> > tmp_neighbors;
                tmp_neighbors.reserve(4);
                std::vector< int * > tmp_neighbors_dir(3);
                // dir = 1
                tmp_neighbors_dir[0] = &tmp_occupation_map[this->coord_1x1(tmp_x+1,tmp_y+0)];
                tmp_neighbors_dir[1] = &tmp_occupation_map[this->coord_2x2(tmp_x+1,tmp_y+0)];
                tmp_neighbors_dir[2] = &tmp_occupation_map[this->coord_2x2(tmp_x+1,tmp_y+1)];
                tmp_neighbors.push_back(tmp_neighbors_dir);
                // dir = 2
                tmp_neighbors_dir[0] = &tmp_occupation_map[this->coord_1x1(tmp_x+0,tmp_y+1)];
                tmp_neighbors_dir[1] = &tmp_occupation_map[this->coord_2x2(tmp_x+0,tmp_y+2)];
                tmp_neighbors_dir[2] = &tmp_occupation_map[this->coord_2x2(tmp_x-1,tmp_y+2)];
                tmp_neighbors.push_back(tmp_neighbors_dir);
                // dir = 3
                tmp_neighbors_dir[0] = &tmp_occupation_map[this->coord_1x1(tmp_x-1,tmp_y+0)];
                tmp_neighbors_dir[1] = &tmp_occupation_map[this->coord_2x2(tmp_x-2,tmp_y+1)];
                tmp_neighbors_dir[2] = &tmp_occupation_map[this->coord_2x2(tmp_x-2,tmp_y+0)];
                tmp_neighbors.push_back(tmp_neighbors_dir);
                // dir = 4
                tmp_neighbors_dir[0] = &tmp_occupation_map[this->coord_1x1(tmp_x+0,tmp_y-1)];
                tmp_neighbors_dir[1] = &tmp_occupation_map[this->coord_2x2(tmp_x-1,tmp_y-1)];
                tmp_neighbors_dir[2] = &tmp_occupation_map[this->coord_2x2(tmp_x+0,tmp_y-1)];
                tmp_neighbors.push_back(tmp_neighbors_dir);
                //
                this->m_neighbors.push_back(tmp_neighbors);
        }

        for(int n = this->m_number_of_sites; n < 2*this->m_number_of_sites; n++) {
                tmp_x = n/this->m_grid_size;
                tmp_y = n%this->m_grid_size;
                std::vector< std::vector< int *> > tmp_neighbors;
                tmp_neighbors.reserve(4);
                std::vector< int * > tmp_neighbors_dir(6);
                // dir = 1
                tmp_neighbors_dir[0] = &tmp_occupation_map[this->coord_2x2(tmp_x+1,tmp_y+0)];
                tmp_neighbors_dir[1] = &tmp_occupation_map[this->coord_1x1(tmp_x+2,tmp_y-1)];
                tmp_neighbors_dir[2] = &tmp_occupation_map[this->coord_1x1(tmp_x+2,tmp_y+0)];
                tmp_neighbors_dir[3] = &tmp_occupation_map[this->coord_2x2(tmp_x+3,tmp_y+1)];
                tmp_neighbors_dir[4] = &tmp_occupation_map[this->coord_2x2(tmp_x+3,tmp_y+0)];
                tmp_neighbors_dir[5] = &tmp_occupation_map[this->coord_2x2(tmp_x+3,tmp_y-1)];
                tmp_neighbors.push_back(tmp_neighbors_dir);
                // dir = 2
                tmp_neighbors_dir[0] = &tmp_occupation_map[this->coord_2x2(tmp_x+0,tmp_y+1)];
                tmp_neighbors_dir[1] = &tmp_occupation_map[this->coord_1x1(tmp_x+1,tmp_y+1)];
                tmp_neighbors_dir[2] = &tmp_occupation_map[this->coord_1x1(tmp_x+0,tmp_y+1)];
                tmp_neighbors_dir[3] = &tmp_occupation_map[this->coord_2x2(tmp_x+1,tmp_y+3)];
                tmp_neighbors_dir[4] = &tmp_occupation_map[this->coord_2x2(tmp_x+0,tmp_y+3)];
                tmp_neighbors_dir[5] = &tmp_occupation_map[this->coord_2x2(tmp_x-1,tmp_y+3)];
                tmp_neighbors.push_back(tmp_neighbors_dir);
                // dir = 3
                tmp_neighbors_dir[0] = &tmp_occupation_map[this->coord_2x2(tmp_x-1,tmp_y+0)];
                tmp_neighbors_dir[1] = &tmp_occupation_map[this->coord_1x1(tmp_x-1,tmp_y+0)];
                tmp_neighbors_dir[2] = &tmp_occupation_map[this->coord_1x1(tmp_x-1,tmp_y-1)];
                tmp_neighbors_dir[3] = &tmp_occupation_map[this->coord_2x2(tmp_x-3,tmp_y+1)];
                tmp_neighbors_dir[4] = &tmp_occupation_map[this->coord_2x2(tmp_x-3,tmp_y+0)];
                tmp_neighbors_dir[5] = &tmp_occupation_map[this->coord_2x2(tmp_x-3,tmp_y-1)];
                tmp_neighbors.push_back(tmp_neighbors_dir);
                // dir = 4
                tmp_neighbors_dir[0] = &tmp_occupation_map[this->coord_2x2(tmp_x+0,tmp_y-1)];
                tmp_neighbors_dir[1] = &tmp_occupation_map[this->coord_1x1(tmp_x+0,tmp_y-2)];
                tmp_neighbors_dir[2] = &tmp_occupation_map[this->coord_1x1(tmp_x+1,tmp_y-2)];
                tmp_neighbors_dir[3] = &tmp_occupation_map[this->coord_2x2(tmp_x-1,tmp_y-3)];
                tmp_neighbors_dir[4] = &tmp_occupation_map[this->coord_2x2(tmp_x+0,tmp_y-3)];
                tmp_neighbors_dir[5] = &tmp_occupation_map[this->coord_2x2(tmp_x+1,tmp_y-3)];
                tmp_neighbors.push_back(tmp_neighbors_dir);
                //
                this->m_neighbors.push_back(tmp_neighbors);
        }
        D( std::cout << "Finishing Neighbor setup " << std::endl );
        D( this->print_neighbors() );
        //D( this->print_sites() );
        //D( std::cout << "Size of m_sites_1x1: " << this->m_sites_1x1.size() << std::endl << "Size of m_sites_2x2: " << this->m_sites_2x2.size() << std::endl);
}

void Lattice::setup_movement_selection_list()
{
        //
        // structure of the selection list vector:
        //  T Y P E 1               T Y P E 2
        // 0 0 0 1 1 1 2 2 2 ... 10 10 11 11 12 12
        // in this case, rate_1/rate_2 = 3/2

        // use the movement selector list to randomly pick particles to move, with probabilities proportional to the ratios of the stepping rates and numbers of particles
        int tmp_idx = 0;
        // The first 1x1 tracer is this->Tracers[n_2x2], the first 2x2 tracer is this->Tracers[0]
        int tmp_idx_1x1 = this->m_number_of_tracers_2x2;
        // tmp_m1 - tmp_m3 : how many times each tracer idx for tracers of types 1-3 appears in the list
        int tmp_m1 = (int)(1.0/this->m_step_rate_2x2);
        int tmp_m2 = (int)(1.0/this->m_step_rate_1x1);
        // tmp_n1 - tmp_n3 : how many list entires do types 1 - 3 get in the final list
        int tmp_n1 = this->m_number_of_tracers_1x1*tmp_m1;
        int tmp_n2 = this->m_number_of_tracers_2x2*tmp_m2;
        //D( std::cout << "Setup movement selector list: " << tmp_m1 << " " << tmp_m2 << " " << tmp_n1 << " " << tmp_n2 << std::endl );
        // fill the movement selector vector according to the numbers of tracers
        if(this->m_number_of_tracers_2x2) {
                for(int n = 0; n < tmp_n2; n++)
                {
                        this->m_movement_selector[tmp_idx] = n/tmp_m2;
                        tmp_idx++;
                }
        }
        if(this->m_number_of_tracers_1x1) {
                for(int n = 0; n < tmp_n1; n++)
                {
                        this->m_movement_selector[tmp_idx] =  tmp_idx_1x1 + n/tmp_m1;
                        tmp_idx++;
                }
        }
        D( std::cout << "Size of m_movement_selector: " << this->m_movement_selector.size() << std::endl << "Step attempts per timestep: " << this->m_step_attempts_per_timestep << std::endl );
        D( print_vector(this->m_movement_selector) );
}

void Lattice::setup_tracers()
{
        int tmp_id = 0;
        int tmp_x;
        int tmp_y;
        int tmp_pos;
        Site * tmp_site;
        std::vector<int> tmp_occupied_sites(this->m_grid_size * this->m_grid_size,0);
        D( std::cout << "Starting Tracer setup" << std::endl );
        if(this->m_number_of_tracers_2x2)
        {
                this->m_tracers_2x2.reserve(this->m_number_of_tracers_2x2);
                // tmp vector to handle start position allocation for 2x2 tracers
                std::vector<int> tmp_start_positions_2x2;
                tmp_start_positions_2x2.reserve(this->m_grid_size/2 * this->m_grid_size/2);
                for (int y=0; y < this->m_grid_size; y+=2)
                {
                        for (int x=0; x < this->m_grid_size; x+=2)
                        {
                                // TODO: add check for occupancy
                                tmp_start_positions_2x2.push_back(this->coord(x,y));
                        }
                }
                //D( std::cout << "Start positions" << std::endl );
                //D( this->print_vector(tmp_start_positions_2x2) );
                // create the 2x2 tracers:
                std::shuffle(tmp_start_positions_2x2.begin(),tmp_start_positions_2x2.end(),generator);
                // shuffle_vector(tmp_start_positions_2x2);
                //D( std::cout << "Start positions (shuffled)" << std::endl );
                //D( this->print_vector(tmp_start_positions_2x2) );
                for (int n = 0; n < this->m_number_of_tracers_2x2; n++)
                {
                        tmp_pos = tmp_start_positions_2x2[n];
                        // x and y location of the starting site
                        tmp_x = this->coord_to_x(tmp_pos);
                        tmp_y = this->coord_to_y(tmp_pos);
                        this->m_tracer_locations[tmp_id] = tmp_pos;
                        this->m_occupation_map[this->coord_2x2(tmp_x,tmp_y)] = tmp_id;
                        // create new tracer at starting site
                        this->m_tracers.push_back(new Tracer_2x2(tmp_id,tmp_pos));
                        this->m_tracers_2x2.push_back(this->m_tracers.back());
                        // increment id counter,
                        tmp_id++;
                        tmp_occupied_sites[coord(tmp_x+0,tmp_y+0)] = 1;
                        tmp_occupied_sites[coord(tmp_x+1,tmp_y+0)] = 1;
                        tmp_occupied_sites[coord(tmp_x+0,tmp_y+1)] = 1;
                        tmp_occupied_sites[coord(tmp_x+1,tmp_y+1)] = 1;
                }
                //D( std::cout << "occupied sites after 2x2 setup: " << std::endl );
                //D( this->print_vector(tmp_occupied_sites) );
                //D( std::cout << "Number of Tracers 2x2: " << this->m_tracers_2x2.size() << std::endl );
        }
        if(this->m_number_of_tracers_1x1)
        {
                this->m_tracers_1x1.reserve(this->m_number_of_tracers_1x1);
                // tmp vector to handle creation of the 1x1 tracers
                std::vector<int> tmp_start_positions_1x1;
                tmp_start_positions_1x1.reserve(this->m_grid_size * this->m_grid_size - 4 * this->m_number_of_tracers_2x2);
                for (int i = 0; i < this->m_grid_size*this->m_grid_size; i++)
                {
                        if(!tmp_occupied_sites[i])
                        {
                                tmp_start_positions_1x1.push_back(i);
                        }
                }
                //D( std::cout << "Start positions" << std::endl );
                //D( this->print_vector(tmp_start_positions_1x1) );
                // shuffle_vector(tmp_start_positions_1x1);
                std::shuffle(tmp_start_positions_1x1.begin(),tmp_start_positions_1x1.end(),generator);
                //D( std::cout << "Start positions (shuffled)" << std::endl );
                //D( this->print_vector(tmp_start_positions_1x1) );
                for (int n = 0; n < this->m_number_of_tracers_1x1; n++)
                {
                        tmp_pos = tmp_start_positions_1x1[n];
                        // x and y location of the starting site
                        tmp_x = this->coord_to_x(tmp_pos);
                        tmp_y = this->coord_to_y(tmp_pos);
                        this->m_tracer_locations[tmp_id] = tmp_pos;
                        this->m_occupation_map[this->coord_1x1(tmp_x,tmp_y)] = tmp_id;
                        // create new tracer at starting site
                        this->m_tracer_locations[tmp_id] = tmp_start_positions_1x1[n];
                        this->m_tracers.push_back(new Tracer(tmp_id,tmp_pos));
                        this->m_tracers_1x1.push_back(this->m_tracers.back());
                        // increment id counter,
                        tmp_id++;
                }
                //D( std::cout << "Number of Tracers 1x1: " << this->m_tracers_1x1.size() << std::endl);
        }
        D( print_vector(this->m_tracer_locations) );
        //D( std::cout << "Number of Tracers: " << this->m_tracers.size() << std::endl);
}

// = = = = = = = = = = = = = = = = = = = = = =
// = = = M E M B E R   F U N C T I O N S = = =
// = = = = = = = = = = = = = = = = = = = = = =

inline int Lattice::coord(int x, int y){
        //
        return ((x+this->m_grid_size)%this->m_grid_size)*this->m_grid_size+((y+this->m_grid_size)%this->m_grid_size);
}

inline int Lattice::coord_1x1(int x, int y){
        //
        return ((x+this->m_grid_size)%this->m_grid_size)*this->m_grid_size+((y+this->m_grid_size)%this->m_grid_size);
}

inline int Lattice::coord_2x2(int x, int y){
        //
        return ((x+this->m_grid_size)%this->m_grid_size)*this->m_grid_size+((y+this->m_grid_size)%this->m_grid_size)+this->m_number_of_sites;
}

inline int Lattice::coord_to_x(int coord){
        return coord/this->m_grid_size;
}
inline int Lattice::coord_to_y(int coord){
        return coord%this->m_grid_size;
}

void Lattice::timestep(){
        int tmp_rnd; // random number between 0 and 4*length of selector list
        int tmp_par; // which particle moves
        int tmp_dir; // which direction it moves in
        for (int n=0; n < this->m_step_attempts_per_timestep; n++)
        {
                tmp_rnd = random_int(0,4*this->m_movement_selector_length-1);
                tmp_par = this->m_movement_selector[tmp_rnd/4];
                tmp_dir = 1+tmp_rnd%4;
                // if( this->m_tracer_locations[tmp_par]->step_is_valid(tmp_dir)) {
                //      this->m_tracer_locations[tmp_par] = this->m_tracer_locations[tmp_par]->get_neighbor_by_dir(tmp_dir);
                this->m_tracers[tmp_par]->step(tmp_dir);
                // }
                // eacth timestep consists of m_step_attempts_per_timestep equal intervals, thus the time at which a given step attempt takes place is t = m_t + n / m_step_attempts_per_timestep
                //this->m_tracers[this->m_movement_selector[tmp_rnd/4]]->step(1+tmp_rnd%4);

        }
}

void Lattice::timestep_warmup(){
        int tmp_rnd; // random number between 0 and 4*length of selector list
        int tmp_par; // which particle moves
        int tmp_dir; // which direction it moves in
        for (long n=0; n < this->m_step_attempts_per_timestep; n++)
        {
                tmp_rnd = random_int(0,4*this->m_movement_selector_length-1);
                tmp_par = this->m_movement_selector[tmp_rnd/4];
                tmp_dir = 1+tmp_rnd%4;
                this->m_tracers[tmp_par]->step_warmup(tmp_dir);
        }
}
// - - - - - - - - - - - - - - - - - - - - - - - -
void Lattice::timestep_no_interaction(){
        int tmp_rnd; // random number between 0 and 4*length of selector list
        int tmp_par; // which particle moves
        int tmp_dir; // which direction it moves in
        for (long n=0; n < this->m_step_attempts_per_timestep; n++)
        {
                tmp_rnd = random_int(0,4*this->m_movement_selector_length);
                tmp_par = this->m_movement_selector[tmp_rnd/4];
                tmp_dir = 1+tmp_rnd%4;
                this->m_tracers[tmp_par]->step_unhindered(tmp_dir);
        }
}

int Lattice::get_number_of_tracers_1x1()
{
        return this->m_number_of_tracers_1x1;
}

int Lattice::get_number_of_tracers_2x2()
{
        return this->m_number_of_tracers_2x2;
}

int Lattice::get_grid_size()
{
        return this->m_grid_size;
}

int Lattice::get_tracer_size(int id)
{
        return this->m_tracers[id]->get_size();
}

int Lattice::get_tracer_position_x(int id)
{
        return this->m_tracers[id]->get_dx();
}

int Lattice::get_tracer_position_y(int id)
{
        return this->m_tracers[id]->get_dy();
}

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

std::vector<int> Lattice::get_tracer_positions()
{
        std::vector<int> tmp_tracer_positions;
        tmp_tracer_positions.reserve(4*this->m_number_of_tracers);
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
// Setup funtions for the neighbor relationships between the lattice sites
// - - - - - - - - - - - - - - - - - - - - - - - -
// for the neighbors of 1x1 type sites
// - - - - - - - - - - - - - - - - - - - - - - - -

void Lattice::set_neighbor_sites(Site_1x1 * s)
{
        int tmp_x = s->get_x();
        int tmp_y = s->get_y();
        std::vector<Site *> tmp_sites;
        tmp_sites.reserve(12);
        tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x+1,tmp_y+0)]);
        tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x+0,tmp_y+1)]);
        tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x-1,tmp_y-0)]);
        tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x-0,tmp_y-1)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+2,tmp_y+1)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+2,tmp_y+0)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+1,tmp_y-1)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+0,tmp_y-1)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x-1,tmp_y-0)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x-1,tmp_y+1)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x-0,tmp_y+2)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+1,tmp_y+2)]);
        s->set_neighbors(tmp_sites);
}

// - - - - - - - - - - - - - - - - - - - - - - - -
// for the neighbors of 2x2 type sites
// - - - - - - - - - - - - - - - - - - - - - - - -

void Lattice::set_neighbor_sites(Site_2x2 * s)
{
        int tmp_x = s->get_x();
        int tmp_y = s->get_y();
        std::vector<Site *> tmp_sites;
        tmp_sites.reserve(24);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+1,tmp_y+0)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+0,tmp_y+1)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x-1,tmp_y-0)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x-0,tmp_y-1)]);
        tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x+1,tmp_y+0)]);
        tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x+1,tmp_y-1)]);
        tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x+0,tmp_y-2)]);
        tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x-1,tmp_y-2)]);
        tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x-2,tmp_y-1)]);
        tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x-2,tmp_y+0)]);
        tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x-1,tmp_y+1)]);
        tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x-0,tmp_y+1)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+2,tmp_y+1)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+2,tmp_y+0)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+2,tmp_y-1)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+1,tmp_y-2)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+0,tmp_y-2)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x-1,tmp_y-2)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x-2,tmp_y-1)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x-2,tmp_y-0)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x-2,tmp_y+1)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x-1,tmp_y+2)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+0,tmp_y+2)]);
        tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+1,tmp_y+2)]);
        s->set_neighbors(tmp_sites);
}

// - - - - - - - - - - - - - - - - - - - - - - - -
// D E B U G G I N G  O U T P U T
// - - - - - - - - - - - - - - - - - - - - - - - -
//
// void Lattice::print_tracer_positions()
// {
//         // Prints Tracer information to std output
//         for(Tracer * tr : this->m_tracers)
//         {
//                 std::cout << "T: " << this->m_t << " ID: " << tr->get_id() << " X: " << tr->get_x() << " Y: " << tr->get_y() << " dx: " << tr->get_dx() << " dy: " << tr->get_dy() << " l^2 " << tr->get_lsquared() << " Size: " << tr->get_size() << std::endl;
//         }
// }

void Lattice::print_sites()
{
        // print all sites (1x1 and 2x2) to std output
        std::cout << "Sites:" << std::endl << "1x1" << std::endl;
        for(Site * s : this->m_sites_1x1) {
                std::cout << s->get_id() << " " << s->get_x() << " " << s->get_y() << std::endl;
        }
        std::cout << "2x2" << std::endl;
        for(Site * s : this->m_sites_2x2) {
                std::cout << s->get_id() << " " << s->get_x() << " " << s->get_y() << std::endl;
        }
}

void Lattice::print_neighbors()
{
        // print all neighbor connections, for debugging purposes
        std::cout << "Neighbors:" << std::endl;
        int dir;
        for(auto n = 0; n < this->m_neighbors.size(); n++)
        {
                auto vec1 = this->m_neighbors[n];
                std::cout << "Site: " << n << std::endl;
                for(auto vec2 : vec1) {
                        for(auto val : vec2)
                        {
                                std::cout << *val << " ";
                        }
                }
                std::cout << std::endl;
        }
}

void Lattice::print_occupation_map()
{
        std::cout << "Occupation map:";
        for(int y = 0; y < this->m_grid_size; y++)
        {
                std::cout << std::endl;
                for(int x = 0; x < this->m_grid_size; x++)
                {

                }
        }
}
