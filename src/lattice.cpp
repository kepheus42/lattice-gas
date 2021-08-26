#include "lattice.hpp"

#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif
// To add debugging messages, use D(std::cerr << "Debugging message 1 2 3!" << std::endl; )
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
Lattice::Lattice(int timesteps,
                 int timesteps_w,
                 int grid_size,
                 int number_of_tracers_1x1,
                 int number_of_tracers_2x2,
                 double step_rate_1x1,
                 double step_rate_2x2,
                 unsigned int rng_seed) :
        m_t(0),
        m_w(0),
        m_timesteps(timesteps),
        m_timesteps_w(timesteps_w),
        m_grid_size(grid_size),
        m_number_of_sites(grid_size*grid_size),
        m_number_of_tracers(number_of_tracers_1x1+number_of_tracers_2x2),
        m_number_of_tracers_1x1(number_of_tracers_1x1),
        m_number_of_tracers_2x2(number_of_tracers_2x2),
        m_step_rate_1x1(this->m_number_of_tracers_1x1 > 0 ? step_rate_1x1 : 1.0),
        m_step_rate_2x2(this->m_number_of_tracers_2x2 > 0 ? step_rate_2x2 : 1.0),

        m_step_attempts_per_timestep((int)(this->m_number_of_tracers_1x1*this->m_step_rate_1x1)+(int)(this->m_number_of_tracers_2x2*this->m_step_rate_2x2)),
        m_movement_selector_length(this->m_number_of_tracers_1x1*(int)(1/this->m_step_rate_2x2)+this->m_number_of_tracers_2x2*(int)(1/this->m_step_rate_1x1)),
        m_movement_selector(this->m_movement_selector_length,0),

        m_dpoints(this->m_timesteps ? 9*(int)std::log10(this->m_timesteps)+this->m_timesteps/(int)std::pow(10,(int) std::log10(this->m_timesteps)) : 0),
        m_dinterval(0),

        m_rng(rng_seed),
        m_random_par(0,this->m_movement_selector_length-1),
        m_random_dir(1,4),
        // precompute these:
        m_one_over_n_1x1(1/(double)this->m_number_of_tracers_1x1),
        m_one_over_four_n_1x1(1/(double)(4*this->m_number_of_tracers_1x1)),
        m_one_over_n_2x2(1/(double)this->m_number_of_tracers_2x2),
        m_one_over_four_n_2x2(1/(double)(4*this->m_number_of_tracers_2x2))
{
        // check if there's enough space on the grid to place the tracers
        if (!this->m_grid_size || (!this->m_number_of_tracers_1x1 && !this->m_number_of_tracers_2x2))
        {
                throw std::invalid_argument("No grid or no tracers!");
        }
        if (this->m_grid_size*this->m_grid_size < this->m_number_of_tracers_1x1 + 4*this->m_number_of_tracers_2x2)
        {
                throw std::invalid_argument("Too many tracers for grid of the given size!");
        }
        else if (this->m_number_of_tracers_2x2&&(grid_size < 3))
        {
                throw std::invalid_argument("Grid is too small for 2x2 Tracer!");
        }
        //
        this->m_avg_rate_1x1.reserve(this->m_dpoints);
        this->m_avg_rate_2x2.reserve(this->m_dpoints);

        this->m_avg_lsq_1x1.reserve(this->m_dpoints);
        this->m_avg_lsq_2x2.reserve(this->m_dpoints);

        this->m_pos_1x1.reserve(this->m_number_of_tracers_1x1 ? this->m_number_of_tracers_1x1 * (this->m_dpoints+1) : 0);
        this->m_pos_2x2.reserve(this->m_number_of_tracers_2x2 ? this->m_number_of_tracers_2x2 * (this->m_dpoints+1) : 0);

        this->m_displacements_1x1.reserve(this->m_number_of_tracers_1x1 ? 2 * this->m_number_of_tracers_1x1 * this->m_dpoints : 0);
        this->m_displacements_2x2.reserve(this->m_number_of_tracers_2x2 ? 2 * this->m_number_of_tracers_2x2 * this->m_dpoints : 0);

        this->m_site_corr_1x1.reserve(this->m_number_of_tracers_2x2 ?  4*this->m_dpoints : 1*this->m_dpoints );
        this->m_site_corr_2x2.reserve(this->m_number_of_tracers_1x1 ?  9*this->m_dpoints : 3*this->m_dpoints );

        //this->m_site_corr_1x1_counter.reserve(this->m_number_of_tracers_1x1 ?  4*this->m_number_of_tracers_1x1*this->m_dpoints : 0 );
        //this->m_site_corr_2x2_counter.reserve(this->m_number_of_tracers_2x2 ?  4*this->m_number_of_tracers_2x2*this->m_dpoints : 0 );

        this->m_sublattice_concentrations.reserve(this->m_number_of_tracers_2x2 ? 5 * this->m_dpoints : 0);

        // only allocate memory for sites if the respective tracer species is present
        this->m_sites_1x1.reserve(this->m_number_of_tracers_1x1 > 0 ? this->m_number_of_sites : 0);
        this->m_sites_2x2.reserve(this->m_number_of_tracers_2x2 > 0 ? this->m_number_of_sites : 0);
        // store sites_2x2 ordered by sublattice (1..4)
        this->m_sites_2x2_by_sublattice.reserve(this->m_number_of_tracers_2x2 > 0 ? 4 : 0);
        this->m_sites_2x2_by_sublattice[0].reserve(this->m_number_of_tracers_2x2 > 0 ? this->m_number_of_sites/4 : 0);
        this->m_sites_2x2_by_sublattice[1].reserve(this->m_number_of_tracers_2x2 > 0 ? this->m_number_of_sites/4 : 0);
        this->m_sites_2x2_by_sublattice[2].reserve(this->m_number_of_tracers_2x2 > 0 ? this->m_number_of_sites/4 : 0);
        this->m_sites_2x2_by_sublattice[3].reserve(this->m_number_of_tracers_2x2 > 0 ? this->m_number_of_sites/4 : 0);
        // this->m_sites.reserve(this->m_sites_1x1.capacity()+this->m_sites_2x2.capacity());
        // allocate sufficient memory for the tracers
        this->m_tracers.reserve(this->m_number_of_tracers);
        this->m_tracers_1x1.reserve(this->m_number_of_tracers_1x1);
        this->m_tracers_2x2.reserve(this->m_number_of_tracers_2x2);
        //D( std::cout << "Capacity of m_sites_1x1: " << this->m_sites_1x1.capacity() << std::endl << "Capacity of m_sites_2x2: " << this->m_sites_2x2.capacity() << std::endl );
        //D( std::cout << "Parameters of lattice: R1 " << this->m_step_rate_1x1 << " R2 " << this->m_step_rate_2x2 << std::endl );
        this->setup_sites();
        this->setup_movement_selection_list();
        this->setup_tracers();

}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
// Setup functions
void Lattice::setup_sites(){
        D( std::cout << "Setting up sites..." << std::endl );
        // init sites for 1x1 tracers
        D( std::cout << "1x1: " << std::endl );
        if(this->m_number_of_tracers_1x1)
        {
                for(int tmp_y = 0; tmp_y < this->m_grid_size; tmp_y++)
                {
                        for(int tmp_x = 0; tmp_x < this->m_grid_size; tmp_x++)
                        {
                                //D( std::cout << "X: " << tmp_y << " Y: " << tmp_x << std::endl );
                                // this->m_sites.push_back(new Site(this->coord(tmp_x,tmp_y), tmp_y, tmp_x, 1));
                                // this->m_sites_1x1.push_back(this->m_sites.back());
                                this->m_sites_1x1.push_back(new Site(this->coord(tmp_x,tmp_y), tmp_x, tmp_y, 1));
                        }
                }
        }
        D( std::cout << "Number of 1x1 sites created: " << this->m_sites_1x1.size() << std::endl );
        // init sites for 2x2 tracers
        D( std::cout << "2x2: " << std::endl );
        if(this->m_number_of_tracers_2x2)
        {
                for(int tmp_y = 0; tmp_y < this->m_grid_size; tmp_y++)
                {
                        for(int tmp_x = 0; tmp_x < this->m_grid_size; tmp_x++)
                        {
                                //D( std::cout << "X: " << tmp_y << " Y: " << tmp_x << std::endl );
                                // this->m_sites.push_back(new Site(this->coord(tmp_x,tmp_y), tmp_y, tmp_x, 2));
                                // this->m_sites_2x2.push_back(this->m_sites.back());
                                this->m_sites_2x2.push_back(new Site(this->coord(tmp_x,tmp_y), tmp_x, tmp_y, 2));
                                this->m_sites_2x2_by_sublattice[(tmp_x%2)+2*(tmp_y%2)].push_back(this->m_sites_2x2.back());
                        }
                }
        }
        D( std::cout << "Number of 2x2 sites created: " << this->m_sites_2x2.size() << std::endl );
        // link sites
        int tmp_y = 0;
        int tmp_x = 0;

        std::vector<Site *> tmp_sites;
        tmp_sites.reserve(4);

        std::vector<std::vector<Site *> > tmp_blocking;
        tmp_blocking.reserve(4);

        for(Site * s : this->m_sites_1x1) {
                tmp_x = s->get_x();
                tmp_y = s->get_y();
                //
                tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x+1,tmp_y+0)]);
                tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x+0,tmp_y+1)]);
                tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x-1,tmp_y+0)]);
                tmp_sites.push_back(this->m_sites_1x1[this->coord(tmp_x+0,tmp_y-1)]);
                //
                s->set_neighbor_sites(tmp_sites);
                tmp_sites.clear();
                if(!this->m_number_of_tracers_2x2) {
                        // dir = 1
                        tmp_blocking.push_back({this->m_sites_1x1[this->coord(tmp_x+1,tmp_y+0)]});
                        // dir = 2
                        tmp_blocking.push_back({this->m_sites_1x1[this->coord(tmp_x+0,tmp_y+1)]});
                        // dir = 3
                        tmp_blocking.push_back({this->m_sites_1x1[this->coord(tmp_x-1,tmp_y+0)]});
                        // dir = 4
                        tmp_blocking.push_back({this->m_sites_1x1[this->coord(tmp_x+0,tmp_y-1)]});
                }
                else {
                        // dir = 1
                        tmp_blocking.push_back({this->m_sites_1x1[this->coord(tmp_x+1,tmp_y+0)],
                                                this->m_sites_2x2[this->coord(tmp_x+2,tmp_y-1)],
                                                this->m_sites_2x2[this->coord(tmp_x+2,tmp_y+0)]});
                        // dir = 2
                        tmp_blocking.push_back({this->m_sites_1x1[this->coord(tmp_x+0,tmp_y+1)],
                                                this->m_sites_2x2[this->coord(tmp_x+1,tmp_y+1)],
                                                this->m_sites_2x2[this->coord(tmp_x+0,tmp_y+1)]});
                        // dir = 3
                        tmp_blocking.push_back({this->m_sites_1x1[this->coord(tmp_x-1,tmp_y+0)],
                                                this->m_sites_2x2[this->coord(tmp_x-1,tmp_y+0)],
                                                this->m_sites_2x2[this->coord(tmp_x-1,tmp_y-1)]});
                        // dir = 4
                        tmp_blocking.push_back({this->m_sites_1x1[this->coord(tmp_x+0,tmp_y-1)],
                                                this->m_sites_2x2[this->coord(tmp_x+0,tmp_y-2)],
                                                this->m_sites_2x2[this->coord(tmp_x+1,tmp_y-2)]});
                }
                s->set_blocking_sites(tmp_blocking);
                tmp_blocking.clear();
        }
        for(Site * s : this->m_sites_2x2) {
                tmp_x = s->get_x();
                tmp_y = s->get_y();
                //
                tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+1,tmp_y+0)]);
                tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+0,tmp_y+1)]);
                tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x-1,tmp_y+0)]);
                tmp_sites.push_back(this->m_sites_2x2[this->coord(tmp_x+0,tmp_y-1)]);
                //
                s->set_neighbor_sites(tmp_sites);
                tmp_sites.clear();
                if(!this->m_number_of_tracers_1x1) {
                        // dir = 1
                        tmp_blocking.push_back({this->m_sites_2x2[this->coord(tmp_x+2,tmp_y-1)],
                                                this->m_sites_2x2[this->coord(tmp_x+2,tmp_y+0)],
                                                this->m_sites_2x2[this->coord(tmp_x+2,tmp_y+1)]});
                        // dir = 2
                        tmp_blocking.push_back({this->m_sites_2x2[this->coord(tmp_x+1,tmp_y+2)],
                                                this->m_sites_2x2[this->coord(tmp_x+0,tmp_y+2)],
                                                this->m_sites_2x2[this->coord(tmp_x-1,tmp_y+2)]});
                        // dir = 3
                        tmp_blocking.push_back({this->m_sites_2x2[this->coord(tmp_x-2,tmp_y+1)],
                                                this->m_sites_2x2[this->coord(tmp_x-2,tmp_y+0)],
                                                this->m_sites_2x2[this->coord(tmp_x-2,tmp_y-1)]});
                        // dir = 4
                        tmp_blocking.push_back({this->m_sites_2x2[this->coord(tmp_x-1,tmp_y-2)],
                                                this->m_sites_2x2[this->coord(tmp_x+0,tmp_y-2)],
                                                this->m_sites_2x2[this->coord(tmp_x+1,tmp_y-2)]});
                }
                else{
                        // dir = 1
                        tmp_blocking.push_back({this->m_sites_1x1[this->coord(tmp_x+1,tmp_y+0)],
                                                this->m_sites_1x1[this->coord(tmp_x+1,tmp_y+1)],
                                                this->m_sites_2x2[this->coord(tmp_x+2,tmp_y-1)],
                                                this->m_sites_2x2[this->coord(tmp_x+2,tmp_y+0)],
                                                this->m_sites_2x2[this->coord(tmp_x+2,tmp_y+1)]});
                        // dir = 2
                        tmp_blocking.push_back({this->m_sites_1x1[this->coord(tmp_x+0,tmp_y+2)],
                                                this->m_sites_1x1[this->coord(tmp_x-1,tmp_y+2)],
                                                this->m_sites_2x2[this->coord(tmp_x+1,tmp_y+2)],
                                                this->m_sites_2x2[this->coord(tmp_x+0,tmp_y+2)],
                                                this->m_sites_2x2[this->coord(tmp_x-1,tmp_y+2)]});
                        // dir = 3
                        tmp_blocking.push_back({this->m_sites_1x1[this->coord(tmp_x-2,tmp_y+1)],
                                                this->m_sites_1x1[this->coord(tmp_x-2,tmp_y+0)],
                                                this->m_sites_2x2[this->coord(tmp_x-2,tmp_y+1)],
                                                this->m_sites_2x2[this->coord(tmp_x-2,tmp_y+0)],
                                                this->m_sites_2x2[this->coord(tmp_x-2,tmp_y-1)]});
                        // dir = 4
                        tmp_blocking.push_back({this->m_sites_1x1[this->coord(tmp_x-1,tmp_y-1)],
                                                this->m_sites_1x1[this->coord(tmp_x+0,tmp_y-1)],
                                                this->m_sites_2x2[this->coord(tmp_x-1,tmp_y-2)],
                                                this->m_sites_2x2[this->coord(tmp_x+0,tmp_y-2)],
                                                this->m_sites_2x2[this->coord(tmp_x+1,tmp_y-2)]});
                }
                s->set_blocking_sites(tmp_blocking);
                tmp_blocking.clear();
        }
        D( std::cout << "Finished Neighbor setup..." << std::endl );
        D( std::cout << "Size of m_sites_1x1: " << this->m_sites_1x1.size() << std::endl << "Size of m_sites_2x2: " << this->m_sites_2x2.size() << std::endl);
        D( this->print_sites() );
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Lattice::setup_movement_selection_list(){
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
        // D( std::cout << "Size of m_movement_selector: " << this->m_movement_selector.size() << std::endl << "Step attempts per timestep: " << this->m_step_attempts_per_timestep << std::endl );
        // D( print_vector(this->m_movement_selector) );
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Lattice::setup_tracers(){
        int tmp_id = 0;
        int tmp_y = 0;
        int tmp_x = 0;
        int tmp_pos = 0;
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
                std::shuffle(tmp_start_positions_2x2.begin(),tmp_start_positions_2x2.end(),this->m_rng);
                // shuffle_vector(tmp_start_positions_2x2);
                //D( std::cout << "Start positions (shuffled)" << std::endl );
                //D( this->print_vector(tmp_start_positions_2x2) );
                for (int n = 0; n < this->m_number_of_tracers_2x2; n++)
                {
                        tmp_pos = tmp_start_positions_2x2[n];
                        tmp_x = this->coord_to_x(tmp_pos);
                        tmp_y = this->coord_to_y(tmp_pos);
                        tmp_site = this->m_sites_2x2[tmp_pos];
                        // create new tracer at starting site
                        this->m_tracers.push_back(new Tracer(tmp_id,2,tmp_site));
                        this->m_tracers_2x2.push_back(this->m_tracers.back());
                        // increment id counter,
                        tmp_id++;
                        tmp_occupied_sites[this->coord(tmp_x+0,tmp_y+0)] = 1;
                        tmp_occupied_sites[this->coord(tmp_x+1,tmp_y+0)] = 1;
                        tmp_occupied_sites[this->coord(tmp_x+0,tmp_y+1)] = 1;
                        tmp_occupied_sites[this->coord(tmp_x+1,tmp_y+1)] = 1;
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
                std::shuffle(tmp_start_positions_1x1.begin(),tmp_start_positions_1x1.end(),this->m_rng);
                //D( std::cout << "Start positions (shuffled)" << std::endl );
                //D( this->print_vector(tmp_start_positions_1x1) );
                for (int n = 0; n < this->m_number_of_tracers_1x1; n++)
                {
                        tmp_pos = tmp_start_positions_1x1[n];
                        tmp_site = this->m_sites_1x1[tmp_pos];
                        // create new tracer at starting site
                        this->m_tracers.push_back(new Tracer(tmp_id,1,tmp_site));
                        this->m_tracers_1x1.push_back(this->m_tracers.back());
                        // increment id counter,
                        tmp_id++;
                }
                //D( std::cout << "Number of Tracers 1x1: " << this->m_tracers_1x1.size() << std::endl);
        }
        // Why did I implement storage of start position?
        // for(Tracer * tr : this->m_tracers_1x1) { this->m_pos_1x1.push_back(tr->get_pos()); }
        // for(Tracer * tr : this->m_tracers_2x2) { this->m_pos_2x2.push_back(tr->get_pos()); }
        // D( std::cout << "Number of Tracers: " << this->m_tracers.size() << std::endl);
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
// Transforms x,y to linear coordinate and back
inline int Lattice::coord(int x, int y){
        // compute linear coordinate from (x,y), with periodic boundary conditions enforced
        return ( ( y + this->m_grid_size ) % this->m_grid_size ) * this->m_grid_size + ( ( x + this->m_grid_size ) % this->m_grid_size );
}
inline int Lattice::coord_to_y(int coord){
        return coord/this->m_grid_size;
}
inline int Lattice::coord_to_x(int coord){
        return coord%this->m_grid_size;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Lattice::warmup(){
        while(this->m_w < this->m_timesteps_w) { this->timestep_warmup(); }
}
void Lattice::evolve(){
        this->store_initial_positions();
        while(this->m_t < this->m_timesteps) { this->timestep(); }
}
void Lattice::evolve_no_interaction(){
        this->store_initial_positions();
        while(this->m_t < this->m_timesteps) {  this->timestep_no_interaction(); }
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
inline void Lattice::timestep(){
        for (int n=0; n < this->m_step_attempts_per_timestep; n++)
        {
                this->m_tracers[this->m_movement_selector[this->m_random_par(this->m_rng)]]->step(this->m_random_dir(this->m_rng));
        }
        this->m_t++;
        D( this->print_tracers() );
        if(this->m_t % this->m_dinterval) { return; }
        this->update_data();
        // D( this->print_tracer_positions() );
        // D( this->print_site_states() );
}
inline void Lattice::timestep_warmup(){
        for (long n=0; n < this->m_step_attempts_per_timestep; n++)
        {
                this->m_tracers[this->m_movement_selector[this->m_random_par(this->m_rng)]]->step_warmup(this->m_random_dir(this->m_rng));
        }
        this->m_w++;
}
inline void Lattice::timestep_no_interaction(){
        for (long n=0; n < this->m_step_attempts_per_timestep; n++)
        {
                this->m_tracers[this->m_movement_selector[this->m_random_par(this->m_rng)]]->step_unhindered(this->m_random_dir(this->m_rng));
        }
        this->m_t++;
        if(this->m_t % this->m_dinterval) { return; }
        this->update_data();
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Lattice::update_data(){
        //
        double tmp_sum_stp_1x1 = 0;
        double tmp_sum_lsq_1x1 = 0;
        double tmp_sum_stp_2x2 = 0;
        double tmp_sum_lsq_2x2 = 0;
        double tmp_one_over_t = 1/(double)this->m_t;
        std::vector<int> tmp_site_corr_1x1;
        std::vector<int> tmp_site_corr_2x2;
        tmp_site_corr_1x1.reserve(4*this->m_number_of_tracers_1x1);
        tmp_site_corr_2x2.reserve(4*this->m_number_of_tracers_2x2);
        //
        this->m_dinterval = (int)std::pow(10,(int)std::log10(this->m_t));
        // double tmp_avg_rate_div = this->m_t > 1 ? (int)std::pow(10,(int)std::log10(this->m_t-1)) : 1.0;
        D( std::cout << this->m_t << " Data Storage Interval: " << this->m_dinterval << std::endl );
        //
        for(Tracer * tr : this->m_tracers_1x1) {
                tmp_sum_stp_1x1 += tr->get_steps_taken();
                tmp_sum_lsq_1x1 += tr->get_lsq();
                this->m_pos_1x1.push_back(tr->get_pos());
                for(int s : tr->get_site_correlation()) { tmp_site_corr_1x1.push_back(s); }
        }
        for(Tracer * tr : this->m_tracers_2x2) {
                tmp_sum_stp_2x2 += tr->get_steps_taken();
                tmp_sum_lsq_2x2 += tr->get_lsq();
                this->m_pos_2x2.push_back(tr->get_pos());
                for(int s : tr->get_site_correlation()) { tmp_site_corr_2x2.push_back(s); }
        }
        //
        if(this->m_number_of_tracers_1x1) {
                this->m_avg_rate_1x1.push_back(tmp_sum_stp_1x1*tmp_one_over_t*this->m_one_over_n_1x1);
                this->m_avg_lsq_1x1.push_back(tmp_sum_lsq_1x1*this->m_one_over_n_1x1);
                if(!this->m_number_of_tracers_2x2)
                {
                        std::vector<unsigned long> tmp_site_corr_1x1_counter;
                        tmp_site_corr_1x1_counter.reserve(1);
                        tmp_site_corr_1x1_counter.push_back(std::count(tmp_site_corr_1x1.begin(),
                                                                       tmp_site_corr_1x1.end(),1));
                        for(unsigned long n : tmp_site_corr_1x1_counter) { this->m_site_corr_1x1.push_back(((double)n)*this->m_one_over_four_n_1x1); }

                }
                else
                {
                        std::vector<unsigned long> tmp_site_corr_1x1_counter;
                        tmp_site_corr_1x1_counter.reserve(3);
                        tmp_site_corr_1x1_counter.push_back(std::count(tmp_site_corr_1x1.begin(),
                                                                       tmp_site_corr_1x1.end(),1));
                        tmp_site_corr_1x1_counter.push_back(std::count(tmp_site_corr_1x1.begin(),
                                                                       tmp_site_corr_1x1.end(),2));
                        tmp_site_corr_1x1_counter.push_back(std::count(tmp_site_corr_1x1.begin(),
                                                                       tmp_site_corr_1x1.end(),4));
                        for(unsigned long n : tmp_site_corr_1x1_counter) { this->m_site_corr_1x1.push_back((double)n*this->m_one_over_four_n_1x1); }

                }
        }
        if(this->m_number_of_tracers_2x2) {
                this->m_avg_rate_2x2.push_back(tmp_sum_stp_2x2*tmp_one_over_t*this->m_one_over_n_2x2);
                this->m_avg_lsq_2x2.push_back(tmp_sum_lsq_2x2*this->m_one_over_n_2x2);
                if(!this->m_number_of_tracers_1x1)
                {
                        std::vector<unsigned long> tmp_site_corr_2x2_counter;
                        tmp_site_corr_2x2_counter.reserve(4);
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),1));
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),2));
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),4));
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),5));
                        for(unsigned long n : tmp_site_corr_2x2_counter) { this->m_site_corr_2x2.push_back((double)n*this->m_one_over_four_n_2x2); }
                }
                else
                {
                        std::vector<unsigned long> tmp_site_corr_2x2_counter;
                        tmp_site_corr_2x2_counter.reserve(9);
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),1));
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),2));
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),4));
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),8));
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),16));
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),3));
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),6));
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),17));
                        tmp_site_corr_2x2_counter.push_back(std::count(tmp_site_corr_2x2.begin(),
                                                                       tmp_site_corr_2x2.end(),20));
                        for(unsigned long n : tmp_site_corr_2x2_counter) { this->m_site_corr_2x2.push_back((double)n*this->m_one_over_four_n_2x2); }
                }
        }
        //
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Lattice::store_initial_positions()
{
        for(Tracer * tr : this->m_tracers_1x1) {
                this->m_pos_1x1.push_back(tr->get_pos());
        }
        for(Tracer * tr : this->m_tracers_2x2) {
                this->m_pos_2x2.push_back(tr->get_pos());
        }
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
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
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
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
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<double> Lattice::get_avg_rate_1x1(){
        return this->m_avg_rate_1x1;
}
std::vector<double> Lattice::get_avg_rate_2x2(){
        return this->m_avg_rate_2x2;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<double> Lattice::get_avg_lsq_1x1(){
        return this->m_avg_lsq_1x1;
}
std::vector<double> Lattice::get_avg_lsq_2x2(){
        return this->m_avg_lsq_2x2;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<double> Lattice::get_avg_corr_1x1(){
        std::vector<unsigned long long> tmp_corr_count(84,0);
        for(Tracer * tr : this->m_tracers_1x1) {
                std::transform(tmp_corr_count.begin(),tmp_corr_count.end(),tr->get_correlations().begin(),tmp_corr_count.begin(),std::plus<>{});
        }
        // calculate norm for 1-3 step sequences
        double tmp_norm_1 = 1/(double)std::accumulate(tmp_corr_count.begin()+0,tmp_corr_count.begin()+4,0);
        double tmp_norm_2 = 1/(double)std::accumulate(tmp_corr_count.begin()+4,tmp_corr_count.begin()+20,0);
        double tmp_norm_3 = 1/(double)std::accumulate(tmp_corr_count.begin()+20,tmp_corr_count.end(),0);
        // cast count vector to double
        std::vector<double> tmp_avg_corr(tmp_corr_count.begin(),tmp_corr_count.end());
        // apply norm
        std::transform(tmp_corr_count.begin()+0,
                       tmp_corr_count.begin()+4,
                       tmp_avg_corr.begin(),
                       std::bind(std::multiplies<double>(), std::placeholders::_1, tmp_norm_1));
        std::transform(tmp_corr_count.begin()+4,
                       tmp_corr_count.begin()+20,
                       tmp_avg_corr.begin()+4,
                       std::bind(std::multiplies<double>(), std::placeholders::_1, tmp_norm_2));
        std::transform(tmp_corr_count.begin()+20,
                       tmp_corr_count.end(),
                       tmp_avg_corr.begin()+20,
                       std::bind(std::multiplies<double>(), std::placeholders::_1, tmp_norm_3));
        // return result
        return tmp_avg_corr;
}
std::vector<double> Lattice::get_avg_corr_2x2(){
        std::vector<unsigned long long> tmp_corr_count(84,0);
        for(Tracer * tr : this->m_tracers_2x2) {
                std::transform(tmp_corr_count.begin(),tmp_corr_count.end(),tr->get_correlations().begin(),tmp_corr_count.begin(),std::plus<>{});
        }
        // calculate norm for 1-3 step sequences
        double tmp_norm_1 = 1/(double)std::accumulate(tmp_corr_count.begin()+0,tmp_corr_count.begin()+4,0);
        double tmp_norm_2 = 1/(double)std::accumulate(tmp_corr_count.begin()+4,tmp_corr_count.begin()+20,0);
        double tmp_norm_3 = 1/(double)std::accumulate(tmp_corr_count.begin()+20,tmp_corr_count.end(),0);
        // cast count vector to double
        std::vector<double> tmp_avg_corr(tmp_corr_count.begin(),tmp_corr_count.end());
        // apply norm
        std::transform(tmp_corr_count.begin()+0,
                       tmp_corr_count.begin()+4,
                       tmp_avg_corr.begin(),
                       std::bind(std::multiplies<double>(), std::placeholders::_1, tmp_norm_1));
        std::transform(tmp_corr_count.begin()+4,
                       tmp_corr_count.begin()+20,
                       tmp_avg_corr.begin()+4,
                       std::bind(std::multiplies<double>(), std::placeholders::_1, tmp_norm_2));
        std::transform(tmp_corr_count.begin()+20,
                       tmp_corr_count.end(),
                       tmp_avg_corr.begin()+20,
                       std::bind(std::multiplies<double>(), std::placeholders::_1, tmp_norm_3));
        // return result
        return tmp_avg_corr;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<int> Lattice::get_pos_1x1(){
        return this->m_pos_1x1;
}
std::vector<int> Lattice::get_pos_2x2(){
        return this->m_pos_2x2;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<double> Lattice::get_site_corr_1x1(){
        return this->m_site_corr_1x1;
}
std::vector<double> Lattice::get_site_corr_2x2(){
        return this->m_site_corr_2x2;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Lattice::print_tracers()
{
        // Prints Tracer information to std output
        for(Tracer * tr : this->m_tracers)
        {
                printf(" %i,%i at (%i,%i) D: %3.2f R: %3.2f\n",tr->get_type(),tr->get_id(),tr->get_x(),tr->get_y(),(double)tr->get_lsq()/(double)this->m_t,(double)tr->get_steps_taken()/(double)this->m_t);
        }
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Lattice::print_tracer_positions()
{
        // Prints Tracer information to std output
        for(Tracer * tr : this->m_tracers)
        {
                printf(" %i,%i at %i,%i\n",tr->get_type(),tr->get_id(),tr->get_x(),tr->get_y());
        }
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Lattice::print_site_1x1_states()
{
        for(Site * s : this->m_sites_1x1)
        {
                printf(" %i,%i state: %i \n",s->get_x(),s->get_y(),s->is_empty());
        }
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Lattice::print_site_2x2_states()
{
        for(Site * s : this->m_sites_2x2)
        {
                printf(" %i,%i state: %i \n",s->get_x(),s->get_y(),s->is_empty());
        }
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Lattice::print_sites()
{
        std::cout << "Site Connections:" << std::endl;
        if(this->m_number_of_tracers_1x1)
        {
                std::cout << "1x1 Neighbors Sites" << std::endl;
                for(Site * s : this->m_sites_1x1)
                {
                        s->print_neighbor_sites();
                }
                std::cout << "1x1 Blocking Sites" << std::endl;
                for(Site * s : this->m_sites_1x1)
                {
                        s->print_blocking_sites();
                }
        }
        if(this->m_number_of_tracers_2x2)
        {
                std::cout << "2x2 Neighbors Sites" << std::endl;
                for(Site * s : this->m_sites_2x2)
                {
                        s->print_neighbor_sites();
                }
                std::cout << "2x2 Blocking Sites" << std::endl;
                for(Site * s : this->m_sites_2x2)
                {
                        s->print_blocking_sites();
                }
        }
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
template <typename T>
void Lattice::db_print_vector(std::vector<T> vec){
        std::cout << "db_print_vector output:" << std::endl;
        for (T e : vec) {
                std::cout << e << ",";
        }
}
