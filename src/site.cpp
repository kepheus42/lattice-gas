#include "site.hpp"

// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif

Site::Site(int id, int x, int y, int t) :
        m_id(id),
        m_x(x),
        m_y(y),
        m_is_empty(true),
        m_type(t),
        m_powers_of_two({1,2,4,8,16})
        // m_blocking_sites_per_dir(b_per_dir)
{
        this->m_neighbor_sites.reserve(4);
        this->m_blocking_sites.reserve(4);
        D(std::cout << "Creating site: " << id << " at " << x << "," << y << std::endl);
}

bool Site::step_is_valid(int dir){
        /*
           ~ returns AND over all site states in the selected direction
           ~ true if all sites are empty
           ~ false otherwise
         */
        return std::all_of(this->m_blocking_sites[dir-1].begin(), this->m_blocking_sites[dir-1].end(), []( Site * s ) -> bool {
                return s->is_empty();
        });
}

bool Site::get_blocking_site_state(int dir, int which){
        return this->m_blocking_sites[dir-1][which]->is_empty();
}

std::vector<std::vector<Site * > > Site::get_blocking_sites(){
        return this->m_blocking_sites;
}

std::vector<int> Site::get_site_correlation(){
        std::vector<int> tmp_site_corr;
        tmp_site_corr.reserve(4);
        for(auto direction : this->m_blocking_sites)
        {
                /* */
                tmp_site_corr.push_back(std::inner_product(direction.begin(),direction.end(),this->m_powers_of_two.begin(),0,std::plus<>(),[](Site * s, int i) -> int {
                        return (!(s->is_empty())) * i;
                }));
        }
        return tmp_site_corr;
}

void Site::set_neighbor_sites(std::vector<Site *> sites){
        //D( std::cout << "Site: " << this->m_id << "(" << this->m_x << "," << this->m_y << ")" << std::endl );
        //D( std::cout << "# of neighbors: "<< sites.size() << std::endl );
        //D( std::cout << "Setting neighbors to: " << std::endl );
        //D( this->db_print_vector(sites) );
        this->m_neighbor_sites = sites;
}

void Site::set_blocking_sites(std::vector<std::vector<Site *> > vec){
        this->m_blocking_sites = vec;
}

Site * Site::jump_in_direction(int dir){
        this->m_is_empty = true;
        this->m_neighbor_sites[dir-1]->set_not_empty();
        return this->m_neighbor_sites[dir-1];
}

void Site::swap_state() {
        this->m_is_empty = !this->m_is_empty;
}

void Site::set_not_empty(){
        this->m_is_empty = false;
}

bool Site::is_empty() {
        return this->m_is_empty;
}

bool * Site::get_state_ptr(){
        return &(this->m_is_empty);
}

Site * Site::get_neighbor_by_dir(int dir){
        return this->m_neighbor_sites[dir-1];
}

int Site::get_id(){
        return this->m_id;
}

int Site::get_x() {
        return this->m_x;
}

int Site::get_y() {
        return this->m_y;
}

void Site::db_print_vector(std::vector<Site *> sites){
        for(Site* s : sites) {
                std::cout << "ID: " << s->get_id() << " (" << s->get_x() << "," << s->get_y() << ")" << std::endl;
        }
}

void Site::db_print_properties(){
        std::cout << "Site: " << this->m_type << " " << this->m_id << " " << this->m_x << " " << this->m_y << std::endl;
}

void Site::print_neighbor_sites(){
        std::cout << "Site: " << this->m_id << " at " << this->m_x << "," << this->m_y << std::endl;
        for(int dir = 1; dir < 5; dir++)
        {
                std::cout << "Dir "<< dir << " : " << this->m_neighbor_sites[dir-1]->get_id() << std::endl;
        }
}

void Site::print_blocking_sites(){
        std::cout << "Site: " << this->m_id << " at " << this->m_x << "," << this->m_y << std::endl;
        for(int dir = 1; dir < 5; dir++)
        {
                std::cout << "Dir "<< dir << " : ";
                for(Site * s : this->m_blocking_sites[dir-1]) {
                        std::cout << s->get_id() << " ";
                }
                std::cout << std::endl;
        }
}
