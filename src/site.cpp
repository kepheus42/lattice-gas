#include "site.hpp"

// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif

Site::Site(int id, int x, int y, int t) : m_id(id), m_x(x), m_y(y), m_is_empty(true), m_type(t)
{
        //this->m_neighbor_sites.reserve(4);
        D(std::cout << "Creating site: " << id << " at " << x << "," << y << std::endl);
}

bool Site::step_is_valid(int dir)
{
        int tmp = 1;
        for(auto bp : this->m_blocking_sites_occ[dir-1]) { tmp *= (*bp); }
        return tmp;
}

void Site::set_neighbors(std::vector<Site *> sites){
        //D( std::cout << "Site: " << this->m_id << "(" << this->m_x << "," << this->m_y << ")" << std::endl );
        //D( std::cout << "# of neighbors: "<< sites.size() << std::endl );
        //D( std::cout << "Setting neighbors to: " << std::endl );
        //D( this->db_print_vector(sites) );
        this->m_neighbor_sites = sites;
}

void Site::set_neighbors_occ(std::vector<bool*> vec){
        this->m_neighbor_is_empty = vec;
}

void Site::set_blocking_sites_occ(std::vector<std::vector<bool*> > vec){
        this->m_blocking_sites_occ = vec;
}

Site * Site::move(int dir){
        this->m_is_empty = true;
        this->m_neighbor_is_empty[dir-1]->set_not_empty();
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
