#include "site.hpp"
#include "global.hpp"
#include <vector>
#include <iostream>


// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif

Site::Site(int id, int x, int y, int t) : m_id(id), m_x(x), m_y(y), m_is_empty(true), m_type(t)
{
        //this->m_neighbors.reserve(4);
        D(std::cout << "Creating site: " << id << " at " << x << "," << y << std::endl);
}

bool Site::step_is_valid(int dir)
{
        int tmp = 1;
        for(auto bp : this->m_neighbor_occupancy[dir-1]) { tmp *= (*bp); }
        return tmp;
}

void Site::set_neighbors(std::vector<Site *> sites){
        //D( std::cout << "Site: " << this->m_id << "(" << this->m_x << "," << this->m_y << ")" << std::endl );
        //D( std::cout << "# of neighbors: "<< sites.size() << std::endl );
        //D( std::cout << "Setting neighbors to: " << std::endl );
        //D( this->db_print_vector(sites) );
        this->m_neighbors = sites;
}

void Site::set_neighbor_occupancy(std::vector<std::vector<bool*> > vec){
        this->m_neighbor_occupancy = vec;
}

Site * Site::move_to(int dir){
        this->swap_state();
        this->m_neighbors[dir-1]->swap_state();
        return this->m_neighbors[dir-1];
}

void Site::swap_state() {
        this->m_is_empty = !this->m_is_empty;
}

bool Site::is_empty() {
        return this->m_is_empty;
}

bool * Site::get_state_ptr(){
        return &(this->m_is_empty);
}

Site * Site::get_neighbor_by_dir(int dir){
        return this->m_neighbors[dir-1];
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

Site_1x1::Site_1x1(int id, int x, int y) : Site(id, x, y)
{
        this->m_type = 1;
        D(std::cout << "Creating site_1x1: " << id << " at " << x << "," << y << std::endl);
        this->m_neighbors.reserve(12);
}

bool Site_1x1::step_is_valid(int dir){
/* = = = = = = = = = = = = = = = = = = =
   blocking sites coordinates
   = = = = = = = = = = = = = = = = = = =
   m_neighbors[0]  @  X+1  Y
   m_neighbors[1]  @  X    Y+1
   m_neighbors[2]  @  X-1  Y
   m_neighbors[3]  @  X    Y-1
   = = = = = = = = = = = = = = = = = = =
   m_neighbors[0]  @  X+2  Y-1
   m_neighbors[1]  @  X+2  Y
   m_neighbors[2]  @  X+1  Y+1
   m_neighbors[3]  @  X    Y+1
   m_neighbors[4]  @  X-1  Y
   m_neighbors[5]  @  X-1  Y-1
   m_neighbors[6]  @  X    Y-2
   m_neighbors[7]  @  X+1  Y-2
   = = = = = = = = = = = = = = = = = = = */

        switch(dir) {
        case 1: return (this->m_neighbors[0]->is_empty() &&
                        this->m_neighbors[4]->is_empty() &&
                        this->m_neighbors[5]->is_empty());
        case 2: return (this->m_neighbors[1]->is_empty() &&
                        this->m_neighbors[6]->is_empty() &&
                        this->m_neighbors[7]->is_empty());
        case 3: return (this->m_neighbors[2]->is_empty() &&
                        this->m_neighbors[8]->is_empty() &&
                        this->m_neighbors[9]->is_empty());
        case 4: return (this->m_neighbors[3]->is_empty() &&
                        this->m_neighbors[10]->is_empty() &&
                        this->m_neighbors[11]->is_empty());
        default: return false;
        }
}

Site_2x2::Site_2x2(int id, int x, int y) : Site(id, x, y)
{
        this->m_type = 2;
        D(std::cout << "Creating site_2x2: " << id << " at " << x << "," << y << std::endl);
        this->m_neighbors.reserve(24);
}

bool Site_2x2::step_is_valid(int dir){
/* = = = = = = = = = = = = = = = = = = =
   blocking sites coordinates
   = = = = = = = = = = = = = = = = = = =
   m_neighbors[0]  @  X+1  Y
   m_neighbors[1]  @  X+1  Y+1
   m_neighbors[2]  @  X    Y+2
   m_neighbors[3]  @  X-1  Y+2
   m_neighbors[4]  @  X-2  Y+1
   m_neighbors[5]  @  X-2  Y
   m_neighbors[6]  @  X-1  Y-1
   m_neighbors[7]  @  X    Y-1
   = = = = = = = = = = = = = = = = = = =
   m_neighbors[0]  @  X+2  Y-1
   m_neighbors[1]  @  X+2  Y
   m_neighbors[2]  @  X+2  Y+1
   m_neighbors[3]  @  X+1  Y+2
   m_neighbors[4]  @  X    Y+2
   m_neighbors[5]  @  X-1  Y+2
   m_neighbors[6]  @  X-2  Y+1
   m_neighbors[7]  @  X-2  Y
   m_neighbors[8]  @  X-2  Y-1
   m_neighbors[9]  @  X-1  Y-2
   m_neighbors[10] @  X    Y-2
   m_neighbors[11] @  X+1  Y-2
   = = = = = = = = = = = = = = = = = = = */
        switch(dir) {
        case 1: return (this->m_neighbors[4]->is_empty() &&
                        this->m_neighbors[5]->is_empty() &&
                        this->m_neighbors[12]->is_empty() &&
                        this->m_neighbors[13]->is_empty() &&
                        this->m_neighbors[14]->is_empty());
        case 2: return (this->m_neighbors[6]->is_empty() &&
                        this->m_neighbors[7]->is_empty() &&
                        this->m_neighbors[15]->is_empty() &&
                        this->m_neighbors[16]->is_empty() &&
                        this->m_neighbors[17]->is_empty());
        case 3: return (this->m_neighbors[8]->is_empty() &&
                        this->m_neighbors[9]->is_empty() &&
                        this->m_neighbors[18]->is_empty() &&
                        this->m_neighbors[19]->is_empty() &&
                        this->m_neighbors[20]->is_empty());
        case 4: return (this->m_neighbors[10]->is_empty() &&
                        this->m_neighbors[11]->is_empty() &&
                        this->m_neighbors[21]->is_empty() &&
                        this->m_neighbors[22]->is_empty() &&
                        this->m_neighbors[23]->is_empty());
        default: return false;
        }
}
