#include "site.hpp"
#include "global.hpp"
#include <vector>

// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif

Site::Site(int id, int x, int y) : m_id(id), m_x(x), m_y(y), m_is_empty(true), m_neighbors(4)
{

}

void Site::set_neighbors(std::vector<Site*> sites){
        this->m_neighbors = sites;
}

void Site::set_neighbors_1x1(std::vector<Site*> sites) {
        this->m_neighbors_1x1 = sites;
}

void Site::set_neighbors_2x2(std::vector<Site*> sites) {
        this->m_neighbors_2x2 = sites;
}

void Site::set_empty() {
        this->m_is_empty = true;
}

void Site::set_occupied() {
        this->m_is_empty = false;
}

bool Site::is_empty() {
        return this->m_is_empty;
}

bool Site::is_blocked() {
        return !this->m_is_empty;
}

Site * Site::get_neighbor_by_dir(int dir){
        return this->m_neighbors[dir-1];
}

virtual int Site::step_is_invalid(dir){

}

int Site::pos_x() {
        return this->m_x;
}

int Site::pos_y() {
        return this->m_y;
}

Site_1x1::Site_1x1(int id, int x, int y) : Site(id, x, y), m_neighbors_1x1(4), m_neighbors_2x2(8)
{
        //...
}

bool Site_1x1::step_is_valid(int dir){
/* = = = = = = = = = = = = = = = = = = =
   blocking sites coordinates
   = = = = = = = = = = = = = = = = = = =
   m_neighbors_1x1[0]  @  X+1  Y
   m_neighbors_1x1[1]  @  X    Y+1
   m_neighbors_1x1[2]  @  X-1  Y
   m_neighbors_1x1[3]  @  X    Y-1
   = = = = = = = = = = = = = = = = = = =
   m_neighbors_2x2[0]  @  X+2  Y-1
   m_neighbors_2x2[1]  @  X+2  Y
   m_neighbors_2x2[2]  @  X+1  Y+1
   m_neighbors_2x2[3]  @  X    Y+1
   m_neighbors_2x2[4]  @  X-1  Y
   m_neighbors_2x2[5]  @  X-1  Y-1
   m_neighbors_2x2[6]  @  X    Y-2
   m_neighbors_2x2[7]  @  X+1  Y-2
   = = = = = = = = = = = = = = = = = = = */

        switch(dir) {
        case 1: return (this->m_neighbors_1x1[0]->is_empty() &&
                        this->m_neighbors_2x2[0]->is_empty() &&
                        this->m_neighbors_2x2[1]->is_empty());
        case 2: return (this->m_neighbors_1x1[1]->is_empty() &&
                        this->m_neighbors_2x2[2]->is_empty() &&
                        this->m_neighbors_2x2[3]->is_empty());
        case 3: return (this->m_neighbors_1x1[2]->is_empty() &&
                        this->m_neighbors_2x2[4]->is_empty() &&
                        this->m_neighbors_2x2[5]->is_empty());
        case 4: return (this->m_neighbors_1x1[3]->is_empty() &&
                        this->m_neighbors_2x2[6]->is_empty() &&
                        this->m_neighbors_2x2[7]->is_empty());
        default: return false;
        }
}

int Site_1x1::step_is_invalid(dir){
        switch(dir)
        {
        case 1: return(1*int(this->m_neighbors_1x1[0]->is_empty())+
                       2*int(this->m_neighbors_2x2[0]->is_empty())+
                       2*int(this->m_neighbors_2x2[1]->is_empty()));
        case 2: return(1*int(this->m_neighbors_1x1[1]->is_empty())+
                       2*int(this->m_neighbors_2x2[2]->is_empty())+
                       2*int(this->m_neighbors_2x2[3]->is_empty()));
        case 3: return(1*int(this->m_neighbors_1x1[2]->is_empty())+
                       2*int(this->m_neighbors_2x2[4]->is_empty())+
                       2*int(this->m_neighbors_2x2[5]->is_empty()));
        case 4: return(1*int(this->m_neighbors_1x1[3]->is_empty())+
                       2*int(this->m_neighbors_2x2[6]->is_empty())+
                       2*int(this->m_neighbors_2x2[7]->is_empty()));
        }
}

Site_2x2::Site_2x2(int id, int x, int y) : Site(id, x, y), m_neighbors_1x1(8), m_neighbors_2x2(12)
{
        // ...
}

bool Site_2x2::step_is_valid(int dir){
// = = = = = = = = = = = = = = = = = = =
// blocking sites coordinates
// = = = = = = = = = = = = = = = = = = =
// m_neighbors_1x1[0]  @  X+1  Y
// m_neighbors_1x1[1]  @  X+1  Y+1
// m_neighbors_1x1[2]  @  X    Y+2
// m_neighbors_1x1[3]  @  X-1  Y+2
// m_neighbors_1x1[4]  @  X-2  Y+1
// m_neighbors_1x1[5]  @  X-2  Y
// m_neighbors_1x1[6]  @  X-1  Y-1
// m_neighbors_1x1[7]  @  X    Y-1
// = = = = = = = = = = = = = = = = = = =
// m_neighbors_2x2[0]  @  X+2  Y-1
// m_neighbors_2x2[1]  @  X+2  Y
// m_neighbors_2x2[2]  @  X+2  Y+1
// m_neighbors_2x2[3]  @  X+1  Y+2
// m_neighbors_2x2[4]  @  X    Y+2
// m_neighbors_2x2[5]  @  X-1  Y+2
// m_neighbors_2x2[6]  @  X-2  Y+1
// m_neighbors_2x2[7]  @  X-2  Y
// m_neighbors_2x2[8]  @  X-2  Y-1
// m_neighbors_2x2[9]  @  X-1  Y-2
// m_neighbors_2x2[10] @  X    Y-2
// m_neighbors_2x2[11] @  X+1  Y-2
// = = = = = = = = = = = = = = = = = = =
        switch(dir) {
        case 1: return (this->m_neighbors_1x1[0]->is_empty() &&
                        this->m_neighbors_1x1[1]->is_empty() &&
                        this->m_neighbors_2x2[0]->is_empty() &&
                        this->m_neighbors_2x2[1]->is_empty() &&
                        this->m_neighbors_2x2[2]->is_empty());
        case 2: return (this->m_neighbors_1x1[2]->is_empty() &&
                        this->m_neighbors_1x1[3]->is_empty() &&
                        this->m_neighbors_2x2[3]->is_empty() &&
                        this->m_neighbors_2x2[4]->is_empty() &&
                        this->m_neighbors_2x2[5]->is_empty());
        case 3: return (this->m_neighbors_1x1[4]->is_empty() &&
                        this->m_neighbors_1x1[5]->is_empty() &&
                        this->m_neighbors_2x2[6]->is_empty() &&
                        this->m_neighbors_2x2[7]->is_empty() &&
                        this->m_neighbors_2x2[8]->is_empty());
        case 4: return (this->m_neighbors_1x1[6]->is_empty() &&
                        this->m_neighbors_1x1[7]->is_empty() &&
                        this->m_neighbors_2x2[9]->is_empty() &&
                        this->m_neighbors_2x2[10]->is_empty() &&
                        this->m_neighbors_2x2[11]->is_empty());
        default: return false;
        }
}

int Site_2x2::step_is_invalid(dir){
        switch(dir)
        {

        case 1: return(1*int(this->m_neighbors_1x1[0]->is_blocked() || this->m_neighbors_1x1[1]->is_blocked())+
                       3*int(this->m_neighbors_1x1[0]->is_blocked() && this->m_neighbors_1x1[1]->is_blocked())+
                       3*int(this->m_neighbors_2x2[0]->is_blocked())+
                       2*int(this->m_neighbors_2x2[1]->is_blocked())+
                       3*int(this->m_neighbors_2x2[2]->is_blocked())+
                       1*int(this->m_neighbors_1x1[0]->is_blocked() && this->m_neighbors_2x2[2]->is_blocked())+
                       1*int(this->m_neighbors_1x1[1]->is_blocked() && this->m_neighbors_2x2[0]->is_blocked()));

        case 2: return(1*int(this->m_neighbors_1x1[2]->is_blocked() || this->m_neighbors_1x1[3]->is_blocked())+
                       3*int(this->m_neighbors_1x1[2]->is_blocked() && this->m_neighbors_1x1[3]->is_blocked())+
                       3*int(this->m_neighbors_2x2[3]->is_blocked())+
                       2*int(this->m_neighbors_2x2[4]->is_blocked())+
                       3*int(this->m_neighbors_2x2[5]->is_blocked())+
                       1*int(this->m_neighbors_1x1[2]->is_blocked() && this->m_neighbors_2x2[5]->is_blocked())+
                       1*int(this->m_neighbors_1x1[3]->is_blocked() && this->m_neighbors_2x2[3]->is_blocked()));

        case 3: return(1*int(this->m_neighbors_1x1[4]->is_blocked() || this->m_neighbors_1x1[5]->is_blocked())+
                       3*int(this->m_neighbors_1x1[4]->is_blocked() && this->m_neighbors_1x1[5]->is_blocked())+
                       3*int(this->m_neighbors_2x2[6]->is_blocked())+
                       2*int(this->m_neighbors_2x2[7]->is_blocked())+
                       3*int(this->m_neighbors_2x2[8]->is_blocked())+
                       1*int(this->m_neighbors_1x1[4]->is_blocked() && this->m_neighbors_2x2[8]->is_blocked())+
                       1*int(this->m_neighbors_1x1[5]->is_blocked() && this->m_neighbors_2x2[6]->is_blocked()));

        case 4: return(1*int(this->m_neighbors_1x1[6]->is_blocked() || this->m_neighbors_1x1[7]->is_blocked())+
                       3*int(this->m_neighbors_1x1[6]->is_blocked() && this->m_neighbors_1x1[7]->is_blocked())+
                       3*int(this->m_neighbors_2x2[9]->is_blocked())+
                       2*int(this->m_neighbors_2x2[10]->is_blocked())+
                       3*int(this->m_neighbors_2x2[11]->is_blocked())+
                       1*int(this->m_neighbors_1x1[6]->is_blocked() && this->m_neighbors_2x2[11]->is_blocked())+
                       1*int(this->m_neighbors_1x1[7]->is_blocked() && this->m_neighbors_2x2[9]->is_blocked()));

        }
}
