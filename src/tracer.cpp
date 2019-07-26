#include "tracer.hpp"
#include "global.hpp"

// constructor
Tracer::Tracer(int id, int x, int y, int grid_size_x, int grid_size_y) : m_id(id), m_x(x), m_y(y), m_dx(0), m_dy(0),  m_grid_size_x(grid_size_x), m_grid_size_y(grid_size_y), m_lsquared(0.0), m_isstuck(false), m_size(1){

}

void Tracer::unhindered_step(){
        // function for stepping without collisions/tracer-tracer interaction
        int dir = random_int(1,4);
        switch(dir) {
        case 1:
        {
                int new_x = (this->m_x + 1)%this->m_grid_size_x;
                this->m_x = new_x;
                this->m_dx++;
                break;
        }
        case 2:
        {
                int new_y = (this->m_y + 1)%this->m_grid_size_y;
                this->m_y = new_y;
                this->m_dy++;
                break;
        }
        case 3:
        {
                int new_x = (this->m_x - 1 + this->m_grid_size_x)%this->m_grid_size_x;
                this->m_x = new_x;
                this->m_dx--;
                break;
        }
        case 4:
        {
                int new_y = (this->m_y - 1 + this->m_grid_size_y)%this->m_grid_size_y;
                this->m_y = new_y;
                this->m_dy--;
                break;
        }
        }
        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
}
// function for random walk stepping
void Tracer::step(std::vector<int> &grid_occupation_vector){
        if(this->m_isstuck) {
                return;
        }
        int dir = random_int(1,4);
        switch(dir) {
        case 1:
        {
                int new_x = (this->m_x + 1)%this->m_grid_size_x;
                int & new_site = grid_occupation_vector[new_x * this->m_grid_size_y + this->m_y];
                if(!new_site) {
                        int & old_site = grid_occupation_vector[this->m_x * this->m_grid_size_y + this->m_y];
                        old_site = 0;
                        this->m_x = new_x;
                        this->m_dx++;
                        new_site = this->m_id;
                }
                break;
        }
        case 2:
        {
                int new_y = (this->m_y + 1)%this->m_grid_size_y;
                int & new_site = grid_occupation_vector[this->m_x * this->m_grid_size_y + new_y];
                if(!new_site) {
                        int & old_site = grid_occupation_vector[this->m_x * this->m_grid_size_y + this->m_y];
                        old_site = 0;
                        this->m_y = new_y;
                        this->m_dy++;
                        new_site = this->m_id;
                }
                break;
        }
        case 3:
        {
                int new_x = (this->m_x - 1 + this->m_grid_size_x)%this->m_grid_size_x;
                int & new_site = grid_occupation_vector[new_x * this->m_grid_size_y + this->m_y];
                if(!new_site) {
                        int & old_site = grid_occupation_vector[this->m_x * this->m_grid_size_y + this->m_y];
                        old_site = 0;
                        this->m_x = new_x;
                        this->m_dx--;
                        new_site = this->m_id;
                }
                break;
        }
        case 4:
        {
                int new_y = (this->m_y - 1 + this->m_grid_size_y)%this->m_grid_size_y;
                int & new_site = grid_occupation_vector[this->m_x * this->m_grid_size_y + new_y];
                if(!new_site) {
                        int & old_site = grid_occupation_vector[this->m_x * this->m_grid_size_y + this->m_y];
                        old_site = 0;
                        this->m_y = new_y;
                        this->m_dy--;
                        new_site = this->m_id;
                }
                break;
        }
        }
        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
}

// getters
int Tracer::get_id()
{
        return this->m_id;
}
int Tracer::get_x()
{
        return this->m_x;
}
int Tracer::get_dx()
{
        return this->m_dx;
}
int Tracer::get_y()
{
        return this->m_y;
}
int Tracer::get_dy()
{
        return this->m_dy;
}
int Tracer::get_size()
{
        return this->m_size;
}
double Tracer::get_lsquared()
{
        return this->m_lsquared;
}
bool Tracer::get_isstuck()
{
        return this->m_isstuck;
}
// change mobility state
void Tracer::stuck()
{
        this->m_isstuck = true;
}
void Tracer::unstuck()
{
        this->m_isstuck = false;
}
//
Tracer_2x2::Tracer_2x2(int id, int x, int y, int grid_size_x, int grid_size_y) : Tracer(id,x,y,grid_size_x,grid_size_y)
{
        // Change size from default to 4
        this->m_size = 4;
}

void Tracer_2x2::step(std::vector<int> &grid_occupation_vector){
        if(this->m_isstuck) {
                return;
        }
        // old coordinates of the tracers cells
        int old_x_1 = this->m_x;
        int old_x_2 = (this->m_x + 1)%this->m_grid_size_x;
        int old_y_1 = this->m_y;
        int old_y_2 = (this->m_y + 1)%this->m_grid_size_y;
        // chose randomly one of the 4 directions
        int dir = random_int(1,4);
        switch(dir) {
        case 1:
        {
                int new_x = (this->m_x + 2)%this->m_grid_size_x;
                // get new positions occupation state by reference
                int &new_site_1 = grid_occupation_vector[new_x * this->m_grid_size_y + old_y_1];
                int &new_site_2 = grid_occupation_vector[new_x * this->m_grid_size_y + old_y_2];
                if(!(new_site_1 || new_site_2)) {
                        // if new sites are free, get old sites by reference
                        int &old_site_1 = grid_occupation_vector[old_x_1 * this->m_grid_size_y + old_y_1];
                        int &old_site_2 = grid_occupation_vector[old_x_1 * this->m_grid_size_y + old_y_2];
                        // set them to free (=0)
                        old_site_1 = 0;
                        old_site_2 = 0;
                        // move the tracer to new position
                        this->m_x = old_x_2;
                        this->m_dx++;
                        // set new sites to occupied (!=0)
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                }
                break;
        }
        case 2:
        {
                int new_y = (this->m_y + 2)%this->m_grid_size_y;
                int &new_site_1 = grid_occupation_vector[old_x_1 * this->m_grid_size_y + new_y];
                int &new_site_2 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + new_y];
                if(!(new_site_1 || new_site_2)) {
                        int &old_site_1 = grid_occupation_vector[old_x_1 * this->m_grid_size_y + old_y_1];
                        int &old_site_2 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + old_y_1];
                        old_site_1 = 0;
                        old_site_2 = 0;
                        this->m_y = old_y_2;
                        this->m_dy++;
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                }
                break;
        }
        case 3:
        {
                int new_x = (this->m_x - 1 + this->m_grid_size_x)%this->m_grid_size_x;
                int &new_site_1 = grid_occupation_vector[new_x * this->m_grid_size_y + old_y_1];
                int &new_site_2 = grid_occupation_vector[new_x * this->m_grid_size_y + old_y_2];
                if(!(new_site_1 || new_site_2)) {
                        int &old_site_1 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + old_y_1];
                        int &old_site_2 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + old_y_2];
                        old_site_1 = 0;
                        old_site_2 = 0;
                        this->m_x = new_x;
                        this->m_dx--;
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                }
                break;
        }
        case 4:
        {
                int new_y = (this->m_y - 1 + this->m_grid_size_y)%this->m_grid_size_y;
                int &new_site_1 = grid_occupation_vector[old_x_1 * this->m_grid_size_y + new_y];
                int &new_site_2 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + new_y];
                if(!(new_site_1 || new_site_2)) {
                        int &old_site_1 = grid_occupation_vector[old_x_1 * this->m_grid_size_y + old_y_2];
                        int &old_site_2 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + old_y_2];
                        old_site_1 = 0;
                        old_site_2 = 0;
                        this->m_y = new_y;
                        this->m_dy--;
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                }
                break;
        }
        }
        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
}
void Tracer_2x2::unhindered_step(){
        // function for stepping without collisions/tracer-tracer interaction
        int dir = random_int(1,4);
        switch(dir) {
        case 1:
        {
                int new_x = (this->m_x + 1)%this->m_grid_size_x;
                this->m_x = new_x;
                this->m_dx++;
                break;
        }
        case 2:
        {
                int new_y = (this->m_y + 1)%this->m_grid_size_y;
                this->m_y = new_y;
                this->m_dy++;
                break;
        }
        case 3:
        {
                int new_x = (this->m_x - 1 + this->m_grid_size_x)%this->m_grid_size_x;
                this->m_x = new_x;
                this->m_dx--;
                break;
        }
        case 4:
        {
                int new_y = (this->m_y - 1 + this->m_grid_size_y)%this->m_grid_size_y;
                this->m_y = new_y;
                this->m_dy--;
                break;
        }
        }
        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
}

Tracer_3x3::Tracer_3x3(int id, int x, int y, int grid_size_x, int grid_size_y) : Tracer(id,x,y,grid_size_x,grid_size_y)
{
        // Change size from default to 9
        this->m_size = 9;
}

void Tracer_3x3::step(std::vector<int> &grid_occupation_vector){
        if(this->m_isstuck) {
                return;
        }
        // old coordinates of the tracers cells
        int old_x_1 = this->m_x;
        int old_x_2 = (this->m_x + 1)%this->m_grid_size_x;
        int old_x_3 = (this->m_x - 1 + this->m_grid_size_x)%this->m_grid_size_x;
        int old_y_1 = this->m_y;
        int old_y_2 = (this->m_y + 1)%this->m_grid_size_y;
        int old_y_3 = (this->m_y - 1 + this->m_grid_size_y)%this->m_grid_size_y;
        // chose randomly one of the 4 directions
        int dir = random_int(1,4);
        switch(dir) {
        case 1:
        {
                // case 1 -> move +1 in x-direction
                int new_x = (this->m_x + 2)%this->m_grid_size_x;
                // get new positions occupation state by reference
                int &new_site_1 = grid_occupation_vector[new_x * this->m_grid_size_y + old_y_1];
                int &new_site_2 = grid_occupation_vector[new_x * this->m_grid_size_y + old_y_2];
                int &new_site_3 = grid_occupation_vector[new_x * this->m_grid_size_y + old_y_3];
                if(!(new_site_1 || new_site_2 || new_site_3)) {
                        // if new sites are free, get old sites by reference
                        int &old_site_1 = grid_occupation_vector[old_x_3 * this->m_grid_size_y + old_y_1];
                        int &old_site_2 = grid_occupation_vector[old_x_3 * this->m_grid_size_y + old_y_2];
                        int &old_site_3 = grid_occupation_vector[old_x_3 * this->m_grid_size_y + old_y_3];
                        // set them to free (=0)
                        old_site_1 = 0;
                        old_site_2 = 0;
                        old_site_3 = 0;
                        // move the tracer to new position
                        this->m_x = old_x_2;
                        this->m_dx++;
                        // set new sites to occupied (!=0)
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                        new_site_3 = this->m_id;
                }
                break;
        }
        case 2:
        {
                // -> +1 in y direction
                int new_y = (this->m_y + 2)%this->m_grid_size_y;
                int &new_site_1 = grid_occupation_vector[old_x_1 * this->m_grid_size_y + new_y];
                int &new_site_2 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + new_y];
                int &new_site_3 = grid_occupation_vector[old_x_3 * this->m_grid_size_y + new_y];
                if(!(new_site_1 || new_site_2 || new_site_3)) {
                        int &old_site_1 = grid_occupation_vector[old_x_1 * this->m_grid_size_y + old_y_3];
                        int &old_site_2 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + old_y_3];
                        int &old_site_3 = grid_occupation_vector[old_x_3 * this->m_grid_size_y + old_y_3];
                        old_site_1 = 0;
                        old_site_2 = 0;
                        old_site_3 = 0;
                        this->m_y = old_y_2;
                        this->m_dy++;
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                        new_site_3 = this->m_id;
                }
                break;
        }
        case 3:
        {
                // -> -1 in x direction
                int new_x = (this->m_x - 2 + this->m_grid_size_x)%this->m_grid_size_x;
                int &new_site_1 = grid_occupation_vector[new_x * this->m_grid_size_y + old_y_1];
                int &new_site_2 = grid_occupation_vector[new_x * this->m_grid_size_y + old_y_2];
                int &new_site_3 = grid_occupation_vector[new_x * this->m_grid_size_y + old_y_3];
                if(!(new_site_1 || new_site_2 || new_site_3)) {
                        int &old_site_1 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + old_y_1];
                        int &old_site_2 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + old_y_2];
                        int &old_site_3 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + old_y_3];
                        old_site_1 = 0;
                        old_site_2 = 0;
                        old_site_3 = 0;
                        this->m_x = new_x;
                        this->m_dx--;
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                        new_site_3 = this->m_id;
                }
                break;
        }
        case 4:
        {
                // -> -1 in y direction
                int new_y = (this->m_y - 2 + this->m_grid_size_y)%this->m_grid_size_y;
                int &new_site_1 = grid_occupation_vector[old_x_1 * this->m_grid_size_y + new_y];
                int &new_site_2 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + new_y];
                int &new_site_3 = grid_occupation_vector[old_x_3 * this->m_grid_size_y + new_y];
                if(!(new_site_1 || new_site_2 || new_site_3)) {
                        int &old_site_1 = grid_occupation_vector[old_x_1 * this->m_grid_size_y + old_y_2];
                        int &old_site_2 = grid_occupation_vector[old_x_2 * this->m_grid_size_y + old_y_2];
                        int &old_site_3 = grid_occupation_vector[old_x_3 * this->m_grid_size_y + old_y_2];
                        old_site_1 = 0;
                        old_site_2 = 0;
                        old_site_3 = 0;
                        this->m_y = new_y;
                        this->m_dy--;
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                        new_site_3 = this->m_id;
                }
                break;
        }
        }
        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
}


void Tracer_3x3::unhindered_step(){
        // function for stepping without collisions/tracer-tracer interaction
        int dir = random_int(1,4);
        switch(dir) {
        case 1:
        {
                int new_x = (this->m_x + 1)%this->m_grid_size_x;
                this->m_x = new_x;
                this->m_dx++;
                break;
        }
        case 2:
        {
                int new_y = (this->m_y + 1)%this->m_grid_size_y;
                this->m_y = new_y;
                this->m_dy++;
                break;
        }
        case 3:
        {
                int new_x = (this->m_x - 1 + this->m_grid_size_x)%this->m_grid_size_x;
                this->m_x = new_x;
                this->m_dx--;
                break;
        }
        case 4:
        {
                int new_y = (this->m_y - 1 + this->m_grid_size_y)%this->m_grid_size_y;
                this->m_y = new_y;
                this->m_dy--;
                break;
        }
        }
        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
}
