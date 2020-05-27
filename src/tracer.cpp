#include "tracer.hpp"
#include "global.hpp"


// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif


// - - - - - - - - - - - - - - - - - - - - - - - - -
// C O N S T R U C T O R
// - - - - - - - - - - - - - - - - - - - - - - - - -
Tracer::Tracer(int id, int x, int y, int grid_size_x, int grid_size_y, double step_rate, int wtd_max, int wtd_res, int step_attempts_per_timestep) :
        m_id(id),
        m_size(1),
        m_x(x),
        m_y(y),
        m_dx(0),
        m_dy(0),
        m_step_rate(step_rate),
        m_last_step(0),
        m_last_step_dir(4,0),
        m_last_step_idx(4,0),
        m_last_step_wtd_idx(0),
        m_wtd_max(wtd_max),
        m_wtd_res(wtd_res),
        m_wtd_max_index(16*wtd_max*wtd_res),
        m_time_of_last_step(0),
        m_time_since_last_step(0),
        m_step_attempts_per_timestep(step_attempts_per_timestep),
        m_steps_taken(0),
        m_grid_size_x(grid_size_x),
        m_grid_size_y(grid_size_y),
        m_lsquared(0),
        m_isstuck(false)
{

}
// - - - - - - - - - - - - - - - - - - - - - - - - -
inline void Tracer::update_last_step(int last_step_dir)
{
        // shift m_last_step_dir by one position
        std::rotate(this->m_last_step_dir.rbegin(),this->m_last_step_dir.rbegin()+1,this->m_last_step_dir.rend());
        // overwrite the move that's now at pos 0 with the most recent one
        this->m_last_step_dir[0] = last_step_dir;
        int tmp_idx = 0;
        int tmp_n = 0;
        for(int step : this->m_last_step_dir)
        {
                if(!step) { break; }
                tmp_idx += (int)std::pow(4,tmp_n)*step;
                this->m_last_step_idx[tmp_n] = tmp_idx;
                tmp_n++;
        }
        if(!this->m_last_step_idx[1]) { return; }
        this->m_last_step_wtd_idx = std::min((int)(16*(((this->m_time_since_last_step-1)*this->m_wtd_res)/this->m_step_attempts_per_timestep)),this->m_wtd_max_index)+this->m_last_step_idx[1]-4;
        // = this->m_last_step_idx[1];
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
// coordinate conversion x,y -> vector idx
// - - - - - - - - - - - - - - - - - - - - - - - - -
inline int Tracer::coord(int x, int y){
        return x*this->m_grid_size_y+y;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
// function for random walk stepping
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer::step(std::vector<int> &grid_occupation_vector, int dir, int current_time, int current_step_attempt_number){
        //
        this->m_last_step = false;
        if(this->m_isstuck) { return; }
        switch(dir) {
        case 1:
        {
                int new_x = (this->m_x + 1)%this->m_grid_size_x;
                int & new_site = grid_occupation_vector[this->coord(new_x,this->m_y)];
                if(!new_site) {
                        int & old_site = grid_occupation_vector[this->coord(this->m_x,this->m_y)];
                        old_site = 0;
                        this->m_x = new_x;
                        this->m_dx++;
                        new_site = this->m_id;

                        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                        this->m_time_of_last_step = current_time;
                        this->m_last_step_attempt_number = current_step_attempt_number;

                        this->m_last_step = true;
                        this->m_steps_taken++;

                        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                        this->update_last_step(dir);
                }
                return;
        }
        case 2:
        {
                int new_y = (this->m_y + 1)%this->m_grid_size_y;
                int & new_site = grid_occupation_vector[this->coord(this->m_x,new_y)];
                if(!new_site) {
                        int & old_site = grid_occupation_vector[this->coord(this->m_x,this->m_y)];
                        old_site = 0;
                        this->m_y = new_y;
                        this->m_dy++;
                        new_site = this->m_id;

                        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                        this->m_time_of_last_step = current_time;
                        this->m_last_step_attempt_number = current_step_attempt_number;

                        this->m_last_step = true;
                        this->m_steps_taken++;

                        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                        this->update_last_step(dir);
                }
                return;
        }
        case 3:
        {
                int new_x = (this->m_x - 1 + this->m_grid_size_x)%this->m_grid_size_x;
                int & new_site = grid_occupation_vector[this->coord(new_x,this->m_y)];
                if(!new_site) {
                        int & old_site = grid_occupation_vector[this->coord(this->m_x,this->m_y)];
                        old_site = 0;
                        this->m_x = new_x;
                        this->m_dx--;
                        new_site = this->m_id;

                        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                        this->m_time_of_last_step = current_time;
                        this->m_last_step_attempt_number = current_step_attempt_number;

                        this->m_last_step = true;
                        this->m_steps_taken++;

                        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                        this->update_last_step(dir);
                }
                return;
        }
        case 4:
        {
                int new_y = (this->m_y - 1 + this->m_grid_size_y)%this->m_grid_size_y;
                int & new_site = grid_occupation_vector[this->coord(this->m_x,new_y)];
                if(!new_site) {
                        int & old_site = grid_occupation_vector[this->coord(this->m_x,this->m_y)];
                        old_site = 0;
                        this->m_y = new_y;
                        this->m_dy--;
                        new_site = this->m_id;
                        // for waiting time distribution
                        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                        this->m_time_of_last_step = current_time;
                        this->m_last_step_attempt_number = current_step_attempt_number;

                        this->m_last_step = true;
                        this->m_steps_taken++;

                        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                        this->update_last_step(dir);
                }
                return;
        }
        default: { return; }
        }
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
// Function for steps during warm_up
// contains only logic for site blocking and trapping
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer::step_warmup(std::vector<int> &grid_occupation_vector, int dir){
        // if the tracer is trapped, do nothing
        if(this->m_isstuck) { return; }
        switch(dir) {
        case 1:
        {
                int new_x = (this->m_x + 1)%this->m_grid_size_x;
                int & new_site = grid_occupation_vector[this->coord(new_x,this->m_y)];
                if(!new_site) {
                        int & old_site = grid_occupation_vector[this->coord(this->m_x,this->m_y)];
                        old_site = 0;
                        this->m_x = new_x;
                        new_site = this->m_id;
                }
                return;
        }
        case 2:
        {
                int new_y = (this->m_y + 1)%this->m_grid_size_y;
                int & new_site = grid_occupation_vector[this->coord(this->m_x,new_y)];
                if(!new_site) {
                        int & old_site = grid_occupation_vector[this->coord(this->m_x,this->m_y)];
                        old_site = 0;
                        this->m_y = new_y;
                        new_site = this->m_id;
                }
                return;
        }
        case 3:
        {
                int new_x = (this->m_x - 1 + this->m_grid_size_x)%this->m_grid_size_x;
                int & new_site = grid_occupation_vector[this->coord(new_x,this->m_y)];
                if(!new_site) {
                        int & old_site = grid_occupation_vector[this->coord(this->m_x,this->m_y)];
                        old_site = 0;
                        this->m_x = new_x;
                        new_site = this->m_id;
                }
                return;
        }
        case 4:
        {
                int new_y = (this->m_y - 1 + this->m_grid_size_y)%this->m_grid_size_y;
                int & new_site = grid_occupation_vector[this->coord(this->m_x,new_y)];
                if(!new_site) {
                        int & old_site = grid_occupation_vector[this->coord(this->m_x,this->m_y)];
                        old_site = 0;
                        this->m_y = new_y;
                        new_site = this->m_id;
                }
                return;
        }
        default: { return; }
        }
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer::step_unhindered(int dir, int current_time, int current_step_attempt_number){
        // function for stepping without collisions/tracer-tracer interaction
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
        default: { return; }
        }
        // since all step attempts are successful, I can do this outside of the switch-case
        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
        this->m_time_of_last_step = current_time;
        this->m_last_step_attempt_number = current_step_attempt_number;
        // last but not least, update lsq, last moves
        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
        this->update_last_step(dir);

}
// - - - - - - - - - - - - - - - - - - - - - - - - -
// G E T T E R S
// - - - - - - - - - - - - - - - - - - - - - - - - -
int Tracer::get_id()
{
        return this->m_id;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
int Tracer::get_x()
{
        return this->m_x;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
int Tracer::get_dx()
{
        return this->m_dx;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
int Tracer::get_y()
{
        return this->m_y;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
int Tracer::get_dy()
{
        return this->m_dy;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
int Tracer::get_size()
{
        return this->m_size;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
double Tracer::get_lsquared()
{
        return this->m_lsquared;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
bool Tracer::get_isstuck()
{
        return this->m_isstuck;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
bool Tracer::get_last_step()
{
        return this->m_last_step;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
int Tracer::get_wtd_idx()
{
        return this->m_last_step_wtd_idx;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector< int> Tracer::get_last_step_idx(){
        return this->m_last_step_idx;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
int Tracer::get_steps_taken()
{
        return this->m_steps_taken;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
// change mobility state (e.g. by trapping mechanism)
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer::stuck()
{
        this->m_isstuck = true;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer::unstuck()
{
        this->m_isstuck = false;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
Tracer_2x2::Tracer_2x2(int id, int x, int y, int grid_size_x, int grid_size_y, double step_rate, int wtd_max, int wtd_res, int step_attempts_per_timestep) : Tracer(id,x,y,grid_size_x,grid_size_y,step_rate,wtd_max,wtd_res,step_attempts_per_timestep)
{
        // Change size from default to 4
        this->m_size = 4;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer_2x2::step(std::vector<int> &grid_occupation_vector, int dir, int current_time, int current_step_attempt_number){
        D(std::cout << "STEP 2X2" << std::endl);
        this->m_last_step = false;
        if(this->m_isstuck) {
                return;
        }
        // old coordinates of the tracers cells
        int old_x_1 = this->m_x;
        int old_x_2 = (this->m_x + 1)%this->m_grid_size_x;
        int old_y_1 = this->m_y;
        int old_y_2 = (this->m_y + 1)%this->m_grid_size_y;
        // chose randomly one of the 4 directions
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

                        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                        this->m_time_of_last_step = current_time;
                        this->m_last_step_attempt_number = current_step_attempt_number;

                        this->m_last_step = true;
                        this->m_steps_taken++;

                        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                        this->update_last_step(dir);
                }
                return;
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

                        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                        this->m_time_of_last_step = current_time;
                        this->m_last_step_attempt_number = current_step_attempt_number;

                        this->m_last_step = true;
                        this->m_steps_taken++;

                        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                        this->update_last_step(dir);
                }
                return;
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

                        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                        this->m_time_of_last_step = current_time;
                        this->m_last_step_attempt_number = current_step_attempt_number;

                        this->m_last_step = true;
                        this->m_steps_taken++;

                        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                        this->update_last_step(dir);
                }
                return;
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

                        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                        this->m_time_of_last_step = current_time;
                        this->m_last_step_attempt_number = current_step_attempt_number;

                        this->m_last_step = true;
                        this->m_steps_taken++;

                        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                        this->update_last_step(dir);
                }
                return;
        }
        default: { return; }
        }
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer_2x2::step_warmup(std::vector<int> &grid_occupation_vector, int dir){
        this->m_last_step = false;
        if(this->m_isstuck) {
                return;
        }
        // old coordinates of the tracers cells
        int old_x_1 = this->m_x;
        int old_x_2 = (this->m_x + 1)%this->m_grid_size_x;
        int old_y_1 = this->m_y;
        int old_y_2 = (this->m_y + 1)%this->m_grid_size_y;
        // chose randomly one of the 4 directions
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
                        // set new sites to occupied (!=0)
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                }
                return;
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
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                }
                return;
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
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                }
                return;
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
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                }
                return;
        }
        default: { return; }
        }
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer_2x2::step_unhindered(int dir, int current_time, int current_step_attempt_number){
        // function for stepping without collisions/tracer-tracer interaction
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
        default: { return; }
        }
        // since all step attempts are successful, I can do this outside of the switch-case
        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
        this->m_time_of_last_step = current_time;
        this->m_last_step_attempt_number = current_step_attempt_number;
        // last but not least, update lsq, last moves
        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
        this->update_last_step(dir);
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
Tracer_3x3::Tracer_3x3(int id, int x, int y, int grid_size_x, int grid_size_y, double step_rate, int wtd_max, int wtd_res, int step_attempts_per_timestep) : Tracer(id,x,y,grid_size_x,grid_size_y,step_rate,wtd_max,wtd_res,step_attempts_per_timestep)
{
        // Change size from default to 9
        this->m_size = 9;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer_3x3::step(std::vector<int> &grid_occupation_vector, int dir, int current_time, int current_step_attempt_number){
        this->m_last_step = false;
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

                        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                        this->m_time_of_last_step = current_time;
                        this->m_last_step_attempt_number = current_step_attempt_number;

                        this->m_last_step = true;
                        this->m_steps_taken++;

                        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                        this->update_last_step(dir);
                }
                return;
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

                        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                        this->m_time_of_last_step = current_time;
                        this->m_last_step_attempt_number = current_step_attempt_number;

                        this->m_last_step = true;
                        this->m_steps_taken++;

                        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                        this->update_last_step(dir);
                }
                return;
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

                        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                        this->m_time_of_last_step = current_time;
                        this->m_last_step_attempt_number = current_step_attempt_number;

                        this->m_last_step = true;
                        this->m_steps_taken++;

                        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                        this->update_last_step(dir);
                }
                return;
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

                        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                        this->m_time_of_last_step = current_time;
                        this->m_last_step_attempt_number = current_step_attempt_number;

                        this->m_last_step = true;
                        this->m_steps_taken++;

                        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                        this->update_last_step(dir);
                }
                return;
        }
        default: { return; }
        }
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer_3x3::step_warmup(std::vector<int> &grid_occupation_vector, int dir){
        this->m_last_step = false;
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
                        // set new sites to occupied (!=0)
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                        new_site_3 = this->m_id;
                }
                return;
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
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                        new_site_3 = this->m_id;
                }
                return;
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
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                        new_site_3 = this->m_id;
                }
                return;
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
                        new_site_1 = this->m_id;
                        new_site_2 = this->m_id;
                        new_site_3 = this->m_id;
                }
                return;
        }
        default: { return; }
        }
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer_3x3::step_unhindered(int dir, int current_time, int current_step_attempt_number){
        // function for stepping without collisions/tracer-tracer interaction
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
        default: { return; }
        }

        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
        this->m_time_of_last_step = current_time;
        this->m_last_step_attempt_number = current_step_attempt_number;

        this->m_last_step = true;
        this->m_steps_taken++;
        this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
        this->update_last_step(dir);
        return;
}
