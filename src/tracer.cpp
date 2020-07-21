#include "tracer.hpp"
#include "global.hpp"


// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif

//
// map<int, int> m_next_sites;
// map<int, std::vector<int>> m_neighbor_coords;
//

// - - - - - - - - - - - - - - - - - - - - - - - - -
// C O N S T R U C T O R
// - - - - - - - - - - - - - - - - - - - - - - - - -
Tracer::Tracer(int id, Site * starting_site, Lattice * parent_lattice, int wtd_max, int wtd_res, int step_attempts_per_timestep) :
        m_id(id),
        m_parent_lattice(parent_lattice),
        m_site(starting_site),
        m_size(1),
        m_dx(0),
        m_dy(0),
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
        m_step_attempt_result(0),
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
// function for random walk stepping
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer::step(int dir){
        //
        this->m_last_step = false;
        if(this->m_isstuck) { return; }
        //
        this->m_step_attempt_result = this->m_site->step_is_invalid(dir);
        if(!this->m_step_attempt_result) {
                // make the old site empty
                this->m_site->set_empty();
                // occupy the new site
                this->m_site = this->m_site->get_neighbor_by_dir(dir);
                this->m_site->set_occupied();
                //
                switch(dir)
                {
                case 1: this->m_dx++;
                case 2: this->m_dy++;
                case 3: this->m_dx--;
                case 4: this->m_dx--;
                }
                //
                int current_time = this->m_parent_lattice->get_t();
                int current_step_attempt = this->m_parent_lattice->get_current_step_attempt();
                //
                this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
                this->m_time_of_last_step = current_time;
                this->m_last_step_attempt_number = current_step_attempt_number;
                //
                this->m_last_step = true;
                this->m_steps_taken++;
                //
                this->m_lsquared = pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
                this->update_last_step(dir);
        }
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
// Function for steps during warm_up
// contains only logic for site blocking and trapping
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer::step_warmup(int dir){
        //
        if(this->m_isstuck) { return; }
        //
        if(this->m_site->step_is_valid(dir)) {
                // make the old site empty
                this->m_site->set_empty();
                // occupy the new site
                this->m_site = this->m_site->get_neighbor_by_dir(dir);
                this->m_site->set_occupied();
        }
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer::step_unhindered(int dir){
        //
        this->m_last_step = false;
        // make the old site empty
        this->m_site->set_empty();
        // occupy the new site
        this->m_site = this->m_site->get_neighbor_by_dir(dir);
        this->m_site->set_occupied();
        //
        switch(dir)
        {
        case 1: this->m_dx++;
        case 2: this->m_dy++;
        case 3: this->m_dx--;
        case 4: this->m_dx--;
        }
        //
        int current_time = this->m_parent_lattice->get_t();
        int current_step_attempt = this->m_parent_lattice->get_current_step_attempt();
        //
        this->m_time_since_last_step = (current_time-this->m_time_of_last_step < this->m_wtd_max ? (current_time-this->m_time_of_last_step)*this->m_step_attempts_per_timestep + current_step_attempt_number - this->m_last_step_attempt_number : this->m_wtd_max);
        this->m_time_of_last_step = current_time;
        this->m_last_step_attempt_number = current_step_attempt_number;
        //
        this->m_last_step = true;
        this->m_steps_taken++;
        //
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
int Tracer::get_step_attempt_result()
{
        return this->m_step_attempt_result;
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
Tracer_2x2::Tracer_2x2(int id, Site * starting_site, Lattice * parent_lattice, int wtd_max, int wtd_res, int step_attempts_per_timestep) : Tracer(id,starting_site,parent_lattice,wtd_max,wtd_res,step_attempts_per_timestep)
{
        // Change size from default to 4
        this->m_size = 4;
}
