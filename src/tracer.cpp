#include "tracer.hpp"

// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif
// To add debugging messages, use D(std::cerr << "Debugging message 1 2 3!" << std::endl; )

// - - - - - - - - - - - - - - - - - - - - - - - - -
// C O N S T R U C T O R
// - - - - - - - - - - - - - - - - - - - - - - - - -
Tracer::Tracer(int id, Site * starting_site) :
        m_id(id),
        m_size(1),
        m_site(starting_site),
        m_dx(0),
        m_dy(0),
        m_lsquared(0),
        m_last_step(false),
        m_last_step_dir(3,0),
        m_correlations(84,0),
        //m_steps_taken(0),
        m_steps_taken(0),
        m_isstuck(false)
{
        //D( std::cout << "Creating 1x1 Tracer " << this->m_id << " at " << "(" << this->m_site->get_x() << "," << this->m_site->get_y() << ")" << std::endl );
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
                this->m_correlations[tmp_idx-1]++;
                //this->m_last_step_idx[tmp_n] = tmp_idx;
                tmp_n++;
        }
        // for(int idx : this->m_last_step_idx) { this->m_correlations[idx]++; }
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
// function for random walk stepping
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer::step(int dir){
        //D( std::cout << "STEP!" << std::endl );
        //D( std::cout << "dx " << this->m_dx << " dy " << this->m_dy << " MSD " << this->m_lsquared << std::endl );
        //D( std::cout << "Dir: " << dir << std::endl );
        // if(this->m_isstuck) { return; }
        // if(!this->m_site->step_is_valid(dir)) { return; }
        // { return; }
        // this->m_site=this->m_site->move_to_neighbor(dir);
        this->m_site = this->m_site->get_neighbor_by_dir(dir);
        this->m_steps_taken++;
        this->update_last_step(dir);
        switch(dir)
        {
        case 1: this->m_dx++; return;
        case 2: this->m_dy++; return;
        case 3: this->m_dx--; return;
        case 4: this->m_dy--; return;
        }
        //this->m_steps_taken++;

        //D( std::cout << "dx " << this->m_dx << " dy " << this->m_dy << " MSD " << this->m_lsquared << std::endl );
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
// Function for steps during warm_up
// contains only logic for site blocking and trapping
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer::step_warmup(int dir){
        // if(this->m_isstuck) { return; }
        if(!(this->m_site->step_is_valid(dir))) { return; }
        this->m_site=this->m_site->move_to_neighbor(dir);
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Tracer::step_unhindered(int dir){
        this->m_site=this->m_site->move_to_neighbor(dir);
        switch(dir)
        {
        case 1: this->m_dx++; break;
        case 2: this->m_dy++; break;
        case 3: this->m_dx--; break;
        case 4: this->m_dy--; break;
        }
        //this->m_steps_taken++;
        this->m_steps_taken++;
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
        return this->m_site->get_x();
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
int Tracer::get_dx()
{
        return this->m_dx;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
int Tracer::get_y()
{
        return this->m_site->get_y();
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
        return pow((double)this->m_dx,2.0)+pow((double)this->m_dy,2.0);
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
std::vector<long> Tracer::get_correlations(){
        return this->m_correlations;
}

int Tracer::get_steps_taken()
{
        int tmp_steps_taken = this->m_steps_taken;
        this->m_steps_taken = 0;
        return tmp_steps_taken;
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
Tracer_2x2::Tracer_2x2(int id, Site * starting_site) : Tracer(id,starting_site)
{
        // change default size
        this->m_size = 4;
        //D( std::cout << "Creating 2x2 Tracer " << this->m_id << " at " << "(" << this->m_site->get_x() << "," << this->m_site->get_y() << ")" << std::endl );
}
