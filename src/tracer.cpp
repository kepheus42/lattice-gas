#include "tracer.hpp"

// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   To add debugging messages, use:
   D(std::cerr << "Debugging message 1 2 3!" << std::endl );
   D(std::cerr << debugging_code() << std::endl );
   etc.
   = = = = = = = = = = = = = = = = = = = = = = = = = = = = */

/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
Tracer::Tracer(int id, int type, Site * site) :
        m_id(id),
        m_type(type),
        m_site(site),
        m_dx(0),
        m_dy(0),
        m_last_step_dir(3,0),
        m_correlations(84,0),
        m_steps_taken(0),
        m_isstuck(false)
{
        this->m_site->swap_state();
        //D( std::cout << "Creating 1x1 Tracer " << this->m_id << " at " << "(" << this->m_site->get_x() << "," << this->m_site->get_y() << ")" << std::endl );
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
inline
void Tracer::update_last_step(int last_step_dir)
{
        // shift m_last_step_dir by one position
        std::rotate(this->m_last_step_dir.rbegin(),this->m_last_step_dir.rbegin()+1,this->m_last_step_dir.rend());
        // overwrite the move that's now at pos 0 with the most recent one
        // this->m_last_step_dir[2] = this->m_last_step_dir[1];
        // this->m_last_step_dir[1] = this->m_last_step_dir[0];
        this->m_last_step_dir[0] = last_step_dir;
        // this->m_last_sequence = last_step_dir;
        // this->m_correlations[] += 1;
        // this->m_correlations[] += 1;
        int tmp_n = 1;
        int tmp_idx = 0;
        /* */
        for(int step : this->m_last_step_dir)
        {
                tmp_idx += tmp_n*step;
                tmp_n *= 4;
                this->m_correlations[tmp_idx-1] += !!step; // adding "not not step" results in incrementation, but only if step is not zero
        }
        /* */
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
// function for random walk stepping
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Tracer::step(int dir){
        D( std::cout << this->m_id << " -> " << dir << std::endl );
        //D( std::cout << "STEP!" << std::endl );
        //D( std::cout << "dx " << this->m_dx << " dy " << this->m_dy << " MSD " << this->m_lsq << std::endl );
        D( std::cout << "Dir: " << dir << std::endl );
        // if(this->m_isstuck) { return; }
        if( this->m_site->step_is_valid(dir) )
        { // return; }
                this->m_site=this->m_site->move(dir);
                this->m_steps_taken++;
                this->update_last_step(dir);
                // Try:
                // this->m_dx -= (dir-2)%2;
                // this->m_dy -= (dir-3)%2;
                //
                switch(dir)
                {
                case 1: this->m_dx++; return;
                case 2: this->m_dy++; return;
                case 3: this->m_dx--; return;
                case 4: this->m_dy--; return;
                }
        }
        //D( std::cout << "dx " << this->m_dx << " dy " << this->m_dy << " MSD " << this->m_lsq << std::endl );
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
// Function for steps during warm_up
// contains only logic for site blocking and trapping
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Tracer::step_warmup(int dir){
        // if(this->m_isstuck) { return; }
        if( this->m_site->step_is_valid(dir) ) {
                // return; }
                this->m_site = this->m_site->move(dir);
        }
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Tracer::step_unhindered(int dir){
        this->m_steps_taken++;
        this->update_last_step(dir);
        switch(dir)
        {
        case 1: this->m_dx++; return;
        case 2: this->m_dy++; return;
        case 3: this->m_dx--; return;
        case 4: this->m_dy--; return;
        }
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
int Tracer::get_id()
{
        return this->m_id;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
int Tracer::get_pos()
{
        return this->m_site->get_id();
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
int Tracer::get_dx()
{
        return this->m_dx;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
int Tracer::get_dy()
{
        return this->m_dy;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
int Tracer::get_x()
{
        return this->m_site->get_x();
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
int Tracer::get_y()
{
        return this->m_site->get_y();
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
int Tracer::get_type()
{
        return this->m_type;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
unsigned long long Tracer::get_lsq()
{
        return (unsigned long long)((this->m_dx*this->m_dx)+(this->m_dy*this->m_dy));
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
bool Tracer::get_isstuck()
{
        return this->m_isstuck;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<unsigned long long> Tracer::get_correlations(){
        return this->m_correlations;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
unsigned long long Tracer::get_steps_taken()
{
        unsigned long long tmp_steps_taken = this->m_steps_taken;
        this->m_steps_taken = 0;
        return tmp_steps_taken;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
unsigned long long Tracer::get_steps_taken_total()
{
        return this->m_steps_taken;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Tracer::stuck()
{
        this->m_isstuck = true;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Tracer::unstuck()
{
        this->m_isstuck = false;
}
