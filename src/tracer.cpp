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
        // m_steps_taken(0),
        m_powers_of_four({1,4,16}),
        m_last_step_dir(3,0),
        m_correlations(84,0),
        m_isstuck(false)
{
        this->m_site->set_not_empty();
        //D( std::cout << "Creating 1x1 Tracer " << this->m_id << " at " << "(" << this->m_site->get_x() << "," << this->m_site->get_y() << ")" << std::endl );
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
// try: inlining
// inline
void Tracer::update_last_step(int last_step_dir)
{
        // shift m_last_step_dir by one position
        std::rotate(this->m_last_step_dir.rbegin(),this->m_last_step_dir.rbegin()+1,this->m_last_step_dir.rend());
        // overwrite the move that's now at pos 0 with the most recent one
        this->m_last_step_dir[0] = last_step_dir;
        // temporary variables for sorting the resulting step sequences into the counting vector
        int tmp_n = 1;
        int tmp_idx = 0;
        /*
           TODO: rewrite with STL function

         */
        // std::accumulate(this->m_last_step_dir.begin(),this->m_last_step_dir.end(),this->m_powers_of_four.begin(),0,[](int dir, int pow) -> int { })
        for(int step : this->m_last_step_dir)
        {
                tmp_idx += tmp_n*step;
                tmp_n   *= 4;
                this->m_correlations[tmp_idx-1] += !!step; // adding "not not step" results in incrementation, but only if step is not zero
        }
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
// function for random walk stepping
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Tracer::step(int dir){
        D( std::cout << "Step attempt:" << std::endl );
        D( std::cout << this->m_id << " -> " << dir << std::endl );
        //D( std::cout << "STEP!" << std::endl );
        //D( std::cout << "dx " << this->m_dx << " dy " << this->m_dy << " MSD " << this->m_lsq << std::endl );
        D( std::cout << "Dir: " << dir << std::endl );
        // if(this->m_isstuck) { return; }
        if( !(this->m_site->step_is_invalid(dir)) )
        {
                // this->m_steps_taken++;
                this->m_site = this->m_site->jump_in_direction(dir);
                this->update_last_step(dir);
        }
        //D( std::cout << "dx " << this->m_dx << " dy " << this->m_dy << " MSD " << this->m_lsq << std::endl );
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
// Function for steps during warm_up
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Tracer::step_warmup(int dir){
        // if(this->m_isstuck) { return; }
        if( !(this->m_site->step_is_invalid(dir)) ) {
                // return; }
                this->m_site = this->m_site->jump_in_direction(dir);
        }
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Tracer::step_unhindered(int dir){
        // this->m_steps_taken++;
        this->update_last_step(dir);
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
long Tracer::get_dx()
{
        return this->m_correlations[0]-this->m_correlations[2];
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
long Tracer::get_dy()
{
        return this->m_correlations[1]-this->m_correlations[3];
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
double Tracer::get_lsq()
{
        long dx = this->m_correlations[0]-this->m_correlations[2];
        long dy = this->m_correlations[1]-this->m_correlations[3];
        return (double)(dx*dx+dy*dy);
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
bool Tracer::get_isstuck()
{
        return this->m_isstuck;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<long> Tracer::get_correlations(){
        return this->m_correlations;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<std::vector<Site*> > Tracer::get_blocking_sites(){
        return this->m_site->get_blocking_sites();
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
//unsigned long long
double Tracer::get_steps_taken()
{
        // return (double)this->m_steps_taken;
        return (double)std::accumulate(this->m_correlations.begin(),this->m_correlations.begin()+4,0);
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<int> Tracer::get_site_correlation()
{
        return this->m_site->get_site_correlation();
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
