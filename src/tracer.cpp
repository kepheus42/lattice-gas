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
        // m_powers_of_four({1,4,16}),
        m_last_step_dir(3), // ,0
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
        this->m_last_step_dir.push_front(last_step_dir); // -1
        // = = = = = = = = = = = =
        // TODO: update to storing of only 2 steps, use function variable last_step_dir for calculations and push_front(last_step_dir) at the end
        // this->m_correlations[last_step_dir - 1]++;  // NEW
        // = = = = = = = = = = = =
        // favored version:
        // uses one STL function call
        // this->m_correlations[std::transform_reduce(this->m_last_step_dir.begin(),this->m_last_step_dir.end(),this->m_powers_of_four.begin(),-1)]++;
        // also workable:
        // this->m_correlations[this->m_last_step_dir[0]+4*this->m_last_step_dir[1]+16*this->m_last_step_dir[2]-1]++;
        // both versions require reworking getter functions, as 1- and 2-step correlations are no longer stored, but have to be calculated from 3-step correlations

        // temporary variables for sorting the resulting step sequences into the counting vector
        int tmp_n = 1;
        // int tmp_n = 4; // NEW
        int tmp_idx = 0;
        // int tmp_idx = last_step_dir;  // NEW
        for(int step : this->m_last_step_dir)
        {
                tmp_idx += tmp_n * step;
                tmp_n   *= 4;
                this->m_correlations[tmp_idx-1] += !!step; // adding "not not step" results in incrementation, but only if step is not zero
        }
        // this->m_last_step_dir.push_front(last_step_dir); // NEW
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
        if( this->m_site->step_is_valid(dir) )
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
        if( this->m_site->step_is_valid(dir) ) {
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
std::vector<double> Tracer::get_probabilities(){
        // copy cast vector for probabilities (double)
        std::vector<double> tmp_probabilities(this->m_correlations.begin(),this->m_correlations.end());
        // calculate norm for 1-3 step sequences
        double tmp_norm_1 = 1/(double)std::accumulate(this->m_correlations.begin()+0,this->m_correlations.begin()+4,0);
        double tmp_norm_2 = 1/(double)std::accumulate(this->m_correlations.begin()+4,this->m_correlations.begin()+20,0);
        double tmp_norm_3 = 1/(double)std::accumulate(this->m_correlations.begin()+20,this->m_correlations.end(),0);
        // cast count vector to double
        // apply norm
        std::transform(tmp_probabilities.begin()+0,
                       tmp_probabilities.begin()+4,
                       tmp_probabilities.begin(),
                       std::bind(std::multiplies<double>(), std::placeholders::_1, tmp_norm_1));
        std::transform(tmp_probabilities.begin()+4,
                       tmp_probabilities.begin()+20,
                       tmp_probabilities.begin()+4,
                       std::bind(std::multiplies<double>(), std::placeholders::_1, tmp_norm_2));
        std::transform(tmp_probabilities.begin()+20,
                       tmp_probabilities.end(),
                       tmp_probabilities.begin()+20,
                       std::bind(std::multiplies<double>(), std::placeholders::_1, tmp_norm_3));
        // return result
        return tmp_probabilities;
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
