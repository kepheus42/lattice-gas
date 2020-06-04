#include "lattice.hpp"
#include "wrapper.hpp"


// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif

// - - - - - - - - - - - - - - - - - - - - - - - -
// S I M U L A T I O N W R A P P E R
// to store multiple lattices and perform mc simulations simultaneously
// - - - - - - - - - - - - - - - - - - - - - - - -
Wrapper::Wrapper(int number_of_lattices,
                 int grid_size_x,
                 int grid_size_y,
                 int number_of_timesteps,
                 int number_of_timesteps_warmup,
                 int number_of_tracers_1x1,
                 int number_of_tracers_2x2,
                 double step_rate_1x1,
                 double step_rate_2x2,
                 int wtd_max,
                 int wtd_res) :
        // sim parameters
        m_number_of_lattices(number_of_lattices),
        m_grid_size_x(grid_size_x),
        m_grid_size_y(grid_size_y),
        m_number_of_timesteps(number_of_timesteps),
        m_t(0),
        m_number_of_tracers_1x1(number_of_tracers_1x1),
        m_number_of_tracers_2x2(number_of_tracers_2x2),
        m_step_rate_1x1(step_rate_1x1),
        m_step_rate_2x2(step_rate_2x2),
        m_wtd_max(wtd_max),
        m_wtd_res(wtd_res),
        // some auxilliary values
        m_number_of_tracers_times_number_of_lattices((number_of_tracers_1x1+number_of_tracers_2x2)*number_of_lattices),
        m_number_of_tracers_1x1_times_number_of_lattices(number_of_tracers_1x1*number_of_lattices),
        m_number_of_tracers_2x2_times_number_of_lattices(number_of_tracers_2x2*number_of_lattices),
        // time series of ensemble avg step rate of 1x1/2x2 tracers
        m_avg_rate_1x1(number_of_tracers_1x1 > 0 ? number_of_timesteps : 0,0),
        m_avg_rate_2x2(number_of_tracers_2x2 > 0 ? number_of_timesteps : 0,0),
        // time series of ensemble avg lsq of 1x1/2x2 tracers
        m_avg_lsquared_1x1(number_of_tracers_1x1 > 0 ? number_of_timesteps : 0,0),
        m_avg_lsquared_2x2(number_of_tracers_2x2 > 0 ? number_of_timesteps : 0,0),
        // time series of the sublattice coverage
        m_sublattice_conc_1x1(number_of_tracers_1x1 > 0 ? 4*number_of_timesteps : 0,0),
        m_sublattice_conc_2x2(number_of_tracers_2x2 > 0 ? 4*number_of_timesteps : 0,0),
        // waiting time distributions
        m_wtd_1x1(number_of_tracers_1x1 > 0 ? 1+16*(wtd_max*wtd_res+1) : 0,0),
        m_wtd_2x2(number_of_tracers_2x2 > 0 ? 1+16*(wtd_max*wtd_res+1) : 0,0),
        // time- and ensemble averaged conditional probabilities for series of steps with length 2/3/4 for 1x1/2x2 tracers
        // 1x1
        m_correlations_1x1(number_of_tracers_1x1 > 0 ?  341 : 0,0),
        m_correlations_2x2(number_of_tracers_2x2 > 0 ?  341 : 0,0),
{
        this->m_lattices.reserve(this->m_number_of_lattices);
        while(this->m_lattices.size() < m_number_of_lattices) {
                this->m_lattices.push_back(new Lattice(this->m_grid_size_x,
                                                       this->m_grid_size_y,
                                                       this->m_number_of_timesteps,
                                                       this->m_number_of_tracers_1x1,
                                                       this->m_number_of_tracers_2x2,
                                                       this->m_step_rate_1x1,
                                                       this->m_step_rate_2x2,
                                                       this->m_wtd_max,
                                                       this->m_wtd_res));
        }
        for(Lattice * l : this->m_lattices)
        {
                for(Tracer * tr : l->get_tracers_1x1)
                {
                        this->m_tracers.push_back(tr);
                        this->m_tracers_1x1.push_back(tr);
                }
                for(Tracer * tr : l->get_tracers_2x2)
                {
                        this->m_tracers.push_back(tr);
                        this->m_tracers_2x2.push_back(tr);
                }
        }
        D(std::cerr << "Initializing Wrapper: Size of m_tracers vectors" << std::endl
                    << "All: " << this->m_tracers.size() << " expected: " << this->m_number_of_tracers_times_number_of_lattices << std::endl
                    << "1x1: " << this->m_tracers_1x1.size() << " expected: " << this->m_number_of_tracers_1x1_times_number_of_lattices << std::endl
                    << "2x2: " << this->m_tracers_2x2.size() << " expected: " << this->m_number_of_tracers_2x2_times_number_of_lattices << std::endl);
}
//
void Wrapper::timestep(){
        for(Lattice * l : this->m_lattices)
        {
                l->timestep();
        }
        // update data
        this->update_data();
        this->m_t++;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
void Wrapper::timestep_warmup(){
        for(Lattice * l : this->m_lattices)
        {
                l->timestep_warmup();
        }
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
int Wrapper::get_t()
{
        return this->m_t;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
inline int Lattice::coord(int x, int y){
        //
        return ((x+this->m_grid_size_x)%this->m_grid_size_x)*this->m_grid_size_y+((y+this->m_grid_size_y)%this->m_grid_size_y);
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
void Wrapper::update_data(){
        // local vars (initialized only if needed)
        if(this->m_number_of_tracers_1x1)
        {
                unsigned long tmp_steps_taken_1x1 = 0;
                double tmp_lsquared_1x1 = 0.0;
        }
        if(this->m_number_of_tracers_2x2)
        {
                unsigned long tmp_steps_taken_2x2 = 0;
                double tmp_lsquared_2x2 = 0.0;
        }
        for(Tracer * tr: this->m_tracers_1x1)
        {
                tmp_steps_taken_1x1 += tr->get_steps_taken();
                tmp_lsquared_1x1    += tr->get_lsquared();
                this->m_step_attempt_results_1x1[tr->get_step_attempt_result()]++;
                if(tr->get_last_step())
                {
                        this->m_wtd_1x1[tr->get_wtd_idx()]++;
                        for(int idx : tr->get_last_step_idx()) { this->m_correlations_1x1[idx]++; }
                }
        }
        for(Tracer * tr: this->m_tracers_2x2)
        {
                tmp_steps_taken_2x2 += tr->get_steps_taken();
                tmp_lsquared_2x2    += tr->get_lsquared();
                this->m_step_attempt_results_2x2[tr->get_step_attempt_result()]++;
                if(tr->get_last_step())
                {
                        this->m_wtd_2x2[tr->get_wtd_idx()]++;
                        for(int idx : tr->get_last_step_idx()) { this->m_correlations_2x2[idx]++; }
                }
        }

        if(this->m_number_of_tracers_1x1)
        {
                this->m_avg_rate_1x1[this->m_t] = (double)tmp_steps_taken_1x1/(this->m_t+1)/this->m_number_of_tracers_1x1_times_number_of_lattices;
                this->m_avg_lsquared_1x1[this->m_t] = tmp_sum_1x1/(this->m_t+1)/this->m_number_of_tracers_1x1_times_number_of_lattices;
        }
        if(this->m_number_of_tracers_2x2)
        {
                this->m_avg_rate_2x2[this->m_t] = (double)tmp_steps_taken_2x2/(this->m_t+1)/this->m_number_of_tracers_2x2_times_number_of_lattices;
                this->m_avg_lsquared_2x2[this->m_t] = tmp_sum_2x2/(this->m_t+1)/this->m_number_of_tracers_2x2_times_number_of_lattices;
        }
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
// Get Simulation Results:
// = = = = = = = = = = = = = = = = = = = = = = = = =
// avg_lsquared(t) for t = [0,...,m_number_of_timesteps]
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> Wrapper::get_result_lsquared_1x1(){
        return this->m_avg_lsquared_1x1;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> Wrapper::get_result_lsquared_2x2(){
        return this->m_avg_lsquared_2x2;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
// avg_rate(t) for t = [0,...,m_number_of_timesteps]
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> Wrapper::get_result_rate_1x1(){
        return this->m_avg_rate_1x1;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> Wrapper::get_result_rate_2x2(){
        return this->m_avg_rate_2x2;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> Wrapper::get_result_step_attempts_1x1(){
        std::vector<double> tmp_step_attempts(3);
        double tmp_div = (double)std::accumulate(this->m_step_attempt_results_1x1.begin(),this->m_step_attempt_results_1x1.end(),0);
        tmp_step_attempts.push_back((double)this->m_step_attempt_results_1x1[0]/tmp_div);
        tmp_step_attempts.push_back((double)this->m_step_attempt_results_1x1[1]/tmp_div);
        tmp_step_attempts.push_back((double)this->m_step_attempt_results_1x1[2]/tmp_div);
        return tmp_step_attempts;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> Wrapper::get_result_step_attempts_2x2(){
        std::vector<double> tmp_step_attempts(7);
        double tmp_div = (double)std::accumulate(this->m_step_attempt_results_2x2.begin(),this->m_step_attempt_results_2x2.end(),0);
        tmp_step_attempts.push_back((double)this->m_step_attempt_results_2x2[0]/tmp_div);
        tmp_step_attempts.push_back((double)(this->m_step_attempt_results_2x2[1]+
                                             this->m_step_attempt_results_2x2[4]+
                                             this->m_step_attempt_results_2x2[5])/tmp_div);
        tmp_step_attempts.push_back((double)this->m_step_attempt_results_2x2[2]/tmp_div);
        tmp_step_attempts.push_back((double)(this->m_step_attempt_results_2x2[3]+
                                             this->m_step_attempt_results_2x2[5]+
                                             this->m_step_attempt_results_2x2[6])/tmp_div);
        tmp_step_attempts.push_back((double)this->m_step_attempt_results_2x2[4]/tmp_div);
        tmp_step_attempts.push_back((double)this->m_step_attempt_results_2x2[5]/tmp_div);
        tmp_step_attempts.push_back((double)this->m_step_attempt_results_2x2[6]/tmp_div);
        return tmp_step_attempts;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> Wrapper::get_result_sublattice_conc_1x1(){
        std::vector<double> tmp_v(this->m_sublattice_conc_1x1.size(),0.0);
        int tmp_div = this->m_number_of_lattices * this->m_number_of_tracers_1x1;
        std::transform(this->m_sublattice_conc_1x1.begin(),
                       this->m_sublattice_conc_1x1.end(),
                       tmp_v.begin(),
                       [tmp_div](int n) -> double {
                return (double)n/tmp_div;
        });
        return tmp_v;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> Wrapper::get_result_sublattice_conc_2x2(){
        std::vector<double> tmp_v(this->m_sublattice_conc_2x2.size(),0.0);
        int tmp_div = this->m_number_of_lattices * this->m_number_of_tracers_2x2;
        std::transform(this->m_sublattice_conc_2x2.begin(),
                       this->m_sublattice_conc_2x2.end(),
                       tmp_v.begin(),
                       [tmp_div](int n) -> double {
                return (double)n/tmp_div;
        });
        D( std::cout << "Sublattice conc result: " << tmp_v.size() << std::endl);
        D( std::cout << tmp_v[0] << " "
                     << tmp_v[1] << " "
                     << tmp_v[2] << " "
                     << tmp_v[3] << " "
                     << tmp_v[4] << " "
                     << tmp_v[5] << std::endl);
        return tmp_v;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
// waiting time distributions for the 4 possible two-step configurations
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<unsigned long> Wrapper::get_result_wtd_1x1(){
        std::vector<unsigned long> tmp_wtd;
        tmp_wtd.reserve(this->m_wtd_1x1.size()-1);
        tmp_wtd.assign(this->m_wtd_1x1.begin()+1,this->m_wtd_1x1.end());
        return tmp_wtd;
}
std::vector<unsigned long> Wrapper::get_result_wtd_2x2(){
        std::vector<unsigned long> tmp_wtd;
        tmp_wtd.reserve(this->m_wtd_2x2.size()-1);
        tmp_wtd.assign(this->m_wtd_2x2.begin()+1,this->m_wtd_2x2.end());
        return tmp_wtd;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<unsigned long> Wrapper::get_result_correlations_1x1(){
        std::vector<unsigned long> tmp_corr;
        tmp_corr.reserve(this->m_correlations_1x1.size()-1);
        tmp_corr.assign(this->m_correlations_1x1.begin()+1,this->m_correlations_1x1.end());
        return tmp_corr;
}
std::vector<unsigned long> Wrapper::get_result_correlations_2x2(){
        std::vector<unsigned long> tmp_corr;
        tmp_corr.reserve(this->m_correlations_2x2.size()-1);
        tmp_corr.assign(this->m_correlations_2x2.begin()+1,this->m_correlations_2x2.end());
        return tmp_corr;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
// TODO: write functions to return the normalized [0.0-1.0], sum=1.0 version of the wtd
std::vector<double> Wrapper::get_result_norm_wtd_1x1(){
        std::vector<double> tmp_wtd(this->m_wtd_1x1.size()-1,0.0);
        return tmp_wtd;
}
std::vector<double> Wrapper::get_result_norm_wtd_2x2(){
        std::vector<double> tmp_wtd(this->m_wtd_2x2.size()-1,0.0);
        return tmp_wtd;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
// TODO: write functions to return the normalized [0.0-1.0], sum=1.0 version of the correlations
// TODO: Note that this requires a different normalization constant per n-step correlation
std::vector<double> Wrapper::get_result_norm_correlations_1x1(){
        std::vector<double> tmp_corr(this->m_correlations_1x1.size()-1,0.0);
        return tmp_corr;
}
std::vector<double> Wrapper::get_result_norm_correlations_2x2(){
        std::vector<double> tmp_corr(this->m_correlations_2x2.size()-1,0.0);
        return tmp_corr;
}
