#include "lattice.hpp"
#include "wrapper.hpp"
// - - - - - - - - - - - - - - - - - - - - - - - -
// S I M U L A T I O N W R A P P E R
// to store multiple lattices and perform mc simulations simultaneously
// - - - - - - - - - - - - - - - - - - - - - - - -
Wrapper::Wrapper(int number_of_lattices, int grid_size_x, int grid_size_y, int number_of_timesteps, int number_of_timesteps_warmup, int number_of_tracers_1x1, int number_of_tracers_2x2, int number_of_tracers_3x3, double step_rate_1x1, double step_rate_2x2, double step_rate_3x3, int wtd_max, int wtd_res) :
        // sim parameters
        m_number_of_lattices(number_of_lattices),
        m_grid_size_x(grid_size_x),
        m_grid_size_y(grid_size_y),
        m_number_of_timesteps(number_of_timesteps),
        m_t(0),
        m_number_of_tracers_1x1(number_of_tracers_1x1),
        m_number_of_tracers_2x2(number_of_tracers_2x2),
        m_number_of_tracers_3x3(number_of_tracers_3x3),
        m_step_rate_1x1(step_rate_1x1),
        m_step_rate_2x2(step_rate_2x2),
        m_step_rate_3x3(step_rate_3x3),
        m_wtd_max(wtd_max),
        m_wtd_res(wtd_res),
        // some auxilliary values
        m_this_number_of_tracers_1x1_times_number_of_lattices((double)number_of_tracers_1x1*(double)numnumber_of_lattices),
        m_this_number_of_tracers_2x2_times_number_of_lattices((double)number_of_tracers_2x2*(double)numnumber_of_lattices),
        m_this_number_of_tracers_3x3_times_number_of_lattices((double)number_of_tracers_3x3*(double)numnumber_of_lattices),
        // time series of ensemble avg step rate of 1x1/2x2/3x3 tracers
        m_avg_rate_1x1(number_of_tracers_1x1 > 0 ? number_of_timesteps : 0,0),
        m_avg_rate_2x2(number_of_tracers_2x2 > 0 ? number_of_timesteps : 0,0),
        m_avg_rate_3x3(number_of_tracers_3x3 > 0 ? number_of_timesteps : 0,0),
        // time series of ensemble avg lsq of 1x1/2x2/3x3 tracers
        m_avg_lsquared_1x1(number_of_tracers_1x1 > 0 ? number_of_timesteps : 0,0),
        m_avg_lsquared_2x2(number_of_tracers_2x2 > 0 ? number_of_timesteps : 0,0),
        m_avg_lsquared_3x3(number_of_tracers_3x3 > 0 ? number_of_timesteps : 0,0),
        // waiting time distributions
        m_wtd_1x1(number_of_tracers_1x1 > 0 ? 16*(wtd_max*wtd_res+1) : 0,0),
        m_wtd_2x2(number_of_tracers_2x2 > 0 ? 16*(wtd_max*wtd_res+1) : 0,0),
        m_wtd_3x3(number_of_tracers_3x3 > 0 ? 16*(wtd_max*wtd_res+1) : 0,0),
        // time- and ensemble averaged conditional probabilities for series of steps with length 2/3/4 for 1x1/2x2/3x3 tracers
        // 1x1
        m_correlations_1x1(number_of_tracers_1x1 > 0 ?  341 : 0,0),
        m_correlations_2x2(number_of_tracers_2x2 > 0 ?  341 : 0,0),
        m_correlations_3x3(number_of_tracers_3x3 > 0 ?  341 : 0,0)
{
        this->m_lattices.reserve(this->m_number_of_lattices);
        while(this->m_lattices.size() < m_number_of_lattices) {
                this->m_lattices.push_back(new Lattice(this->m_grid_size_x,this->m_grid_size_y,this->m_number_of_tracers_1x1,this->m_number_of_tracers_2x2,this->m_number_of_tracers_3x3,this->m_step_rate_1x1,this->m_step_rate_2x2,this->m_step_rate_3x3,this->m_wtd_max));
        }
}
//
void Wrapper::timestep(){
        for(Lattice * l : this->m_lattices)
        {
                l->timestep();
        }
        // update data
        this->update_data();
        // increment time
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
inline void Wrapper::update_data(){
        if(this->m_number_of_tracers_1x1) {
                this->m_avg_rate_1x1[this->m_t] = this->get_avg_rate_1x1();
                this->m_avg_lsquared_1x1[this->m_t] = this->get_avg_lsquared_1x1();
                this->update_correlations_1x1();
        }
        if(this->m_number_of_tracers_2x2) {
                this->m_avg_rate_2x2[this->m_t] = this->get_avg_rate_2x2();
                this->m_avg_lsquared_2x2[this->m_t] = this->get_avg_lsquared_2x2();
                this->update_correlations_2x2();
        }
        if(this->m_number_of_tracers_3x3) {
                this->m_avg_rate_3x3[this->m_t] = this->get_avg_rate_3x3();
                this->m_avg_lsquared_3x3[this->m_t] = this->get_avg_lsquared_3x3();
                this->update_correlations_3x3();
        }
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
double Wrapper::get_avg_rate_1x1(){
        double tmp_sum_steps_taken = 0;
        for(Lattice * l : this->m_lattices)
        {
                for(Tracer * tr : l->get_tracers_1x1())
                {
                        tmp_sum_steps_taken += (double)tr->get_last_step();
                }
        }
        return tmp_sum/this->m_this_number_of_tracers_1x1_times_number_of_lattices;
}
double Wrapper::get_avg_rate_2x2(){
        double tmp_sum_steps_taken = 0;
        for(Lattice * l : this->m_lattices)
        {
                for(Tracer * tr : l->get_tracers_2x2())
                {
                        tmp_sum_steps_taken += (double)tr->get_last_step();
                }
        }
        return tmp_sum/this->m_this_number_of_tracers_2x2_times_number_of_lattices;
}
double Wrapper::get_avg_rate_3x3(){
        double tmp_sum_steps_taken = 0;
        for(Lattice * l : this->m_lattices)
        {
                for(Tracer * tr : l->get_tracers_3x3())
                {
                        tmp_sum_steps_taken += (double)tr->get_last_step();
                }
        }
        return tmp_sum/this->m_this_number_of_tracers_3x3_times_number_of_lattices;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
double Wrapper::get_avg_lsquared_1x1(){
        double tmp_sum_lsquared = 0;
        for(Lattice * l : this->m_lattices)
        {
                for(Tracer * tr : l->get_tracers_1x1())
                {
                        tmp_sum_lsquared += tr->get_lsquared();
                }
        }
        return tmp_sum_lsquared/this->m_this_number_of_tracers_1x1_times_number_of_lattices;
}
double Wrapper::get_avg_lsquared_2x2(){
        double tmp_sum_lsquared = 0;
        for(Lattice * l : this->m_lattices)
        {
                for(Tracer * tr : l->get_tracers_2x2())
                {
                        tmp_sum_lsquared += tr->get_lsquared();
                }
        }
        return tmp_sum_lsquared/this->m_this_number_of_tracers_2x2_times_number_of_lattices;
}
double Wrapper::get_avg_lsquared_3x3(){
        double tmp_sum_lsquared = 0;
        for(Lattice * l : this->m_lattices)
        {
                for(Tracer * tr : l->get_tracers_3x3())
                {
                        tmp_sum_lsquared += tr->get_lsquared();
                }
        }
        return tmp_sum_lsquared/this->m_this_number_of_tracers_3x3_times_number_of_lattices;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
void Wrapper::update_correlations_1x1(){
        // iterate over all lattices
        for(Lattice * l : this->lattices)
        {
                // iterate over all tracers (of type X) on current lattice
                for(Tracer * tr : l->get_tracers_1x1())
                {
                        // check if current tracer has a contribution to correlation in this timestep
                        if(tr->get_last_step())
                        {
                                // update the 2-step waiting time distribution
                                this->m_wtd_1x1[tr->get_wtd_index()]++;
                                // count how often the different indexes occured
                                for(int idx : tr->get_last_step_idx()) {
                                        this->m_correlations_1x1[idx]++;
                                }
                        }

                }
        }
}
void Wrapper::update_correlations_2x2(){
        // iterate over all lattices
        for(Lattice * l : this->lattices)
        {
                // iterate over all tracers (of type X) on current lattice
                for(Tracer * tr : l->get_tracers_2x2())
                {
                        // check if current tracer has a contribution to correlation in this timestep
                        if(tr->get_last_step())
                        {
                                // update the 2-step waiting time distribution
                                this->m_wtd_2x2[tr->get_wtd_index()]++;
                                // count how often the different indexes occured
                                for(int idx : tr->get_last_step_idx()) {
                                        this->m_correlations_2x2[idx]++;
                                }
                        }

                }
        }
}
void Wrapper::update_correlations_3x3(){
        // iterate over all lattices
        for(Lattice * l : this->lattices)
        {
                // iterate over all tracers (of type X) on current lattice
                for(Tracer * tr : l->get_tracers_3x3())
                {
                        // check if current tracer has a contribution to correlation in this timestep
                        if(tr->get_last_step())
                        {
                                // update the 2-step waiting time distribution
                                this->m_wtd_3x3[tr->get_wtd_index()]++;
                                // count how often the different indexes occured
                                for(int idx : tr->get_last_step_idx()) {
                                        this->m_correlations_3x3[idx]++;
                                }
                        }

                }
        }
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<int> Wrapper::get_wtd_1x1()
{
        if(this->m_number_of_tracers_1x1 == 0)
        {
                std::vector<int> tmp_wtd = {0};
                return tmp_wtd;
        }
        std::vector<int> tmp_wtd;
        std::vector<int> tmp_current_wtd;
        tmp_wtd.reserve(4*(this->m_wtd_max+1));
        tmp_current_wtd.reserve(4*(this->m_wtd_max+1));

        for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
        {
                tmp_wtd.push_back(0);
        }
        for(Lattice * l : this->m_lattices)
        {
                tmp_current_wtd = l->get_wtd_1x1();
                for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
                {
                        tmp_wtd[n] += tmp_current_wtd[n];
                }

        }
        return tmp_wtd;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> Wrapper::get_wtd_2x2()
{
        if(this->m_number_of_tracers_2x2 == 0)
        {
                std::vector<int> tmp_wtd = {0};
                return tmp_wtd;
        }
        std::vector<int> tmp_wtd;
        std::vector<int> tmp_current_wtd;
        tmp_wtd.reserve(4*(this->m_wtd_max+1));
        tmp_current_wtd.reserve(4*(this->m_wtd_max+1));

        for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
        {
                tmp_wtd.push_back(0);
        }
        tmp_wtd.reserve((this->m_wtd_max+1));
        for(Lattice * l : this->m_lattices)
        {
                tmp_current_wtd = l->get_wtd_2x2();
                for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
                {
                        tmp_wtd[n] += tmp_current_wtd[n];
                }
        }
        return tmp_wtd;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<int> Wrapper::get_wtd_3x3()
{
        if(this->m_number_of_tracers_3x3 == 0)
        {
                std::vector<int> tmp_wtd = {0};
                return tmp_wtd;
        }
        std::vector<int> tmp_wtd;
        std::vector<int> tmp_current_wtd;
        tmp_wtd.reserve(4*(this->m_wtd_max+1));
        tmp_current_wtd.reserve(4*(this->m_wtd_max+1));

        for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
        {
                tmp_wtd.push_back(0);
        }
        for(Lattice * l : this->m_lattices)
        {
                tmp_current_wtd = l->get_wtd_3x3();
                for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
                {
                        tmp_wtd[n] += tmp_current_wtd[n];
                }
        }
        return tmp_wtd;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Wrapper::get_norm_wtd_1x1()
{
        if(this->m_number_of_tracers_1x1 == 0)
        {
                std::vector<double> tmp_wtd = {0};
                return tmp_wtd;
        }
        std::vector<double> tmp_wtd;
        std::vector<double> tmp_current_wtd;
        tmp_current_wtd.reserve(4*(this->m_wtd_max+1));
        tmp_wtd.reserve(4*(this->m_wtd_max+1));
        for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
        {
                tmp_wtd.push_back(0);
        }
        for(Lattice * l : this->m_lattices)
        {
                tmp_current_wtd = l->get_norm_wtd_1x1();
                for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
                {
                        tmp_wtd[n] += (double)tmp_current_wtd[n];
                }
        }
        for(double & d : tmp_wtd)
        {
                d /= this->m_number_of_lattices;
        }
        return tmp_wtd;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Wrapper::get_norm_wtd_2x2()
{
        if(this->m_number_of_tracers_2x2 == 0)
        {
                std::vector<double> tmp_wtd = {0};
                return tmp_wtd;
        }
        std::vector<double> tmp_wtd;
        std::vector<double> tmp_current_wtd;
        tmp_current_wtd.reserve(4*(this->m_wtd_max+1));
        tmp_wtd.reserve(4*(this->m_wtd_max+1));
        for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
        {
                tmp_wtd.push_back(0);
        }
        for(Lattice * l : this->m_lattices)
        {
                tmp_current_wtd = l->get_norm_wtd_2x2();
                for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
                {
                        tmp_wtd[n] += (double)tmp_current_wtd[n];
                }
        }
        for(double & d : tmp_wtd)
        {
                d /= this->m_number_of_lattices;
        }
        return tmp_wtd;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Wrapper::get_norm_wtd_3x3()
{
        if(this->m_number_of_tracers_3x3 == 0)
        {
                std::vector<double> tmp_wtd = {0};
                return tmp_wtd;
        }
        std::vector<double> tmp_wtd;
        std::vector<double> tmp_current_wtd;
        tmp_current_wtd.reserve(4*(this->m_wtd_max+1));
        tmp_wtd.reserve(4*(this->m_wtd_max+1));
        for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
        {
                tmp_wtd.push_back(0);
        }
        for(Lattice * l : this->m_lattices)
        {
                tmp_current_wtd = l->get_norm_wtd_3x3();
                for(int n = 0; n < 4*(this->m_wtd_max+1); n++)
                {
                        tmp_wtd[n] += (double)tmp_current_wtd[n];
                }
        }
        for(double & d : tmp_wtd)
        {
                d /= this->m_number_of_lattices;
        }
        return tmp_wtd;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Wrapper::get_avg_two_step_correlation_1x1()
{
        std::vector<double> tmp_vec(4,0);
        for(Lattice * l : this->m_lattices)
        {
                int idx = 0;
                for(double tmp_dbl : l->get_avg_two_step_correlation_1x1())
                {
                        tmp_vec[idx] += tmp_dbl;
                        idx++;
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= this->m_number_of_lattices;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Wrapper::get_avg_two_step_correlation_2x2()
{
        std::vector<double> tmp_vec(4,0);
        for(Lattice * l : this->m_lattices)
        {
                int idx = 0;
                for(double tmp_dbl : l->get_avg_two_step_correlation_2x2())
                {
                        tmp_vec[idx] += tmp_dbl;
                        idx++;
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= this->m_number_of_lattices;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Wrapper::get_avg_two_step_correlation_3x3()
{
        std::vector<double> tmp_vec(4,0);
        for(Lattice * l : this->m_lattices)
        {
                int idx = 0;
                for(double tmp_dbl : l->get_avg_two_step_correlation_3x3())
                {
                        tmp_vec[idx] += tmp_dbl;
                        idx++;
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= this->m_number_of_lattices;
        }
        return tmp_vec;
}
// average three step correlations (sums of the directional three step waiting time distributions)
std::vector<double> Wrapper::get_avg_three_step_correlation_1x1()
{
        std::vector<double> tmp_vec(16,0);
        for(Lattice * l : this->m_lattices)
        {
                int idx = 0;
                for(double tmp_dbl : l->get_avg_three_step_correlation_1x1())
                {
                        tmp_vec[idx] += tmp_dbl;
                        idx++;
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= this->m_number_of_lattices;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Wrapper::get_avg_three_step_correlation_2x2()
{
        std::vector<double> tmp_vec(16,0);
        for(Lattice * l : this->m_lattices)
        {
                int idx = 0;
                for(double tmp_dbl : l->get_avg_three_step_correlation_2x2())
                {
                        tmp_vec[idx] += tmp_dbl;
                        idx++;
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= this->m_number_of_lattices;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Wrapper::get_avg_three_step_correlation_3x3()
{
        std::vector<double> tmp_vec(16,0);
        for(Lattice * l : this->m_lattices)
        {
                int idx = 0;
                for(double tmp_dbl : l->get_avg_three_step_correlation_3x3())
                {
                        tmp_vec[idx] += tmp_dbl;
                        idx++;
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= this->m_number_of_lattices;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Wrapper::get_avg_four_step_correlation_1x1()
{
        std::vector<double> tmp_vec(64,0);
        for(Lattice * l : this->m_lattices)
        {
                int idx = 0;
                for(double tmp_dbl : l->get_avg_four_step_correlation_1x1())
                {
                        tmp_vec[idx] += tmp_dbl;
                        idx++;
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= this->m_number_of_lattices;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Wrapper::get_avg_four_step_correlation_2x2()
{
        std::vector<double> tmp_vec(64,0);
        for(Lattice * l : this->m_lattices)
        {
                int idx = 0;
                for(double tmp_dbl : l->get_avg_four_step_correlation_2x2())
                {
                        tmp_vec[idx] += tmp_dbl;
                        idx++;
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= this->m_number_of_lattices;
        }
        return tmp_vec;
}
// - - - - - - - - - - - - - - - - - - - - - - - - -
std::vector<double> Wrapper::get_avg_four_step_correlation_3x3()
{
        std::vector<double> tmp_vec(64,0);
        for(Lattice * l : this->m_lattices)
        {
                int idx = 0;
                for(double tmp_dbl : l->get_avg_four_step_correlation_3x3())
                {
                        tmp_vec[idx] += tmp_dbl;
                        idx++;
                }
        }
        for(double & tmp_dbl : tmp_vec)
        {
                tmp_dbl /= this->m_number_of_lattices;
        }
        return tmp_vec;
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
std::vector<double> Wrapper::get_result_lsquared_3x3(){
        return this->m_avg_lsquared_3x3;
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
std::vector<double> Wrapper::get_result_rate_3x3(){
        return this->m_avg_rate_3x3;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
// waiting time distributions for the 4 possible two-step configurations
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> get_result_wtd_1x1(){
        return this->m_wtd_1x1;
        // TODO: normalize the vector, such that the sum over all the wtds is 1.0
        //std::vector<double> tmp_wtd; // (this->m_wtd_1x1.size(),0.0);
        //tmp_wtd.reserve(this->m_wtd_1x1.size());
        //int tmp_norm = std::accumulate(this->m_wtd_1x1.begin(), this->m_wtd_1x1.end(), 0);
        //std::transform(this->m_wtd_1x1.begin(), this->m_wtd_1x1.end(), std::back_inserter(tmp_wtd), );
        //return tmp_wtd;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> get_result_wtd_2x2(){
        return this->m_wtd_2x2;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> get_result_wtd_3x3(){
        return this->m_wtd_3x3;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> get_result_correlations_1x1(){
        return this->m_correlations_1x1;
        //std::vector<double> tmp_correlations(this->m_correlations_1x1.size(),0.0);
        //return tmp_correlations;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> get_result_correlations_2x2(){
        return this->m_correlations_2x2;

        //std::vector<double> tmp_correlations(this->m_correlations_2x2.size(),0.0);
        //return tmp_correlations;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> get_result_correlations_3x3(){
        return this->m_correlations_3x3;

        //std::vector<double> tmp_correlations(this->m_correlations_3x3.size(),0.0);
        //return tmp_correlations;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
