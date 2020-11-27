#include "wrapper.hpp"

// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif
// To add debugging messages, use D(std::cerr << "Debugging message 1 2 3!" << std::endl; )

// - - - - - - - - - - - - - - - - - - - - - - - -
// S I M U L A T I O N W R A P P E R
// to store multiple lattices and perform mc simulations simultaneously
// - - - - - - - - - - - - - - - - - - - - - - - -
Wrapper::Wrapper(int number_of_lattices,
                 int grid_size,
                 int number_of_timesteps,
                 int number_of_timesteps_warmup,
                 int number_of_tracers_1x1,
                 int number_of_tracers_2x2,
                 double step_rate_1x1,
                 double step_rate_2x2) :
        // sim parameters
        m_number_of_lattices(number_of_lattices),
        m_grid_size(grid_size),
        m_number_of_timesteps(number_of_timesteps),
        m_t(0),
        m_number_of_tracers_1x1(number_of_tracers_1x1),
        m_number_of_tracers_2x2(number_of_tracers_2x2),
        m_step_rate_1x1(step_rate_1x1),
        m_step_rate_2x2(step_rate_2x2),
        m_data_points(9*(int)std::log10(number_of_timesteps)+number_of_timesteps/(int)std::pow(10,(int) std::log10(number_of_timesteps))),
        m_data_storage_intervals(this->m_data_points,0),
        m_next_data_storage_interval(0),
        // some auxilliary values
        m_number_of_tracers_total_times_number_of_lattices((number_of_tracers_1x1+number_of_tracers_2x2)*number_of_lattices),
        m_number_of_tracers_1x1_times_number_of_lattices(number_of_tracers_1x1 > 0 ? number_of_tracers_1x1*number_of_lattices : 1),
        m_number_of_tracers_2x2_times_number_of_lattices(number_of_tracers_2x2 > 0 ? number_of_tracers_2x2*number_of_lattices : 1),
        // time- and ensemble averaged conditional probabilities for series of steps with length 2/3/4 for 1x1/2x2 tracers
        m_correlations_1x1(number_of_tracers_1x1 > 0 ?  84 : 0,0),
        m_correlations_2x2(number_of_tracers_2x2 > 0 ?  84 : 0,0)
{
        //
        this->m_avg_rate_1x1.reserve(this->m_data_points);
        this->m_avg_rate_2x2.reserve(this->m_data_points);
        this->m_avg_lsquared_1x1.reserve(this->m_data_points);
        this->m_avg_lsquared_2x2.reserve(this->m_data_points);
        this->m_positions_1x1.reserve(number_of_tracers_1x1 > 0 ? this->m_number_of_tracers_1x1_times_number_of_lattices * this->m_data_points : 0);
        this->m_positions_2x2.reserve(number_of_tracers_2x2 > 0 ? this->m_number_of_tracers_2x2_times_number_of_lattices * this->m_data_points : 0);
        //
        this->m_lattices.reserve(this->m_number_of_lattices);
        this->m_tracers.reserve(this->m_number_of_tracers_total_times_number_of_lattices);
        this->m_tracers_1x1.reserve(this->m_number_of_tracers_1x1_times_number_of_lattices);
        this->m_tracers_2x2.reserve(this->m_number_of_tracers_2x2_times_number_of_lattices);
        //D( std::cout << "Wrapper: " << this->m_number_of_tracers_1x1_times_number_of_lattices << " " << this->m_number_of_tracers_2x2_times_number_of_lattices << std::endl );
        //
        while(this->m_lattices.size() < m_number_of_lattices) {
                this->m_lattices.push_back(new Lattice(this->m_grid_size,
                                                       this->m_number_of_tracers_1x1,
                                                       this->m_number_of_tracers_2x2,
                                                       this->m_step_rate_1x1,
                                                       this->m_step_rate_2x2));
        }
        for(Lattice * l : this->m_lattices)
        {
                for(Tracer * tr : l->get_tracers_1x1())
                {
                        this->m_tracers.push_back(tr);
                        this->m_tracers_1x1.push_back(tr);
                }
                for(Tracer * tr : l->get_tracers_2x2())
                {
                        this->m_tracers.push_back(tr);
                        this->m_tracers_2x2.push_back(tr);
                }
        }
        int tmp_idx = 0;
        int tmp_dsi = 0;
        int exponent = 0;
        while(tmp_dsi <= this->m_number_of_timesteps) {
                for(int n = 1; n < 10; n++) {
                        tmp_dsi = n*(int)std::pow(10,exponent);
                        if(tmp_dsi > this->m_number_of_timesteps) { break; }
                        this->m_data_storage_intervals[tmp_idx] = tmp_dsi;
                        tmp_idx++;
                }
                exponent++;
        }
        D( std::cout << "Data Storage Intervals:" << std::endl );
        D( print_vector(this->m_data_storage_intervals) );
}


inline int Wrapper::coord(int x, int y){
        //
        return ((x+this->m_grid_size)%this->m_grid_size)*this->m_grid_size+((y+this->m_grid_size)%this->m_grid_size);
}

//
void Wrapper::timestep(){
        for(Lattice * l : this->m_lattices)
        {
                l->timestep();
        }
        // update data
        this->m_t++;
        if(!(this->m_t == this->m_data_storage_intervals[this->m_next_data_storage_interval])) { return; }
        this->update_data();
        this->m_next_data_storage_interval++;
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
        int tmp_sum_steps_taken_1x1 = 0;
        int tmp_sum_steps_taken_2x2 = 0;
        double tmp_sum_msd_1x1 = 0;
        double tmp_sum_msd_2x2 = 0;
        for(Tracer * tr : this->m_tracers_1x1) {
                tmp_sum_msd_1x1 += tr->get_lsquared();
                tmp_sum_steps_taken_1x1 += tr->get_steps_taken();
                this->m_positions_1x1.push_back(this->coord(tr->get_x(),tr->get_y()));
        }
        for(Tracer * tr : this->m_tracers_2x2) {
                tmp_sum_msd_2x2 += tr->get_lsquared();
                tmp_sum_steps_taken_2x2 += tr->get_steps_taken();
                this->m_positions_2x2.push_back(this->coord(tr->get_x(),tr->get_y()));
        }
        this->m_avg_lsquared_1x1.push_back((double)this->m_t);
        this->m_avg_lsquared_1x1.push_back(tmp_sum_msd_1x1/(double)this->m_number_of_tracers_1x1_times_number_of_lattices);
        //D( std::cout << "Update Data" << std::endl );
        //D( std::cout << tmp_sum_steps_taken_1x1 << " " << (double)tmp_sum_steps_taken_1x1/this->m_number_of_tracers_1x1_times_number_of_lattices << std::endl );
        this->m_avg_rate_1x1.push_back((double)this->m_t);
        this->m_avg_rate_1x1.push_back((double)tmp_sum_steps_taken_1x1/this->m_number_of_tracers_1x1_times_number_of_lattices);

        this->m_avg_lsquared_2x2.push_back((double)this->m_t);
        this->m_avg_lsquared_2x2.push_back(tmp_sum_msd_2x2/(double)this->m_number_of_tracers_2x2_times_number_of_lattices);

        this->m_avg_rate_2x2.push_back((double)this->m_t);
        this->m_avg_rate_2x2.push_back((double)tmp_sum_steps_taken_2x2/this->m_number_of_tracers_2x2_times_number_of_lattices);
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
std::vector<unsigned long> Wrapper::get_result_correlations_1x1(){
        std::vector<unsigned long> tmp_corr(84,0);
        for(Tracer * tr : this->m_tracers_1x1) {
                std::transform(tmp_corr.begin(),tmp_corr.end(),tr->get_correlations().begin()+1,tmp_corr.begin(),std::plus<>{});
        }
        return tmp_corr;
}
std::vector<unsigned long> Wrapper::get_result_correlations_2x2(){
        std::vector<unsigned long> tmp_corr(84,0);
        for(Tracer * tr : this->m_tracers_2x2) {
                std::transform(tmp_corr.begin(),tmp_corr.end(),tr->get_correlations().begin()+1,tmp_corr.begin(),std::plus<>{});
        }
        return tmp_corr;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<int> Wrapper::get_result_positions_1x1(){
        return this->m_positions_1x1;
}
std::vector<int> Wrapper::get_result_positions_2x2(){
        return this->m_positions_2x2;
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
