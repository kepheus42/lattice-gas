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
                 int number_of_timesteps_w,
                 int number_of_tracers_1x1,
                 int number_of_tracers_2x2,
                 double step_rate_1x1,
                 double step_rate_2x2) :
        // sim parameters
        m_number_of_lattices(number_of_lattices),
        m_grid_size(grid_size),
        m_timesteps(number_of_timesteps),
        m_timesteps_w(number_of_timesteps_w),
        m_number_of_tracers_1x1(number_of_tracers_1x1),
        m_number_of_tracers_2x2(number_of_tracers_2x2),
        m_step_rate_1x1(step_rate_1x1),
        m_step_rate_2x2(step_rate_2x2),
        m_dpoints(this->m_timesteps ? 9*(int)std::log10(this->m_timesteps)+this->m_timesteps/(int)std::pow(10,(int) std::log10(this->m_timesteps)) : 0) //
{
        this->m_lattices.reserve(this->m_number_of_lattices);
        // this->m_tracers.reserve(this->m_number_of_tracers_total_times_number_of_lattices);
        // this->m_tracers_1x1.reserve(this->m_number_of_tracers_1x1_times_number_of_lattices);
        // this->m_tracers_2x2.reserve(this->m_number_of_tracers_2x2_times_number_of_lattices);
        //D( std::cout << "Wrapper: " << this->m_number_of_tracers_1x1_times_number_of_lattices << " " << this->m_number_of_tracers_2x2_times_number_of_lattices << std::endl );
        //
        while(this->m_lattices.size() < m_number_of_lattices) {
                this->m_lattices.push_back(new Lattice(this->m_timesteps,
                                                       this->m_timesteps_w,
                                                       this->m_grid_size,
                                                       this->m_number_of_tracers_1x1,
                                                       this->m_number_of_tracers_2x2,
                                                       this->m_step_rate_1x1,
                                                       this->m_step_rate_2x2));
        }
}

void Wrapper::warmup(){
        D( std::cerr << "Starting Warmup..." << std::endl );
        for(Lattice * l : this->m_lattices) { l->warmup(); }
        D( std::cerr << "Done!" << std::endl );
}

void Wrapper::evolve(){
        D( std::cerr << "Starting Evolution..." << std::endl );
        for(Lattice * l : this->m_lattices) { l->evolve(); }
        D( std::cerr << "Done!" << std::endl );
}

void Wrapper::evolve_no_interaction(){
        D( std::cerr << "Starting Evolution (no_interaction)..." << std::endl );
        for(Lattice * l : this->m_lattices) { l->evolve_no_interaction(); }
        D( std::cerr << "Done!" << std::endl );
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
// Get Simulation Results:
// = = = = = = = = = = = = = = = = = = = = = = = = =
// avg_lsq(t) for t = [0,...,m_timesteps]
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> Wrapper::get_result_lsq_1x1(){
        D( std::cout << "  Wrapper::get_result_lsq_1x1() ");
        std::vector<double> tmp_avg_lsq(this->m_dpoints,0.0);
        for(Lattice * l : this->m_lattices) {
                std::vector<double> tmp_current_avg_lsq = l->get_avg_lsq_1x1();
                std::transform(tmp_current_avg_lsq.begin(),tmp_current_avg_lsq.end(),tmp_avg_lsq.begin(),tmp_avg_lsq.begin(),std::plus<>());
        }
        std::transform(tmp_avg_lsq.begin(),tmp_avg_lsq.end(),tmp_avg_lsq.begin(),std::bind(std::divides<double>(), std::placeholders::_1, (double)this->m_number_of_lattices));
        D( std::cout << "Done!" << std::endl );
        return tmp_avg_lsq;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> Wrapper::get_result_lsq_2x2(){
        D( std::cout << "  Wrapper::get_result_lsq_2x2() ");
        std::vector<double> tmp_avg_lsq(this->m_dpoints,0.0);
        for(Lattice * l : this->m_lattices) {
                std::vector<double> tmp_current_avg_lsq = l->get_avg_lsq_2x2();
                std::transform(tmp_current_avg_lsq.begin(),tmp_current_avg_lsq.end(),tmp_avg_lsq.begin(),tmp_avg_lsq.begin(),std::plus<double>());
        }
        std::transform(tmp_avg_lsq.begin(),tmp_avg_lsq.end(),tmp_avg_lsq.begin(),std::bind(std::divides<double>(), std::placeholders::_1, (double)this->m_number_of_lattices));
        D( std::cout << "Done!" << std::endl );
        return tmp_avg_lsq;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
// avg_rate(t) for t = [0,...,m_timesteps]
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> Wrapper::get_result_rate_1x1(){
        D( std::cout << "  Wrapper::get_result_rate_1x1() ");
        std::vector<double> tmp_avg_rate(this->m_dpoints,0.0);
        for(Lattice * l : this->m_lattices) {
                std::vector<double> tmp_current_avg_rate = l->get_avg_rate_1x1();
                std::transform(tmp_current_avg_rate.begin(),tmp_current_avg_rate.end(),tmp_avg_rate.begin(),tmp_avg_rate.begin(),std::plus<double>());
        }
        std::transform(tmp_avg_rate.begin(),tmp_avg_rate.end(),tmp_avg_rate.begin(),std::bind(std::divides<double>(), std::placeholders::_1, (double)this->m_number_of_lattices));
        D( std::cout << "Done!" << std::endl );
        return tmp_avg_rate;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
std::vector<double> Wrapper::get_result_rate_2x2(){
        D( std::cout << "  Wrapper::get_result_rate_2x2() ");
        std::vector<double> tmp_avg_rate(this->m_dpoints,0.0);
        for(Lattice * l : this->m_lattices) {
                std::vector<double> tmp_current_avg_rate = l->get_avg_rate_2x2();
                std::transform(tmp_current_avg_rate.begin(),tmp_current_avg_rate.end(),tmp_avg_rate.begin(),tmp_avg_rate.begin(),std::plus<double>());
        }
        std::transform(tmp_avg_rate.begin(),tmp_avg_rate.end(),tmp_avg_rate.begin(),std::bind(std::divides<double>(), std::placeholders::_1, (double)this->m_number_of_lattices));
        D( std::cout << "Done!" << std::endl );
        return tmp_avg_rate;
}
// = = = = = = = = = = = = = = = = = = = = = = = = =
/*

   std::vector<unsigned long long> Wrapper::get_result_correlations_1x1(){
        std::vector<unsigned long> tmp_corr(84,0);
        for(Tracer * tr : this->m_tracers_1x1) {
                std::transform(tmp_corr.begin(),tmp_corr.end(),tr->get_correlations().begin()+1,tmp_corr.begin(),std::plus<>{});
        }
        return tmp_corr;
   }
   std::vector<unsigned long long> Wrapper::get_result_correlations_2x2(){
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
 */
