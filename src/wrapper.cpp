#include "wrapper.hpp"

// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif

/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
Wrapper::Wrapper(int rng_seed,
                 int number_of_lattices,
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
        m_dpoints(this->m_timesteps ? 9*(int)std::log10(this->m_timesteps)+this->m_timesteps/(int)std::pow(10,(int) std::log10(this->m_timesteps)) : 0),
        m_one_over_n_lattices(1/(double)this->m_number_of_lattices)
{
        this->m_lattices.reserve(this->m_number_of_lattices);
        // this->m_tracers.reserve(this->m_number_of_tracers_total_times_number_of_lattices);
        // this->m_tracers_1x1.reserve(this->m_number_of_tracers_1x1_times_number_of_lattices);
        // this->m_tracers_2x2.reserve(this->m_number_of_tracers_2x2_times_number_of_lattices);
        //D( std::cout << "Wrapper: " << this->m_number_of_tracers_1x1_times_number_of_lattices << " " << this->m_number_of_tracers_2x2_times_number_of_lattices << std::endl );
        //
        int seed = rng_seed;
        while(this->m_lattices.size() < m_number_of_lattices) {
                // #pragma omp parallel
                this->m_lattices.push_back(new Lattice(this->m_timesteps,
                                                       this->m_timesteps_w,
                                                       this->m_grid_size,
                                                       this->m_number_of_tracers_1x1,
                                                       this->m_number_of_tracers_2x2,
                                                       this->m_step_rate_1x1,
                                                       this->m_step_rate_2x2,
                                                       rng_seed ? seed : this->m_rd()));
                seed++;
        }
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
void Wrapper::warmup(){
        D( std::cerr << "Starting Warmup..." << std::endl );
        // single thread serial
        // for(Lattice * l : this->m_lattices) { l->warmup(); }
        #pragma omp parallel for schedule(static,1)
        for(auto it = this->m_lattices.begin(); it < this->m_lattices.end(); it++) { (*it)->warmup(); }
        D( std::cerr << "Done!" << std::endl );
}
void Wrapper::evolve(){
        D( std::cerr << "Starting Evolution..." << std::endl );
        // single thread serial
        // for(Lattice * l : this->m_lattices) { l->evolve(); }
        // parallelize with omp
        #pragma omp parallel for schedule(static,1)
        for(auto it = this->m_lattices.begin(); it < this->m_lattices.end(); it++) { (*it)->evolve(); }
        D( std::cerr << "Done!" << std::endl );
}
void Wrapper::evolve_no_interaction(){
        D( std::cerr << "Starting Evolution (no_interaction)..." << std::endl );
        for(Lattice * l : this->m_lattices) { l->evolve_no_interaction(); }
        // #pragma omp parallel for
        // for(auto it = this->m_lattices.begin(); it < this->m_lattices.end(); it++) {
        //         (*it)->evolve_no_interaction();
        // }
        D( std::cerr << "Done!" << std::endl );
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<double> Wrapper::get_result_lsq_1x1(){
        D( std::cout << "  Wrapper::get_result_lsq_1x1() ");
        std::vector<double> tmp_avg_lsq(this->m_dpoints,0.0);
        for(Lattice * l : this->m_lattices) {
                std::vector<double> tmp_current = l->get_avg_lsq_1x1();
                std::transform(tmp_avg_lsq.begin(),tmp_avg_lsq.end(),tmp_current.begin(),tmp_avg_lsq.begin(),std::plus<double>());
        }
        std::transform(tmp_avg_lsq.begin(),tmp_avg_lsq.end(),tmp_avg_lsq.begin(),std::bind(std::divides<double>(), std::placeholders::_1, (double)this->m_number_of_lattices));
        D( std::cout << "Done!" << std::endl );
        return tmp_avg_lsq;
}
std::vector<double> Wrapper::get_result_lsq_2x2(){
        D( std::cout << "  Wrapper::get_result_lsq_2x2() ");
        std::vector<double> tmp_avg_lsq(this->m_dpoints,0.0);
        for(Lattice * l : this->m_lattices) {
                std::transform(tmp_avg_lsq.begin(),tmp_avg_lsq.end(),l->get_avg_lsq_2x2().begin(),tmp_avg_lsq.begin(),std::plus<double>());
        }
        std::transform(tmp_avg_lsq.begin(),tmp_avg_lsq.end(),tmp_avg_lsq.begin(),std::bind(std::divides<double>(), std::placeholders::_1, (double)this->m_number_of_lattices));
        D( std::cout << "Done!" << std::endl );
        return tmp_avg_lsq;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
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
std::vector<double> Wrapper::get_result_rate_2x2(){
        D( std::cout << "  Wrapper::get_result_rate_2x2() ");
        std::vector<double> tmp_avg_rate(this->m_dpoints,0);
        for(Lattice * l : this->m_lattices) {
                std::vector<double> tmp_current_avg_rate = l->get_avg_rate_2x2();
                std::transform(tmp_current_avg_rate.begin(),tmp_current_avg_rate.end(),tmp_avg_rate.begin(),tmp_avg_rate.begin(),std::plus<double>());
        }
        std::transform(tmp_avg_rate.begin(),tmp_avg_rate.end(),tmp_avg_rate.begin(),std::bind(std::divides<double>(), std::placeholders::_1, (double)this->m_number_of_lattices));
        D( std::cout << "Done!" << std::endl );
        return tmp_avg_rate;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<double> Wrapper::get_result_diff_1x1(){
        D( std::cout << "  Wrapper::get_result_diff_1x1() ");
        std::vector<double> tmp_avg_diff(this->m_dpoints,0.0);
        for(Lattice * l : this->m_lattices) {
                std::vector<double> tmp_current_avg_diff = l->get_avg_diff_1x1();
                std::transform(tmp_current_avg_diff.begin(),tmp_current_avg_diff.end(),tmp_avg_diff.begin(),tmp_avg_diff.begin(),std::plus<double>());
        }
        std::transform(tmp_avg_diff.begin(),tmp_avg_diff.end(),tmp_avg_diff.begin(),std::bind(std::divides<double>(), std::placeholders::_1, (double)this->m_number_of_lattices));
        D( std::cout << "Done!" << std::endl );
        return tmp_avg_diff;
}
std::vector<double> Wrapper::get_result_diff_2x2(){
        D( std::cout << "  Wrapper::get_result_diff_2x2() ");
        std::vector<double> tmp_avg_diff(this->m_dpoints,0);
        for(Lattice * l : this->m_lattices) {
                std::vector<double> tmp_current_avg_diff = l->get_avg_diff_2x2();
                std::transform(tmp_current_avg_diff.begin(),tmp_current_avg_diff.end(),tmp_avg_diff.begin(),tmp_avg_diff.begin(),std::plus<double>());
        }
        std::transform(tmp_avg_diff.begin(),tmp_avg_diff.end(),tmp_avg_diff.begin(),std::bind(std::divides<double>(), std::placeholders::_1, (double)this->m_number_of_lattices));
        D( std::cout << "Done!" << std::endl );
        return tmp_avg_diff;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<double> Wrapper::get_result_step_corr_1x1(){
        std::vector<double> tmp_corr(84,0);
        double tmp_norm = 1/(double)this->m_number_of_lattices;
        for(Lattice * l : this->m_lattices) {
                std::transform(tmp_corr.begin(),tmp_corr.end(),l->get_avg_corr_1x1().begin(),tmp_corr.begin(),std::plus<>{});
        }
        std::transform(tmp_corr.begin(),tmp_corr.end(),tmp_corr.begin(),std::bind(std::multiplies<double>(), std::placeholders::_1, tmp_norm));
        return tmp_corr;
}
std::vector<double> Wrapper::get_result_step_corr_2x2(){
        std::vector<double> tmp_corr(84,0);
        double tmp_norm = 1/(double)this->m_number_of_lattices;
        for(Lattice * l : this->m_lattices) {
                std::transform(tmp_corr.begin(),tmp_corr.end(),l->get_avg_corr_2x2().begin(),tmp_corr.begin(),std::plus<>{});
        }
        std::transform(tmp_corr.begin(),tmp_corr.end(),tmp_corr.begin(),std::bind(std::multiplies<double>(), std::placeholders::_1, tmp_norm));
        return tmp_corr;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<double> Wrapper::get_result_site_corr_1x1(){
        std::vector<double> tmp_site_corr;
        tmp_site_corr.reserve(this->m_lattices[0]->get_site_corr_1x1().size() * this->m_number_of_lattices);
        for(Lattice * l : this->m_lattices)
        {
                for(double d : l->get_site_corr_1x1())
                {
                        tmp_site_corr.push_back(d);
                }
        }
        return tmp_site_corr;
}
std::vector<double> Wrapper::get_result_site_corr_2x2(){
        std::vector<double> tmp_site_corr;
        tmp_site_corr.reserve(this->m_lattices[0]->get_site_corr_2x2().size() * this->m_number_of_lattices);
        for(Lattice * l : this->m_lattices)
        {
                for(double d : l->get_site_corr_2x2())
                {
                        tmp_site_corr.push_back(d);
                }
        }
        return tmp_site_corr;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
std::vector<int> Wrapper::get_result_pos_1x1(){
        std::vector<int> tmp_pos;
        tmp_pos.reserve(this->m_number_of_lattices * this->m_number_of_tracers_1x1 * this->m_dpoints);
        for(Lattice * l : this->m_lattices) {
                for(int pos : l->get_pos_1x1()) {
                        tmp_pos.push_back(pos);
                }
        }
        return tmp_pos;
}
std::vector<int> Wrapper::get_result_pos_2x2(){
        std::vector<int> tmp_pos;
        tmp_pos.reserve(this->m_number_of_lattices * this->m_number_of_tracers_2x2 * this->m_dpoints);
        for(Lattice * l : this->m_lattices) {
                for(int pos : l->get_pos_2x2()) {
                        tmp_pos.push_back(pos);
                }
        }
        return tmp_pos;
}
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
