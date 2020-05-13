#include "sim_parameters.h"

sim_param::sim_param(unsigned int start_pop_size, //!OCLINT
                     unsigned int exp_new_pop_size,
                     double min_dist,
                     int grid_side,
                     double diff_coeff,
                     double init_food,
                     double mutation_prob,
                     double mutation_step,
                     double base_disp_prob,
                     double spore_advantage,
                     double reproduction_prob,
                     double metab_degrad_rate,
                     int n_cycles,
                     int cycle_duration):
    m_env_param{grid_side,
                diff_coeff,
                init_food,
                metab_degrad_rate},
    m_meta_param{n_cycles,
                 cycle_duration},
    m_pop_param{start_pop_size,
                exp_new_pop_size,
                min_dist,
                mutation_prob,
                mutation_step,
                base_disp_prob,
                spore_advantage,
                reproduction_prob}
{

}

sim_param::sim_param(env_param e, meta_param m, pop_param p):
    m_env_param{e},
    m_meta_param{m},
    m_pop_param{p}
{

}

  void test_sim_param() noexcept //!OCLINT
{

    ///A simulation can be initialized by sets of parameters
    /// for population, environment and meta
    /// It suffices for this test to pass to make the code compile and execute
    {
        pop_param p;
        env_param e;
        meta_param m;
        sim_param{e, m, p};
    }

    //A sim_paramaeters can be intitialized from a file
    {
    }
}
