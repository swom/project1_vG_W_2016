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
                     double metab_degrad_rate):
    m_pop_param{start_pop_size,
                exp_new_pop_size,
                min_dist,
                mutation_prob,
                mutation_step,
                base_disp_prob,
                spore_advantage,
                reproduction_prob},
    m_env_param{grid_side,
                diff_coeff,
                init_food,
                metab_degrad_rate}
{

}


void test_sim_parameters() noexcept //!OCLINT
{

    //A sim_paramaeters can be intitialized from a file
    {
    }
}
