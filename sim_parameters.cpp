#include "sim_parameters.h"

sim_param::sim_param(unsigned int start_pop_size,//!OCLINT
                     int exp_new_pop_size,
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
    m_base_disp_prob{base_disp_prob},
    m_diff_coeff{diff_coeff},
    m_exp_new_pop_size{exp_new_pop_size},
    m_grid_side{grid_side},
    m_init_food{init_food},
    m_min_init_dist_btw_inds{min_dist},
    m_metab_degradation_rate{metab_degrad_rate},
    m_mutation_prob{mutation_prob},
    m_mutation_step{mutation_step},
    m_repr_prob{reproduction_prob},
    m_spore_advantage{spore_advantage},
    m_start_pop_size{start_pop_size}
{

}


void test_sim_parameters() noexcept //!OCLINT
{

    //A sim_paramaeters can be intitialized from a file
    {
        sim_param s_p;
    }
}
