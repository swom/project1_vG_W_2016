#ifndef SIM_PARAMETERS_H
#define SIM_PARAMETERS_H
#include"env_param.h"
#include"pop_param.h"

class sim_param
{
public:
    sim_param(unsigned int start_pop_size = 1,
              unsigned int exp_new_pop_size = 1,
              double min_dist = 0.1,
              int grid_side = 1,
              double diff_coeff = 0.1,
              double init_food = 1.0,
              double mutation_prob = 0.01,
              double mutation_step = 0.1,
              double base_disp_prob = 0.01,
              double spore_advantage = 10.0,
              double reproduction_prob = 0.1,
              double metab_degrad_rate = 0.01);


    ///Gets const reference to population parameter
    const pop_param& get_pop_param() const noexcept {return m_pop_param;}

    ///Gets const reference to population parameter
    const env_param& get_env_param() const noexcept {return m_env_param;}

private:

    ///Parameters concerning the population of individuals
    pop_param m_pop_param;

    ///Parameters for the environment
    env_param m_env_param;

};

void test_sim_parameters() noexcept;

#endif // SIM_PARAMETERS_H
