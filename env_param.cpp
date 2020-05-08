#include "env_param.h"
#include <cassert>

env_param::env_param(int grid_side,
                     double diff_coeff,
                     double init_food,
                     double metab_degrad_rate):
    m_diff_coeff{diff_coeff},
    m_grid_side{grid_side},
    m_init_food{init_food},
    m_metab_degradation_rate{metab_degrad_rate}
{

    assert(m_diff_coeff > -0.000000000001 &&
           m_diff_coeff < 1.000000000001);
    assert(m_grid_side >= 0);
    assert(m_init_food > -0.000000000001);
    assert(m_metab_degradation_rate > -0.000000001 &&
           m_metab_degradation_rate < 1.000000001);
}

void test_env_param() noexcept //!OCLINT

{

}
