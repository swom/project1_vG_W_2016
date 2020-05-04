#include "env_param.h"

env_param::env_param(int grid_side,
                     double diff_coeff,
                     double init_food,
                     double metab_degrad_rate):
    m_diff_coeff{diff_coeff},
    m_grid_side{grid_side},
    m_init_food{init_food},
    m_metab_degradation_rate{metab_degrad_rate}
{

}
