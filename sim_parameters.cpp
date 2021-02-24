#include "sim_parameters.h"
#include <cassert>

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

sim_param::sim_param(env_param e, ind_param i, meta_param m, pop_param p):
    m_env_param{e},
    m_ind_param{i},
    m_meta_param{m},
    m_pop_param{p}
{

}

std::ostream& operator<<(std::ostream& os, const sim_param& p)
{
    os << p.get_env_param()
       << p.get_meta_param()
       << p.get_pop_param();
    return os;
}

bool operator==(const sim_param& lhs, const sim_param& rhs) noexcept
{
    return
            lhs.get_env_param() == rhs.get_env_param()
            && lhs.get_meta_param() == rhs.get_meta_param()
            &&lhs.get_pop_param() == lhs.get_pop_param()
            ;

}

void save_sim_parameters(
        const sim_param& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    f << p.get_env_param()  << " , " <<
         p.get_ind_param() << " , " <<
         p.get_meta_param() << " , " <<
         p.get_pop_param();
}

sim_param load_sim_parameters(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    if(!f.is_open())
    {
        std::cout << "Could not find specified sim_par.csv file. \n";
        abort();
    }

    pop_param p;
    meta_param m;
    ind_param i;
    env_param e;

    std::string dummy; // To remove the annotation (" , ") in the file

    f >> e >> dummy;
    f >> i >> dummy;
    f >> m >> dummy;
    f >> p;

    sim_param s{e, i, m, p};
    return  s;

}

