#include "sim_parameters.h"
#include <cassert>

sim_param::sim_param(env_changer e, ind_param i, meta_param m, pop_param p):
    m_env_changer{e},
    m_ind_param{i},
    m_meta_param{m},
    m_pop_param{p}
{

}

std::ostream& operator<<(std::ostream& os, const sim_param& p)
{
    os << p.get_env_changer()
       << p.get_meta_param()
       << p.get_pop_param();
    return os;
}

bool operator==(const sim_param& lhs, const sim_param& rhs) noexcept
{
    return
            lhs.get_env_changer() == rhs.get_env_changer()
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
    f << p.get_env_changer()  << " , " <<
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
    env_changer e;

    std::string dummy; // To remove the annotation (" , ") in the file

    f >> e >> dummy;
    f >> i >> dummy;
    f >> m >> dummy;
    f >> p;

    sim_param s{e, i, m, p};
    return  s;

}

void save_sim_parameters_json(
        const sim_param& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    nlohmann::json json_out;
    json_out = p;
    f << p;
}

sim_param load_sim_parameters_json(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    nlohmann::json json_in;
    f >> json_in;
    sim_param s = json_in;
    return  s;
}

