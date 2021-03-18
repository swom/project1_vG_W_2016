#include "demographic_cycle.h"
#include<cassert>

demographic_cycle::demographic_cycle(int n_actives,
                                     int n_spores,
                                     int n_sporulating,
                                     int n_ticks,
                                     env_param env_p,
                                     ind_param ind_p):
    m_n_actives{n_actives},
    m_n_spores{n_spores},
    m_n_sporulating{n_sporulating},
    m_n_ticks{n_ticks},
    m_env_param{env_p},
    m_ind_param{ind_p}
{

}

std::ostream& operator<<(std::ostream& os, const demographic_cycle& d_c)
{
    os << d_c.get_n_actives() << " , "
       << d_c.get_n_spores() << " , "
       << d_c.get_n_sporulating() << " , "
       << d_c.get_n_ticks() << " , ";

    os << d_c.get_env_param() << " , ";
    os << d_c.get_ind_param() << std::endl;
    return os;
}

std::ifstream& operator>>(std::ifstream& is, demographic_cycle& d_c)
{
    env_param e;
    ind_param i;
    int n_spores;
    int n_sporulating;
    int n_actives;
    int n_ticks;
    std::string dummy; // To remove the annotation in the file
    is >>
            n_actives >> dummy >>
            n_spores >> dummy >>
            n_sporulating >> dummy >>
            n_ticks >> dummy;

    is >> e >> dummy;
    is >> i;
    d_c = demographic_cycle {n_actives,
            n_spores,
            n_sporulating,
            n_ticks,
            e,
            i
};
    return is;
}


bool operator==(const demographic_cycle& lhs, const demographic_cycle& rhs) noexcept
{
    return
            lhs.get_n_spores() == rhs.get_n_spores()
            && lhs.get_n_actives() == rhs.get_n_actives()
            && lhs.get_n_sporulating() == rhs.get_n_sporulating()
            && lhs.get_n_ticks() == rhs.get_n_ticks()
            && lhs.get_env_param() == rhs.get_env_param()
            && lhs.get_ind_param() == rhs.get_ind_param()
            ;
}

bool operator!=(const demographic_cycle& lhs, const demographic_cycle& rhs) noexcept
{
    return !(lhs == rhs);
}

demographic_cycle load_demographic_cycle(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    nlohmann::json json_in;
    demographic_cycle d_c;
    f >> json_in;
    return d_c = json_in;
}

void save_demographic_cycle(
        const demographic_cycle& d_c,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    nlohmann::json json_out;
    json_out = d_c;
    f << json_out;
}

