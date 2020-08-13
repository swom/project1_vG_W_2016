#include "demographic_cycle.h"
#include<cassert>

demographic_cycle::demographic_cycle(int n_actives,
                                     int n_spores,
                                     int n_sporulating,
                                     env_param env_p,
                                     ind_param ind_p):
    m_n_actives{n_actives},
    m_n_spores{n_spores},
    m_n_sporulating{n_sporulating},
    m_env_param{env_p},
    m_ind_param{ind_p}
{

}

std::ostream& operator<<(std::ostream& os, const demographic_cycle& d_c)
{
    os << d_c.get_n_actives() << " , "
       << d_c.get_n_spores() << " , "
       << d_c.get_n_sporulating() << " , ";
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
    std::string dummy; // To remove the annotation in the file
    is >>
            n_actives >> dummy >>
            n_spores >> dummy >>
            n_sporulating >> dummy
            ;
    is >> e >> dummy;
    is >> i;
    d_c = demographic_cycle {n_actives,
            n_spores,
            n_sporulating,
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
            && lhs.get_env_param() == rhs.get_env_param()
            && lhs.get_ind_param() == rhs.get_ind_param()
            ;
}

bool operator!=(const demographic_cycle& lhs, const demographic_cycle& rhs) noexcept
{
    return !(lhs == rhs);
}

int count_actives(const population& pop)
{
    const auto& p = pop.get_v_ind();
    return static_cast<int>(std::count_if(p.begin(),p.end(), [](const individual& ind) { 
        return ind.get_phen() == phenotype::active;}
    ));
}

int count_spores(const population& pop)
{
    const auto& p = pop.get_v_ind();
    return static_cast<int>(std::count_if(p.begin(),p.end(), [](const individual& ind) { 
        return ind.get_phen() == phenotype::spore;}
    ));
}

int count_sporulating(const population& pop)
{
    const auto& p = pop.get_v_ind();
    return static_cast<int>(std::count_if(p.begin(),p.end(), [](const individual& ind) { 
        return ind.get_phen() == phenotype::sporulating;}
    ));
}

demographic_cycle demographics(const population &p, const env_param &e) noexcept
{
    return demographic_cycle{count_actives(p),
                count_spores(p),
                count_sporulating(p),
                e,
                p.get_param().get_ind_param()};
}


demographic_cycle load_demographic_cycle(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    demographic_cycle d_c{0,0,0,env_param{},ind_param{}};
    f >> d_c;
    return d_c;
}

void save_demographic_cycle(
        const demographic_cycle& d_c,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    f << d_c;
}

