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
    demographic_cycle d_c{0,0,0,0,env_param{},ind_param{}};
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


void test_demographic_cycle() noexcept
{
    //The demograhic of a cycle can be initiated to a given value
    //And with given env and ind param
    {
        env_param e;
        ind_param i;
        int n_spores = 2;
        int n_sporulating = 3;
        int n_actives = 4;
        int n_ticks = 0;
        demographic_cycle d_c{n_actives,
                    n_spores,
                    n_sporulating,
                    n_ticks,
                    e,
                    i};
        assert(d_c.get_n_spores() == n_spores);
        assert(d_c.get_n_sporulating() == n_sporulating);
        assert(d_c.get_n_actives() == n_actives);
        assert(d_c.get_env_param() == e);
        assert(d_c.get_ind_param() == i);
    }


    //demographic_cycle object can be loaded and saved to a given file name
    {


        int n_spores = 2;
        int n_sporulating = 3;
        int n_actives = 4;
        int n_ticks = 5;
        env_param e;
        ind_param i;
        demographic_cycle p{n_actives,
                    n_spores,
                    n_sporulating,
                    n_ticks,
                    e,
                    i};

        //Test load and save
        const std::string filename = "demographic_cycle.csv";
        save_demographic_cycle(p, filename);
        const demographic_cycle q = load_demographic_cycle(filename);
        assert(p == q);
        //Test >> operator overload
        std::ifstream f(filename);

        demographic_cycle s{45,
                            46,
                            58,
                            68,
                            env_param{42,
                                      0.424,
                                      42,
                                      0.42},
                                            ind_param{42,
                                                      42,
                                                      42,
                                                      42}};
        assert(s != p);
        f >> s;
        assert(s == p);
    }

    {
        const std::string filename = "demographic_cycle2.csv";
        env_param e;
        ind_param i;
        int n_spores = 2;
        int n_sporulating = 3;
        int n_actives = 4;
        int n_ticks = 5;
        demographic_cycle p{n_actives,
                    n_spores,
                    n_sporulating,
                    n_ticks,
                    e,
                    i};
        //Make two to check that it writes them in
        // 2 different lines
        demographic_cycle p1{n_actives + 1,
                    n_spores + 1,
                    n_sporulating + 1,
                    n_ticks + 1,
                    env_param{42,
                              0.424,
                              42,
                              0.42},
                             ind_param{42,
                                       42,
                                       42,
                                       42}};
        std::ofstream s(filename);
        s << p << p1;
        const demographic_cycle q = load_demographic_cycle(filename);
        assert(p == q);
        assert(p1 != q);

    }
}
