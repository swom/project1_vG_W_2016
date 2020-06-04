#include "demographic_cycle.h"
#include<cassert>

demographic_cycle::demographic_cycle(int n_actives,
                                     int n_spores,
                                     int n_sporulating):
    m_n_actives{n_actives},
    m_n_spores{n_spores},
    m_n_sporulating{n_sporulating}
{

}

std::ostream& operator<<(std::ostream& os, const demographic_cycle& d_c)
{
    os << d_c.get_n_actives() << " , "
       << d_c.get_n_spores() << " , "
       << d_c.get_n_sporulating() << std::endl;
    ;
    return os;
}

std::ifstream& operator>>(std::ifstream& is, demographic_cycle& d_c)
{

    int n_spores;
    int n_sporulating;
    int n_actives;
    std::string dummy; // To remove the annotation in the file
    is >>
            n_actives >> dummy >>
            n_spores >> dummy >>
            n_sporulating;

    d_c = demographic_cycle {n_actives,
            n_spores,
            n_sporulating
};
    return is;
}


bool operator==(const demographic_cycle& lhs, const demographic_cycle& rhs) noexcept
{
    return
            lhs.get_n_spores() == rhs.get_n_spores()
            && lhs.get_n_actives() == rhs.get_n_actives()
            && lhs.get_n_sporulating() == rhs.get_n_sporulating()
            ;
}

demographic_cycle demographics(const population& p) noexcept
{
    int n_spores = std::count_if(p.get_v_ind().begin(),p.get_v_ind().end(),
                                 [](const individual& ind){ return ind.get_phen() == phenotype::spore;});
    int n_sporulating = std::count_if(p.get_v_ind().begin(),p.get_v_ind().end(),
                                      [](const individual& ind){ return ind.get_phen() == phenotype::sporulating;});
    int n_actives = std::count_if(p.get_v_ind().begin(),p.get_v_ind().end(),
                                  [](const individual& ind){ return ind.get_phen() == phenotype::active;});

    return demographic_cycle{n_actives, n_spores, n_sporulating};
}


demographic_cycle load_demographic_cycle(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    demographic_cycle d_c{0,0,0};
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
    {
        int n_spores = 2;
        int n_sporulating = 3;
        int n_actives = 4;
        demographic_cycle d_c{n_actives, n_spores, n_sporulating};
        assert(d_c.get_n_spores() == n_spores);
        assert(d_c.get_n_sporulating() == n_sporulating);
        assert(d_c.get_n_actives() == n_actives);
    }

    //It is possible to extract the demographic state of a population
    {
        population p{0};
        assert(p.get_pop_size() == 0);

        int n_spores = 2;
        int n_sporulating = 3;
        int n_actives = 4;
        individual ind;

        ind.set_phen(phenotype::spore);
        for(int i = 0; i != n_spores; i++)
        {
            p.get_v_ind().push_back(ind);
        }

        ind.set_phen(phenotype::sporulating);
        for(int i = 0; i != n_sporulating; i++)
        {
            p.get_v_ind().push_back(ind);
        }

        ind.set_phen(phenotype::active);
        for(int i = 0; i != n_actives; i++)
        {
            p.get_v_ind().push_back(ind);
        }

        assert(p.get_pop_size() == n_spores + n_sporulating + n_actives);

        auto d_c = demographics(p);

        assert(d_c.get_n_spores() == n_spores);
        assert(d_c.get_n_sporulating() == n_sporulating);
        assert(d_c.get_n_actives() == n_actives);
    }

    //demographic_cycle object can be loaded and saved to a given file name
    {


        int n_spores = 2;
        int n_sporulating = 3;
        int n_actives = 4;
        demographic_cycle p{n_actives,
                    n_spores,
                    n_sporulating};

        //Test load and save
        const std::string filename = "demographic_cycle.csv";
        save_demographic_cycle(p, filename);
        const demographic_cycle q = load_demographic_cycle(filename);
        assert(p == q);
        //Test >> operator overload
        std::ifstream f(filename);
        demographic_cycle s{0,0,0};
        f >> s;
        assert(s == p);
    }
    {
        const std::string filename = "demographic_cycle2.csv";
        int n_spores = 2;
        int n_sporulating = 3;
        int n_actives = 4;
        demographic_cycle p{n_actives,
                    n_spores,
                    n_sporulating};
        //Make two to check that it writes them in
        // 2 different lines
        demographic_cycle p1{n_actives + 1,
                    n_spores + 1,
                    n_sporulating + 1};
        std::ofstream s(filename);
        s << p << p1;
        const demographic_cycle q = load_demographic_cycle(filename);
        assert(p == q);

    }
}
