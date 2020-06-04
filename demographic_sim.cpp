#include "demographic_sim.h"
#include<cassert>

demographic_sim::demographic_sim()
{

}

bool operator==(const demographic_sim& lhs, const demographic_sim& rhs) noexcept
{
    return
            std::equal(lhs.get_demo_cycles().begin(),
                       lhs.get_demo_cycles().end(),
                       rhs.get_demo_cycles().begin());
}

demographic_sim load_demographic_sim(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    demographic_sim d_s;
    demographic_cycle d_c{0,0,0};
    std::string dummy;
    while(f >> dummy)//Skips the cycle number
    {
        f >> dummy; //skips the comma
        f >> d_c;
        d_s.get_demo_cycles().push_back(d_c);
    }

    return d_s;
}

void save_demographic_sim(
        const demographic_sim& d_s,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    int cycle_number = 0;
    for(const auto& demo_cycle : d_s.get_demo_cycles())
    {
        f << cycle_number << " , " << demo_cycle;
        cycle_number++;
    }
}

void test_demographic_sim()//!OCLINT
{
    //demographic_sim has contains a vector
    //of demographic_cycle objects
    {
        demographic_sim d_s;
        assert(d_s.get_demo_cycles().size() >= 0);
    }

    //demographic_sim object can be loaded and saved to a given file name
    {


        int n_spores = 2;
        int n_sporulating = 3;
        int n_actives = 4;
        demographic_cycle ex{n_actives,
                    n_spores,
                    n_sporulating};

        demographic_sim p;
        int lenght_demo_sim = 3;
        for(int i = 0; i != lenght_demo_sim; i++)
        {
            p.get_demo_cycles().push_back(ex);
        }

        //Test load and save
        const std::string filename = "demographic_sim.csv";
        save_demographic_sim(p, filename);
        const demographic_sim q = load_demographic_sim(filename);
        assert(p == q);
    }

}
