#include "tests.h"

void test_load_save_json()
{

    env_param e;
    ind_param ind;
    int n_spores = 2;
    int n_sporulating = 3;
    int n_actives = 4;
    int n_ticks = 5;

    demographic_cycle ex{n_actives,
                n_spores,
                n_sporulating,
                n_ticks,
                e,
                ind};

    demographic_sim p;
    int lenght_demo_sim = 3;
    for(int i = 0; i != lenght_demo_sim; i++)
    {
        p.get_demo_cycles().push_back(ex);
    }

    //Test load and save
    const std::string filename = "demographic_sim.json";
    save_demographic_sim_json(p, filename);
    const demographic_sim q = load_demographic_sim_json(filename);
    assert(p == q);
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

        env_param e;
        ind_param ind;
        int n_spores = 2;
        int n_sporulating = 3;
        int n_actives = 4;
        int n_ticks = 5;

        demographic_cycle ex{n_actives,
                    n_spores,
                    n_sporulating,
                    n_ticks,
                    e,
                    ind};

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
