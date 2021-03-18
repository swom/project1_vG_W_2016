#include "test_demographic_cycle.h"

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
    }
}
