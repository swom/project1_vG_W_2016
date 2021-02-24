#include"tests.h"

void test_meta_param() noexcept //!OCLINT
{

    //A simulation's meta parameters object can be initialized with a number of cycles
    {
        int n_cycles = 10;
        meta_param m{n_cycles};
        assert(m.get_n_cycles() == n_cycles);
    }

    //A simulation's meta parameters object can be initialized
    //with a number of ticks per cycle
    {
        int cycle_duration = 10;
        meta_param m{1, cycle_duration};
        assert(m.get_cycle_duration() == cycle_duration);
    }
    //A simulation's meta parameters object can be initialized
    //with a seed number
    {
        int seed = 1;
        meta_param m{1, 1, seed};
        assert(m.get_seed() == seed);
    }

    //Metaparameters can be initialized with a change_frequency value
    {
        int change_frequency = 3;
        meta_param m{1,1,1, change_frequency};
        assert(m.get_change_freq() == change_frequency);
    }

    //Meta parameters can be saved and loaded correctly
    {
        int n_cycle = 34736;
        int cycle_duration = 3267;
        int seed = 3287;
        meta_param p{ n_cycle, cycle_duration, seed};
        const std::string filename = "meta_param.csv";
        save_meta_parameters(p, filename);
        const meta_param q = load_meta_parameters(filename);
        assert(p == q);
        //Test >> operator overload
        std::ifstream f(filename);
        meta_param s;
        f >> s;
        assert(s == p);
    }
    {
        const meta_param p;
        std::ostringstream s;
        s << p;
    }



}
