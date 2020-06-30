#include "meta_param.h"
#include <cassert>


meta_param::meta_param(int n_cycles,
                       int cycle_duration,
                       int seed):
    m_cycle_duration{cycle_duration},
    m_n_cycles{n_cycles},
    m_seed{seed}
{
    assert(m_n_cycles > 0);
    assert(m_cycle_duration > 0);
}

std::ifstream& operator>> (std::ifstream& is, meta_param& p)
{
    int n_cycles;
    int cycle_duration;
    int seed;
    std::string dummy; // To remove the annotation in the file

    is
            >> n_cycles >> dummy
            >> cycle_duration >> dummy
            >> seed;

    p = meta_param {n_cycles,
            cycle_duration,
            seed};

    return is;
}

std::ostream& operator<<(std::ostream& os, const meta_param& p)
{
    os << p.get_n_cycles() << " , "
       << p.get_cycle_duration() << " , "
       << p.get_seed();

    return os;
}

bool operator==(const meta_param& lhs, const meta_param& rhs) noexcept
{
    return
            lhs.get_n_cycles() == rhs.get_n_cycles()
            && lhs.get_cycle_duration() == rhs.get_cycle_duration()
            && lhs.get_seed() == rhs.get_seed();
}

void save_meta_parameters(
        const meta_param& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    f << p;
}

meta_param load_meta_parameters(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    meta_param p;
    f >> p;

    return p;
}

void test_meta_param() noexcept
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
