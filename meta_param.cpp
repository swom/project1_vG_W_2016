#include "meta_param.h"
#include <cassert>


meta_param::meta_param(int n_cycles,
                       int cycle_duration,
                       int seed,
                       int change_frequency,
                       int pop_max,
                       int collision_check_interval):
    m_change_frequency{change_frequency},
    m_collision_check_interval{collision_check_interval},
    m_cycle_duration{cycle_duration},
    m_n_cycles{n_cycles},
    m_pop_max{pop_max},
    m_seed{seed}
{
    assert(m_n_cycles > 0);
    assert(m_cycle_duration > 0);
    assert(m_pop_max > 0);
    assert(m_collision_check_interval >= 0);
}

std::ifstream& operator>> (std::ifstream& is, meta_param& p)
{
    int n_cycles;
    int cycle_duration;
    int seed;
    int change_frequency;
    std::string dummy; // To remove the annotation in the file

    is
            >> n_cycles >> dummy
            >> cycle_duration >> dummy
            >> seed >> dummy
            >> change_frequency;

    p = meta_param {n_cycles,
            cycle_duration,
            seed,
            change_frequency};

    return is;
}

std::ostream& operator<<(std::ostream& os, const meta_param& p)
{
    os << p.get_n_cycles() << " , "
       << p.get_cycle_duration() << " , "
       << p.get_seed() << " , "
       << p.get_change_freq();

    return os;
}

bool operator==(const meta_param& lhs, const meta_param& rhs) noexcept
{
    return
            lhs.get_n_cycles() == rhs.get_n_cycles()
            && lhs.get_cycle_duration() == rhs.get_cycle_duration()
            && lhs.get_seed() == rhs.get_seed()
            && lhs.get_change_freq() == rhs.get_change_freq();
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
