#include "meta_param.h"
#include <cassert>
meta_param::meta_param( int n_cycles, int cycle_duration):
    m_cycle_duration{cycle_duration},
    m_n_cycles{n_cycles}
{
    assert(m_n_cycles > 0);
    assert(m_cycle_duration > 0);
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

}
