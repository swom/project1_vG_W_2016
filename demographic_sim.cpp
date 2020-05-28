#include "demographic_sim.h"
#include<cassert>

demographic_sim::demographic_sim()
{

}

void test_demographic_sim()
{
    //demographic_sim has contains a vector
    //of demographic_cycle objects
    {
        demographic_sim d_s;
        assert(d_s.get_demo_cycles().size() >= 0);
    }
}
