#ifndef DEMOGRAPHIC_SIM_H
#define DEMOGRAPHIC_SIM_H
#include "demographic_cycle.h"
#include <vector>

class demographic_sim
{
public:
    demographic_sim();

    ///Returns const ref to vector m_demo_cycles
    const std::vector<demographic_cycle>& get_demo_cycles() const noexcept {return m_demo_cycles;}
private:
    std::vector<demographic_cycle> m_demo_cycles;
};

void test_demographic_sim();
#endif // DEMOGRAPHIC_SIM_H
