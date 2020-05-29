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

    ///Returns ref to vector m_demo_cycles
    std::vector<demographic_cycle>& get_demo_cycles()  noexcept {return m_demo_cycles;}

    ///Instantiates m_demo_cycles to a new vector
    void set_demo_cycles(const std::vector<demographic_cycle>& v_d_c) noexcept {m_demo_cycles = v_d_c;}

private:
    std::vector<demographic_cycle> m_demo_cycles;
};

///Saves the object demographic_sim to a given filename
void save_demographic_sim(const demographic_sim& p, const std::string& filename);

void test_demographic_sim();
#endif // DEMOGRAPHIC_SIM_H
