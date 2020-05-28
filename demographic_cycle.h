#ifndef DEMOGRAPHIC_CYCLE_H
#define DEMOGRAPHIC_CYCLE_H
#include "population.h"

class demographic_cycle
{
public:
    demographic_cycle(int n_actives,
                      int n_spores,
                      int n_sporulating);

    ///Returns number of spores
    int get_n_spores() const noexcept {return m_n_spores;}
    ///Returns number of sporulating
    int get_n_sporulating() const noexcept {return m_n_sporulating;}
    ///Returns number of spores
    int get_n_actives() const noexcept {return m_n_actives;}

private:

    ///number of active individuals in the pop at that moment
    int m_n_actives;
    ///number of spores in the pop at that moment
    int m_n_spores;
    ///number of sporulating individuals in the pop at that moment
    int m_n_sporulating;
};

///Returns a demographic cycle object storing data about a population
demographic_cycle demographics(const population& p) noexcept;

void test_demographic_cycle() noexcept;

#endif // DEMOGRAPHIC_CYCLE_H
