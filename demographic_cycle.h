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
///Prints a demographic_cycle object to an output stream
std::ostream& operator<<(std::ostream& os, const demographic_cycle& p);

///Instantiates a demographic object from a file stream (no internal checks for input correctness)
std::ifstream& operator>>(std::ifstream& is, demographic_cycle& p);

///Compares 2 demographic cycle objects to see if all their memebers are equal
bool operator==(const demographic_cycle& lhs, const demographic_cycle& rhs) noexcept;

///Returns a demographic cycle object storing data about a population
demographic_cycle demographics(const population& p) noexcept;

///Instantiate a demographic_cycle object from a given file name
demographic_cycle load_demographic_cycle(const std::string& filename);

///Saves the demographic of one cycle
void save_demographic_cycle(const demographic_cycle& p, const std::string& filename);

void test_demographic_cycle() noexcept;

#endif // DEMOGRAPHIC_CYCLE_H
