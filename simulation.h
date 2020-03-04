#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include"individual.h"

class simulation
{
public:
    simulation(int pop_size = 1);

    ///Gets the vector containing all the individuals of the population
    const std::vector<individual>& get_v_individuals() const noexcept {return m_vector_of_individuals;}

    ///Gets an inidividual at a certain index in the vector
    const individual& get_individual(int i) const noexcept {return m_vector_of_individuals[i];}

    ///Gets an inidividual at a certain index in the vector and allows to change it
    individual& get_individual(int i) noexcept {return m_vector_of_individuals[i];}

    ///Gets the size of the population
    int get_pop_size() const noexcept {return static_cast<int>(m_vector_of_individuals.size());}

    ///The individuals in the vector are copied at the end of the pop vector
    void individuals_divide(std::vector<individual> dividing_individuals);

private:
    std::vector<individual> m_vector_of_individuals;
};

void test_simulation();

#endif // SIMULATION_H
