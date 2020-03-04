#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include"individual.h"

class simulation
{
public:
    simulation(int pop_size = 1);

    ///finds the individuals ready to divide
    // to be changed to vectorof references
    std::vector<individual> get_dividing_individuals();

    ///Gets an inidividual at a certain index in the vector
    const individual& get_individual(int i) const noexcept {return population[i];}

    ///Gets an inidividual at a certain index in the vector and allows to change it
    individual& get_individual(int i) noexcept {return population[i];}

    ///Gets the size of the population
    int get_pop_size() const noexcept {return static_cast<int>(population.size());}

    ///Gets the vector containing all the individuals of the population
    const std::vector<individual>& get_pop() const noexcept {return population;}

    ///The individuals in the vector are copied at the end of the pop vector
    // To be chaged to take vec of references as an argument
    void individuals_divide(std::vector<individual> dividing_individuals);

private:
    std::vector<individual> population;
};

void test_simulation();

#endif // SIMULATION_H
