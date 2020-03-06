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
    std::vector<int> get_dividing_individuals();

    ///Gets an inidividual at a certain index in the vector
    const individual get_individual(int i) const {return population[i];}

    ///Gets an inidividual at a certain index in the vector and allows to change it
    individual& get_individual(int i) {return population[i];}

    ///Gets an individual's energy
    double get_ind_en(int i) {return population[i].get_energy();}

    ///Gets an individual's treshold energy
    double get_ind_tr_en(int i) {return population[i].get_treshold_energy();}

    ///Gets the size of the population
    int get_pop_size() const noexcept {return static_cast<int>(population.size());}

    ///Gets the vector containing all the individuals of the population
    const std::vector<individual>& get_pop() const noexcept {return population;}

    ///The individuals in the vector are copied at the end of the pop vector
    // To be chaged to take vec of references as an argument
    void division(std::vector<int> dividing_individuals);

    ///Individuals with enough energy divide and redistribute energy
    void do_reprduction() noexcept;

    ///Places starting cells in hexagonal packing pattern
    void place_start_cells() noexcept;

    ///Sets and individual's energy
    void set_ind_en(int i, double en) {population[i].set_energy(en);}


private:
    std::vector<individual> population;
    double m_min_init_dist_btw_cells;
};

bool has_collision(simulation s);
void test_simulation();

#endif // SIMULATION_H
