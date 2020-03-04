#include "simulation.h"
#include <cassert>

simulation::simulation(int pop_size):
    population(pop_size,individual(0,0))
{

}

std::vector<individual> simulation::get_dividing_individuals(){
    std::vector<individual> dividing_individuals;

    for(auto ind : population){
        if(ind.get_energy()>= ind.get_treshold_energy()){dividing_individuals.push_back(ind);}
    }
    return dividing_individuals;
}

void simulation::individuals_divide(std::vector<individual> dividing_individuals){

    for(auto ind : dividing_individuals){
        population.push_back(ind);
    }
}




void test_simulation()
{
    //Simulation is initialized with a certain number of individuals
    {
        int pop_size = 100;
        simulation s(pop_size);
        // The value 1234567890 is irrelevant: just get this to compile
        for(unsigned int i = 0; i < s.get_pop().size(); ++i)
        {
        assert( s.get_individual(i).get_x() > -1234567890 );
        }
    }

    //The size of a population is equal to the size of the vector containing its individuals
    {
        simulation s;
        assert( s.get_pop_size() == static_cast<int>(s.get_pop().size()));
    }

    //when an individual reproduces it adds a copy of itself to the population/individual vector(i.e divides)
    {
        simulation s;
        double lhs = s.get_pop_size();
        //let's allow all individuals in the population to reproduce, to facilitate the testing conditions
        s.individuals_divide(s.get_pop());
        double rhs = s.get_pop_size();
        assert(lhs * 2 - rhs < 0.0000001);

    }

    //Individuals that have more than the energy treshold value are put in the dividing individuals vector

    //NOT SURE HOW TO BE MORE PRECISE, should add ID to individuals and create a == operator for individual class
    //for now I will just check that the reproducing_individuals vector has the same size as the number of reproducing individuals

    //No individuals are dividing at the start of the simulation
    {
      std::vector<individual> dividing_individuals;
      std::vector<int> future_dividing_ind_indexes{0,1};
      simulation s(3);

      int lhs = dividing_individuals.size();
      assert(static_cast<int>(s.get_dividing_individuals().size())  - lhs < 0.00000000001);

    }
}
