#include "simulation.h"
#include <cassert>
#include <numeric>

simulation::simulation(int pop_size):
    population(pop_size,individual(0,0))
{

}

std::vector<int> simulation::get_dividing_individuals()
{

    std::vector<int> dividing_individuals;
    for(int i = 0; i != get_pop_size(); i++)
    {
        if(population[i].get_energy()>= population[i].get_treshold_energy())
        {
            dividing_individuals.push_back(i);
        }
    }
    return dividing_individuals;
}

void simulation::division(std::vector<int> dividing_individuals)
{

    for(size_t i = 0; i != dividing_individuals.size(); i++)
    {
        int div_ind = dividing_individuals[i];
        double offspring_initial_energy = population[div_ind].split_excess_energy();
        population[div_ind].set_energy(offspring_initial_energy);
        population.push_back(population[div_ind]);
    }
}

void simulation::do_reprduction() noexcept
{
    division(get_dividing_individuals());
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
        std::vector<int> dividing_pop(s.get_pop_size()) ;
        std::iota (std::begin(dividing_pop), std::end(dividing_pop), 0);
        s.division(dividing_pop);

        double rhs = s.get_pop_size();
        assert(lhs * 2 - rhs < 0.0000001);

    }

    //No individuals are dividing at the start of the simulation
    {
      //initiate empty vector and leave it empty
      std::vector<individual> dividing_individuals;
      simulation s(3);
      assert(static_cast<int>(s.get_dividing_individuals().size()) < 0.00000000001);

    }

    //Only individuals with energy >= than the treshold will divide
    {
        simulation s(3);
        s.get_individual(0).set_energy(s.get_individual(0).get_treshold_energy()-1);
        s.get_individual(1).set_energy(s.get_individual(1).get_treshold_energy());
        s.get_individual(2).set_energy(s.get_individual(2).get_treshold_energy()+1);

        std::vector<int> div_ind = s.get_dividing_individuals();

        assert(div_ind.size() > 0);
        assert(s.get_pop()[div_ind[0]].get_energy() >= s.get_individual(1).get_treshold_energy());
        assert(s.get_pop()[div_ind[1]].get_energy() >= s.get_individual(2).get_treshold_energy());

    }

    //Individuals that have more than the energy treshold value are put in the dividing individuals vector
    //NOT SURE HOW TO BE MORE PRECISE, should add ID to individuals and create a == operator for individual class
    //for now I will just check that the reproducing_individuals vector has the same size as the number of reproducing individuals

    //I do not like this test, i do not know how to test it otherwise though
    {
        std::vector<int> dividing_individuals;
        simulation s(3);
        //At the beginning no individual reproduce
        assert(static_cast<int>(s.get_dividing_individuals().size()) < 0.00000000001);

        for (int i = 0; i != s.get_pop_size(); i++)
        {
            double en_tresh_ind = s.get_individual(i).get_treshold_energy();
            s.get_individual(i).set_energy(en_tresh_ind);
            dividing_individuals = s.get_dividing_individuals();
            assert(static_cast<int>(dividing_individuals.size()) == i+1);
        }
    }


    //After a reproduction round individuals with enough energy will divide
    //(already tested)
    //and redistribute their remaining energy to their daughter cells

    //I do not like this test, i do not know how to test it otherwise though
    {
        simulation s(2);
        //This individual will not reproduce
        s.get_individual(0).set_energy(s.get_individual(0).get_treshold_energy()-1);
        //This individual will reproduce
        s.get_individual(1).set_energy(s.get_individual(1).get_treshold_energy());

        int original_pop = s.get_pop_size();
        int new_ind = s.get_dividing_individuals().size();
        double excess_energy_of_reproducing_ind = s.get_individual(1).get_energy() - s.get_individual(1).get_treshold_energy();

        s.do_reprduction();
        assert(s.get_pop_size() - (original_pop + new_ind) < 0.000000001);
        assert(s.get_individual(1).get_energy() - excess_energy_of_reproducing_ind/2 < 0.00000001);
        assert(s.get_individual(2).get_energy() - excess_energy_of_reproducing_ind/2 < 0.00000001);

    }

}

















