#include "simulation.h"
#include <cassert>

simulation::simulation(int pop_size):
    m_vector_of_individuals(pop_size,individual(0,0))
{

}

void simulation::individuals_divide(std::vector<individual> dividing_individuals){

    for(auto ind : dividing_individuals){
        m_vector_of_individuals.push_back(ind);
    }
}



void test_simulation()
{
    //Simulation is initialized with a certain number of individuals
    {
        int pop_size = 100;
        simulation s(pop_size);
        // The value 1234567890 is irrelevant: just get this to compile
        for(unsigned int i = 0; i < s.get_v_individuals().size(); ++i)
        {
        assert( s.get_individual(i).get_x() > -1234567890 );
        }
    }

    //The size of a population is equal to the size of the vector containing its individuals
    {
        simulation s;
        assert( s.get_pop_size() == static_cast<int>(s.get_v_individuals().size()));
    }

    //when an individual reproduces it adds a copy of itself to the population/individual vector(i.e divides)
    {
        simulation s;
        double lhs = s.get_pop_size();
        s.individuals_divide(s.get_v_individuals());
        double rhs = s.get_pop_size();
        assert(lhs * 2 - rhs < 0.0000001);

    }

}
