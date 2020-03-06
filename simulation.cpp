#include "simulation.h"
#include <cassert>
#include <numeric>
#include <algorithm>
#include <math.h>

simulation::simulation(int pop_size, double min_dist):
    population(pop_size,individual(0,0)),
    m_min_init_dist_btw_cells(min_dist)
{
    place_start_cells();
}

std::vector<int> simulation::get_dividing_individuals() const noexcept
{

    std::vector<int> dividing_individuals;
    for(auto i = 0; i != get_pop_size(); i++)
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

std::vector<double> simulation::get_excess_energies(std::vector<int>& v_ind) const noexcept
{
    std::vector<double> excess_energy;
    for(auto i=0; i != static_cast<int>(v_ind.size()); i++)
    {
        int ind_index = v_ind[i];
        excess_energy.push_back(
                    get_ind_en(ind_index) - get_ind_tr_en(ind_index)
                    );
    }
    return excess_energy;
}

std::vector<std::pair<int,int>> simulation::get_off_indexes() const noexcept
{
    std::vector<std::pair<int,int>> off_indexes;
    int off_dist;
    auto div_ind = get_dividing_individuals();
    std::pair<int,int> daughters;
    for (int ind_index = 0; ind_index < get_pop_size(); ind_index++)
    {
        off_dist = get_pop_size() - div_ind[ind_index];
        daughters.first = ind_index;
        daughters.second = ind_index + off_dist;
        off_indexes.push_back(daughters);
    }
    return off_indexes;
}

void simulation::place_start_cells() noexcept
{
    int n = count_hex_layers(get_pop_size());
    int placed_ind = 0;

    // d is the distance between 2 individuals's centers
    float d = 2 * population[0].get_size() + m_min_init_dist_btw_cells;

    for(int i=0; i != n; i++)
    {
        double y = (sqrt(3)*i*d)/2.0;

        for (int j = 0; j < (2*n-1-i); j++)
        {
            double x = (-(2*n-i-2)*d)/2.0 + j*d;

            population[placed_ind].set_x(x);
            population[placed_ind].set_y(y);
            placed_ind++;
            if(placed_ind == get_pop_size()) return;

            if (y!=0)
            {
                population[placed_ind].set_x(x);
                population[placed_ind].set_y(-y);
                placed_ind++;
                if(placed_ind == get_pop_size()) return;

            }
        }

    }
}

bool has_collision(simulation s)
{
    const auto n_ind = s.get_pop().size();
    for (unsigned int i = 0; i < n_ind - 1; ++i)
    {
        for (unsigned int j = i + 1; j < n_ind; ++j)
        {
            if (are_colliding(s.get_individual(i), s.get_individual(j)))
            {
                return true;
            }
        }
    }
    return false;

}

int count_hex_layers(int pop_size)  noexcept
{
    int n = 1;
    if(pop_size>0){while(3*n*(n-1)+1 < pop_size) n++;}
    else {return 0;};
    return n;
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
            double x = s.get_individual(i).get_x();
            assert( x > -1234567890 );
        }
    }


    //The size of a population is equal to the size of the vector containing its individuals
    {
        simulation s;
        assert( s.get_pop_size() == static_cast<int>(s.get_pop().size()));
    }


    //No individuals are dividing at the start of the simulation
    {
        //initiate empty vector and leave it empty
        std::vector<individual> dividing_individuals;
        simulation s(3);
        assert(static_cast<int>(s.get_dividing_individuals().size()) < 0.00000000001);
    }


    // The number n of hex_layers required for x cells is
    // 3n(n-1)+1 <= x < 3(n+1)((n+1)-1)+1
    {
        std::vector<int> x{0,1,7,19,22};
        std::vector<int> n{0,1,2,3,4};

        for(auto i = 0; i < static_cast<int>(x.size()); i++ )
        {
            assert(count_hex_layers(x[i]) == n[i]);
        }
    }


    //Individuals' centers are placed in an hexagonal pattern
    {
        simulation s(7);
        //The angle formed by any 3 individual_centers is a multiple of PI/6*number_of_hex_layers
        //This is true in an hex_patt but i do not know if it is sufficient

        for(int i = 0; i != s.get_pop_size()-2; i++)
            for(int j = i+1; j != s.get_pop_size(); j++)
                for(int k = j+1; k != s.get_pop_size(); k++)
                {
                    individual P1 = s.get_individual(i);
                    individual P2 = s.get_individual(j);
                    individual P3 = s.get_individual(k);
                    double angle =
                            std::atan2(P3.get_y() - P1.get_y(), P3.get_x() - P1.get_x()) -
                            atan2(P2.get_y() - P1.get_y(), P2.get_x() - P1.get_x());
                    double modulus = abs(fmod(angle,M_PI/12));
                    assert( modulus < 0.0000000001 || modulus > M_PI/12 - 0.000001);
                }
    }

    //No individuals are overlapping at the start of the simulation
    {
        simulation s(50);
        assert(!has_collision(s));
    }

    //when an individual divides it adds a copy of itself to the population/individual vector(i.e divides)
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


    //Only individuals with energy >= than the treshold will divide
    {
        simulation s(3);
        //This ind will not reproduce
        s.get_individual(0).set_energy(s.get_individual(0).get_treshold_energy()-1);
        //This ind will reproduce with no extra energy
        s.get_individual(1).set_energy(s.get_individual(1).get_treshold_energy());
        //This ind will reproduce with extra energy
        s.get_individual(2).set_energy(s.get_individual(2).get_treshold_energy()+1);

        std::vector<int> div_ind = s.get_dividing_individuals();
        assert(div_ind.size() > 0);

        for(int i = 0; i != static_cast<int>(div_ind.size()); i++)
            assert(s.get_individual(div_ind[i]).get_energy() >= s.get_individual(div_ind[i]).get_treshold_energy());

    }

    //The excess energy of dividing individuals is equal to the diference between their energy and their treshold energy
    {
        simulation s(2);
        double energy = s.get_ind_tr_en(0);
        s.set_ind_en(0,energy);
        s.set_ind_en(1,energy*2);
        auto div_ind = s.get_dividing_individuals();
        auto v_ex_en = s.get_excess_energies(div_ind);
        for(int i =0; i != static_cast<int>(v_ex_en.size()); i++){
            assert(v_ex_en[i] - (s.get_ind_en(i) - s.get_ind_tr_en(i)) < 0.00000001);
        }
    }


    //The distance between the two offspring from same mother in population vector
    //is the distance betwen the mother index and the end of the population vector

    {
        simulation s(1);
        //First individual reproduces
        s.set_ind_en(0,s.get_ind_tr_en(0)*2);
        auto daughters = s.get_off_indexes();
        //the distance between the two offspring will be 1 element of the vector in next gen
        for(auto sisters : daughters)
        {
        assert( sisters.second - sisters.first - 1 < 0.000001);
        }
    }


    //After a reproduction round individuals with enough energy will divide
    //and redistribute their remaining energy to their daughter cells
    //I do not like this test, i do not know how to test it otherwise though
//    {
//        simulation s(3);
//        //This individual will not reproduce
//        s.set_ind_en(0,s.get_ind_tr_en(0) - 1);
//        //This individual will reproduce with no extra energy
//        s.set_ind_en(1,s.get_ind_tr_en(1));
//        //This ind will reproduce with extra energy
//        s.set_ind_en(2,s.get_ind_tr_en(2) + 1);


//        //Horrible code to find both daughter cells and check that their energy corresponds
//        //To half the energy of the mother cell



//        std::vector<double> excess_energy = s.get_excess_energies(divided_ind);

//        int original_pop_size = s.get_pop_size();
//        int par_off_dist = s.get_off_dist();
//        s.do_reprduction();
//        for(int i = 0; i != original_pop_size; i++)
//        {
//            if(std::count(divided_ind.begin(),divided_ind.end(),i)){

//                auto en_ind = std::distance(
//                            divided_ind.begin(),
//                            std::find(divided_ind.begin(),divided_ind.end(),i)
//                            );

//                assert(s.get_ind_en(i) - excess_energy[en_ind]/2 < 0.00000001);
//                assert(s.get_ind_en(i + par_off_dist) - excess_energy[en_ind]/2 < 0.00000001);

//            }
//        }
//    }



}

















