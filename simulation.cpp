#include "simulation.h"
#include <cassert>
#include <numeric>
#include <algorithm>
#include <math.h>

simulation::simulation(int pop_size):
  population(pop_size,individual(0,0))
{

}

std::vector<int> simulation::get_dividing_individuals()
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

void simulation::place_start_cells() noexcept
{
  //dispose cells in a hexagonal packing patter
  std::vector<double> angles{60,120,180,240,300,360};
  for(int i = 0; i != get_pop_size(); i++)
    {
      double ref_x = population[i].get_x();
      double ref_y = population[i].get_y();
      for(int j = 0 ; j!= 6; j++)//place the next six individuals around the focal
        {
          int s_ind = i*6+1+j;
          if(s_ind > get_pop_size()){return;}

          double s_x = std::cos(angles[j])
              * (population[i].get_size()
                 + population[s_ind].get_size()
              + m_min_init_dist_btw_cells
              );
          double s_y = std::sin(angles[j])
              * (population[i].get_size()
                 + population[s_ind].get_size()
              + m_min_init_dist_btw_cells
              );
          individual s_i (ref_x + s_x, ref_y + s_y);
          if(std::find(population.begin(),population.end(), s_i) != population.end())
            {
              population[s_ind] = s_i;
            }
        }
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

  //No individuals are dividing at the start of the simulation
  {
    //initiate empty vector and leave it empty
    std::vector<individual> dividing_individuals;
    simulation s(3);
    assert(static_cast<int>(s.get_dividing_individuals().size()) < 0.00000000001);

  }

  //Individuals' centers are placed in an hexagonal pattern
  {
    simulation s(7);
    s.place_start_cells();
    //The angle formed by any 3 individual_centers is a multiple of PI/12
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

            double modulus = fmod(angle,M_PI/3);
            assert( modulus < 0.0000000001 || modulus > M_PI/3 - 0.000001);
          }
  }
  //  //No individuals are overlapping at the start of the simulation
  //  {
  //    simulation s(2);
  //    assert(!has_collision(s));
  //  }

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
    assert(s.get_pop()[div_ind[0]].get_energy() >= s.get_individual(1).get_treshold_energy());
    assert(s.get_pop()[div_ind[1]].get_energy() >= s.get_individual(2).get_treshold_energy());

  }

  //Individuals that have more than the energy treshold value are put in the dividing individuals vector

  //NOT SURE HOW TO BE MORE PRECISE, should add ID to individuals and create a == operator
  //for individual class. For now I will just check that the reproducing_individuals vector
  //has the same size as the number of reproducing individuals

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
  //and redistribute their remaining energy to their daughter cells

  //I do not like this test, i do not know how to test it otherwise though
  {
    simulation s(3);
    //This individual will not reproduce
    s.set_ind_en(0,s.get_ind_tr_en(0) - 1);
    //This individual will reproduce with no extra energy
    s.set_ind_en(1,s.get_ind_tr_en(1));
    //This ind will reproduce with extra energy
    s.set_ind_en(2,s.get_ind_tr_en(2) + 1);


    //Horrible code to find both daughter cells and check that their energy corresponds
    //To half the energy of the mother cell

    int original_pop = s.get_pop_size();
    std::vector<int> divided_ind = s.get_dividing_individuals();
    int par_off_dist = original_pop - divided_ind[0];
    std::vector<double> excess_energy;
    for(auto i=0; i != static_cast<int>(divided_ind.size()); i++)
      {
        int ind_index = divided_ind[i];
        excess_energy.push_back(
              s.get_ind_en(ind_index) - s.get_ind_tr_en(ind_index)
              );
      }

    s.do_reprduction();
    for(int i = 0; i != original_pop; i++)
      {
        if(std::count(divided_ind.begin(),divided_ind.end(),i)){

            auto en_ind = std::distance(
                  divided_ind.begin(),
                  std::find(divided_ind.begin(),divided_ind.end(),i)
                  );

            assert(s.get_ind_en(i) - excess_energy[en_ind]/2 < 0.00000001);
            assert(s.get_ind_en(i + par_off_dist) - excess_energy[en_ind]/2 < 0.00000001);

          }
      }
  }



}

















