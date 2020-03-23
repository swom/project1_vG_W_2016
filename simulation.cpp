#include "simulation.h"
#include <cassert>
#include <numeric>
#include <algorithm>
#include <math.h>

simulation::simulation(int pop_size,
                       int grid_side,
                       double min_dist,
                       double starting_food):
  population(static_cast<unsigned int>(pop_size),individual(0,0)),
  m_min_init_dist_btw_cells(min_dist),
  m_e(grid_side,0.1,starting_food)
{
  place_start_cells();
}

std::vector<int> simulation::get_dividing_individuals() const noexcept
{

  std::vector<int> dividing_individuals;
  for(unsigned int i = 0; i != get_pop().size(); i++)
    {
      if(population[i].get_energy() >= population[i].get_treshold_energy()
         && population[i].get_type() == individual_type::living)
        {
          dividing_individuals.push_back(static_cast<int>(i));
        }
    }
  return dividing_individuals;
}
void simulation::divide_ind(individual& i)
{
  double offs_init_en = i.split_excess_energy();
  i.set_energy(offs_init_en);
  individual d2 = i;
  population.push_back(d2);
  set_ind_pos(population.back(),get_d2_pos(d2));
}

void simulation::do_division(std::vector<int> dividing_individuals)
{

  for(size_t i = 0; i != dividing_individuals.size(); i++)
    {
      int div_ind = dividing_individuals[i];
      divide_ind(get_ind(div_ind));
    }
}

void simulation::do_reprduction() noexcept
{
  do_division(get_dividing_individuals());
}

void feeding(simulation& s)
{
  for(auto& ind : s.get_pop())
    {
      auto index_grid = find_grid_index(ind,s.get_env().get_grid_side());
      if(index_grid == -100 ||
         ind.get_type() != individual_type::living)
        {continue;}
      feed(ind,s.get_env().get_cell(index_grid));
    }
}

void exec(simulation& s, int n_tick) noexcept
{
  while(s.get_tick() != n_tick){tick(s);}
}

std::vector<double> simulation::get_excess_energies(std::vector<int> v_ind) const noexcept
{
  std::vector<double> excess_energy;
  for(unsigned int i=0; i != v_ind.size(); i++)
    {
      int ind_index = v_ind[i];
      excess_energy.push_back(
            get_ind_en(ind_index) - get_ind_tr_en(ind_index)
            );
    }
  return excess_energy;
}

std::vector<double> simulation::get_excess_energies(std::vector<int>& v_ind) const noexcept
{
  std::vector<double> excess_energy;
  for(unsigned int i=0; i != v_ind.size(); i++)
    {
      int ind_index = v_ind[i];
      excess_energy.push_back(
            get_ind_en(ind_index) - get_ind_tr_en(ind_index)
            );
    }
  return excess_energy;
}

const std::pair<double,double> simulation::get_ind_pos(int i)
{
  std::pair<double, double> pos = get_pos(get_ind(i));
  return pos;
}

std::vector<std::pair<int,int>> simulation::get_sisters_index_offset() const noexcept
{
  std::vector<std::pair<int,int>> sis_indexes;
  std::pair<int,int> daughters;
  int sis_dist;
  auto div_ind = get_dividing_individuals();

  for (int ind_index = 0; ind_index < static_cast<int>(div_ind.size()); ind_index++)
    {
      sis_dist = get_pop_size() - div_ind[static_cast<unsigned int>(ind_index)];
      daughters.first = ind_index;
      daughters.second = ind_index + sis_dist;
      sis_indexes.push_back(daughters);
    }
  return sis_indexes;
}


void metabolism_pop(simulation& s)
{
  for(auto& ind : s.get_pop())
    {
      if(ind.get_type() != individual_type::spore)
      metabolism(ind);
    }
}
void simulation::set_ind_pos(individual& i, double x, double y)
{
  i.set_x(x);
  i.set_y(y);
}

void simulation::set_ind_pos(individual& i, const std::pair<double, double> pos)
{
  set_pos(i, pos);
}


void tick(simulation& s)
{
  feeding(s);
  metabolism_pop(s);
  s.do_reprduction();
  manage_static_collisions(s);
  diffusion(s.get_env());
  s.tick_up();
}
void simulation::place_start_cells() noexcept
{
  int n = count_hex_layers(get_pop_size());
  unsigned int placed_ind = 0;

  // d is the distance between 2 individuals's centers
  double d = 2 * (population[0].get_size() + m_min_init_dist_btw_cells);

  for(int i=0; i != n; i++)
    {
      double y = (sqrt(3) * i * d) / 2.0;

      for (int j = 0; j < (2 * n - 1 - i); j++)
        {
          double x = (-(2 * n - i - 2) * d) / 2.0 + j * d;

          population[placed_ind].set_x(x);
          population[placed_ind].set_y(y);
          placed_ind++;
          if(placed_ind == get_pop().size()) return;

          if (y > 0.000001 || y < -0.000001)
            {
              population[placed_ind].set_x(x);
              population[placed_ind].set_y(-y);
              placed_ind++;
              if(placed_ind == get_pop().size()) return;

            }
        }

    }
}

bool has_collision(const simulation& s)
{
  int n_ind = s.get_pop_size();
  for ( int i = 0; i < n_ind; ++i)
    {
      for ( int j = i ; j < n_ind; ++j)
        {
          if(i == j)
            {
              continue;
            }
          else if (are_colliding(s.get_ind(i), s.get_ind(j)))
            {
              return true;
            }
        }
    }
  return false;
}

void manage_static_collisions(simulation& s)
{
  while(has_collision(s))
    {
      const auto n_ind = s.get_pop_size();
      for ( int i = 0; i < n_ind; ++i)
        {
          for ( int j = 0 ; j < n_ind; ++j)
            {
              if (i == j)
                {
                  continue;
                }
              else if(are_colliding(s.get_ind(i), s.get_ind(j)))
                {
                  displace(s.get_ind(i), s.get_ind(j));
                }
            }
        }
    }
}

int count_hex_layers(int pop_size)  noexcept
{
  int n = 1;
  if(pop_size>0){while(3*n*(n-1)+1 < pop_size) n++;}
  else {return 0;}
  return n;
}

std::vector<double> modulus_of_btw_ind_angles(simulation& s, double ang_rad)
{
  std::vector<double> v_modulus;
  for(int i = 0; i != s.get_pop_size()-2; i++)
    for(int j = i+1; j != s.get_pop_size(); j++)
      for(int k = j+1; k != s.get_pop_size(); k++)
        {
          individual P1 = s.get_ind(i);
          individual P2 = s.get_ind(j);
          individual P3 = s.get_ind(k);
          double angle =
              std::atan2(P3.get_y() - P1.get_y(), P3.get_x() - P1.get_x()) -
              atan2(P2.get_y() - P1.get_y(), P2.get_x() - P1.get_x());
          v_modulus.push_back((abs(fmod(angle,ang_rad))));
        }
  return v_modulus;
}

void test_simulation()//!OCLINT tests may be many
{
  //Simulation is initialized with a certain number of individuals
  // The value 1234567890 is irrelevant: just get this to compile

  {
    int pop_size = 100;
    simulation s(pop_size);
    for(int i = 0; i < s.get_pop_size(); ++i)
      {
        double x = s.get_ind(i).get_x();
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
    simulation s(3);
    assert(static_cast<int>(s.get_dividing_individuals().size()) < 0.00000000001);
  }


  // The number n of hex_layers required for x individual is
  // 3n(n-1)+1 <= x < 3(n+1)((n+1)-1)+1
  {
    std::vector<int> x{0,1,7,19,22};
    std::vector<int> n{0,1,2,3,4};

    for(unsigned int i = 0; i < x.size(); i++ )
      {
        assert(count_hex_layers(x[i]) == n[i]);
      }
  }

  //Individuals' centers are placed in an hexagonal pattern
  {
    simulation s(7);

    //The angle formed by any 3 individual_centers is a multiple of PI/6*number_of_hex_layers
    //This is true in an hex_patt but i do not know if it is sufficient

    auto v_modulus = modulus_of_btw_ind_angles(s,M_PI/12);
    for(auto ind_modulus : v_modulus)
      assert( ind_modulus < 0.0000000001 || ind_modulus > M_PI/12 - 0.000001);

  }


  //No individuals are overlapping at the start of the simulation
  {
    simulation s(2);
    assert(!has_collision(s));
  }


  //An individual position can be retrieved as a pair object x,y
  {
    simulation s;
    std::pair<double, double> v{0,0};//default coordinates of individuals
    std::pair<double, double> v2{1,1};//different coordinates from default

    assert(s.get_ind_pos(0) == v);
    assert(s.get_ind_pos(0) != v2);
  }


  // An individaul can be placed to some given coordinates after initialization
  {
    simulation s(2);
    for(int i = 0; i < s.get_pop_size(); i++)
      {
        s.set_ind_pos(s.get_ind(i),i,i);
      }

    for (int i = 0; i < s.get_pop_size() - 1; i++)
      {
        for(int j = i+1 ; j < s.get_pop_size(); j++)
          {
            assert(s.get_ind_pos(i) != s.get_ind_pos(j));
          }
      }

  }


  //When an individual divides it adds a copy of itself
  //to the population/individual vector(i.e divides)
  {
    simulation s;
    double lhs = s.get_pop_size();
    //let's allow all individuals in the population to reproduce,
    //to facilitate the testing conditions
    std::vector<int> dividing_pop(s.get_pop().size()) ;
    std::iota (std::begin(dividing_pop), std::end(dividing_pop), 0);

    s.do_division(dividing_pop);

    double rhs = s.get_pop_size();
    assert(lhs * 2 - rhs < 0.0000001);

  }


  //Only individuals with energy >= than the treshold will divide
  {
    simulation s(3);
    //This ind will not reproduce
    s.get_ind(0).set_energy(s.get_ind(0).get_treshold_energy()-1);
    //This ind will reproduce with no extra energy
    s.get_ind(1).set_energy(s.get_ind(1).get_treshold_energy());
    //This ind will reproduce with extra energy
    s.get_ind(2).set_energy(s.get_ind(2).get_treshold_energy()+1);

    std::vector<int> div_ind = s.get_dividing_individuals();
    assert(div_ind.size() > 0);

    for(unsigned int i = 0; i != div_ind.size(); i++)
      assert(s.get_ind(div_ind[i]).get_energy()
             >= s.get_ind(div_ind[i]).get_treshold_energy());

  }

  //The excess energy of dividing individuals is equal
  //to the diference between their energy and their treshold energy
  {
    simulation s(2);
    double energy = s.get_ind_tr_en(0);
    s.set_ind_en(0,energy);
    s.set_ind_en(1,energy*2);
    auto v_ex_en = s.get_excess_energies(s.get_dividing_individuals());
    for(int i =0; i != static_cast<int>(v_ex_en.size()); i++){
        assert(v_ex_en[static_cast<unsigned int>(i)]
            - (s.get_ind_en(i) - s.get_ind_tr_en(i)) < 0.00000001);
      }
  }


  //The distance between the two offspring from same mother in population vector
  //is the distance betwen the mother index and the end of the population vector

  {
    simulation s(1);
    //First individual reproduces
    s.set_ind_en(0,s.get_ind_tr_en(0)*2);
    auto daughters = s.get_sisters_index_offset();
    // In this case the distance between the two offspring
    //will be 1 element of the vector in next gen
    auto sister_distances = 1;
    for(auto sisters : daughters)
      {
        assert( sisters.second - sisters.first - sister_distances < 0.000001);
      }
  }


  //After a reproduction round individuals with enough energy will divide
  //and redistribute their remaining energy to their daughter cells
  {
    simulation s(3);
    //This individual will not reproduce
    s.set_ind_en(0,s.get_ind_tr_en(0) - 1);
    //This individual will reproduce with no extra energy
    s.set_ind_en(1,s.get_ind_tr_en(1));
    //This ind will reproduce with extra energy
    s.set_ind_en(2,s.get_ind_tr_en(2) + 1);

    auto mother_excess_energy = s.get_excess_energies(s.get_dividing_individuals());
    auto sister_cells_vector = s.get_sisters_index_offset();
    s.do_reprduction();
    for(unsigned int i = 0; i < s.get_sisters_index_offset().size(); i++)
      {
        assert(s.get_ind_en(sister_cells_vector[i].first)
               - mother_excess_energy[i]/2 < 0.00001);
        assert(s.get_ind_en(sister_cells_vector[i].second)
               - mother_excess_energy[i]/2 < 0.00001);
        assert(s.get_ind_en(sister_cells_vector[i].first)
               - s.get_ind_en(sister_cells_vector[i].second) < 0.00001);
      }
  }


  //After reproduction the first daughter individual takes the position of the mother
  {
    simulation s;
    auto parent_pop = s.get_pop();
    //setting energy high enough for the individual
    //to reproduce without offspring having negative enrgies
    s.set_ind_en(0,s.get_ind_tr_en(0));

    s.divide_ind(s.get_ind(0));
    //The first daughter cell is at the same index of the mother cell
    assert(s.get_ind_pos(0) == get_pos(parent_pop[0]) );
  }

  //After reproduction the second daughter individual
  //is placed just outside the position of the first daughter cell
  {
    simulation s;
    //setting energy high enough for the individual to reproduce
    //without offspring having negative enrgies
    s.set_ind_en(0,s.get_ind_tr_en(0));
    s.divide_ind(s.get_ind(0));
    //The first daughter cell is at the same index of the mother cell
    assert(!has_collision(s));
  }

  //If there are collisions individuals are displaced
  //until there are no more collisions
  {
    simulation s(2);
    assert(!has_collision(s));
    s.set_ind_pos(s.get_ind(1),s.get_ind_pos(0));
    assert(has_collision(s));
    manage_static_collisions(s);
    assert(!has_collision(s));

  }
  //After reproduction new collisions caused by new individuals
  //being placed where other individuals already are managed
  {
    simulation s(7);
    assert(!has_collision(s));
    //add 1 individual overlapping with central individual
    s.get_pop().push_back(individual(0,0));
    assert(has_collision(s));

    while(has_collision(s))
      {
        manage_static_collisions(s);
      }
    assert(!has_collision(s));
  }

  //A simulation has an environment
  // The value -1234567890 is irrelevant: just get this to compile

  {
    simulation s;
    assert(s.get_env().get_grid_side() > -1234567890);
  }

  //A simulation can be initialized with a certain amount
  //of food in all the cells of its environment grid
  //1 by default
  {
    simulation s;
    for( auto& grid_cell : s.get_env().get_grid())
      {
        assert(grid_cell.get_food() - 1 < 0.000001
               && grid_cell.get_food() - 1 > -0.0000001);
      }

    double starting_food = 3.14;
    s = simulation(1, 1, 0.1, starting_food);
    for( auto& grid_cell : s.get_env().get_grid())
      {
        assert(grid_cell.get_food() - starting_food < 0.000001
               && grid_cell.get_food() - starting_food > -0.0000001);
      }
  }
  //Individuals can take up energy from the environment
  {
    simulation s;
    double total_food_init = std::accumulate
        (
          s.get_env().get_grid().begin(),
          s.get_env().get_grid().end(),0,
          [](double sum,const env_grid_cell& c){return sum + c.get_food();}
    );

    double total_en_init = std::accumulate
        (
          s.get_pop().begin(),
          s.get_pop().end(),0,
          [](double sum,const individual& i){return sum + i.get_energy();}
    );

    feeding(s);

    double total_food_after = std::accumulate
        (
          s.get_env().get_grid().begin(),
          s.get_env().get_grid().end(),0.0,
          [](double sum, const env_grid_cell& c){return sum + c.get_food();}
    );

    double total_en_after = std::accumulate
        (
          s.get_pop().begin(),
          s.get_pop().end(),0.0,
          [](double sum,const individual& i){
        return sum + i.get_energy();}
    );

    assert(total_food_init > total_food_after);
    assert(total_en_init < total_en_after);
    assert(total_en_after - s.get_pop_size() * s.get_ind(0).get_uptake_rate()
           - total_en_init < 0.00000001
           &&
           total_en_after - s.get_pop_size() * s.get_ind(0).get_uptake_rate()
           - total_en_init > -0.00000001);
  }
  //Individuals outside the grid do not feed
  {
    simulation s;

    double total_food_init = std::accumulate
        (
          s.get_env().get_grid().begin(),
          s.get_env().get_grid().end(),0.0,
          [](double sum,const env_grid_cell& c){return sum + c.get_food();}
    );

    double total_en_init = std::accumulate
        (
          s.get_pop().begin(),
          s.get_pop().end(),0.0,
          [](double sum,const individual& i){return sum + i.get_energy();}
    );

    set_pos(s.get_ind(0),std::pair<double,double>(-42,42));

    feeding(s);

    double total_food_after = std::accumulate
        (
          s.get_env().get_grid().begin(),
          s.get_env().get_grid().end(),0.0,
          [](double sum, const env_grid_cell& c){return sum + c.get_food();}
    );

    double total_en_after = std::accumulate
        (
          s.get_pop().begin(),
          s.get_pop().end(),0.0,
          [](double sum,const individual& i){
        return sum + i.get_energy();}
    );

    assert(total_food_init - total_food_after < 0.00001
           && total_food_init - total_food_after > -0.00001);
    assert(total_en_init - total_en_after < 0.000001
           && total_en_init - total_en_after > -0.000001);
  }

  //Individuals lose energy through metabolism
  {
    simulation s;
    feeding(s);
    double init_en_tot = std::accumulate
        (
          s.get_pop().begin(),
          s.get_pop().end(),0.0,
          [](double sum,const individual& i){return sum + i.get_energy();}
    );

    metabolism_pop(s);
    double after_en_tot = std::accumulate
        (
          s.get_pop().begin(),
          s.get_pop().end(),0.0,
          [](double sum,const individual& i){return sum + i.get_energy();}
    );
    assert(after_en_tot < init_en_tot);
  }

  //In one tick/timestep individuals first feed,
  //than reproduce, than substances diffuse
  {
    simulation s(1,3);

    //The single individual in this population
    //after a tick should reproduce
    auto init_pop_size = s.get_pop().size();
    s.get_ind(0).set_energy(s.get_ind_tr_en(0)
                            + s.get_ind(0).get_metab_rate() + 0.01
                            - s.get_ind(0).get_uptake_rate());

    //and the grid_cell where it is should recieve
    //food nutrients
    auto grid_index_ind =  find_grid_index(s.get_ind(0),s.get_env().get_grid_side());
    double food_after_feeding =
        s.get_env().get_cell(grid_index_ind).get_food() - s.get_ind(0).get_uptake_rate();

    tick(s);

    assert(!(food_after_feeding - s.get_env().get_cell(grid_index_ind).get_food() < 0.00001
             && food_after_feeding -s.get_env().get_cell(grid_index_ind).get_food() > -0.0001));
    assert(s.get_pop().size() == 2 * init_pop_size);
    assert(!has_collision(s));
  }
  //every timestep/tick collisions are handled
  {
    simulation s(7,3);
    //The central individual in this population
    //after a tick should reproduce
    auto init_pop_size = s.get_pop().size();
    s.get_ind(1).set_energy(s.get_ind_tr_en(1)
                            + s.get_ind(1).get_metab_rate() + 0.01
                            - s.get_ind(1).get_uptake_rate());
    tick(s);

    assert(s.get_pop().size() == init_pop_size + 1);
    assert(!has_collision(s));
  }
  //A simulation is initialized with a m_tick = 0;
  {
    simulation s;
    assert(s.get_tick() == 0);
  }
  //After each tick the simulation updates its m_tick
  //by one
  {
    simulation s;
    for(int i = 0; i != 3; i++)
      {
        assert(s.get_tick() == i);
        tick(s);
      }
  }

  //A simulation can be run for a certain amount of ticks
  {
    simulation s;
    auto n_ticks = 100;
    exec(s, n_ticks);
    assert(s.get_tick() - n_ticks == 0);
  }

  //Spores do not get detected when looking for dividing individual
  {
    simulation s;
    s.set_ind_en(0,s.get_ind_tr_en(0));
    assert(s.get_dividing_individuals()[0] == 0);
    s.get_ind(0).set_type(individual_type::spore);
    assert(s.get_dividing_individuals().empty());
  }

  //Spores do not reproduce
  {
    simulation s;
    s.get_ind(0).set_type(individual_type::spore);
    s.set_ind_en(0,s.get_ind_tr_en(0));
    auto init_pop_size = s.get_pop_size();
    s.set_ind_en(0,s.get_ind_tr_en(0));
    tick(s);
    assert(init_pop_size == s.get_pop_size());
  }

  //Spores do not feed
  {
    simulation s;
    auto init_food = s.get_env().get_cell(0).get_food();
    s.get_ind(0).set_type(individual_type::spore);
    feeding(s);
    assert(init_food - s.get_env().get_cell(0).get_food() < 0.000001
           && init_food - s.get_env().get_cell(0).get_food() > -0.000001);
  }
  //Spores do not lose energy
  {
    simulation s;
    s.get_ind(0).set_type(individual_type::spore);
    s.set_ind_en(0,s.get_ind_tr_en(0));
    auto init_en_ind0 = s.get_ind_en(0);
    tick(s);
    assert(s.get_ind_en(0) - init_en_ind0 < 0.000001
           && s.get_ind_en(0) - init_en_ind0 > -0.000001);
  }

  //Sporulating individuals do not get detected when looking for dividing individual
  {
    simulation s;
    s.set_ind_en(0,s.get_ind_tr_en(0));
    assert(s.get_dividing_individuals()[0] == 0);
    s.get_ind(0).set_type(individual_type::sporulating);
    assert(s.get_dividing_individuals().empty());
  }

  //Sporulating individuals cannot reproduce
  {
    simulation s;
    s.get_ind(0).set_type(individual_type::sporulating);
    s.set_ind_en(0,s.get_ind_tr_en(0));
    auto init_pop_size = s.get_pop_size();
    s.set_ind_en(0,s.get_ind_tr_en(0));
    tick(s);
    assert(init_pop_size == s.get_pop_size());
  }

  //Sporulating individuals do not feed but they lose energy
  {
    simulation s;
    s.get_ind(0).set_type(individual_type::sporulating);
    s.set_ind_en(0,1);
    auto init_en_ind0 = s.get_ind_en(0);
    auto init_food = s.get_env().get_cell(0).get_food();
    tick(s);
    assert(init_en_ind0 - s.get_ind_en(0) - s.get_ind(0).get_metab_rate() < 0.000001
           && init_en_ind0 - s.get_ind_en(0) - s.get_ind(0).get_metab_rate() > -0.000001);
    assert(init_food - s.get_env().get_cell(0).get_food() < 0.000001
           && init_food - s.get_env().get_cell(0).get_food() > -0.000001);
  }


}

















