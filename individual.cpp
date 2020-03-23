#include "individual.h"
#include <cassert>
#include <math.h>

individual::individual(double x_pos, double y_pos,
                       double size, double energy,
                       double treshold_energy, double uptake_rate, double metabolic_rate,
                       individual_type individual_type):
  m_x(x_pos),
  m_y(y_pos),
  m_size(size),
  m_energy(energy),
  m_treshold_energy(treshold_energy),
  m_uptake_rate(uptake_rate),
  m_metab_rate(metabolic_rate),
  m_individual_type(individual_type)
{

}

bool operator == (const individual& lhs, const individual& rhs) noexcept {
  //check that two individuals' centers are at the same coordinates
  return lhs.get_x() - rhs.get_x() < 0.00001
      && lhs.get_y() - rhs.get_y() < 0.000001;
}

bool are_colliding(const individual &lhs, const individual &rhs) noexcept
{
  const double dx = lhs.get_x() - rhs.get_x();
  const double dy = lhs.get_y() - rhs.get_y();
  const double actual_distance = (dx * dx) + (dy * dy) ;
  const double collision_distance = lhs.get_size() + rhs.get_size();
  return actual_distance < collision_distance;
}

double distance(const individual& lhs, const individual& rhs) noexcept
{
  return sqrt((lhs.get_x() - rhs.get_x()) * (lhs.get_x() - rhs.get_x())
              + (lhs.get_y() - rhs.get_y()) * (lhs.get_y() - rhs.get_y()));
}

void feed(individual& i, env_grid_cell& c) noexcept
{
  if(c.get_food() > i.get_uptake_rate())
    {
      c.increment_food(- i.get_uptake_rate());
      i.change_en(i.get_uptake_rate());
    }
  else
    {
      i.change_en(c.get_food());
      c.set_food(0);
    }
}


int find_grid_index( individual& i, double grid_side)
{
  auto x_offset = i.get_x() + grid_side/2;
  auto y_offset = i.get_y() + grid_side/2;
  //cast correctly negative numbers so to keep them outside the grid
  int x_index_offset = x_offset < 0 ?
        static_cast<int>(x_offset - 1) : static_cast<int>(x_offset);
  int y_index_offset = y_offset < 0 ?
        static_cast<int>(y_offset - 1) : static_cast<int>(y_offset);

  if(x_offset >= grid_side || x_offset < 0
     || y_offset >= grid_side || y_offset < 0)
    {return  -100;}
  return x_index_offset + y_index_offset * static_cast<int>(grid_side);
}

const std::pair<double,double> get_d2_pos(individual& i) noexcept
{
  std::pair<double, double> pos;
  pos.first = i.get_x() + cos(0)*2*i.get_size();
  pos.second += i.get_y() + sin(0)*2*i.get_size();
  return pos;
}


const std::pair<double, double> get_pos(individual& i)  noexcept
{
  std::pair<double, double> pos;
  pos.first = i.get_x();
  pos.second = i.get_y();
  return pos;
}



double overlap(const individual& lhs, const individual& rhs) noexcept
{
  return (distance(lhs,rhs) - lhs.get_size() - rhs.get_size())/2;
}

void set_pos(individual& i, std::pair<double, double> pos)  noexcept
{
  i.set_x(pos.first);
  i.set_y(pos.second);
}

std::pair<double, double> get_displacement(const individual& lhs, const individual& rhs) noexcept
{
  std::pair<double, double> displ_x_y;
  displ_x_y.first = overlap(lhs,rhs) *
      (distance(lhs, rhs) > 0 ? (lhs.get_x() - rhs.get_x()) / distance(lhs, rhs) :
                                (lhs.get_size() + rhs.get_size())/(2*sqrt(2)));
  displ_x_y.second = overlap(lhs,rhs) *
      (distance(lhs, rhs) > 0 ? (lhs.get_y() - rhs.get_y()) / distance(lhs, rhs) :
                                (lhs.get_size() + rhs.get_size())/(2*sqrt(2)));
  return  displ_x_y;
}

void displace(individual& lhs, individual& rhs) noexcept
{
  auto displacement = get_displacement(lhs,rhs);
  lhs.change_x(displacement.first);
  lhs.change_y(displacement.second);
  rhs.change_x(-displacement.first);
  rhs.change_y(-displacement.second);
}


void metabolism(individual& i) noexcept
{
  if(i.get_energy() >= i.get_metab_rate()){i.change_en(-i.get_metab_rate());}
  else{i.set_energy(0);}
}

void test_individual()//!OCLINT tests may be many
{

  //An individual should be initialized with the defined starting size
  {
    double starting_size = 14.0;
    individual i(0,0,starting_size);
    assert(i.get_size() - starting_size < 0.0000001);
  }

  //An individual should be initialized at a certain position
  {
    double x = 100;
    double y = 100;
    individual i(x,y);
    assert(i.get_x() - x < 0.0000001);
    assert(i.get_y() - y < 0.0000001);
  }

  //An individual is initialized with an m_food_uptake value
  //0.1 by default
  {
    individual i;
    assert(i.get_uptake_rate() - 0.1 < 0.000001);

    double uptake_rate = 0.3;
    i = individual(0, 0, 0, 0, 0, uptake_rate);
    assert(i.get_uptake_rate() - uptake_rate < 0.000001);
  }

  //An individual is initialized with a m_metab_rate
  //0.01 by default
  {
    individual i;
    assert(i.get_metab_rate() - 0.01 < 0.000001
           && i.get_metab_rate() - 0.01 > -0.000001);
    double metabolic_rate = 2;
    i = individual(0,0,0,0,0,0,metabolic_rate);
    assert(i.get_metab_rate() - metabolic_rate < 0.000001
           && i.get_metab_rate() - metabolic_rate > -0.000001);

  }
  //Individuals should be initialized with 0 internal energy
  {
    individual i;
    assert(i.get_energy() - 0.0 < 0.000001);
  }

  //An individual is initialized with a treshold level of energy
  {
    double treshold_energy = 3;
    individual i(0,0,0,0,treshold_energy);
    assert(i.get_treshold_energy() - treshold_energy < 0.00000001);
  }

  // an individual's energy can be set
  {
    individual i(0,0);
    double lhs = i.get_energy();
    double new_energy = 3 + i.get_energy();
    i.set_energy(new_energy);
    double rhs = i.get_energy();
    assert(abs(rhs - lhs)>0.0000001);
  }

  //an individual energy can be changed
  {
    individual i;
    double init_en = i.get_energy();
    double en_change = 3.14;
    i.change_en(en_change);
    assert(i.get_energy() - en_change - init_en < 0.000001);
  }

  //Energy after reproduction should be half of the excess of energy
  {
    individual i(0,0,0,4,2);//this individual should have energy in excess = 2
    //after division
    double excess_energy = i.get_energy()-i.get_treshold_energy();
    assert(i.split_excess_energy() - excess_energy/2 < 0.000000001);

  }
  // An individual position can be extracted as a pair of doubles x,y
  {
    individual i(0,0);
    assert(get_pos(i).first - i.get_x() < 0.00001);
    assert(get_pos(i).second - i.get_y() < 0.00001);
  }

  //The distance between two individuals can be calculated
  {
    std::vector<individual> pop(2, individual(0,0));
    //place individuals along x axis at 1 distance from each other
    for (auto i = 0; i != static_cast<int>(pop.size()); i++)
      {
        pop[static_cast<unsigned int>(i)].set_x(i);
      }
    for (unsigned int i = 0; i != pop.size(); i++)
      {
        for (unsigned int j = 0; j != pop.size(); j++)
          {
            auto distance_i_j = distance(pop[i],pop[j]);
            if(i == j)
              {
                assert(distance_i_j < 0.000001 && distance_i_j > -0.00001);
              }
            else if(i != j)
              {
                assert(distance_i_j - 1 < 0.000001 && distance_i_j - 1 > -0.00001);
              }
          }
      }
  }

  //The position of the second daughter cell is just outside the mother cell
  {
    individual mom(0,0);
    auto daughter_pos = get_d2_pos(mom);
    individual daughter2(daughter_pos.first, daughter_pos.second);
    assert(!are_colliding(mom,daughter2));
  }

  //The overlap of two individuals can be calculated
  {
    //Two individuals overlap by half their radius
    individual lhs(0,0);
    individual rhs(0,lhs.get_size());
    assert(overlap(lhs, rhs) - lhs.get_size() < 0.000001);
  }

  //The necessary displacement along the x and y axis can be calculated
  //for two individuals to not overlap
  //This is used to manage static_collisions in simulation.cpp
  {
    individual lhs(0,0);
    individual rhs(0,0);
    assert(get_displacement(lhs,rhs).first -
           (lhs.get_size() + rhs.get_size())/2 < 0.00001);
    assert(get_displacement(lhs,rhs).second -
           (lhs.get_size() + rhs.get_size())/2 < 0.00001);
  }

  //Two cells can be displaced so not to overlap anymore
  {
    individual lhs(0,0);
    individual rhs(0,0);
    displace(lhs, rhs);
    assert(distance(lhs,rhs) - (lhs.get_size() + rhs.get_size()) < 0.00001
           && distance(lhs,rhs) - (lhs.get_size() + rhs.get_size()) > -0.00001);
  }

  //The grid cell where an individual should be
  //can be deduced from individual's coordinates
  {
    individual i(0,0);
    int grid_side = 1;
    assert(find_grid_index(i, grid_side) == 0);
    grid_side = 3;
    assert(find_grid_index(i, grid_side) == 4);
    auto pos = std::pair<double,double>(-0.2, 1.2);
    set_pos(i,pos);
    assert(find_grid_index(i, grid_side) == 7);
  }

  //If an individual is past the grid side
  //find_grid_index() will return -100
  {
    individual i(0,0);
    int grid_side = 3;
    auto pos = std::pair<double,double>(1.7, 0);
    set_pos(i,pos);
    assert( find_grid_index(i, grid_side) == -100);

    pos = std::pair<double,double>(-1.7, 0);
    set_pos(i,pos);
    assert(find_grid_index(i, grid_side) == -100);
  }

  //If an individual is above or below the grid limit
  //find_grid_index() will return -100
  {
    individual i(0,0);
    int grid_side = 3;
    auto pos = std::pair<double,double>(0, -1.7);
    set_pos(i,pos);
    assert(find_grid_index(i, grid_side) == -100);

    pos = std::pair<double,double>(0, 1.7);
    set_pos(i,pos);
    assert(find_grid_index(i, grid_side) == -100);
  }

  //An individual can increase its energy by eating food
  //It also depletes the food
  {
    individual i;
    double init_en = i.get_energy();
    double init_food = 0.1;
    env_grid_cell c(0,init_food);
    feed(i,c);
    assert(i.get_energy() > init_en);
    assert(init_food - i.get_uptake_rate() < 0.00001);
    assert((i.get_energy() + i.get_uptake_rate() + c.get_food()) -
           (init_en + init_food) > 0.000001);

  }
  //if there is less food than a individual
  //can normally eat, it eats what is available
  {
    individual i;
    double init_en = i.get_energy();
    double init_food = 0.01;
    env_grid_cell c(0,init_food);
    assert(i.get_uptake_rate() > c.get_food());

    feed(i, c);
    assert(c.get_food() < 0.000001 && c.get_food() > -0.000001);
    assert(i.get_energy() - init_food - init_en < 0.000001);
  }

  //individuals lose energy due to their metabolism
  {
    individual i(0,0,0,1);
    auto init_en = i.get_energy();
    metabolism(i);
    auto en_after = i.get_energy();
    assert(init_en > en_after);
  }

  //An individual's energy cannot go below 0
  {
    individual i(0,0,0,0,0,0);
    assert(i.get_energy() - i.get_metab_rate() < 0);
    metabolism(i);
    assert(i.get_energy() < 0.000001
           && i.get_energy() > -0.0000001);
  }

  //An individual is initialized with a individual_type(living, sporulating, spore)
  //By default it is initialized with  individual_type::living
  {
    individual i;
    assert(to_str(i.get_type()) == "living");
    i = individual (0,0,0,0,0,0,0,individual_type::spore);
    assert(to_str(i.get_type()) == "spore");
    i = individual (0,0,0,0,0,0,0,individual_type::sporulating);
    assert(to_str(i.get_type()) == "sporulating");
  }

  //An individuals type can be changed after initialization
  {
    individual i;
    assert(to_str(i.get_type()) == "living");
    i.set_type(individual_type::spore);
    assert(to_str(i.get_type()) == "spore");
    assert(to_str(i.get_type()) != "living");


  }

}
