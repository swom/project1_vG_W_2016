#include "individual.h"
#include <cassert>
#include <math.h>

individual::individual(double x_pos, double y_pos,
                       double size, double energy,
                       double treshold_energy):
  m_x(x_pos),
  m_y(y_pos),
  m_size(size),
  m_energy(energy),
  m_treshold_energy(treshold_energy)
{

}

bool operator == (const individual& lhs, const individual& rhs) noexcept {
  //check that two individuals' centers are at the same coordinates
  return lhs.get_x() - rhs.get_x() < 0.00001
      && lhs.get_y() - rhs.get_y() < 0.000001;
}

bool are_colliding(const individual &lhs, const individual &rhs) noexcept
{
  const double dx = std::abs(lhs.get_x() - rhs.get_x());
  const double dy = std::abs(lhs.get_y() - rhs.get_y());
  const double actual_distance = ((dx * dx) + (dy * dy)) ;
  const double collision_distance =
      (lhs.get_size() + rhs.get_size()) * (lhs.get_size() + rhs.get_size()) ;
  return actual_distance < collision_distance;
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

double distance(const individual& lhs, const individual& rhs) noexcept
{
  return sqrt((lhs.get_x() - rhs.get_x()) * (lhs.get_x() - rhs.get_x())
              + (lhs.get_y() - rhs.get_y()) * (lhs.get_y() - rhs.get_y()));
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

  //Individuals should be initialized with 0 internal energy
  {
    individual i(0,0);
    assert(i.get_energy() - 0.0 < 0.000001);
  }

  //An individual is initialized with a treshold level of energy
  {
    double treshold_energy = 3;
    individual i(0,0,0,0,treshold_energy);
    assert(i.get_treshold_energy() - treshold_energy < 0.00000001);
  }

  // an individual's energy can be changed
  {
    individual i(0,0);
    double lhs = i.get_energy();
    double new_energy = 3 + i.get_energy();
    i.set_energy(new_energy);
    double rhs = i.get_energy();
    assert(abs(rhs - lhs)>0.0000001);
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

  //An individual can be displaced given a certain overlap and distance value
  //This is used to manage static_collisions in simulation.cpp
//  {
//    individual lhs(0,0);
//    individual rhs(0,0);
//    displace(lhs,);
//    assert()
//  }

}
