#include "env_grid_cell.h"
#include <cassert>
#include <cmath>
#include <numeric>


env_grid_cell::env_grid_cell(double metabolite, double food, double max_food, double max_metab):
  m_metabolite{metabolite},
  m_food{food},
  m_max_food{max_food},
  m_max_metab{max_metab}
{

}

bool operator == (const env_grid_cell& lhs, const env_grid_cell& rhs) noexcept
{
  return lhs.get_food() - rhs.get_food() < 0.00001 &&
      lhs.get_food() - rhs.get_food() > -0.00001 &&
      lhs.get_metab() - rhs.get_metab() < 0.00001 &&
      lhs.get_metab() - rhs.get_metab() > -0.00001;
}

bool operator != (const env_grid_cell& lhs, const env_grid_cell& rhs) noexcept
{
  return !(lhs == rhs);
}

double calc_food_flux( const env_grid_cell& cell, double av_food_flux,
                          double diffusion_coeff) noexcept
{
  auto in_out_flux = cell.get_v_neighbors().size() * av_food_flux * diffusion_coeff;
  if(cell.get_food() + in_out_flux < 0)
    {
      return - cell.get_food();
    }
  return in_out_flux;
}


double calc_metab_flux( const env_grid_cell& cell, double av_metab_flux,
                                double diffusion_coeff) noexcept
{
  auto in_out_flux = cell.get_v_neighbors().size() * av_metab_flux * diffusion_coeff;
  if(cell.get_metab() + in_out_flux < 0)
    {
      return - cell.get_metab();
    }
  return in_out_flux;
}

double food_flux(const env_grid_cell &lhs, const env_grid_cell &rhs) noexcept
{
  return  rhs.get_food() - lhs.get_food();
}

double metab_flux(const env_grid_cell &lhs, const env_grid_cell &rhs) noexcept
{
  return rhs.get_metab()- lhs.get_metab();
}

void metabolite_degrades(env_grid_cell& g, double degrad_rate) noexcept
{
   auto new_metab = g.get_metab() - degrad_rate;
   if(new_metab < 0)
   {
       g.set_metab(0);
       return;
   }
   g.set_metab(new_metab);
}
void test_env_grid_cell()//!OCLINT tests may be many
{
#ifndef NDEBUG

  //A grid cell has a metabolite value,
  //The value -123456789 is irrelevant
  //is just to make sure that get_metabolite()
  //is called properly
  {
    env_grid_cell g;
    assert( g.get_metab() > -123456789);
  }

  //Concentration of metabolite can be set
  {
    env_grid_cell g;
    double new_metabolite_conc = 3;
    assert(g.get_metab() - new_metabolite_conc > 0.0000001
           || g.get_metab() - new_metabolite_conc < -0.0000001);
    g.set_metab(new_metabolite_conc);
    assert(g.get_metab() - new_metabolite_conc < 0.00001);
  }

  //Concentration of metabolite can be changed
  {
    env_grid_cell g;
    env_grid_cell g1 = g;
    double change_in_metabolite = 3;
    assert(g.get_metab() - change_in_metabolite > 0.0000001
           || g.get_metab() - change_in_metabolite < -0.0000001);
    g1.set_metab_change(change_in_metabolite);
    g1.increment_metabolite();
    assert(g1.get_metab_change() < 0.00001 && g1.get_metab_change() > -0.0001);

    assert(g.get_metab() - g1.get_metab() > 0.000000001
           || g.get_metab() - g1.get_metab() < -0.000000001);
    assert(std::abs(g.get_metab() - g1.get_metab()) -
           change_in_metabolite < 0.000001);

  }

  //Metabolite concentration cannot be less than 0
  {
    env_grid_cell g;
    g.set_metab_change(-10000);
    g.increment_metabolite();
    assert(g.get_metab_change() < 0.00001 && g.get_metab_change() > -0.0001);
    assert(g.get_metab() >= 0);
  }

  //A grid cell has a food value,
  //The value -123456789 is irrelevant
  //is just to make sure that get_food()
  //is called properly
  {
    env_grid_cell g;
    assert( g.get_food() > -123456789);
  }

  //Concentration of food can be set
  {
    env_grid_cell g;
    double new_food_conc = 3;
    assert(g.get_food() - new_food_conc > 0.00000001
           || g.get_food() - new_food_conc < -0.00000001);
    g.set_food(new_food_conc);
    assert(g.get_food() - new_food_conc < 0.00001);
  }

  //Concentration of food can be changed
  {
    env_grid_cell g;
    env_grid_cell g1 = g;
    double change_in_food = 3;
    assert(g.get_food() - change_in_food > 0.0000001
           || g.get_food() - change_in_food < -0.0000001);
    g1.set_food_change(change_in_food);
    g1.increment_food();
    assert(g1.get_food_change() < 0.00001 && g1.get_food_change() > -0.0001);
    assert(g.get_food() - g1.get_food() > 0.00000001
           || g.get_food() - g1.get_food() < -0.00000001);
    assert(std::abs(g.get_food() - g1.get_food()) - change_in_food < 0.000001);

  }

  //Food concentration cannot be less than 0
  {
    env_grid_cell g;
    g.set_food_change(-10000);
    g.increment_food();
    assert(g.get_food_change() < 0.00001 && g.get_food_change() > -0.0001);
    assert(g.get_food() >= 0);
  }

  //Metabolite degrades at a constant rate
    {
        double init_metab = 3.0;
        double degrad_rate = init_metab;
        env_grid_cell g;
        g.set_metab(init_metab);
        metabolite_degrades(g, degrad_rate);
        assert(init_metab > g.get_metab());
        //Metabolite con should be 0
        assert(g.get_metab() < 0.000001 && g.get_metab() > -0.000001);
        //The metabolite would go below 0 so it stays 0
        metabolite_degrades(g, degrad_rate);
        assert(g.get_metab() < 0.000001 && g.get_metab() > -0.000001);
    }

  //A grid_cell has a vector that can store
  //the indexes of its neighbouring cells if placed in a grid
  {
    env_grid_cell g;
    assert(g.get_v_neighbors().size() == 0);
  }

  //A grid cell vector of neighbors can be set
  {
    //all neighors will have the same address for the purpose of this test
    int neighbors_address = 42;
    std::vector<int> list_of_neighbors(42,neighbors_address);
    env_grid_cell g;
    g.set_v_neighbors(list_of_neighbors);
    assert(g.get_v_neighbors().size() - list_of_neighbors.size() < 0.0000001);
    for(auto neighbor : g.get_v_neighbors())
      {
        assert(neighbor == neighbors_address);
      }
  }

  //The flux between two cells can be found
  //by subtracting the amount of substance of the focal cell
  //from the neighbor,
  {
    auto food_lhs = 0.0;
    auto food_rhs = 1.0;
    auto metab_lhs = 0.0;
    auto metab_rhs = 1.0;
    env_grid_cell lhs(metab_lhs, food_lhs);
    env_grid_cell rhs(metab_rhs, food_rhs);

    auto flux_food = food_flux(lhs,rhs) - (food_rhs - food_lhs) ;
    assert(flux_food < 0.0000001
           && flux_food > -0.0000001);
    auto flux_metab = metab_flux(lhs,rhs) - (metab_rhs - metab_lhs) ;

    assert(flux_metab < 0.0000001
           && flux_metab > -0.0000001);
  }


  //Env_grid_cell has a boolean operator
  {
    env_grid_cell c;
    auto c1 = c;
    assert(c == c1);
    c1.set_food(42);
    assert(c != c1);
  }

#endif
}

