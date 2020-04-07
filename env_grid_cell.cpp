#include "env_grid_cell.h"
#include <cassert>
#include <cmath>
#include <numeric>


env_grid_cell::env_grid_cell(double metabolite, double food, double max_food):
  m_metabolite(metabolite),
  m_food(food),
  m_max_food(max_food)
{

}

bool operator == (const env_grid_cell& lhs, const env_grid_cell& rhs) noexcept
{
  return lhs.get_food() - rhs.get_food() < 0.00001 &&
      lhs.get_food() - rhs.get_food() > -0.00001 &&
      lhs.get_metabolite() - rhs.get_metabolite() < 0.00001 &&
      lhs.get_metabolite() - rhs.get_metabolite() > -0.00001;
}

bool operator != (const env_grid_cell& lhs, const env_grid_cell& rhs) noexcept
{
  return !(lhs == rhs);
}

double calc_in_out_flux( const env_grid_cell& cell, double av_food_flux,
                          double diffusion_coeff) noexcept
{
  auto in_out_flux = cell.get_v_neighbors().size() * av_food_flux * diffusion_coeff;
  if(cell.get_food() + in_out_flux < 0)
    {
      return - cell.get_food();
    }
  return in_out_flux;
}


double calc_exiting_metabolite(const env_grid_cell& cell, double tot_metab_delta,
                               double diffusion_coeff) noexcept
{
  return tot_metab_delta * diffusion_coeff >= 1 ?
        cell.get_metabolite() :
        std::max(cell.get_metabolite() * tot_metab_delta * diffusion_coeff, 0.0);
}

double food_flux(const env_grid_cell &lhs, const env_grid_cell &rhs) noexcept
{
  return  rhs.get_food() - lhs.get_food();
}

double metab_difference(const env_grid_cell &lhs, const env_grid_cell &rhs) noexcept
{
  return lhs.get_metabolite() - rhs.get_metabolite() > 0 ?
        lhs.get_metabolite() - rhs.get_metabolite() : 0;
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
    assert( g.get_metabolite() > -123456789);
  }

  //Concentration of metabolite can be set
  {
    env_grid_cell g;
    double new_metabolite_conc = 3;
    assert(g.get_metabolite() - new_metabolite_conc > 0.0000001
           || g.get_metabolite() - new_metabolite_conc < -0.0000001);
    g.set_metabolite(new_metabolite_conc);
    assert(g.get_metabolite() - new_metabolite_conc < 0.00001);
  }

  //Concentration of metabolite can be changed
  {
    env_grid_cell g;
    env_grid_cell g1 = g;
    double change_in_metabolite = 3;
    assert(g.get_metabolite() - change_in_metabolite > 0.0000001
           || g.get_metabolite() - change_in_metabolite < -0.0000001);
    g1.increment_metabolite(change_in_metabolite);
    assert(g.get_metabolite() - g1.get_metabolite() > 0.000000001
           || g.get_metabolite() - g1.get_metabolite() < -0.000000001);
    assert(std::abs(g.get_metabolite() - g1.get_metabolite()) -
           change_in_metabolite < 0.000001);

  }

  //Metabolite concentration cannot be less than 0
  {
    env_grid_cell g;
    g.increment_metabolite(-10000);
    assert(g.get_metabolite() >= 0);
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
    g1.increment_food(change_in_food);
    assert(g.get_food() - g1.get_food() > 0.00000001
           || g.get_food() - g1.get_food() < -0.00000001);
    assert(std::abs(g.get_food() - g1.get_food()) - change_in_food < 0.000001);

  }

  //Food concentration cannot be less than 0
  {
    env_grid_cell g;
    g.increment_food(-10000);
    assert(g.get_food() >= 0);
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

  //The substance difference between two cells can be found
  //by subtracting the amount of substance of the neighbor cell
  //from the focal,
  //if the difference is negative, it is considered 0
  {
    auto food_lhs = 0.0;
    auto food_rhs = 1.0;
    auto metab_lhs = 0.0;
    auto metab_rhs = 1.0;
    env_grid_cell lhs(metab_lhs, food_lhs);
    env_grid_cell rhs(metab_rhs, food_rhs);

    auto a = food_flux(lhs,rhs) - (food_rhs - food_lhs) ;
    assert(a < 0.0000001
           && a > -0.0000001);
    assert(metab_difference(lhs,rhs)< 0.0000001
           && metab_difference(lhs,rhs) > -0.0000001);
  }

  //A cell will lose due to diffusion an proportion of his metabolite equal to
  //The sum of differences in metabolite with all its neighbors * the  diffusion coefficient
  //If this proportion is bigger than 1 it will lose all its food
  {
    //Check for case in which exiting metabolite >= cell_metabolite
    //And therefore exiting metabolite = cell_metabolite
    double diffusion_coeff = 1;
    double tot_delta = 3;//Three neighbours each with 1 metabolite less than focal cell
    env_grid_cell c(1,0);
    assert(calc_exiting_metabolite(c,tot_delta,diffusion_coeff) - c.get_metabolite() < 0.000001 &&
           calc_exiting_metabolite(c,tot_delta,diffusion_coeff) - c.get_metabolite() > -0.000001);
    //Check for case in which exiting metabolite < cell_metabolite
    //-> exiting food = cell_metabolite * tot_difference * diffusion_metabolite
    diffusion_coeff = 0.1;
    assert(calc_exiting_metabolite(c,tot_delta,diffusion_coeff)
           - tot_delta * diffusion_coeff < 0.000001 &&
           calc_exiting_metabolite(c,tot_delta,diffusion_coeff)
           - tot_delta * diffusion_coeff > -0.000001);
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

