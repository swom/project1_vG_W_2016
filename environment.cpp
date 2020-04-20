#include "environment.h"
#include <algorithm>
#include <cassert>
#include <numeric>

environment::environment(int grid_side, double diff_coeff, double init_food):
  m_side(grid_side),
  m_diffusion_coefficient(diff_coeff),
  m_init_food(init_food),
  m_grid(static_cast<unsigned int>(m_side * m_side),env_grid_cell(0,m_init_food))
{
  find_neighbors_all_grid(*this);
}


bool operator == (const environment& lhs, const environment& rhs) noexcept
{

  return lhs.get_grid_size() == rhs.get_grid_size() &&
      lhs.get_diff_coeff() - rhs.get_diff_coeff() < 0.00001 &&
      lhs.get_diff_coeff() - rhs.get_diff_coeff() > -0.00001 &&
      std::equal(lhs.get_grid().begin(),lhs.get_grid().end(),
                 rhs.get_grid().begin(),rhs.get_grid().end(),
                 [](const env_grid_cell& l, const env_grid_cell& r) {return l == r;});
}

bool operator != (const environment& lhs, const environment& rhs) noexcept
{
  return !(lhs == rhs);
}

void calc_change_food(environment& e, int index_focal_cell) noexcept
{
  auto& focal_gridcell = e.get_cell(index_focal_cell);
  auto v_food_fluxes = get_neighbors_food_fluxes(focal_gridcell, e);
  if(std::none_of(v_food_fluxes.begin(), v_food_fluxes.end(),
                 [](double f){return  f > 0.000001 || f < -0.000001;}))
    {
      assert(focal_gridcell.get_food_change() < 0.000001 &&
             focal_gridcell.get_food_change() > -0.000001);
          return;
    }
  auto av_food_flux =
      std::accumulate(v_food_fluxes.begin(),v_food_fluxes.end(),0.0) /
      (v_food_fluxes.empty() ? 1 : v_food_fluxes.size());

  auto in_out_flux_food = calc_food_flux(focal_gridcell, av_food_flux, e.get_diff_coeff());
  focal_gridcell.set_food_change(in_out_flux_food);
}

void calc_change_metab(environment& e, int index_focal_cell) noexcept
{
  auto& focal_gridcell = e.get_cell(index_focal_cell);
  auto v_metab_fluxes = get_neighbors_metab_fluxes(focal_gridcell, e);
  if(std::none_of(v_metab_fluxes.begin(), v_metab_fluxes.end(),
                 [](double f){return  f < 0.000001 || f > -0.000001;}))
    {
      assert(focal_gridcell.get_metab_change() < 0.00000001 &&
             focal_gridcell.get_metab_change() > -0.00000001);
      return;
    }
  auto av_metab_flux =
      std::accumulate(v_metab_fluxes.begin(),v_metab_fluxes.end(),0.0) /
      (v_metab_fluxes.empty() ? 1 : v_metab_fluxes.size());

  auto in_out_flux_metab = calc_metab_flux(focal_gridcell, av_metab_flux, e.get_diff_coeff());
  focal_gridcell.set_metab_change(in_out_flux_metab);
}

void calc_diffusion(environment& e) noexcept
{
  for(int i = 0; i != e.get_grid_size(); i++)
    {
      calc_change_food(e, i);
      calc_change_metab(e, i);
    }
}

void calc_diffusion_food(environment& e) noexcept
{
  for(size_t i = 0; i != e.get_grid().size(); i++)
    {
      auto& focal_gridcell = e.get_grid()[i];
      auto v_food_fluxes = get_neighbors_food_fluxes(focal_gridcell, e);
      auto av_food_flux =
          std::accumulate(v_food_fluxes.begin(),v_food_fluxes.end(),0.0) /
          (v_food_fluxes.empty() ? 1 : v_food_fluxes.size());

      auto in_out_flux_food = calc_food_flux(focal_gridcell, av_food_flux, e.get_diff_coeff());
      focal_gridcell.set_food_change(in_out_flux_food);
    }
}

void calc_diffusion_metab(environment &e) noexcept
{
  for(size_t i = 0; i != e.get_grid().size(); i++)
    {
      auto v_metab_fluxes = get_neighbors_metab_fluxes(e.get_grid()[i], e);
      auto av_metab_flux =
          std::accumulate(v_metab_fluxes.begin(),v_metab_fluxes.end(),0.0) /
          (v_metab_fluxes.empty() ? 1 : v_metab_fluxes.size());

      auto in_out_flux_metab = calc_metab_flux(e.get_grid()[i], av_metab_flux, e.get_diff_coeff());
      e.get_grid()[i].set_metab_change(in_out_flux_metab);
    }
}

void diffusion(environment& e) noexcept
{
  calc_diffusion(e);
  apply_diffusion(e);
}

void apply_diffusion(environment& e) noexcept
{
  for(auto & cell : e.get_grid())
    {
      cell.increment_food();
      cell.increment_metabolite();
    }
}

void diff_food(environment& e) noexcept
{
  calc_diffusion_food(e);
  for(auto & cell : e.get_grid())
    {
      cell.increment_food();
    }
}

void diff_metab(environment& e) noexcept
{
  calc_diffusion_metab(e);
  for(auto & cell : e.get_grid())
    {
      cell.increment_metabolite();
    }
}


std::vector<int> find_neighbors(int grid_size, int grid_side, int index) noexcept
{
  std::vector<int> n_indexes;
  for(int column = -1; column != 2; column++)
    {
      if(is_over_sides(index, grid_side, column)){continue;}
      for(int row = -1; row != 2; row++ )
        {
          if(is_past_limit(index, grid_side, grid_size, row)){continue;}
          else if(is_same_cell(column, row)){continue;}

          int neighbor_index = index + grid_side * row + column;
          n_indexes.push_back(neighbor_index);

        }
    }
  return n_indexes;
}

void find_neighbors_all_grid(environment& e) noexcept
{
  for (int i = 0; i != e.get_grid_size(); i++)
    {
      e.get_cell(i).set_v_neighbors(find_neighbors(e.get_grid_size(), e.get_grid_side(), i));
    }
}

double food_balance(const std::vector<env_grid_cell>& lhs, const std::vector<env_grid_cell>& rhs)
{
  auto grid_food = std::accumulate(
        lhs.begin(),lhs.end(),0.0,
        [](double sum, const env_grid_cell& c) {return sum + c.get_food();});
  auto grid_new_food = std::accumulate(
        rhs.begin(), rhs.end(),0.0,
        [](double sum, const env_grid_cell& c) {return sum + c.get_food();}
  );

  return grid_food - grid_new_food;
}

std::vector<double> get_neighbors_food_fluxes(const env_grid_cell& c, const environment &e) noexcept
{
  std::vector<double> v_food_deltas;
  for(auto neighbor : c.get_v_neighbors())
    {
      v_food_deltas.push_back(food_flux(c, e.get_cell(neighbor)));
    }
  return v_food_deltas;
}

std::vector<double> get_neighbors_metab_fluxes(const env_grid_cell& c,
                                               const environment &e) noexcept
{
  std::vector<double> v_metab_deltas;
  for(auto neighbor : c.get_v_neighbors())
    {
      v_metab_deltas.push_back(metab_flux(c, e.get_cell(neighbor)));
    }
  return v_metab_deltas;
}

bool is_over_sides(int index, int grid_side, int column) noexcept
{
  return (index % grid_side == 0 && column == -1) ||
      (index % grid_side == grid_side - 1 && column == 1);
}

bool is_past_limit(int index, int grid_side, int grid_size, int row) noexcept
{
  return (index - grid_side < 0 && row == -1) ||(index + grid_side >= grid_size && row == 1);
}

bool is_same_cell(int column, int row) noexcept
{
  return column == 0 && row == 0;
}


void reset_env(environment& e)
{
  e = environment(e.get_grid_side(), e.get_diff_coeff(), e.get_init_food());
}

void reset_env(environment& e, int grid_side, double diff_coeff, double food)
{
  e = environment(grid_side, diff_coeff, food);
}

void test_environment()//!OCLINT tests may be many
{
#ifndef NDEBUG

  //An enevironment has a grid of cells
  //The .size() call and the -123456789
  //value are irrelevant is just to see if the command
  //is called properly
  {
    environment e;
    assert(static_cast<int>(e.get_grid().size()) > -123456789);
  }


  //An environment has a grid with  n*n elements
  //where n is the size of the side of the square
  //constituting the environment
  {
    int n = 2;
    environment e(n);
    assert(e.get_grid_size() == 2 * 2);
  }

  //One can get the gride side's size
  {
    int side = 42;
    environment e(side);
    assert(e.get_grid_side() == side );
  }

  //An environment is initialized by default with
  //all his gridcells with the same amount of food
  {
    environment g;
    for (auto& cell : g.get_grid())
      {
        for(auto& comparison_cell :g.get_grid())
          {
            assert(cell.get_food() - comparison_cell.get_food() < 0.00001);
          }
      }
  }

  //An environment is initialized by default with
  //all his gridcells with the 0 metabolite
  {
    environment g;
    for (auto& cell : g.get_grid())
      {
        assert(cell.get_metab() < 0.000001
               || cell.get_metab() > -0.000001);
      }
  }

  //An environment is initialized by default with
  //all his gridcells with 1 food
  {
    environment g;
    for (auto& cell : g.get_grid())
      {
        assert(cell.get_food() - 1 < 0.000001
               || cell.get_food() -1 > -0.000001);
      }
  }

  //An environment is initialized with a specific diffusion coefficient
  //same for food and metabolite (this might change in the future)
  //0 by default
  {
    environment e;
    assert(e.get_diff_coeff() < 0.000001
           || e.get_diff_coeff() > -0.000001);

    double diff_coeff = 3.14;
    environment e1(1, diff_coeff);
    assert(e1.get_diff_coeff() - diff_coeff < 0.0001);
  }

  //On initializtion each grid_cells has its neighbor vector
  //updated with the indexes of its neighbors
  {
    auto env_side = 3;
    auto env_size = 3 * 3;
    environment e3x3(3);//3x3 env just to have some actual neighbours
    for(int i = 0; i != e3x3.get_grid_size(); i++)
      {
        assert(e3x3.get_cell(i).get_v_neighbors()
               == find_neighbors(env_size, env_side, i));
      }
  }
  //The environment can get a specific a cell's address
  //The -123456789 value is irrelevant, just to assure the function is called properly
  {
    environment e;
    assert(e.get_cell(0).get_food() > -1234546789);
  }

  //The environment can modify information for
  //cell in the grid
  {
    environment e;
    double food = 0.42;
    double metabolite = 3.14;
    int target_cell = 0;

    e.get_cell(target_cell).set_food(food);
    assert(e.get_cell(target_cell).get_food() - food < 0.000001);

    e.get_cell(target_cell).set_metabolite(metabolite);
    assert(e.get_cell(target_cell).get_metab() - metabolite < 0.000001);
  }


  //find_neigbours returns a vector of indexes of other grid_cells
  {
    environment e;
    //the focal cell in this case is the first and only cell of the grid
    auto focal_cell_index = 0;
    assert(find_neighbors(e.get_grid_size(),e.get_grid_side(),focal_cell_index).empty());

    //the central grid_cell in a 3x3 grid (index = 4) should have 8 neighbors
    environment e3x3(3);
    focal_cell_index = 4;
    assert(
          find_neighbors
          (
            e3x3.get_grid_size(),
            e3x3.get_grid_side(),
            focal_cell_index
            ).size()
          == 8
          );
  }

  //A cell will lose due to diffusion a proportion of his food equal to
  //The average delta_food with its neighbors * the  diffusion coefficient *
  //the number of neighbors
  //If this proportion is bigger than 1 it will lose all its food
  {
    //Check for case in which exiting food >= cell_food
    //And therefore exiting food = cell_food
    double diffusion_coeff = 1;
    double av_diff = - 1;//Three neighbours each with 1 food less than focal cell
    environment e(2,diffusion_coeff);
    auto c = e.get_cell(0);
    //Check for case in which exiting food < cell_food
    //-> exiting food = cell_food * av_difference * diffusion_coeff * neighbors_
    diffusion_coeff = 0.1;
    auto predicted_flux = av_diff * diffusion_coeff * c.get_food() * c.get_v_neighbors().size();
    assert(calc_food_flux(c,av_diff,diffusion_coeff) -
           predicted_flux < 0.000001 &&
           calc_food_flux(c,av_diff,diffusion_coeff) -
           predicted_flux > -0.000001);
  }

  //A cell will lose due to diffusion a proportion of his metab equal to
  //The average delta_metab with its neighbors * the  diffusion coefficient *
  //the number of neighbors
  //If this proportion is bigger than 1 it will lose all its metab
  {
    //Check for case in which exiting metab >= cell_metab
    //And therefore exiting metab = cell_metab
    double diffusion_coeff = 1;
    double av_diff = - 1;//Three neighbours each with 1 metab less than focal cell
    environment e(2,diffusion_coeff);
    auto c = e.get_cell(0);
    //Check for case in which exiting metab < cell_metab
    //-> exiting metab = cell_metab * av_difference * diffusion_coeff * neighbors_
    diffusion_coeff = 0.1;
    auto predicted_flux = av_diff * diffusion_coeff * c.get_metab() * c.get_v_neighbors().size();
    assert(calc_metab_flux(c,av_diff,diffusion_coeff) -
           predicted_flux < 0.000001 &&
           calc_metab_flux(c,av_diff,diffusion_coeff) -
           predicted_flux > -0.000001);
  }

  //A cell with more substance than its neighbours diffuses food to them
  {
    environment e(3, 0.1, 0);
    auto focal_cell_index = 4;
    e.get_cell(focal_cell_index).set_food(1);
    e.get_cell(focal_cell_index).set_metabolite(1);

    //register strting values
    std::vector<double> v_orig_food_val;
    for(auto cell : e.get_grid()) {v_orig_food_val.push_back(cell.get_food());}
    std::vector<double> v_orig_metab_val;
    for(const auto& cell : e.get_grid())
      {v_orig_metab_val.push_back(cell.get_metab());}

    diffusion(e);

    //register new values
    std::vector<double> v_new_food_values;
    for(const auto& cell : e.get_grid()) {v_new_food_values.push_back(cell.get_food());}
    std::vector<double> v_new_metab_val;
    for(const auto& cell : e.get_grid()) {v_new_metab_val.push_back(cell.get_metab());}

    //Check they are different
    assert(v_orig_food_val != v_new_food_values);
    assert(v_orig_metab_val != v_new_metab_val);

    //Check they diminish or augment as expected
    for(unsigned int i = 0; i != e.get_grid().size(); i++)
      {
        if(static_cast<int>(i) == focal_cell_index)
          {assert(v_new_food_values[i] < v_orig_food_val[i]);}
        else {assert(v_new_food_values[i] > v_orig_food_val[i]);}
      }
    for(unsigned int i = 0; i != e.get_grid().size(); i++)
      {
        if(static_cast<int>(i) == focal_cell_index)
          {assert(v_new_metab_val[i] < v_orig_metab_val[i]);}
        else {assert(v_new_metab_val[i] > v_orig_metab_val[i]);}
      }
  }

  //A cell with less food than its neighbours  receives food
  {
    environment e(3, 0.1, 4);
    auto focal_cell_index = 4;
    e.get_cell(focal_cell_index).set_food(1);

    //register starting values
    std::vector<double> v_orig_food_val;
    for(auto cell : e.get_grid()) {v_orig_food_val.push_back(cell.get_food());}
    auto tot_food_init = std::accumulate(v_orig_food_val.begin(), v_orig_food_val.end(), 0.0);

    diff_food(e);
    auto new_grid = e.get_grid();

    //register new values
    std::vector<double> v_new_food_values;
    v_new_food_values.reserve(new_grid.size());
    for(auto cell : new_grid) {v_new_food_values.push_back(cell.get_food());}
    auto tot_food_new = std::accumulate(v_new_food_values.begin(), v_new_food_values.end(), 0.0);

    //Check they are different
    assert(v_orig_food_val != v_new_food_values);
    //Check the total amount of food does not change
    assert(tot_food_new -  tot_food_init < 0.00001 &&
           tot_food_new -  tot_food_init > -0.00001);

    //Check they diminish or augment as expected
    for(unsigned int i = 0; i != new_grid.size(); i++)
      {
        if(static_cast<int>(i) == focal_cell_index)
          {assert(v_new_food_values[i] > v_orig_food_val[i]);}
        else {assert(v_new_food_values[i] < v_orig_food_val[i]);}
      }
  }

  //A cell with less metab than its neighbours  receives metab
  {
    environment e(3, 0.1);
    //Let's add metabolite to all grid_cells
    double metab_value = 4;
    for(auto& cell: e.get_grid()){cell.set_metabolite(metab_value);}
    //Let's have one cell with less metabolite
    auto focal_cell_index = 4;
    e.get_cell(focal_cell_index).set_metabolite(1);

    //register starting values
    std::vector<double> v_orig_metab_val;
    for(const auto& cell : e.get_grid()) {v_orig_metab_val.push_back(cell.get_metab());}
    auto tot_metab_init = std::accumulate(v_orig_metab_val.begin(), v_orig_metab_val.end(), 0.0);

    diff_metab(e);
    auto new_grid = e.get_grid();

    //register new values
    std::vector<double> v_new_metab_values;
    v_new_metab_values.reserve(new_grid.size());
    for(const auto& cell : new_grid) {v_new_metab_values.push_back(cell.get_metab());}
    auto tot_metab_new = std::accumulate(v_new_metab_values.begin(), v_new_metab_values.end(), 0.0);

    //Check they are different
    assert(v_orig_metab_val != v_new_metab_values);
    //Check the total amount of metab does not change
    assert(tot_metab_new -  tot_metab_init < 0.00001 &&
           tot_metab_new -  tot_metab_init > -0.00001);

    //Check they diminish or augment as expected
    for(unsigned int i = 0; i != new_grid.size(); i++)
      {
        if(static_cast<int>(i) == focal_cell_index)
          {assert(v_new_metab_values[i] > v_orig_metab_val[i]);}
        else {assert(v_new_metab_values[i] < v_orig_metab_val[i]);}
      }
  }

  //Diffusion does not change the total amount of substance present in a system
  {
    //for food
    environment e(3, 0.1, 4);
    auto focal_cell_index = 4;
    e.get_cell(focal_cell_index).set_food(1);
    e.get_cell(focal_cell_index).set_metabolite(1);

    auto tot_food_bef = std::accumulate(
          e.get_grid().begin(), e.get_grid().end(), 0.0,
          [](double sum, const env_grid_cell& c){return sum + c.get_food();}
    );

    auto tot_metab_bef = std::accumulate(
          e.get_grid().begin(), e.get_grid().end(), 0.0,
          [](double sum, const env_grid_cell& c){return sum + c.get_metab();}
    );
    diffusion(e);

    auto tot_food_after = std::accumulate(
          e.get_grid().begin(), e.get_grid().end(), 0.0,
          [](double sum, const env_grid_cell& c){return sum + c.get_food();}
    );

    auto tot_metab_after = std::accumulate(
          e.get_grid().begin(), e.get_grid().end(), 0.0,
          [](double sum, const env_grid_cell& c){return sum + c.get_metab();}
    );

    assert(tot_food_bef - tot_food_after < 0.0001 &&
           tot_food_bef - tot_food_after > -0.0001);
    assert(tot_metab_bef - tot_metab_after < 0.0001 &&
           tot_metab_bef - tot_metab_after > -0.0001);
  }

  //environment has a boolean operator (it checks if two grids have all
  //the same cells)
  {
    environment e;
    auto e1 = e;
    environment e2(10);
    environment e3(1,0.2);
    assert(e == e1 );
    assert(e != e2 && e != e3 && e2 != e3);
    e1.get_cell(0).set_food(42);
    assert(e != e1);
  }

  // The environment can be reset to its initial state
  {
    environment e;
    double food_quantity = 3.14;
    double metabolite_quantity = 3.14;
    for(auto& grid_cell : e.get_grid())
      {
        grid_cell.set_food(food_quantity);
        grid_cell.set_metabolite(metabolite_quantity);
      }
    auto e2 = e;
    reset_env(e2);
    assert(e2 != e);
    assert(e2.get_grid_side() == e.get_grid_side());
    assert(e2.get_diff_coeff() - e.get_diff_coeff() < 0.00001 &&
           e2.get_diff_coeff() - e.get_diff_coeff() > -0.000001);
    for (const auto& grid_cell : e2.get_grid())
      {
        assert(grid_cell.get_food() - e.get_init_food() < 0.00001 &&
               grid_cell.get_food() - e.get_init_food() > -0.00001 );
        assert(grid_cell.get_metab() < 0.000001 && grid_cell.get_metab() > -0.000001);
      }
  }

  //The environment can be reset to a new env ( no metabolite)
  //with set size, diffusion coefficent and food
  {
    auto init_grid_side = 2;
    auto init_diff_coeff = 2;
    auto init_food = 2;
    environment e(init_grid_side, init_diff_coeff, init_food);
    environment e2 = e;

    auto grid_side = 42;
    auto diff_coeff = 42;
    auto food = 42;
    reset_env(e2, grid_side, diff_coeff, food);
    assert( e != e2);
    assert(e2.get_grid_side() == grid_side);
    assert(e2.get_diff_coeff() - diff_coeff < 0.00001 &&
           e2.get_diff_coeff() - diff_coeff > -0.000001);
    assert(e2.get_init_food() - food < 0.00001 &&
           e2.get_init_food() - food > -0.000001);
  }

#endif
}
