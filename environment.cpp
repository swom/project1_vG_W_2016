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

void diffusion(environment& e) noexcept
{
  calc_diffusion_food(e);
  calc_diffusion_metab(e);
  redistribute_substances(e);
}

void calc_diffusion_food(environment& e) noexcept
{
  for(auto& cell : e.get_grid())
    {
      auto v_food_deltas = get_neighbors_food_deltas(cell, e);
      auto exiting_food = calc_exiting_food(cell, v_food_deltas, e.get_diff_coeff());
      cell.increment_food_change(- exiting_food);

      std::vector<double> v_recieved_foods;
      for(unsigned int i = 0; i != cell.get_v_neighbors().size(); i++)
        {
          auto total_diff = cell.get_tot_food_delta();
          auto recieved_food = exiting_food * ( total_diff > 0 ? v_food_deltas[i]/total_diff : 0);
          v_recieved_foods.push_back(recieved_food);
          e.get_cell(cell.get_v_neighbors()[i]).increment_food_change(recieved_food);
        }
      auto tot_recieved_food = std::accumulate(v_recieved_foods.begin(), v_recieved_foods.end(), 0.0);
      assert(tot_recieved_food - exiting_food < 0.000001 &&
             tot_recieved_food - exiting_food > -0.000001);
    }
}

void calc_diffusion_metab(environment& e) noexcept
{
  for(auto& cell : e.get_grid())
    {
      std::vector<double> v_metabolite_diff;
      for(auto neighbor : cell.get_v_neighbors())
        {
          v_metabolite_diff.push_back(metab_difference(cell, e.get_cell(neighbor)));
        }
      auto total_diff = std::accumulate(v_metabolite_diff.begin(), v_metabolite_diff.end(), 0);
      auto exiting_metabolite = total_diff * e.get_diff_coeff() >= 1 ?
            cell.get_metabolite() :
            cell.get_metabolite() * total_diff * e.get_diff_coeff();

      cell.increment_metabolite_change(exiting_metabolite > 0 ? -exiting_metabolite : 0);

      for(unsigned int i = 0; i != v_metabolite_diff.size(); i++)
        {
          auto recieved_metabolite = exiting_metabolite *
              (total_diff > 0 ? v_metabolite_diff[i]/total_diff : 0);
          e.get_cell(cell.get_v_neighbors()[i]).increment_metabolite_change(recieved_metabolite);
        }
    }
}

void find_neighbors_all_grid(environment& e) noexcept
{
  for (int i = 0; i != e.get_grid_size(); i++)
    {
      e.get_cell(i).set_v_neighbors(find_neighbors(e.get_grid_size(), e.get_grid_side(), i));
    }
}

std::vector<double> get_neighbors_food_deltas(const env_grid_cell& c, const environment &e) noexcept
{
  std::vector<double> v_food_deltas;
  for(auto neighbor : c.get_v_neighbors())
    {
      v_food_deltas.push_back(food_difference(c, e.get_cell(neighbor)));
    }
  return v_food_deltas;
}

std::vector<double> get_neighbors_metab_deltas(const env_grid_cell& c, const environment &e) noexcept
{
  std::vector<double> v_metab_deltas;
  for(auto neighbor : c.get_v_neighbors())
    {
      v_metab_deltas.push_back(metab_difference(c, e.get_cell(neighbor)));
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

const std::vector<int> find_neighbors(int grid_size, int grid_side, int index) noexcept
{
  std::vector<int> n_indexes;
  for(int column = -1; column != 2; column++)
    {
      if(is_over_sides(index, grid_side, column)){continue;}
      else {
          for(int row = -1; row != 2; row++ )
            {
              if(is_past_limit(index, grid_side, grid_size, row)){continue;}
              else if(is_same_cell(column, row)){continue;}
              else
                {
                  int neighbor_index = index + grid_side * row + column;
                  n_indexes.push_back(neighbor_index);
                }
            }
        }
    }
  return n_indexes;
}

void redistribute_substances(environment& e) noexcept
{
  for(auto& cell : e.get_grid())
    {
      cell.increment_food(cell.get_food_change());
      cell.reset_food_change();

      cell.increment_metabolite(cell.get_metabolite_change());
      cell.reset_metabolite_change();
    }
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
        assert(cell.get_metabolite() < 0.000001
               || cell.get_metabolite() > -0.000001);
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
    assert(e.get_cell(target_cell).get_metabolite() - metabolite < 0.000001);
  }


  //find_neigbours returns a vector of indexes of other grid_cells
  {
    environment e;
    //the focal cell in this case is the first and only cell of the grid
    auto focal_cell_index = 0;
    assert(find_neighbors(e.get_grid_size(),e.get_grid_side(),focal_cell_index).size() == 0);

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


  //It is possible to return a vector containing the difference in food between a cell and its neighbors
  //If difference is negative than it is equaled to 0
  {

    //Case where all neighbors have more food
    double init_food = 1;
    auto cell_food = 0;
    environment e{3, 0, init_food};
    auto c = &e.get_cell(0);
    c->set_food(cell_food);
    auto v_food_deltas = get_neighbors_food_deltas(e.get_cell(0), e);
    for(size_t i = 0; i != v_food_deltas.size(); i++)
      assert(v_food_deltas[i] < 0.00001 &&
             v_food_deltas[i] > -0.00001);

    //Case where all neighbors have less food
    init_food = 0;
    cell_food = 1;
    e = environment{3, 0, init_food};
    c = &e.get_cell(0);
    c->set_food(cell_food);
    v_food_deltas = get_neighbors_food_deltas(*c, e);
    for(size_t i = 0; i != v_food_deltas.size(); i++)
      assert(v_food_deltas[i] - food_difference(*c, e.get_cell(c->get_v_neighbors()[i])) < 0.00001 &&
             v_food_deltas[i] - food_difference(*c, e.get_cell(c->get_v_neighbors()[i])) > -0.00001);

    //Case where some neighbors have less food and some more
    double more_food = 2;
    e = environment{3, 0, init_food};
    c = &e.get_cell(0);
    c->set_food(cell_food);
    auto c_more_food = c + 1;
    c_more_food->set_food(more_food);
    assert(c->get_food() < c_more_food->get_food());
    v_food_deltas = get_neighbors_food_deltas(*c, e);
    for(size_t i = 0; i != v_food_deltas.size(); i++)
      {
        auto food_diff = food_difference(*c,e.get_cell(c->get_v_neighbors()[i]));
        assert(v_food_deltas[i] - food_diff < 0.00001 &&
               v_food_deltas[i] - food_diff > -0.00001);

      }
  }

  //It is possible to return a vector containing the difference in metabolite between a cell and its neighbors
  //If difference is negative than it is equaled to 0
  {

    //Case where all neighbors have more metabolite
    double init_metabolite = 1;
    auto cell_metabolite = 0;
    environment e{3, init_metabolite, 0};
    auto c = &e.get_cell(0);
    c->set_metabolite(cell_metabolite);
    auto v_metabolite_deltas = get_neighbors_metab_deltas(e.get_cell(0), e);
    for(size_t i = 0; i != v_metabolite_deltas.size(); i++)
      assert(v_metabolite_deltas[i] < 0.00001 &&
             v_metabolite_deltas[i] > -0.00001);

    //Case where all neighbors have less metabolite
    init_metabolite = 0;
    cell_metabolite = 1;
    e = environment{3, init_metabolite, 0};
    c = &e.get_cell(0);
    c->set_metabolite(cell_metabolite);
    v_metabolite_deltas = get_neighbors_metab_deltas(*c, e);
    for(size_t i = 0; i != v_metabolite_deltas.size(); i++)
      assert(v_metabolite_deltas[i] - metab_difference(*c, e.get_cell(c->get_v_neighbors()[i])) < 0.00001 &&
             v_metabolite_deltas[i] - metab_difference(*c, e.get_cell(c->get_v_neighbors()[i])) > -0.00001);

    //Case where some neighbors have less metabolite and some more
    double more_metabolite = 2;
    e = environment{3, init_metabolite, 0};
    c = &e.get_cell(0);
    c->set_metabolite(cell_metabolite);
    auto c_more_metabolite = c + 1;
    c_more_metabolite->set_metabolite(more_metabolite);
    assert(c->get_metabolite() < c_more_metabolite->get_metabolite());
    v_metabolite_deltas = get_neighbors_metab_deltas(*c, e);
    for(size_t i = 0; i != v_metabolite_deltas.size(); i++)
      {
        auto metabolite_diff = metab_difference(*c,e.get_cell(c->get_v_neighbors()[i]));
        assert(v_metabolite_deltas[i] - metabolite_diff < 0.00001 &&
               v_metabolite_deltas[i] - metabolite_diff > -0.00001);

      }
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
    for(auto cell : e.get_grid())
      {v_orig_metab_val.push_back(cell.get_metabolite());}

    diffusion(e);

    //register new values
    std::vector<double> v_new_food_values;
    for(auto cell : e.get_grid()) {v_new_food_values.push_back(cell.get_food());}
    std::vector<double> v_new_metab_val;
    for(auto cell : e.get_grid()) {v_new_metab_val.push_back(cell.get_metabolite());}

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

  //A cell with less substance than its neighbours  receives substances
  {
    environment e(3, 0.1, 4);
    auto focal_cell_index = 4;
    e.get_cell(focal_cell_index).set_food(1);

    //register strting values
    std::vector<double> v_orig_food_val;
    for(auto cell : e.get_grid()) {v_orig_food_val.push_back(cell.get_food());}

    calc_diffusion_food(e);
    redistribute_substances(e);

    //register new values
    std::vector<double> v_new_food_values;
    for(auto cell : e.get_grid()) {v_new_food_values.push_back(cell.get_food());}

    //Check they are different
    assert(v_orig_food_val != v_new_food_values);

    //Check they diminish or augment as expected
    for(unsigned int i = 0; i != e.get_grid().size(); i++)
      {
        if(static_cast<int>(i) == focal_cell_index)
          {assert(v_new_food_values[i] > v_orig_food_val[i]);}
        else {assert(v_new_food_values[i] < v_orig_food_val[i]);}
      }
  }

  //Diffusion does not change the total amount of substanec present in a system
  {
    environment e(3, 0.1, 4);
    auto focal_cell_index = 4;
    e.get_cell(focal_cell_index).set_food(1);
    auto before_tot_food = std::accumulate(e.get_grid().begin(),e.get_grid().end(),0.0,
                                           [](int sum, const env_grid_cell & g){return sum + g.get_food();});
    calc_diffusion_food(e);

    auto after_tot_food = std::accumulate(e.get_grid().begin(),e.get_grid().end(),0.0,
                                          [](int sum, const env_grid_cell & g){return sum + g.get_food();});
    assert(before_tot_food - after_tot_food < 0.000001 &&
           before_tot_food - after_tot_food > 0.0000001);
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
        assert(grid_cell.get_metabolite() < 0.000001 && grid_cell.get_metabolite() > -0.000001);
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
