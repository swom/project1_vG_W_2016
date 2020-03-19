#include "environment.h"
#include <cassert>
#include <numeric>

environment::environment(int grid_side, double diff_coeff, double food):
  m_grid(static_cast<unsigned int>(grid_side * grid_side),env_grid_cell(0,food)),
  m_side(grid_side),
  m_diffusion_coefficient(diff_coeff)

{
  find_neighbors_all_grid(*this);
}

void diffusion(environment& e) noexcept
{
  diffusion_food(e);
  diffusion_metabolite(e);
}

void diffusion_food(environment& e) noexcept
{
  for(auto& cell : e.get_grid())
    {
      std::vector<double> v_food_diff;
      for(const auto& neighbor : cell.get_v_neighbors())
        {
          auto food_diff = cell.get_food() - e.get_cell(neighbor).get_food();
          food_diff > 0 ? v_food_diff.push_back(food_diff) : v_food_diff.push_back(0) ;
        }
      auto total_diff = std::accumulate(v_food_diff.begin(), v_food_diff.end(), 0);
      auto max_diff = cell.get_food() * cell.get_v_neighbors().size();

      auto exiting_food = cell.get_food() *
          (max_diff > 0 ? total_diff / max_diff : 0) *
          e.get_diff_coeff();

      cell.increment_food_change(exiting_food > 0 ? -exiting_food : 0);

      for(unsigned int i = 0; i != v_food_diff.size(); i++)
        {
          auto recieved_food = exiting_food * (total_diff > 0 ? v_food_diff[i]/total_diff : 0);
          e.get_cell(cell.get_v_neighbors()[i]).increment_food_change(recieved_food);
        }
    }

  for(auto& cell : e.get_grid())
    {
      cell.increment_food(cell.get_food_change());
      cell.reset_food_change();
    }

}

void diffusion_metabolite(environment& e) noexcept
{
  for(auto& cell : e.get_grid())
    {
      std::vector<double> v_metabolite_diff;
      for(auto neighbor : cell.get_v_neighbors())
        {
          auto metabolite_diff = cell.get_metabolite() - e.get_cell(neighbor).get_metabolite();
          metabolite_diff > 0 ? v_metabolite_diff.push_back(metabolite_diff)
                              : v_metabolite_diff.push_back(0) ;
        }
      auto total_diff = std::accumulate(v_metabolite_diff.begin(), v_metabolite_diff.end(), 0);
      auto max_diff = cell.get_metabolite() * cell.get_v_neighbors().size();

      auto exiting_metabolite = cell.get_metabolite() *
          (max_diff > 0 ? total_diff / max_diff : 0) *
          e.get_diff_coeff();

      cell.increment_metabolite_change(exiting_metabolite > 0 ? -exiting_metabolite : 0);

      for(unsigned int i = 0; i != v_metabolite_diff.size(); i++)
        {
          auto recieved_metabolite = exiting_metabolite *
              (total_diff > 0 ? v_metabolite_diff[i]/total_diff : 0);
          e.get_cell(cell.get_v_neighbors()[i]).increment_metabolite_change(recieved_metabolite);
        }
    }

  for(auto& cell : e.get_grid())
    {
      cell.increment_metabolite(cell.get_metabolite_change());
      cell.reset_metabolite_change();
    }

}

void find_neighbors_all_grid(environment& e) noexcept
{
  for (int i = 0; i != e.get_env_size(); i++)
    {
      e.get_cell(i).set_v_neighbors(find_neighbors(e.get_env_size(), e.get_grid_side(), i));
    }
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

void test_environment()//!OCLINT tests may be many
{

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
    assert(e.get_env_size() == 2 * 2);
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
    for(int i = 0; i != e3x3.get_env_size(); i++)
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
    assert(find_neighbors(e.get_env_size(),e.get_grid_side(),focal_cell_index).size() == 0);

    //the central grid_cell in a 3x3 grid (index = 4) should have 8 neighbors
    environment e3x3(3);
    focal_cell_index = 4;
    assert(
          find_neighbors
          (
            e3x3.get_env_size(),
            e3x3.get_grid_side(),
            focal_cell_index
            ).size()
          == 8
          );
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

}
