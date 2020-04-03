#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H
#include "vector"
#include "env_grid_cell.h"

class environment
{
public:
    environment(int grid_side = 1, double diff_coeff = 0, double init_food = 1);

    ///Gets const reference of cell at a certain index
    const env_grid_cell& get_cell(int i) const noexcept
    {
      return m_grid[static_cast<unsigned int>(i)];
    }

    ///Gets reference of cell at a certain index
    env_grid_cell& get_cell(int i) noexcept
    {
      return m_grid[static_cast<unsigned int>(i)];
    }

    ///Get the diffusion coefficent of the environment
    const double& get_diff_coeff() const noexcept {return m_diffusion_coefficient;}

    ///Gets the size of the grid of cells(number of grid cells)
    int get_grid_size() const noexcept {return static_cast<int>(m_grid.size());}

    ///Gets the size of the grid of cells(number of grid cells)
    int get_grid_side() const noexcept {return m_side;}

    ///Gets the entire grid by const reference
    const std::vector<env_grid_cell>& get_grid() const noexcept {return m_grid;}

    ///Gets the entire grid by reference
    std::vector<env_grid_cell>& get_grid() noexcept {return m_grid;}

    ///Gets the amount of initial food in the grid
    double get_init_food() const noexcept {return  m_init_food;}


private:
    int m_side;
    //Always needs to be between 0 and 1!!!
    double m_diffusion_coefficient;
    double m_init_food;
    std::vector<env_grid_cell> m_grid;
};
///Checks if two environment have the same grid_size and the the same amount of food and metabolite in each cell
bool operator == (const environment& lhs, const environment& rhs) noexcept;

///Checks if two environment have the same grid_size and the the same amount of food and metabolite in each cell
bool operator != (const environment& lhs, const environment& rhs) noexcept;

///Finds the indexes of the neighboring grid_cell of a grid_cell at a certain index
///Assumes that the grid is a square
const std::vector<int> find_neighbors(int grid_size, int grid_side, int index) noexcept;

///finds neighbors for all cells in grid
void find_neighbors_all_grid(environment& e) noexcept;

///Diffuses substances in environment
void diffusion(environment& e) noexcept;

///diffuses metabolite in grid
void diffusion_metabolite(environment& e) noexcept;

///diffuses food in grid
void diffusion_food(environment& e) noexcept;

///Checks that cells do not check for neighbors over the side of the grid
bool is_over_sides(int index, int grid_side, int column)  noexcept;

///Checks that cells do not check for neighbors past the limits of the grid
bool is_past_limit(int index, int grid_side, int grid_size, int row) noexcept;

///Checks cells do not look for neighbors at their own coordinates
bool is_same_cell(int column, int row) noexcept;

///Resets the environment to new parameters
void reset_env(environment& e, int grid_side, double diff_coeff, double food);

///Resets the environment to new parameters
void reset_env(environment& e);

void test_environment();

#endif // ENVIRONMENT_H
