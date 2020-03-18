#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H
#include "vector"
#include "env_grid_cell.h"

class environment
{
public:
    environment(int grid_side = 1, double diff_coeff = 0);

    ///Gets the size of the grid of cells(number of grid cells)
    int get_env_size() const noexcept {return static_cast<int>(m_grid.size());}

    ///Gets the size of the grid of cells(number of grid cells)
    int get_grid_side() const noexcept {return m_side;}

    ///Gets the entire grid by const reference
    const std::vector<env_grid_cell>& get_grid() const noexcept {return m_grid;}

    ///Gets the entire grid by reference
    std::vector<env_grid_cell>& get_grid() noexcept {return m_grid;}

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

private:
    std::vector<env_grid_cell> m_grid;
    int m_side;
    double m_diffusion_coefficient;
};

///Finds the indexes of the neighboring grid_cell of a grid_cell at a certain index
///Assumes that the grid is a square
const std::vector<int> find_neighbors(int grid_size, int grid_side, int index) noexcept;

///finds neighbors for all cells in grid
void find_neighbors_all_grid(environment& e) noexcept;

///diffuses metabolite in grid
void diffusion_metabolite(environment& e) noexcept;

///diffuses food in grid
void diffusion_food(environment& e) noexcept;

void test_environment();

#endif // ENVIRONMENT_H
