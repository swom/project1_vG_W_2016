#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H
#include <vector>
#include "env_grid_cell.h"

class environment
{
public:
    environment(int n = 1);

    ///Gets the size of the grid of cells(number of grid cells)
    int get_env_size() const noexcept {return static_cast<int>(m_grid.size());}

    ///Gets the size of the grid of cells(number of grid cells)
    int get_grid_side() const noexcept {return m_side;}

    ///Gets the entire grid by const reference
    const std::vector<env_grid_cell>& get_grid() noexcept {return m_grid;}

    ///Gets const reference of cell at a certain index
    const env_grid_cell& get_cell(int i) const noexcept {return m_grid[i];}

    ///Gets reference of cell at a certain index
    env_grid_cell& get_cell(int i) noexcept {return m_grid[i];}

private:
    std::vector<env_grid_cell> m_grid;
    int m_side;
};

///Finds the indexes of the neighboring grid_cell of a grid_cell at a certain index
const std::vector<int> find_neighbours(int grid_size, int grid_side, int index) noexcept;

void diffusion(environment& e) noexcept;

void test_environment();

#endif // ENVIRONMENT_H
