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

    ///Gets the entire grid by const reference
    const std::vector<env_grid_cell>& get_grid() noexcept {return m_grid;}

    ///Gets const reference of cell at a certain index
    const env_grid_cell& get_cell(int i) const noexcept {return m_grid[i];}

    ///Gets reference of cell at a certain index
    env_grid_cell& get_cell(int i) noexcept {return m_grid[i];}

private:
    std::vector<env_grid_cell> m_grid;
};

void diffuse(environment& e) noexcept;

void test_environment();

#endif // ENVIRONMENT_H
