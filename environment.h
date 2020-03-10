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

private:
    std::vector<env_grid_cell> m_grid;
};

void test_environment();

#endif // ENVIRONMENT_H
