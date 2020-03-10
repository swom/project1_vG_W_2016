#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H
#include <vector>
#include "env_grid_cell.h"

class environment
{
public:
    environment(int n = 1);

private:
    std::vector<env_grid_cell> m_grid;
};

void test_environment();

#endif // ENVIRONMENT_H
