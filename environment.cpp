#include "environment.h"
#include <cassert>

environment::environment(int n):
    m_grid(n*n)
{

}

void test_environment()
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
            assert(cell.get_metabolite() == 0);
        }
    }


}
