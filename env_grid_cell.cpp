#include "env_grid_cell.h"
#include <assert.h>

env_grid_cell::env_grid_cell(double metabolite):
    m_metabolite(metabolite)
{

}

void test_env_grid_cell()
{
    //A grid cell has a signal value,
    //The value -123456789 is irrelevant
    //is just to make sure that get_signal()
    //is called properly
    {
        env_grid_cell g;
        assert( g.get_metabolite() > -123456789);
    }

    //A grid cell concentration of metabolite can be set
    {
        env_grid_cell g;
        double new_metabolite_conc = 3;
        assert(g.get_metabolite() != new_metabolite_conc);
        g.set_metabolite(new_metabolite_conc);
        assert(g.get_metabolite() - new_metabolite_conc < 0.00001);
    }

}
