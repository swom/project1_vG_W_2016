#include "env_grid_cell.h"
#include <cassert>

env_grid_cell::env_grid_cell(double metabolite, double food):
    m_metabolite(metabolite),
    m_food(food)
{

}

void test_env_grid_cell()
{
    //A grid cell has a metabolite value,
    //The value -123456789 is irrelevant
    //is just to make sure that get_metabolite()
    //is called properly
    {
        env_grid_cell g;
        assert( g.get_metabolite() > -123456789);
    }

    //Concentration of metabolite can be set
    {
        env_grid_cell g;
        double new_metabolite_conc = 3;
        assert(g.get_metabolite() != new_metabolite_conc);
        g.set_metabolite(new_metabolite_conc);
        assert(g.get_metabolite() - new_metabolite_conc < 0.00001);
    }

    //Concentration of metabolite can be changed
    {
        env_grid_cell g;
        env_grid_cell g1 = g;
        double change_in_metabolite = 3;
        assert(g.get_metabolite() != change_in_metabolite);
        g1.change_metabolite(change_in_metabolite);
        assert(g.get_metabolite() != g1.get_metabolite());
        assert(abs(g.get_metabolite() - g1.get_metabolite()) - change_in_metabolite < 0.000001);

    }

    //Metabolite concentration cannot be less than 0
    {
        env_grid_cell g;
        g.change_metabolite(-10000);
        assert(g.get_metabolite() >= 0);
    }

    //A grid cell has a food value,
    //The value -123456789 is irrelevant
    //is just to make sure that get_food()
    //is called properly
    {
        env_grid_cell g;
        assert( g.get_food() > -123456789);
    }

    //Concentration of food can be set
    {
        env_grid_cell g;
        double new_food_conc = 3;
        assert(g.get_food() != new_food_conc);
        g.set_food(new_food_conc);
        assert(g.get_food() - new_food_conc < 0.00001);
    }

    //Concentration of food can be changed
    {
        env_grid_cell g;
        env_grid_cell g1 = g;
        double change_in_food = 3;
        assert(g.get_food() != change_in_food);
        g1.change_food(change_in_food);
        assert(g.get_food() != g1.get_food());
        assert(abs(g.get_food() - g1.get_food()) - change_in_food < 0.000001);

    }

    //Food concentration cannot be less than 0
    {
        env_grid_cell g;
        g.change_food(-10000);
        assert(g.get_food() >= 0);
    }
}
