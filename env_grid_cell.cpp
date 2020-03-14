#include "env_grid_cell.h"
#include <cassert>
#include <cmath>
env_grid_cell::env_grid_cell(double metabolite, double food):
    m_metabolite(metabolite),
    m_metabolite_change(0),//this will always be 0 at initialization
    m_food(food),
    m_food_change(0)//this will always be 0 at initialization

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
        assert(g.get_metabolite() - new_metabolite_conc > 0.0000001
               || g.get_metabolite() - new_metabolite_conc < -0.0000001);
        g.set_metabolite(new_metabolite_conc);
        assert(g.get_metabolite() - new_metabolite_conc < 0.00001);
    }

    //Concentration of metabolite can be changed
    {
        env_grid_cell g;
        env_grid_cell g1 = g;
        double change_in_metabolite = 3;
        assert(g.get_metabolite() - change_in_metabolite > 0.0000001
               || g.get_metabolite() - change_in_metabolite < -0.0000001);
        g1.increment_metabolite(change_in_metabolite);
        assert(g.get_metabolite() - g1.get_metabolite() > 0.000000001
               || g.get_metabolite() - g1.get_metabolite() < -0.000000001);
        assert(abs(g.get_metabolite() - g1.get_metabolite()) - change_in_metabolite < 0.000001);

    }

    //Metabolite concentration cannot be less than 0
    {
        env_grid_cell g;
        g.increment_metabolite(-10000);
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
        assert(g.get_food() - new_food_conc > 0.00000001
               || g.get_food() - new_food_conc < -0.00000001);
        g.set_food(new_food_conc);
        assert(g.get_food() - new_food_conc < 0.00001);
    }

    //Concentration of food can be changed
    {
        env_grid_cell g;
        env_grid_cell g1 = g;
        double change_in_food = 3;
        assert(g.get_food() - change_in_food > 0.0000001
               || g.get_food() - change_in_food < -0.0000001);
        g1.increment_food(change_in_food);
        assert(g.get_food() - g1.get_food() > 0.00000001
               || g.get_food() - g1.get_food() < -0.00000001);
        assert(abs(g.get_food() - g1.get_food()) - change_in_food < 0.000001);

    }

    //Food concentration cannot be less than 0
    {
        env_grid_cell g;
        g.increment_food(-10000);
        assert(g.get_food() >= 0);
    }

    //A grid_cell has a vector that can store
    //the indexes of its neighbouring cells if placed in a grid
    {
        env_grid_cell g;
        assert(g.get_v_neighbors().size() == 0);
    }

    //A grid cell vector of neighbors can be set
    {
        //all neighors will have the same address for the purpose of this test
        int neighbors_address = 42;
        std::vector<int> list_of_neighbors(42,neighbors_address);
        env_grid_cell g;
        g.set_v_neighbors(list_of_neighbors);
        assert(g.get_v_neighbors().size() - list_of_neighbors.size() < 0.0000001);
        for(auto neighbor : g.get_v_neighbors())
        {
        assert(neighbor == neighbors_address);
        }
    }

    //A grid_cell has a member variable that stores
    //the amount of change in food at the next tick
    {
        env_grid_cell g;
        assert(g.get_food_change() < 0.0000001 || g.get_food_change() > -0.0000001);
    }

    //A grid_cell has a member variable that stores
    //the amount of change in metabolite at the next tick
    {
        env_grid_cell g;
        assert(g.get_metabolite_change() < 0.0000001 || g.get_metabolite_change() > -0.0000001);
    }
}
