#include "environment.h"
#include <cassert>

environment::environment(int grid_side):
    m_grid(grid_side*grid_side),
    m_side(grid_side)

{

}

const std::vector<int> find_neighbors(int grid_size, int grid_side, int index) noexcept
{

    //vector where indexes of neighbours will be stored
    std::vector<int> n_indexes;

    //find index of neighbor grid_cells
    //considering env is a square of side = grid_side

        //checks columns to left/same/right starting from left
        for(int column = -1; column != 2; column++)
        {
            //if on the left-most cell of the row
            //do not check cells on the left
            if(index % grid_side == 0 && column == -1)
            {
                    continue;
            }

            //if on the right most column of the row
            //do not check cells on the right
            if(index % grid_side == grid_side - 1 && column == 1)
            {
                    continue;
            }

                //checks row up/same/bottom starting from up
                for(int row = -1; row != 2; row++ )
                {
                    //cells does not count itself as a neighbor
                    if(row == 0 && column == 0)
                    {
                        continue;
                    }
                    //if on bottom-most or top-most row do not check rows below or above
                    if(index - grid_side < 0 || index + grid_side >= grid_size )
                    {
                        continue;
                    }
                    else
                    {
                        int neighbor_index = index + grid_side * row + column;
                        n_indexes.push_back(neighbor_index);
                    }
                }
        }

        return n_indexes;
}


//void diffusion(environment& e) noexcept
//{

//}



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

    //One can get the gride side's size
    {
        int side = 42;
        environment e(side);
        assert(e.get_grid_side() == side );
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


    //The environment can get a specific a cell's address
    //The -123456789 value is irrelevant, just to assure the function is called properly
    {
        environment e;
        assert(e.get_cell(0).get_food() > -1234546789);
    }

    //The environment can modify information for
    //cell in the grid
    {
        environment e;
        double food = 0.42;
        double metabolite = 3.14;
        int target_cell = 0;

        e.get_cell(target_cell).set_food(food);
        assert(e.get_cell(target_cell).get_food() - food < 0.000001);

        e.get_cell(target_cell).set_metabolite(metabolite);
        assert(e.get_cell(target_cell).get_metabolite() - metabolite < 0.000001);
    }


    //find_neigbours returns a vector of indexes of other grid_cells
    {
        environment e;
        //the focal cell in this case is the first and only cell of the grid
        auto focal_cell_index = 0;
        assert(find_neighbors(e.get_env_size(),e.get_grid_side(), focal_cell_index).size() == 0);

        //the central grid_cell in a 3x3 grid (index = 4) should have 8 neighbors
        environment e3x3(3);
        focal_cell_index = 4;
        assert(find_neighbors(e3x3.get_env_size(),e3x3.get_grid_side(), focal_cell_index).size() == 8);
    }


    //Diffusion: a grid_cell transfers part of his food
    //to neighbouring grid_cells with less food

//    {
//        environment e(9);
//        //central cell has more food
//        int central_cell = 4;
//        double food_increment = 8;
//        e.get_cell(central_cell).increment_food(food_increment);
//        e.get_cell(central_cell).get_food();

//        std::vector<double> original_grid_food;
//        for(auto& cell : e.get_grid())
//        {
//            original_grid_food.push_back(cell.get_food());
//        }

//        //Maybe put this in test above
//        for(int i = 0; i != e.get_env_size(); i++)
//        {
//            for(int j = 0; j != e.get_env_size(); j++)
//            {
//                if(i != central_cell)
//                {
//                    if(j != central_cell)
//                    {
//                        assert(e.get_grid()[i].get_food() - e.get_grid()[j].get_food() < 0.00001);
//                    }
//                    else
//                    {
//                        assert(e.get_grid()[i].get_food() - e.get_grid()[j].get_food() < 0);
//                    }
//                }
//                else
//                {
//                    if(j != central_cell)
//                    {
//                        assert(e.get_grid()[j].get_food() - e.get_grid()[i].get_food() < 0 );
//                    }
//                    else
//                    {
//                        assert(e.get_grid()[i].get_food() - e.get_grid()[j].get_food() < 0.000001);
//                    }
//                }

//            }
//        }

//        diffusion(e);

//        for(int i = 0; i != e.get_env_size(); i++)
//        {
//            if(i != central_cell)
//            {
//                assert(e.get_cell(i).get_food() > original_grid_food[i]);
//            }
//            else
//            {
//                assert(e.get_cell(i).get_food() < original_grid_food[i]);
//            }
//        }
//    }


    //A grid_cell transfers part of his metabolite
    //to neighbouring grid_cells with less metabolite


//    {
//        environment e(9);
//        //central cell has more food
//        int central_cell = 4;
//        double food_increment = 8;
//        e.get_cell(central_cell).increment_metabolite(food_increment);
//        e.get_cell(central_cell).get_metabolite();

//        std::vector<double> original_grid_metabolite;
//        for(auto& cell : e.get_grid())
//        {
//            original_grid_metabolite.push_back(cell.get_metabolite());
//        }

//        //Maybe put this in test above
//        for(int i = 0; i != e.get_env_size(); i++)
//        {
//            for(int j = 0; j != e.get_env_size(); j++)
//            {
//                if(i != central_cell)
//                {
//                    if(j != central_cell)
//                    {
//                        assert(e.get_grid()[i].get_metabolite() - e.get_grid()[j].get_metabolite() < 0.00001);
//                    }
//                    else
//                    {
//                        assert(e.get_grid()[i].get_metabolite() - e.get_grid()[j].get_metabolite() < 0);
//                    }
//                }
//                else
//                {
//                    if(j != central_cell)
//                    {
//                        assert(e.get_grid()[j].get_metabolite() - e.get_grid()[i].get_metabolite() < 0);
//                    }
//                    else
//                    {
//                        assert(e.get_grid()[i].get_metabolite() - e.get_grid()[j].get_metabolite() < 0.000001);
//                    }
//                }

//            }
//        }

//        diffusion(e);

//        for(int i = 0; i != e.get_env_size(); i++)
//        {
//            if(i != central_cell)
//            {
//                assert(e.get_cell(i).get_metabolite() > original_grid_metabolite[i]);
//            }
//            else
//            {
//                assert(e.get_cell(i).get_metabolite() < original_grid_metabolite[i]);
//            }
//        }
//    }

}
