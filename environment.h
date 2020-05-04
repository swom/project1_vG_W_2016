#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H
#include "vector"
#include "env_grid_cell.h"
#include "env_param.h"

class environment
{
public:

    environment(int grid_side = 1, double diff_coeff = 0, double init_food = 1, double metab_deg_rate = 0.001);

    environment(env_param env_parameters);

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

    ///Gets the size of the grid of cells(number of grid cells)
    int get_grid_size() const noexcept {return static_cast<int>(m_grid.size());}

    ///Gets the entire grid by const reference
    const std::vector<env_grid_cell>& get_grid() const noexcept {return m_grid;}

    ///Gets the entire grid by reference
    std::vector<env_grid_cell>& get_grid() noexcept {return m_grid;}

    ///Returns constant reference to the parameters
    const env_param get_param() const noexcept {return m_env_param;}

private:
    env_param m_env_param;
    std::vector<env_grid_cell> m_grid;
};
///Checks if two environment have the same grid_size and the the same amount of food and metabolite in each cell
bool operator == (const environment& lhs, const environment& rhs) noexcept;

///Checks if two environment have the same grid_size and the the same amount of food and metabolite in each cell
bool operator != (const environment& lhs, const environment& rhs) noexcept;

///Applies the changes in substance concntration in all cells do to diffusion, the claculations are
/// done by calc_diffusion
void apply_diffusion(environment& e) noexcept;

///Calculates the change of food for a given gridcell at a certain index in a given environment
void calc_change_food(environment& e, int index_focal_cell) noexcept;

///Calculates the change of metab for a given gridcell at a certain index in a given environment
void calc_change_metab(environment& e, int index_focal_cell) noexcept;

///Calculates the amount of metabolite that needs to diffuse
void calc_diffusion_metab(environment &e) noexcept;

///Calculates the amoutn of substances that diffuse in or out of each cell
void calc_diffusion(environment& e) noexcept;

///Calculates the amount of food that needs to diffuse
void calc_diffusion_food(environment &e) noexcept;

///Metabolite degrades in each grid_cell
void degradation_metabolite(environment& e) noexcept;

///Diffuses substances in environment
void diffusion(environment& e) noexcept;

///Diffuses food
void diff_food(environment& e) noexcept;

///Diffuses metabolite
void diff_metab(environment& e) noexcept;

///Finds the indexes of the neighboring grid_cell of a grid_cell at a certain index
///Assumes that the grid is a square
std::vector<int> find_neighbors(int grid_size, int grid_side, int index) noexcept;

///finds neighbors for all cells in grid
void find_neighbors_all_grid(environment& e) noexcept;

///Finds the difference in food between one focal cell and other cells(normally intended to be its neighbors)
std::vector<double> get_neighbors_food_fluxes(const env_grid_cell& c, const environment& e) noexcept;

///Finds the difference in metabolite between one focal cell and other cells(normally intended to be its neighbors)
std::vector<double> get_neighbors_metab_fluxes(const env_grid_cell& c, const environment& e) noexcept;

///Checks that cells do not check for neighbors over the side of the grid
bool is_over_sides(int index, int grid_side, int column)  noexcept;

///Checks that cells do not check for neighbors past the limits of the grid
bool is_past_limit(int index, int grid_side, int grid_size, int row) noexcept;

///Checks cells do not look for neighbors at their own coordinates
bool is_same_cell(int column, int row) noexcept;

///Resets the environment to new parameters
void reset_env(environment& e, int grid_side, double diff_coeff, double food);

///Resets the environment to new parameters
void reset_env(environment& e);

void test_environment();

#endif // ENVIRONMENT_H

