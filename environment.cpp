#include "environment.h"
#include <algorithm>
#include <cassert>
#include <numeric>

environment::environment(env_param env_parameters):
    m_env_param(env_parameters),
    m_grid(
        static_cast<unsigned int>(m_env_param.get_grid_side() * m_env_param.get_grid_side()),
        env_grid_cell(0,m_env_param.get_init_food())
        )
{
    find_neighbors_all_grid(*this);
}


bool operator == (const environment& lhs, const environment& rhs) noexcept
{

    return  lhs.get_param() == rhs.get_param()
            && lhs.get_grid() == rhs.get_grid();
}

bool operator != (const environment& lhs, const environment& rhs) noexcept
{
    return  !(lhs == rhs);
}

void calc_change_food(environment& e, int index_focal_cell) noexcept
{
    auto& focal_gridcell = e.get_cell(index_focal_cell);
    auto v_food_deltas = get_neighbors_food_deltas(focal_gridcell, e);

    double av_food_delta;
    if(v_food_deltas.empty())
        av_food_delta = 0.0;
    else
        av_food_delta = std::accumulate(v_food_deltas.begin(),v_food_deltas.end(),0.0) /
                v_food_deltas.size();

    auto in_out_flux_food =
            calc_food_flux(focal_gridcell,
                           av_food_delta,
                           e.get_param().get_diff_coeff());

    focal_gridcell.set_food_change(in_out_flux_food);
}

void calc_change_metab(environment& e, int index_focal_cell) noexcept
{
    auto& focal_gridcell = e.get_cell(index_focal_cell);
    auto v_metab_fluxes = get_neighbors_metab_fluxes(focal_gridcell, e);
    if(std::none_of(v_metab_fluxes.begin(), v_metab_fluxes.end(),
                    [](double f){return  f < 0.000001 || f > -0.000001;}))
    {
        assert(focal_gridcell.get_metab_change() < 0.00000001 &&
               focal_gridcell.get_metab_change() > -0.00000001);
        return;
    }
    auto av_metab_flux =
            std::accumulate(v_metab_fluxes.begin(),v_metab_fluxes.end(),0.0) /
            (v_metab_fluxes.empty() ? 1 : v_metab_fluxes.size());

    auto in_out_flux_metab =
            calc_metab_flux(focal_gridcell,
                            av_metab_flux,
                            e.get_param().get_diff_coeff());

    focal_gridcell.set_metab_change(in_out_flux_metab);
}

void calc_diffusion(environment& e) noexcept
{
    for(int i = 0; i != e.get_grid_size(); i++)
    {
        calc_change_food(e, i);
        calc_change_metab(e, i);
    }
}

void calc_diffusion_food(environment& e) noexcept
{
    for(size_t i = 0; i != e.get_grid().size(); i++)
    {
        auto& focal_gridcell = e.get_grid()[i];
        auto v_food_fluxes = get_neighbors_food_deltas(focal_gridcell, e);
        auto av_food_flux =
                std::accumulate(v_food_fluxes.begin(),v_food_fluxes.end(),0.0) /
                (v_food_fluxes.empty() ? 1 : v_food_fluxes.size());

        auto in_out_flux_food =
                calc_food_flux(focal_gridcell,
                               av_food_flux,
                               e.get_param().get_diff_coeff());

        focal_gridcell.set_food_change(in_out_flux_food);
    }
}

void calc_diffusion_metab(environment &e) noexcept
{
    for(size_t i = 0; i != e.get_grid().size(); i++)
    {
        auto v_metab_fluxes = get_neighbors_metab_fluxes(e.get_grid()[i], e);
        auto av_metab_flux =
                std::accumulate(v_metab_fluxes.begin(),v_metab_fluxes.end(),0.0) /
                (v_metab_fluxes.empty() ? 1 : v_metab_fluxes.size());

        auto in_out_flux_metab = calc_metab_flux(e.get_grid()[i],
                                                 av_metab_flux,
                                                 e.get_param().get_diff_coeff());

        e.get_grid()[i].set_metab_change(in_out_flux_metab);
    }
}

void degradation_metabolite(environment& e) noexcept
{
    for(auto& grid_cell : e.get_grid())
    {
        metabolite_degrades(grid_cell, e.get_param().get_degr_rate());
    }
}

void diffusion(environment& e) noexcept
{
    calc_diffusion(e);
    apply_diffusion(e);
}

void apply_diffusion(environment& e) noexcept
{
    for(auto & cell : e.get_grid())
    {
        cell.increment_food();
        cell.increment_metabolite();
    }
}

void diff_food(environment& e) noexcept
{
    calc_diffusion_food(e);
    for(auto & cell : e.get_grid())
    {
        cell.increment_food();
    }
}

void diff_metab(environment& e) noexcept
{
    calc_diffusion_metab(e);
    for(auto & cell : e.get_grid())
    {
        cell.increment_metabolite();
    }
}


std::vector<int> find_neighbors(int grid_size, int grid_side, int index) noexcept
{
    std::vector<int> n_indexes;
    for(int column = -1; column != 2; column++)
    {
        if(is_over_sides(index, grid_side, column)){continue;}
        for(int row = -1; row != 2; row++ )
        {
            if(is_past_limit(index, grid_side, grid_size, row)){continue;}
            else if(is_same_cell(column, row)){continue;}

            int neighbor_index = index + grid_side * row + column;
            n_indexes.push_back(neighbor_index);

        }
    }
    return n_indexes;
}

void find_neighbors_all_grid(environment& e) noexcept
{
    for (int i = 0; i != e.get_grid_size(); i++)
    {
        e.get_cell(i).set_v_neighbors
                (
                    find_neighbors
                    (e.get_grid_size(),
                     e.get_param().get_grid_side(),
                     i
                     )
                    );
    }
}

double food_balance(const std::vector<env_grid_cell>& lhs, const std::vector<env_grid_cell>& rhs)
{
    auto grid_food = std::accumulate(
                lhs.begin(),lhs.end(),0.0,
                [](double sum, const env_grid_cell& c) {return sum + c.get_food();});
    auto grid_new_food = std::accumulate(
                rhs.begin(), rhs.end(),0.0,
                [](double sum, const env_grid_cell& c) {return sum + c.get_food();}
    );

    return grid_food - grid_new_food;
}

std::vector<double> get_neighbors_food_deltas(const env_grid_cell& c, const environment &e) noexcept
{
    std::vector<double> v_food_deltas;
    for(auto neighbor : c.get_v_neighbors())
    {
        v_food_deltas.push_back(food_flux(c, e.get_cell(neighbor)));
    }
    return v_food_deltas;
}

std::vector<double> get_neighbors_metab_fluxes(const env_grid_cell& c,
                                               const environment &e) noexcept
{
    std::vector<double> v_metab_deltas;
    for(auto neighbor : c.get_v_neighbors())
    {
        v_metab_deltas.push_back(metab_flux(c, e.get_cell(neighbor)));
    }
    return v_metab_deltas;
}

bool is_over_sides(int index, int grid_side, int column) noexcept
{
    return (index % grid_side == 0 && column == -1) ||
            (index % grid_side == grid_side - 1 && column == 1);
}

bool is_past_limit(int index, int grid_side, int grid_size, int row) noexcept
{
    return (index - grid_side < 0 && row == -1) ||(index + grid_side >= grid_size && row == 1);
}

bool is_same_cell(int column, int row) noexcept
{
    return column == 0 && row == 0;
}


void reset_env(environment& e)
{
    e.get_grid() =
            std::vector<env_grid_cell>(
                static_cast<unsigned int>(e.get_grid_size()),
                env_grid_cell{0,e.get_param().get_init_food()});
}

void set_all_cell_food(environment& e, double f)
{
    std::for_each(e.get_grid().begin(),
                  e.get_grid().end(),
                  [=](env_grid_cell& g)
    {g.set_food(f);});
}

void reset_env(environment& e, int grid_side, double diff_coeff, double food)
{
    e = environment(env_param{grid_side,diff_coeff,food});
}
