#include "env_grid_cell.h"
#include <cassert>
#include <cmath>
#include <numeric>


env_grid_cell::env_grid_cell(double metabolite, double food, double max_food, double max_metab):
    m_metabolite{metabolite},
    m_food{food},
    m_max_food{max_food},
    m_max_metab{max_metab}
{

}

bool operator == (const env_grid_cell& lhs, const env_grid_cell& rhs) noexcept
{
    return lhs.get_food() - rhs.get_food() < 0.00001 &&
            lhs.get_food() - rhs.get_food() > -0.00001 &&
            lhs.get_metab() - rhs.get_metab() < 0.00001 &&
            lhs.get_metab() - rhs.get_metab() > -0.00001;
}

bool operator != (const env_grid_cell& lhs, const env_grid_cell& rhs) noexcept
{
    return !(lhs == rhs);
}

double calc_food_flux( const env_grid_cell& cell, double av_food_delta,
                       double diffusion_coeff) noexcept
{
    auto in_out_flux = cell.get_v_neighbors().size() * av_food_delta * diffusion_coeff;
    if(cell.get_food() + in_out_flux < 0)
    {
        return - cell.get_food();
    }
    return in_out_flux;
}


double calc_metab_flux( const env_grid_cell& cell, double av_metab_flux,
                        double diffusion_coeff) noexcept
{
    auto in_out_flux = cell.get_v_neighbors().size() * av_metab_flux * diffusion_coeff;
    if(cell.get_metab() + in_out_flux < 0)
    {
        return - cell.get_metab();
    }
    return in_out_flux;
}

void env_grid_cell::increment_metabolite() noexcept
{
    if(m_metabolite + m_metabolite_change > -0.0000000001)
        m_metabolite += m_metabolite_change;
    else if(m_metabolite + m_metabolite_change <= -0.0000000001)
        m_metabolite = 0;
    else if(m_metabolite + m_metabolite_change > m_max_metab)
        m_metabolite = m_max_metab;

    m_metabolite_change = 0;
}

double food_flux(const env_grid_cell &lhs, const env_grid_cell &rhs) noexcept
{
    return  rhs.get_food() - lhs.get_food();
}

double metab_flux(const env_grid_cell &lhs, const env_grid_cell &rhs) noexcept
{
    return rhs.get_metab()- lhs.get_metab();
}

void metabolite_degrades(env_grid_cell& g, double degrad_rate) noexcept
{
    auto new_metab = g.get_metab() - degrad_rate;
    if(new_metab < 0)
    {
        g.set_metab(0);
        return;
    }
    g.set_metab(new_metab);
}

