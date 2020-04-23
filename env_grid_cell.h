#ifndef ENVIRONMENT_GRID_CELL_H
#define ENVIRONMENT_GRID_CELL_H
#include <vector>
#include <numeric>

class env_grid_cell
{
public:
    env_grid_cell(double metabolite = 0, double food = 0, double max_food = 20, double max_metab = 20);

    ///Gets the amount of metabolite in the cell
    const double& get_metab() const noexcept {return m_metabolite;}

    ///Gets the amount of food in the cell
    const double& get_food() const noexcept {return m_food;}

    ///Sets the future change in metabolite concentration
    double get_food_change() noexcept {return m_food_change;}

    ///Sets the future change in metabolite concentration
    double get_metab_change() noexcept {return m_metabolite_change;}

    ///Gets the reference to the amount of food in the cell
    double& get_food() noexcept {return m_food;}

    ///Gets the reference to the amount of food in the cell
    const double& get_max_food() const noexcept {return m_max_food;}

    ///Gets the reference to the amount of food in the cell
    const double& get_max_metab() const noexcept {return m_max_metab;}

    ///Sets the amount of food
    void increment_food() noexcept {m_food + m_food_change > -0.0000001 ? m_food += m_food_change : m_food = 0;
                                    m_food_change = 0;}

    ///Sets the amount of metabolite
    void increment_metabolite() noexcept {m_metabolite + m_metabolite_change > -0.0000000001 ? m_metabolite += m_metabolite_change : m_metabolite = 0;
                                         m_metabolite_change = 0;}

    ///Sets the amount of food
    void set_food(double f) noexcept {m_food = f;}

    ///Sets the amount of metabolite
    void set_metab(double m) noexcept {m_metabolite = m;}

    ///Sets the future change in metabolite concentration
    void set_food_change(double f_change) noexcept {m_food_change = f_change;}

    ///Sets the future change in metabolite concentration
    void set_metab_change(double m_change) noexcept {m_metabolite_change = m_change;}

    ///Gets the const reference to the vector of the neighbors
    const std::vector<int>& get_v_neighbors() const noexcept {return m_v_neighbors;}

    ///Sets the vector of neigbors to a vector of integers
    void set_v_neighbors(std::vector<int> neighbors) noexcept {m_v_neighbors = std::move(neighbors);}


private:
    std::vector<int> m_v_neighbors;
    double m_metabolite;
    double m_metabolite_change;
    double m_food;
    double m_food_change;
    double m_max_food;
    double m_max_metab;
};

bool operator == (const env_grid_cell& lhs, const env_grid_cell& rhs) noexcept;

bool operator != (const env_grid_cell& lhs, const env_grid_cell& rhs) noexcept;

///Gets the flux in food between two cells, negative if food goes out from lhs positive it goes in
double food_flux(const env_grid_cell& lhs, const env_grid_cell& rhs)  noexcept;

///Given the vector of food differences with the neighbors and the diffusion coefficient
/// Calculates how much food the cell will give away
double calc_food_flux(const env_grid_cell& cell, double av_food_flux, double diffusion_coeff) noexcept;

///Given the vector of food differences with the neighbors and the diffusion coefficient
/// Calculates how much metabolite the cell will give away
double calc_metab_flux(const env_grid_cell &cell, double av_metab_flux, double diffusion_coeff) noexcept;

///Gets the difference in metabolite between two cells
double metab_flux(const env_grid_cell &lhs, const env_grid_cell &rhs) noexcept;

///Degrades the amount of metabolite in the grid_cell by a value equal to the degradation rate
void metabolite_degrades(env_grid_cell& g, double degrad_rate) noexcept;

void test_env_grid_cell();


#endif // ENVIRONMENT_GRID_CELL_H
