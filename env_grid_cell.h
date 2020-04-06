#ifndef ENVIRONMENT_GRID_CELL_H
#define ENVIRONMENT_GRID_CELL_H
#include "vector"

class env_grid_cell
{
public:
    env_grid_cell(double metabolite = 0, double food = 0, double max_food = 20);

    ///Gets the amount of metabolite in the cell
    const double& get_metabolite() const noexcept {return m_metabolite;}

    ///Sets the amount of metabolite
    void set_metabolite(double m) noexcept {m_metabolite = m;}

    ///Sets the amount of metabolite
    void increment_metabolite(double m) noexcept {m_metabolite + m > 0 ? m_metabolite += m : 0;}

    ///Gets the amount of food in the cell
    const double& get_food() const noexcept {return m_food;}

    ///Gets the reference to the amount of food in the cell
    double& get_food() noexcept {return m_food;}

    ///Gets the reference to the amount of food in the cell
    const double& get_max_food() const noexcept {return m_max_food;}
    ///Sets the amount of food
    void set_food(double f) noexcept {m_food = f;}

    ///Sets the amount of food
    void increment_food(double f) noexcept {m_food + f > 0 ? m_food += f : 0;}

    ///Gets the const reference to the vector of the neighbors
    const std::vector<int>& get_v_neighbors() const noexcept {return m_v_neighbors;}

    ///Gets the change in the amount of food in the next tick ofthe simulation
    const double& get_food_change() const noexcept {return m_food_change;}

    ///Gets the change in the amount of metabolite in the next tick of the simulation
    const double& get_metabolite_change() const noexcept {return m_metabolite_change;}

    ///Sets the change in the amount of food in the next tick of the simulation
    void increment_food_change(double change) noexcept {m_food_change += change;}

    ///Sets the change in the amount of food in the next tick of the simulation
    void increment_metabolite_change(double change) noexcept {m_metabolite_change += change;}

    ///resets the change in the amount of food in the next tick of the simulation
    void reset_food_change() noexcept {m_food_change = 0;}

    ///resets the change in the amount of metabolite in the next tick of the simulation
    void reset_metabolite_change() noexcept {m_metabolite_change = 0;}

    ///Sets the vector of neigbors to a vector of integers
    void set_v_neighbors(std::vector<int> neighbors) noexcept {m_v_neighbors = neighbors;}


private:
    std::vector<int> m_v_neighbors;
    double m_metabolite;
    double m_metabolite_change;
    double m_food;
    double m_food_change;
    double m_max_food;
};

bool operator == (const env_grid_cell& lhs, const env_grid_cell& rhs) noexcept;

bool operator != (const env_grid_cell& lhs, const env_grid_cell& rhs) noexcept;

///Gets the difference in food between two cells
double food_difference(const env_grid_cell& lhs, const env_grid_cell& rhs)  noexcept;

///Given the vector of food differences with the neighbors and the diffusion coefficient
/// Calculates how much food the cell will give away
double calc_exiting_food(env_grid_cell& cell, std::vector<double> v_food_diff, double diffusion_coeff) noexcept;

///Gets the difference in metabolite between two cells
double metab_difference(const env_grid_cell &lhs, const env_grid_cell &rhs) noexcept;

///Given the vector of food differences with the neighbors and the diffusion coefficient
/// Calculates how much metabolite the cell will give away
double calc_exiting_metabolite(env_grid_cell& cell, std::vector<double> v_metabolite_diff, double diffusion_coeff) noexcept;

void test_env_grid_cell();


#endif // ENVIRONMENT_GRID_CELL_H
