#ifndef ENVIRONMENT_GRID_CELL_H
#define ENVIRONMENT_GRID_CELL_H
#include <vector>

class env_grid_cell
{
public:
    env_grid_cell(double metabolite = 0, double food = 0);

    ///Gets the amount of metabolite in the cell
    const double& get_metabolite() const noexcept {return m_metabolite;}

    ///Sets the amount of metabolite
    void set_metabolite(double m) noexcept {m_metabolite = m;}

    ///Sets the amount of metabolite
    void increment_metabolite(double m) noexcept {m_metabolite + m > 0 ? m_metabolite += m : 0;}

    ///Gets the amount of metabolite in the cell
    const double& get_food() const noexcept {return m_food;}

    ///Sets the amount of metabolite
    void set_food(double f) noexcept {m_food = f;}

    ///Sets the amount of metabolite
    void increment_food(double f) noexcept {m_food + f > 0 ? m_food += f : 0;}

    ///Gets the const reference to the vector of the neighbors
    const std::vector<int>& get_v_neighbors() const noexcept {return m_v_neighbors;}

    ///Sets the vector of neigbors to a vector of integers
    void set_v_neighbors(std::vector<int> neighbors) noexcept {m_v_neighbors = neighbors;}

private:
    std::vector<int> m_v_neighbors;
    double m_metabolite;
    double m_food;
};

void test_env_grid_cell();

#endif // ENVIRONMENT_GRID_CELL_H
