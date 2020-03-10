#ifndef ENVIRONMENT_GRID_CELL_H
#define ENVIRONMENT_GRID_CELL_H


class env_grid_cell
{
public:
    env_grid_cell(double metabolite = 0);

    ///Gets the amount of metabolite in the cell
    const double& get_metabolite()  noexcept {return m_metabolite;}

    ///Sets the amount of metabolite
    void set_metabolite(double m) noexcept {m_metabolite = m;}

private:
    double m_metabolite;
};

void test_env_grid_cell();

#endif // ENVIRONMENT_GRID_CELL_H
