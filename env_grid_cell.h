#ifndef ENVIRONMENT_GRID_CELL_H
#define ENVIRONMENT_GRID_CELL_H


class env_grid_cell
{
public:
    env_grid_cell(double metabolite = 0);

    ///Gets the amount of metabolite in the cell
    double get_metabolite() const noexcept {return m_metabolite;}

private:
    double m_metabolite;
};

void test_env_grid_cell();

#endif // ENVIRONMENT_GRID_CELL_H
