#ifndef ENV_PARAM_H
#define ENV_PARAM_H
#include <iostream>
#include <fstream>
#include <sstream>

class env_param
{
public:
    env_param(int grid_side = 1,
              double diff_coeff = 0.1,
              double init_food = 10.0,
              double metab_degrad_rate = 0.01);


    ///Returns the diffusion coefficient
    double get_diff_coeff() const noexcept {return  m_diff_coeff;}

    ///Returns the side of the grid, used to determine the total size of the grid(square)
    int get_grid_side() const noexcept {return m_grid_side;}

    ///Gets the initial food that will be provided in each grid_cell of the environment
    double get_init_food() const noexcept {return m_init_food;}

    ///Gets the degradation rate of the metabolite
    double get_metab_degr_rate() const noexcept {return m_metab_degradation_rate;}

    ///Gets the degradation coefficient
    const double& get_degr_coeff() const noexcept {return m_metab_degradation_rate;}

    ///Sets the diffusion coefficient

private:

    /// The diffusion coefficient of substances in the grid
    double m_diff_coeff;

    ///The side of the grid
    int m_grid_side;

    ///The initial amount of food in each grid_cell at the start of a pop cycle
    double m_init_food;

    ///The rate at which metabolite degrades
    double m_metab_degradation_rate;
};

bool operator==(const env_param& lhs, const env_param& rhs) noexcept;

std::ostream& operator<<(std::ostream& os, const env_param& p);

std::ifstream& operator>>(std::ifstream& is, env_param& p);

env_param load_env_parameters( const std::string& filename);

void save_env_parameters( const env_param& p, const std::string& filename);

void test_env_param() noexcept;

#endif // ENV_PARAM_H
