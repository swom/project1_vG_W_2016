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
              double metab_degrad_rate = 0.01,
              double min_step_env_change = 0,
              double range_env_change = 0);


    ///Returns the diffusion coefficient
    double get_diff_coeff() const noexcept {return  m_diff_coeff;}

    ///Returns the side of the grid, used to determine the total size of the grid(square)
    int get_grid_side() const noexcept {return m_grid_side;}

    ///Gets the initial food that will be provided in each grid_cell of the environment
    double get_init_food() const noexcept {return m_init_food;}

    ///Gets the degradation coefficient
    const double& get_degr_rate() const noexcept {return m_metab_degradation_rate;}

    ///Returns the range of environmental change
    double get_range_env_change() const noexcept {return m_range_env_change;}

    ///Returns the magnitude of environmental change
    double get_min_step_env_change() const noexcept {return m_min_step_env_change;}

private:

    /// The diffusion coefficient of substances in the grid
    double m_diff_coeff;

    ///The side of the grid
    int m_grid_side;

    ///The initial amount of food in each grid_cell at the start of a pop cycle
    double m_init_food;

    ///The rate at which metabolite degrades
    double m_metab_degradation_rate;

    ///The minimum step of change from previous value for
    /// changing environmental paramters
    double m_min_step_env_change;

    ///The maximum distance from the original
    ///initialization values to which env parameters
    /// can vary
    double m_range_env_change;
};

bool operator==(const env_param& lhs, const env_param& rhs) noexcept;

bool operator!=( const env_param& lhs, const env_param& rhs) noexcept;

std::ostream& operator<<(std::ostream& os, const env_param& p);

std::ifstream& operator>>(std::ifstream& is, env_param& p);

///Returns a new env_param that is a changed version of the
/// variable env_param modified following the magnitude and range values
env_param change_env_param_incr(const env_param& e) noexcept;

env_param load_env_parameters( const std::string& filename);

void save_env_parameters( const env_param& p, const std::string& filename);

void test_env_param() noexcept;

#endif // ENV_PARAM_H
