#ifndef ENV_PARAM_H
#define ENV_PARAM_H
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>

class env_param
{
public:
    env_param(int grid_side = 1,
              double diff_coeff = 0.1,
              double init_food = 10.0,
              double metab_degrad_rate = 0.1,
              double min_change_fraction = 0, // the nth fraction of the range which is going to be the minimum step
              double range_env_change = 0.05);


    ///Returns the diffusion coefficient
    double get_diff_coeff() const noexcept {return  m_diff_coeff;}

    ///Returns the side of the grid, used to determine the total size of the grid(square)
    int get_grid_side() const noexcept {return m_grid_side;}

    ///Gets the initial food that will be provided in each grid_cell of the environment
    double get_init_food() const noexcept {return m_init_food;}

    ///Gets the degradation coefficient
    const double& get_degr_rate() const noexcept {return m_metab_degradation_rate;}

    ///Returns const ref to the range of diff_coeff values
    const std::uniform_real_distribution<double>& get_range_diff_coeff_change() const noexcept
    {return m_range_diff_coeff_change;}

    ///Returns ref to the range of diff_coeff values
    std::uniform_real_distribution<double>& get_range_diff_coeff_change() noexcept
    {return m_range_diff_coeff_change;}

    ///Returns const ref to the range of metabolite degradation values
    const std::uniform_real_distribution<double>& get_range_metab_degr_change() const noexcept
    {return m_range_metab_degr_change;}

    ///Returns const ref to the range of metabolite degradation values
    std::uniform_real_distribution<double>& get_range_metab_degr_change() noexcept
    {return m_range_metab_degr_change;}

    ///Returns the magnitude of diff change
    double get_min_step_diff_change() const noexcept {return m_step_min_diff_change;}

    ///Returns the magnitude of degr change
    double get_min_step_degr_change() const noexcept {return m_step_min_degr_change;}

    ///Sets the value of the diffusion coefficent
    void set_diff_coeff(double diff_coeff) noexcept {m_diff_coeff = diff_coeff;}

    ///Sets the value of the diffusion coefficent
    void set_metab_degr(double metab_degr) noexcept {m_metab_degradation_rate = metab_degr;}

private:

    /// The diffusion coefficient of substances in the grid
    double m_diff_coeff;

    ///The side of the grid
    int m_grid_side;

    ///The initial amount of food in each grid_cell at the start of a pop cycle
    double m_init_food;

    ///The rate at which metabolite degrades
    double m_metab_degradation_rate;

    ///The distribution from which the value of m_diff_coeff is taken
    /// time the environment changes
    std::uniform_real_distribution<double> m_range_diff_coeff_change;

    ///The distribution from which the value of m_metab_degradation_rate is taken
    /// time the environment changes
    std::uniform_real_distribution<double> m_range_metab_degr_change;

    ///The minimum step of change from previous value for
    /// changing diffusion paramters
    double m_step_min_diff_change;

    ///The minimum step of change from previous value for
    /// changing degradation paramters
    double m_step_min_degr_change;
};

bool operator==(const env_param& lhs, const env_param& rhs) noexcept;

bool operator!=( const env_param& lhs, const env_param& rhs) noexcept;

std::ostream& operator<<(std::ostream& os, const env_param& p);

std::ifstream& operator>>(std::ifstream& is, env_param& p);

///Returns a new env_param that is a changed version of the
/// variable env_param modified following the magnitude and range values
env_param change_env_param_unif(const env_param& e, std::minstd_rand &rng) noexcept;

env_param load_env_parameters( const std::string& filename);

void save_env_parameters( const env_param& p, const std::string& filename);

void test_env_param() noexcept;

#endif // ENV_PARAM_H
