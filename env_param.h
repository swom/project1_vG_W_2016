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
              double mean_diff_coeff = 0.1,
              double mean_degr_rate = 0.1,
              double var_diff_coeff = 0.02,
              double var_degr_coeff = 0.02);


    ///Returns the diffusion coefficient
    double get_diff_coeff() const noexcept {return  m_diff_coeff;}

    ///Returns the side of the grid, used to determine the total size of the grid(square)
    int get_grid_side() const noexcept {return m_grid_side;}

    ///Gets the initial food that will be provided in each grid_cell of the environment
    double get_init_food() const noexcept {return m_init_food;}

    ///Gets the degradation coefficient
    double get_degr_rate() const noexcept {return m_metab_degradation_rate;}

    ///Returns the mean of the distribution of possible diffusion coefficent values
    double get_mean_diff_coeff() const noexcept {return m_mean_diff_coeff;}

    ///Returns the mean of the distribution of possible degradation rates
    double get_mean_degr_rate() const noexcept {return m_mean_degr_rate;}

    ///Returns the variance of the distribution of possible diffusion coefficent values
    double get_var_diff_coeff() const noexcept {return m_var_diff_coeff;}

    ///Returns the variance of the distribution of possible degradation coefficent values
    double get_var_degr_rate() const noexcept {return m_var_degr_rate;}

    ///Sets the value of the diffusion coefficent
    void set_diff_coeff(double diff_coeff) noexcept {m_diff_coeff = diff_coeff;}

    ///Sets the value of the diffusion coefficent
    void set_metab_degr(double metab_degr) noexcept {m_metab_degradation_rate = metab_degr;}

private:

    /// The current diffusion coefficient of substances in the grid
    double m_diff_coeff;

    ///The side of the grid
    int m_grid_side;

    ///The initial amount of food in each grid_cell at the start of a pop cycle
    double m_init_food;

    ///The rate at which metabolite degrades
    double m_metab_degradation_rate;

    ///The mean value of the distribution(normal) of possible values for
    /// the diffusion coefficient
     double m_mean_diff_coeff;

    ///The mean value of the distribution(normal)
    /// degradation rate of the metabolite
     double m_mean_degr_rate;

    ///The variance of the distribution(normal) of possible values
    /// for the diffusion coefficient
     double m_var_diff_coeff;

    ///The variance of the distribution(normal) of possible values
    /// for the degradation rate
     double m_var_degr_rate;
};

bool operator==(const env_param& lhs, const env_param& rhs) noexcept;

bool operator!=( const env_param& lhs, const env_param& rhs) noexcept;

std::ostream& operator<<(std::ostream& os, const env_param& p);

std::ifstream& operator>>(std::ifstream& is, env_param& p);

///Returns a new env_param that is a changed version of the
/// given env_param
/// whith new values drawn from a normal distribution
/// with mean and variance as indicated by the m_mean* and m_var* members
env_param change_env_param_norm(const env_param& e, std::minstd_rand &rng) noexcept;

///Returns a new env_param that is a changed version of the
/// given env_param
/// whith new values drawn from a uniform distribution
/// with mean as indicated by the m_mean* members
/// and range = m_mean* -/+ 3 * m_var* members
env_param change_env_param_unif(const env_param& e, std::minstd_rand& rng) noexcept;


env_param load_env_parameters( const std::string& filename);

void save_env_parameters( const env_param& p, const std::string& filename);

void test_env_param() noexcept;

#endif // ENV_PARAM_H
