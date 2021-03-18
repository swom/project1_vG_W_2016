#ifndef ENV_PARAM_H
#define ENV_PARAM_H
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include "json.hpp"

class env_param
{
public:
    env_param(int grid_side = 1,
              double diff_coeff = 0.1,
              double init_food = 10.0,
              double metab_degrad_rate = 0.1);

    ///For to_json and from_json
    NLOHMANN_DEFINE_TYPE_INTRUSIVE(env_param,
                                   m_grid_side,
                                   m_diff_coeff,
                                   m_init_food,
                                   m_metab_degradation_rate)

    ///Returns the diffusion coefficient
    double get_diff_coeff() const noexcept {return  m_diff_coeff;}

    ///Returns the side of the grid, used to determine the total size of the grid(square)
    int get_grid_side() const noexcept {return m_grid_side;}

    ///Gets the initial food that will be provided in each grid_cell of the environment
    double get_init_food() const noexcept {return m_init_food;}

    ///Gets the degradation coefficient
    double get_degr_rate() const noexcept {return m_metab_degradation_rate;}

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


};

bool operator==(const env_param& lhs, const env_param& rhs) noexcept;

bool operator!=( const env_param& lhs, const env_param& rhs) noexcept;

std::ostream& operator<<(std::ostream& os, const env_param& p);

std::ifstream& operator>>(std::ifstream& is, env_param& p);

env_param load_env_parameters(const std::string& filename);

env_param load_env_parameters_json(const std::string& filename);

void save_env_parameters( const env_param& p, const std::string& filename);

void save_env_parameters_json(const env_param& p, const std::string& filename);

#endif // ENV_PARAM_H
