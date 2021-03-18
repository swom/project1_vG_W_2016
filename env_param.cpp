#include "env_param.h"
#include "json.hpp"
#include <cassert>

env_param::env_param(int grid_side,
                     double diff_coeff,
                     double init_food,
                     double metab_degrad_rate):
    m_diff_coeff{diff_coeff},
    m_grid_side{grid_side},
    m_init_food{init_food},
    m_metab_degradation_rate{metab_degrad_rate}
{
    assert(m_diff_coeff > -0.000000000001 &&
           m_diff_coeff < 1.000000000001);
    assert(m_grid_side >= 0);
    assert(m_init_food > -0.000000000001);
    assert(m_metab_degradation_rate > -0.000000001 &&
           m_metab_degradation_rate < 1.000000001);
}



std::ostream& operator<<(std::ostream& os, const env_param& e)
{

    os << e.get_grid_side() << " , "
       << e.get_init_food() << " , "
       << e.get_diff_coeff()  << " , "
       << e.get_degr_rate();
    return os;
}

std::ifstream& operator >>(std::ifstream& is, env_param& e)
{

    int grid_side;
    double diff_coeff;
    double init_food;
    double metab_degrad_rate;
    double mean_diff_coeff;
    double mean_degr_rate;
    double var_diff_coeff;
    double var_degr_coeff;

    std::string dummy; // To remove the annotation in the file
    is >>
            grid_side >> dummy >>
            init_food >> dummy >>
            diff_coeff >>dummy >>
            metab_degrad_rate >> dummy >>
            mean_diff_coeff >> dummy >>
            mean_degr_rate >> dummy >>
            var_diff_coeff >> dummy >>
            var_degr_coeff;

    e = env_param {
            grid_side,
            diff_coeff,
            init_food,
            metab_degrad_rate}
            ;

    return is;
}

bool operator==(const env_param& lhs, const env_param& rhs) noexcept
{
    auto grid = (lhs.get_grid_side() == rhs.get_grid_side());
    auto food = (lhs.get_init_food() - rhs.get_init_food() < 0.0001
                 && lhs.get_init_food() - rhs.get_init_food() > -0.0001);
    auto degr = (lhs.get_degr_rate() - rhs.get_degr_rate() < 0.0001
                 && lhs.get_degr_rate() - rhs.get_degr_rate() > -0.0001);
    auto diff =  (lhs.get_diff_coeff() - rhs.get_diff_coeff() < 0.0001
                  && lhs.get_diff_coeff() - rhs.get_diff_coeff() > -0.0001);

    return grid && food && degr && diff
            ;
}

bool operator!=( const env_param& lhs, const env_param& rhs) noexcept
{
    return !(lhs == rhs);
}

env_param change_env_param_norm(env_param e, std::minstd_rand& rng) noexcept
{

    auto new_diff_coeff = std::normal_distribution<double>{e.get_mean_diff_coeff(),
            e.get_var_diff_coeff()}(rng);

    e.set_diff_coeff(new_diff_coeff);

    auto new_degr_rate = std::normal_distribution<double>{e.get_mean_degr_rate(),
            e.get_var_degr_rate()}(rng);

    e.set_metab_degr(new_degr_rate);

    return e;
}

env_param change_env_param_unif(env_param e, std::minstd_rand& rng) noexcept
{
    auto new_diff_coeff = std::uniform_real_distribution<double>{e.get_mean_diff_coeff() - 3 * e.get_var_diff_coeff(),
            e.get_mean_diff_coeff() + 3 * e.get_var_diff_coeff()}(rng);

    e.set_diff_coeff(new_diff_coeff);

    auto new_degr_rate = std::uniform_real_distribution<double>{e.get_mean_degr_rate() - 3 * e.get_var_degr_rate(),
            e.get_mean_degr_rate() + 3 * e.get_var_degr_rate()}(rng);

    e.set_metab_degr(new_degr_rate);

    return e;
}

env_param change_range_env_param(const env_param& e, double amplitude)
{
    env_param env{e.get_grid_side(),
                e.get_diff_coeff(),
                e.get_init_food(),
                e.get_degr_rate()}
    ;
    return env;
}

void save_env_parameters(
        const env_param& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    f << p;
}

void save_env_parameters_json(
        const env_param& p,
        const std::string& filename
        )
{
    std::ofstream os(filename);
    nlohmann::json json_out;
    json_out = p;
    os << json_out;
}

env_param load_env_parameters_json(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    env_param e;
    nlohmann::json json_in;
    f >> json_in;

    return e = json_in;
}

env_param load_env_parameters(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    env_param e;
    f >> e;

    return e;
}

