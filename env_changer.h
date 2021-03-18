#ifndef ENV_CHANGER_H
#define ENV_CHANGER_H
#include "env_param.h"
#include "utilities.h"

class env_changer
{
public:
    env_changer(env_param e = env_param{},
                double amplitude = 1,
                int seed  = 0,
                double var_degr_rate = 0.02,
                double var_diff_coeff = 0.02);

    NLOHMANN_DEFINE_TYPE_INTRUSIVE(env_changer,
                                   m_amplitude,
                                   m_mean_par,
                                   m_seed,
                                   m_var_degr_rate,
                                   m_var_diff_coeff)

    double get_amplitude() const noexcept {return m_amplitude;}

    const env_param& get_mean_params() const noexcept {return m_mean_par;}

    std::minstd_rand& get_rng() noexcept {return m_rng;}

    double get_var_degr() const noexcept {return m_var_degr_rate;}

    double get_var_diff() const noexcept {return m_var_diff_coeff;}

    int get_seed() const noexcept {return m_seed;}

private:

    ///The general amplitude of changes around the mean
    /// (variance is adjusted accoridingly for each parameter)
    double m_amplitude;

    ///A parameter file containing the means of the parameters values
    env_param m_mean_par;

    ///Random number generator to create different environment
    std::minstd_rand m_rng;

    ///The variance of the degradation rate of metabolite
    double m_var_degr_rate;

    ///The variance of the diffusion coefficient of food and of metabolite
    double m_var_diff_coeff;

    ///The seed of the rng
    int m_seed;

};

std::ostream& operator<<(std::ostream& os, const env_changer& p);

std::ifstream& operator>>(std::ifstream& is, env_changer& p);

bool operator==(const env_changer& lhs, const env_changer& rhs) noexcept;

bool operator!=(const env_changer& lhs, const env_changer& rhs) noexcept;

///Returns a new env_param that is a changed version of the
/// given env_param
/// whith new values drawn from a normal distribution
/// with mean and variance as indicated by the m_mean* and m_var* members
env_param change_env_param_norm(env_changer &e) noexcept;

///Returns a new env_param that is a changed version of the
/// given env_param
/// whith new values drawn from a uniform distribution
/// with mean as indicated by the m_mean* members
/// and range = m_mean* -/+ 3 * m_var* members
env_param change_env_param_unif(env_changer& e) noexcept;

///Changes param as in change_env_param_unif but only keeps extremes of the distribution
env_param change_env_param_unif_extreme(env_changer& e) noexcept;

env_changer load_env_changer(const std::string& filename);

void save_env_changer(const env_changer& p, const std::string& filename);

env_changer load_env_changer_json(const std::string& filename);

void save_env_changer_json(const env_changer& p, const std::string& filename);

#endif // ENV_CHANGER_H
