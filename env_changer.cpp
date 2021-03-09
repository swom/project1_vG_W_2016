#include "env_changer.h"

env_changer::env_changer(env_param e, double amplitude, int seed, double var_degr_rate, double var_diff_coeff):
    m_amplitude{amplitude},
    m_mean_par{e},
    m_rng{seed},
    m_var_degr_rate{var_degr_rate},
    m_var_diff_coeff{var_diff_coeff},
    m_seed{seed}
{
}

std::ostream& operator<<(std::ostream& os, const env_changer& ec)
{
    os << ec.get_mean_params() << " , "
<< ec.get_amplitude() << " , "
<< ec.get_seed() << " , "
<< ec.get_var_degr() << " , "
<< ec.get_var_diff();
    return os;
}

std::ifstream& operator>>(std::ifstream& is, env_changer& ec)
{
    std::string dummy;
    env_param e;
    double amplitude;
    int seed;
    double var_degr;
    double var_diff;

    is >> e >> dummy

            >> amplitude >> dummy

            >> seed >> dummy

            >> var_degr >> dummy

            >> var_diff;

    ec = env_changer{e,
            amplitude,
            seed,
            var_degr,
            var_diff};

    return is;
}

bool operator==(const env_changer& lhs, const env_changer& rhs) noexcept
{
    return lhs.get_mean_params() == rhs.get_mean_params() &&
            lhs.get_amplitude() == rhs.get_amplitude() &&
            lhs.get_var_degr() == rhs.get_var_degr() &&
            lhs.get_var_diff() == rhs.get_var_diff() &&
            lhs.get_seed() == rhs.get_seed();
}

bool operator!=(const env_changer& lhs, const env_changer& rhs) noexcept
{
    return !(lhs == rhs);
}

void save_env_changer(
        const env_changer& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    f << p;
}

env_changer load_env_changer(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    env_changer e;
    f >> e;

    return e;
}

void save_env_changer_json(
        const env_changer& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    nlohmann::json json_out;
    json_out = p;
    f << json_out;
}

env_changer load_env_changer_json(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    env_changer e;
    nlohmann::json json_in;
    f >> json_in;
    e = json_in;

    return e;
}
