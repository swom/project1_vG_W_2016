#ifndef ENV_CHANGER_H
#define ENV_CHANGER_H
#include "env_param.h"

class env_changer
{
public:
    env_changer(env_param e = env_param{},
                double amplitude = 0);

    double get_amplitude() const noexcept {return m_amplitude;}

    const env_param& get_mean_params() const noexcept {return m_mean_par;}
private:

    ///The general amplitude of changes around the mean
    /// (variance is adjusted accoridingly for each parameter)
    double m_amplitude;

    ///A parameter file containing the means of the parameters values
    env_param m_mean_par;
};

std::ostream& operator<<(std::ostream& os, const env_changer& p);

std::ifstream& operator>>(std::ifstream& is, env_changer& p);

bool operator==(const env_changer& lhs, const env_changer& rhs) noexcept;

bool operator!=(const env_changer& lhs, const env_changer& rhs) noexcept;


#endif // ENV_CHANGER_H
