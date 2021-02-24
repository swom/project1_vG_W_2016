#ifndef ENV_CHANGER_H
#define ENV_CHANGER_H
#include "env_param.h"

class env_changer
{
public:
    env_changer(env_param e, double amplitude);

    double get_amplitude() const noexcept {return m_amplitude;}

    const env_param& get_mean_params() const noexcept {return m_mean_par;}
private:

    ///The general amplitude of changes around the mean
    /// (variance is adjusted accoridingly for each parameter)
    double m_amplitude;

    ///A parameter file containing the means of the parameters values
    env_param m_mean_par;
};

#endif // ENV_CHANGER_H
