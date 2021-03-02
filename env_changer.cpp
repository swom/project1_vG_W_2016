#include "env_changer.h"

env_changer::env_changer(env_param e, double amplitude):
    m_amplitude{amplitude},
    m_mean_par{e}
{
}

std::ostream& operator<<(std::ostream& os, const env_changer& ec)
{
    os << ec.get_mean_params() << " , " << ec.get_amplitude();
    return os;
}

std::ifstream& operator>>(std::ifstream& is, env_changer& ec)
{
 std::string dummy;
 env_param e;
 is >> e;

 is >> dummy;

 double amplitude;
 is >> amplitude;

 ec = env_changer{e, amplitude};

 return is;
}

bool operator==(const env_changer& lhs, const env_changer& rhs) noexcept
{
 return lhs.get_mean_params() == rhs.get_mean_params() &&
         lhs.get_amplitude() == rhs.get_amplitude();
}

bool operator!=(const env_changer& lhs, const env_changer& rhs) noexcept
{
return !(lhs == rhs);
}
