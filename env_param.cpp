#include "env_param.h"
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

std::ostream& operator<<(std::ostream& os, const env_param& p)
{
    os << p.get_grid_side() << " , "
       << p.get_init_food() << " , "
       << p.get_diff_coeff()  << " , "
       << p.get_degr_coeff()
          ;
    return os;
}

std::ifstream& operator>>(std::ifstream& is, env_param& p)
{
    int grid_side;
    double diff_coeff;
    double init_food;
    double metab_degrad_rate;
    std::string dummy; // To remove the annotation in the file
    is >>
            grid_side >> dummy >>
            init_food >> dummy >>
            diff_coeff >>dummy >>
            metab_degrad_rate;

    p = env_param {grid_side,
                diff_coeff,
                init_food,
                metab_degrad_rate

    };
    return is;
}

bool operator==(const env_param& lhs, const env_param& rhs) noexcept
{
    return
            lhs.get_grid_side() == rhs.get_grid_side()
            && lhs.get_init_food() == rhs.get_init_food()
            && lhs.get_degr_coeff() == rhs.get_degr_coeff()
            && lhs.get_diff_coeff() == rhs.get_diff_coeff()
            ;
}

void    save_env_parameters(
        const env_param& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    f << p.get_grid_side() << " , "
      << p.get_init_food() << " , "
      << p.get_diff_coeff() << " , "
      << p.get_metab_degr_rate();
}

env_param load_env_parameters(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    int grid_side;
    double diff_coeff;
    double init_food;
    double metab_degrad_rate;
    std::string dummy; // To remove the annotation in the file
    f >>
            grid_side >> dummy >>
            init_food >> dummy >>
            diff_coeff >>dummy >>
            metab_degrad_rate;

    env_param p{grid_side,
                diff_coeff,
                init_food,
                metab_degrad_rate

    };
    return p;
}

void test_env_param() noexcept //!OCLINT
{
    //env_param object can be loaded and saved to a given file name
    {
        int grid_side = 666;
        double diff_coeff = 0.12;
        double init_food = 14.0;
        double metab_degrad_rate = 0.066;

        env_param p{grid_side,
                    diff_coeff,
                    init_food,
                    metab_degrad_rate};

        //Test load and save
        const std::string filename = "env_param.csv";
        save_env_parameters(p, filename);
        const env_param q = load_env_parameters(filename);
        assert(p == q);
        //Test >> operator overload
        std::ifstream f(filename);
        env_param s;
        f >> s;
        assert(s == p);
    }
    {
        const env_param p;
        std::ostringstream s;
        s << p;
    }
}
