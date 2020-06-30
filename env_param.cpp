#include "env_param.h"
#include <cassert>

env_param::env_param(int grid_side,
                     double diff_coeff,
                     double init_food,
                     double metab_degrad_rate,
                     double min_step_env_change,
                     double range_env_change):
    m_diff_coeff{diff_coeff},
    m_grid_side{grid_side},
    m_init_food{init_food},
    m_metab_degradation_rate{metab_degrad_rate},
    m_min_step_env_change{min_step_env_change},
    m_range_env_change{range_env_change}
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
       << p.get_degr_rate() << " , "
       << p.get_min_step_env_change() << " , "
       << p.get_range_env_change()
          ;
    return os;
}

std::ifstream& operator>>(std::ifstream& is, env_param& p)
{
    int grid_side;
    double diff_coeff;
    double init_food;
    double metab_degrad_rate;
    double min_step_env_change;
    double range_env_change;

    std::string dummy; // To remove the annotation in the file
    is >>
            grid_side >> dummy >>
            init_food >> dummy >>
            diff_coeff >>dummy >>
            metab_degrad_rate >> dummy >>
            min_step_env_change >> dummy >>
            range_env_change;

    p = env_param {grid_side,
            diff_coeff,
            init_food,
            metab_degrad_rate,
            min_step_env_change,
            range_env_change
};

    return is;
}

bool operator==(const env_param& lhs, const env_param& rhs) noexcept
{
    return
            (lhs.get_grid_side() - rhs.get_grid_side() < 0.0001
             && lhs.get_grid_side() - rhs.get_grid_side() > -0.0001)
            && (lhs.get_init_food() - rhs.get_init_food() < 0.0001
                && lhs.get_init_food() - rhs.get_init_food() > -0.0001)
            && (lhs.get_degr_rate() - rhs.get_degr_rate() < 0.0001
                && lhs.get_degr_rate() - rhs.get_degr_rate() > -0.0001)
            && (lhs.get_diff_coeff() - rhs.get_diff_coeff() < 0.0001
                && lhs.get_diff_coeff() - rhs.get_diff_coeff() > -0.0001)
            && (lhs.get_min_step_env_change() - rhs.get_min_step_env_change() < 0.0001
                && lhs.get_min_step_env_change() - rhs.get_min_step_env_change() > -0.0001)
            && (lhs.get_range_env_change() - rhs.get_range_env_change() < 0.0001
                && lhs.get_range_env_change() - rhs.get_range_env_change() > -0.0001)
            ;
}

bool operator!=( const env_param& lhs, const env_param& rhs) noexcept
{
    return !(lhs == rhs);
}

env_param change_env_param_incr(const env_param& e) noexcept
{

}

void save_env_parameters(
        const env_param& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    f << p;
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

void test_env_param() noexcept //!OCLINT
{
    //env_param object can be loaded and saved to a given file name
    {
        int grid_side = 666;
        double diff_coeff = 0.12;
        double init_food = 14.0;
        double metab_degrad_rate = 0.066;
        double min_step_env_change = 0.36;
        double env_change_range = 0.3;

        env_param p{grid_side,
                    diff_coeff,
                    init_food,
                    metab_degrad_rate,
                    min_step_env_change,
                    env_change_range};

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

    //Meta parameters are intialized with:
    //Range of change parameter: dictating the range of env chenge from one
    //Generation to the other
    //Magnitude of change: dictating the minimum amount of change from one generaiton
    //To the other
    {
        auto range_of_env_change = 0.1;
        auto magnitude_of_env_change = range_of_env_change / 10;
        env_param e {
            1,
            1,
            1,
            1,
            magnitude_of_env_change,
                    range_of_env_change,
        };

        assert(e.get_range_env_change() == range_of_env_change &&
               e.get_min_step_env_change() == magnitude_of_env_change);
    }


    //Environmental parameters can be changed based on
    //Two values: magnitude of change and absolute range of change
    {
        env_param e;
        env_param e1;

        double magnitude_of_change = 20;
        double range_of_change = magnitude_of_change * 3;

        double prev_deg_rate = e.get_degr_rate();
        double prev_diff_coeff = e.get_diff_coeff();

        int repeats = 1000;

        for( int i = 0;  i != repeats; i++)
        {
            e1 = change_env_param_incr(e);
            auto deg_rate = e1.get_degr_rate();
            auto diff_coeff = e1.get_diff_coeff();

            assert(e != e1);
            assert(e1.get_degr_rate() < e.get_degr_rate() + range_of_change &&
                   e1.get_degr_rate() > e.get_degr_rate() - range_of_change &&
                   e1.get_diff_coeff() < e.get_diff_coeff() + range_of_change &&
                   e1.get_diff_coeff() > e.get_diff_coeff() - range_of_change);

            assert( (deg_rate >= prev_deg_rate + magnitude_of_change ||
                     deg_rate <= prev_deg_rate - magnitude_of_change)
                    && (diff_coeff >= prev_diff_coeff + magnitude_of_change ||
                        diff_coeff <= prev_diff_coeff - magnitude_of_change));

            prev_deg_rate = deg_rate;
            prev_diff_coeff = diff_coeff;
        }
    }
}
