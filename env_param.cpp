#include "env_param.h"
#include <cassert>

env_param::env_param(int grid_side,
                     double diff_coeff,
                     double init_food,
                     double metab_degrad_rate,
                     double min_change_fraction,
                     double range_env_change):
    m_diff_coeff{diff_coeff},
    m_grid_side{grid_side},
    m_init_food{init_food},
    m_metab_degradation_rate{metab_degrad_rate},
    m_range_diff_coeff_change{(diff_coeff - range_env_change),
                              (diff_coeff + range_env_change)},
    m_range_metab_degr_change{(metab_degrad_rate - range_env_change),
                              (metab_degrad_rate + range_env_change)}
{
    assert(m_diff_coeff > -0.000000000001 &&
           m_diff_coeff < 1.000000000001);
    assert(m_grid_side >= 0);
    assert(m_init_food > -0.000000000001);
    assert(m_metab_degradation_rate > -0.000000001 &&
           m_metab_degradation_rate < 1.000000001);
    assert(min_change_fraction > -0.0000000000001);
    assert(m_range_metab_degr_change.min() > -0.0000000001);
    assert(m_range_diff_coeff_change.min() > -0.0000000001);

    if(min_change_fraction == 0)
    {
        m_step_min_degr_change = 0;
        m_step_min_diff_change = 0;
    }
    else
    m_step_min_diff_change = (m_range_diff_coeff_change.max() - m_range_diff_coeff_change.min()) / min_change_fraction;
    m_step_min_degr_change = (m_range_metab_degr_change.max() - m_range_metab_degr_change.min()) / min_change_fraction;
}

std::ostream& operator<<(std::ostream& os, const env_param& p)
{
    os << p.get_grid_side() << " , "
       << p.get_init_food() << " , "
       << p.get_diff_coeff()  << " , "
       << p.get_degr_rate() << " , "
       << (p.get_range_diff_coeff_change().max() - p.get_range_diff_coeff_change().min()) / p.get_min_step_diff_change()
       << " , "
       << (p.get_range_diff_coeff_change().max() - p.get_range_diff_coeff_change().min()) / 2.0
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
            && (lhs.get_min_step_degr_change() - rhs.get_min_step_degr_change() < 0.0001
                && lhs.get_min_step_degr_change() - rhs.get_min_step_degr_change() > -0.0001)
            && (lhs.get_min_step_diff_change() - rhs.get_min_step_diff_change() < 0.0001
                && lhs.get_min_step_diff_change() - rhs.get_min_step_diff_change() > -0.0001)
            && lhs.get_range_diff_coeff_change() == rhs.get_range_diff_coeff_change()
            && lhs.get_range_metab_degr_change() == rhs.get_range_metab_degr_change()
            ;
}

bool operator!=( const env_param& lhs, const env_param& rhs) noexcept
{
    return !(lhs == rhs);
}

env_param change_env_param_incr(const env_param& e, std::minstd_rand& rng) noexcept
{
    auto e1 = e;
    auto new_diff_coeff = e1.get_diff_coeff();
    auto old_diff_coeff = e.get_diff_coeff();
    auto new_degr_rate = e1.get_degr_rate();
    auto old_degr_rate = e.get_degr_rate();

    while(new_diff_coeff < old_diff_coeff + e.get_min_step_diff_change()
          && new_diff_coeff > old_diff_coeff - e.get_min_step_diff_change())
    {
        new_diff_coeff = e1.get_range_diff_coeff_change()(rng);
    }
    e1.set_diff_coeff(new_diff_coeff);

    while(new_degr_rate < old_degr_rate + e.get_min_step_degr_change()
          && new_degr_rate > old_degr_rate - e.get_min_step_degr_change())
    {
        new_degr_rate = e1.get_range_metab_degr_change()(rng);
    }
    e1.set_metab_degr(new_degr_rate);

    return e1;
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
        double env_change_range = 0.003;

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

    //Environment parameters are intialized with:
    //Range of change parameter: dictating the range of env chenge from one
    //Generation to the other
    //Magnitude of change:  the fraction of the possible values range
    //dictating the minimum amount of change from one generaiton
    //To the other
    {
        auto range_of_env_change = 0.1;
        auto magnitude_of_env_change = 10.0;
        env_param e {
            1,
            1,
            1,
            1,
            magnitude_of_env_change,
                    range_of_env_change,
        };

        auto range_diff_change = (e.get_range_diff_coeff_change().max() - e.get_range_diff_coeff_change().min()) / 2.0;
        auto range_metab_degr_change = (e.get_range_metab_degr_change().max() - e.get_range_metab_degr_change().min()) / 2.0;
        assert(range_diff_change - range_of_env_change < 0.0001
               && range_diff_change - range_of_env_change > -0.0001
               && range_metab_degr_change - range_of_env_change < 0.0001
               && range_metab_degr_change - range_of_env_change > -0.0001
               && range_diff_change * 2 / e.get_min_step_diff_change() == magnitude_of_env_change);
    }


    //Environmental parameters can be changed based on
    //Two values: magnitude of change and absolute range of change
    {
        std::random_device r;
        std::minstd_rand rng{r()};

        double magnitude_of_change = 3;
        double range_of_change =  0.5;
        env_param e{1, 1, 1, 1,
                    magnitude_of_change,
                            range_of_change
                   };

        double prev_deg_rate = e.get_degr_rate();
        double prev_diff_coeff = e.get_diff_coeff();

        int repeats = 1000;

        for( int i = 0;  i != repeats; i++)
        {
            auto e_prev = e;
            e = change_env_param_incr(e, rng);
            auto deg_rate = e.get_degr_rate();
            auto diff_coeff = e.get_diff_coeff();

            assert(e != e_prev);

            assert (deg_rate >= prev_deg_rate + e.get_min_step_degr_change() ||
                    deg_rate <= prev_deg_rate - e.get_min_step_degr_change());

            assert (diff_coeff >= prev_diff_coeff + e.get_min_step_diff_change() ||
                    diff_coeff <= prev_diff_coeff - e.get_min_step_diff_change());

            prev_deg_rate = deg_rate;
            prev_diff_coeff = diff_coeff;
        }
    }
}
