#include "env_param.h"
#include <cassert>

env_param::env_param(int grid_side,
                     double diff_coeff,
                     double init_food,
                     double metab_degrad_rate,
                     double mean_diff_coeff,
                     double mean_degr_rate,
                     double var_diff_coeff,
                     double var_degr_coeff):
    m_diff_coeff{diff_coeff},
    m_grid_side{grid_side},
    m_init_food{init_food},
    m_metab_degradation_rate{metab_degrad_rate},
    m_mean_diff_coeff{mean_diff_coeff},
    m_mean_degr_rate{mean_degr_rate},
    m_var_diff_coeff{var_diff_coeff},
    m_var_degr_rate{var_degr_coeff}
{
    assert(m_diff_coeff > -0.000000000001 &&
           m_diff_coeff < 1.000000000001);
    assert(m_grid_side >= 0);
    assert(m_init_food > -0.000000000001);
    assert(m_metab_degradation_rate > -0.000000001 &&
           m_metab_degradation_rate < 1.000000001);
    assert(m_var_degr_rate > -0.00000001);
    assert(m_var_diff_coeff > -0.00000001);
    assert(m_mean_degr_rate > -0.00000001);
    assert(m_mean_diff_coeff > -0.00000001);
    assert(m_mean_diff_coeff - m_var_diff_coeff * 3 > -0.0000000000001);
    assert(m_mean_degr_rate - m_var_degr_rate * 3 > -0.0000000001);
}

std::ostream& operator<<(std::ostream& os, const env_param& e)
{
    os << e.get_grid_side() << " , "
       << e.get_init_food() << " , "
       << e.get_diff_coeff()  << " , "
       << e.get_degr_rate() << " , "
       << e.get_mean_diff_coeff() << " , "
       << e.get_mean_degr_rate() << " , "
       << e.get_var_diff_coeff() << " , "
       << e.get_var_degr_rate();

    return os;
}

std::ifstream& operator>>(std::ifstream& is, env_param& e)
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
            metab_degrad_rate,
            mean_diff_coeff,
            mean_degr_rate,
            var_diff_coeff,
            var_degr_coeff}
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

env_param change_env_param_norm(const env_param& e, std::minstd_rand& rng) noexcept
{
    auto e1 = e;

    auto new_diff_coeff = std::normal_distribution<double>{e.get_mean_diff_coeff(),
            e.get_var_diff_coeff()}(rng);

    e1.set_diff_coeff(new_diff_coeff);

    auto new_degr_rate = std::normal_distribution<double>{e.get_mean_degr_rate(),
            e.get_var_degr_rate()}(rng);

    e1.set_metab_degr(new_degr_rate);

    return e1;
}

env_param change_env_param_unif(const env_param& e, std::minstd_rand& rng) noexcept
{
    auto e1 = e;

    auto new_diff_coeff = std::uniform_real_distribution<double>{e.get_mean_diff_coeff() - 3 * e.get_var_diff_coeff(),
            e.get_mean_diff_coeff() + 3 * e.get_var_diff_coeff()}(rng);

    e1.set_diff_coeff(new_diff_coeff);

    auto new_degr_rate = std::normal_distribution<double>{e.get_mean_degr_rate() - 3 * e.get_var_degr_rate(),
            e.get_mean_degr_rate() + 3 * e.get_var_degr_rate()}(rng);

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
        double mean_diff_coeff = 0.223;
        double mean_degr_rate = 0.456;
        double var_diff_coeff = 0.05;
        double var_degr_rate = 0.012;

        env_param p{grid_side,
                    diff_coeff,
                    init_food,
                    metab_degrad_rate,
                    mean_diff_coeff,
                    mean_degr_rate,
                    var_diff_coeff,
                    var_degr_rate};

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

    ///Environmental parameters can be changed
    /// according to a normal distribution
    /// or an uniform distribution
    /// by modifying the m_metab_degradation_rate and m_diffusion_coefficient
    /// For the normal distribution
    /// the mean and variance  are specified by:
    /// m_mean_diff_coeff & m_var_diff_coeff
    /// m_mean_degr_rate & m_var_degr_rate
    /// For the uniform dostribution
    /// the range is [m_mean* - 3 * m_var, m_mean* + 3 * m_var]
    {
        std::random_device r;
        std::minstd_rand rng{r()};

        auto mean_diff_coeff = 0.1;
        auto mean_degr_coeff = 0.1;
        auto var_degr_rate = mean_degr_coeff / 3 - 0.01;
        auto var_diff_coeff = mean_diff_coeff / 3 - 0.01;
        env_param e{1, 0.1, 0.1, 0.1,
                    mean_diff_coeff,
                            mean_degr_coeff,
                            var_diff_coeff,
                            var_degr_rate};

        int repeats = 1000;

        double norm_mean_of_diff_coefficients = e.get_diff_coeff();
        double norm_mean_of_degr_rates = e.get_degr_rate();
        std::vector<double> variance_vector_degr_rate{norm_mean_of_degr_rates};
        std::vector<double> variance_vector_diff_coeff{norm_mean_of_diff_coefficients};

        double unif_mean_of_diff_coefficients = e.get_diff_coeff();
        double unif_mean_of_degr_rates = e.get_degr_rate();

        for( int i = 0;  i != repeats; i++)
        {
            auto e_prev = e;
            e = change_env_param_norm(e, rng);
            norm_mean_of_diff_coefficients += e.get_diff_coeff();
            norm_mean_of_degr_rates += e.get_degr_rate();
            variance_vector_diff_coeff.push_back( e.get_diff_coeff());
            variance_vector_degr_rate.push_back(e.get_degr_rate());
            assert(e != e_prev);
            e = e_prev;
            change_env_param_unif(e, rng);
            unif_mean_of_degr_rates += e.get_degr_rate();
            unif_mean_of_diff_coefficients += e.get_diff_coeff();
        }

        norm_mean_of_diff_coefficients /= repeats;
        norm_mean_of_degr_rates /= repeats;
        unif_mean_of_degr_rates /= repeats;
        unif_mean_of_diff_coefficients /= repeats;

        assert(mean_diff_coeff - norm_mean_of_diff_coefficients < 0.01
               && mean_diff_coeff - norm_mean_of_diff_coefficients > -0.01);
        assert(mean_degr_coeff - norm_mean_of_degr_rates < 0.01
               && mean_degr_coeff - norm_mean_of_degr_rates > -0.01);
        assert(mean_diff_coeff - unif_mean_of_diff_coefficients < 0.01
               && mean_diff_coeff - unif_mean_of_diff_coefficients > -0.01);
        assert(mean_degr_coeff - unif_mean_of_degr_rates < 0.01
               && mean_degr_coeff - unif_mean_of_degr_rates > -0.01);
    }
}
