#include "tests.h"

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
    ///
    /// For the normal distribution the mean and variance  are specified by:
    /// m_mean_diff_coeff & m_var_diff_coeff, m_mean_degr_rate & m_var_degr_rate
    ///
    /// For the uniform dostribution the range is:
    /// [m_mean* - 3 * m_var, m_mean + 3 * m_var]
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

    /// It is possible to create env_ param object
    /// starting from initial ones,
    /// that have a wider range of possible values
    {
        env_param e;
        double amplitude = 1.50;

        auto e2 = change_range_env_param(e, amplitude);

        assert(e.get_var_degr_rate() - (e2.get_var_degr_rate() / amplitude) < 0.00001
               && e.get_var_degr_rate() - (e2.get_var_degr_rate() / amplitude) > -0.00001);
        assert(e.get_var_diff_coeff() - (e2.get_var_diff_coeff() / amplitude) < 0.00001
               && e.get_var_diff_coeff() - (e2.get_var_diff_coeff() / amplitude) > -0.00001);
    }

    ///Saving to json is possible
    {
       std::string name = "json_env_test.json";
       env_param e;
       save_env_parameters_json(e, name);
       auto e1 = load_env_parameters_json(name);
       assert(e == e1);
    }
}
