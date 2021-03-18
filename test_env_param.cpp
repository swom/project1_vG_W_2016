#include "tests.h"

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

    ///Saving to json is possible
    {
       std::string name = "env_param.json";
       env_param e;
       save_env_parameters_json(e, name);
       auto e1 = load_env_parameters_json(name);
       assert(e == e1);
    }
}
