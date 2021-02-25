#include"tests.h"
void test_sim_param() noexcept //!OCLINT
{

    ///A simulation can be initialized by sets of parameters
    /// for population, environment and meta
    /// It suffices for this test to pass to make the code compile and execute
    {
        pop_param p;
        env_param e;
        ind_param i;
        meta_param m;
        sim_param{e, i, m, p};
    }

    //A sim_parameters can be loaded from and saved to a file
    {
        {
            //Thanks to the order the tests are run all these parameters file will be alreday
            //be created by other tests
            env_param env = load_env_parameters("env_param.csv");
            meta_param meta = load_meta_parameters("meta_param.csv");
            pop_param pop = load_pop_parameters("pop_param.csv");
            ind_param ind = load_ind_parameters("ind_param.csv");
            //ind_param ind = load_ind_parameters("ind_param.csv");
            sim_param s{env, ind, meta, pop};

            const std::string filename = "sim_param.csv";
            save_sim_parameters(s, filename);
            const sim_param q = load_sim_parameters(filename);
            assert(s == q);
        }
        {
            const meta_param p;
            std::ostringstream s;
            s << p;
        }
    }
}