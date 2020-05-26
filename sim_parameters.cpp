#include "sim_parameters.h"
#include <cassert>

sim_param::sim_param(unsigned int start_pop_size, //!OCLINT
                     unsigned int exp_new_pop_size,
                     double min_dist,
                     int grid_side,
                     double diff_coeff,
                     double init_food,
                     double mutation_prob,
                     double mutation_step,
                     double base_disp_prob,
                     double spore_advantage,
                     double reproduction_prob,
                     double metab_degrad_rate,
                     int n_cycles,
                     int cycle_duration):
    m_env_param{grid_side,
                diff_coeff,
                init_food,
                metab_degrad_rate},
    m_meta_param{n_cycles,
                 cycle_duration},
    m_pop_param{start_pop_size,
                exp_new_pop_size,
                min_dist,
                mutation_prob,
                mutation_step,
                base_disp_prob,
                spore_advantage,
                reproduction_prob}
{

}

sim_param::sim_param(env_param e, meta_param m, pop_param p):
    m_env_param{e},
    m_meta_param{m},
    m_pop_param{p}
{

}

std::ostream& operator<<(std::ostream& os, const sim_param& p)
{
    os << p.get_env_param()
       << p.get_meta_param()
       << p.get_pop_param();
    return os;
}

bool operator==(const sim_param& lhs, const sim_param& rhs) noexcept
{
    return
            lhs.get_env_param() == rhs.get_env_param()
            && lhs.get_meta_param() == rhs.get_meta_param()
            &&lhs.get_pop_param() == lhs.get_pop_param()
            ;

}

void save_sim_parameters(
        const sim_param& p,
        const std::string& filename
        )
{
    std::ofstream f(filename);
    f << p.get_env_param()  << " , " <<
         p.get_meta_param() << " , " <<
         p.get_pop_param();
}

sim_param load_sim_parameters(
        const std::string& filename
        )
{
    std::ifstream f(filename);
    pop_param p;
    meta_param m;
    env_param e;

    std::string dummy; // To remove the annotation (" , ") in the file

    f >> e >> dummy;
    f >> m >> dummy;
    f >> p;

    sim_param s{e,m,p};
    return  s;

}

void test_sim_param() noexcept //!OCLINT
{

    ///A simulation can be initialized by sets of parameters
    /// for population, environment and meta
    /// It suffices for this test to pass to make the code compile and execute
    {
        pop_param p;
        env_param e;
        meta_param m;
        sim_param{e, m, p};
    }

    //A sim_parameters can be loaded from and saved to a file
    {
        {
            //Thanks to the order the tests are run all these parameters file will be alreday
            //be created by other tests
            env_param env = load_env_parameters("env_param.csv");
            meta_param meta = load_meta_parameters("meta_param.csv");
            pop_param pop = load_pop_parameters("pop_param.csv");
            //ind_param ind = load_ind_parameters("ind_param.csv");
            sim_param s{env, meta, pop};

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
