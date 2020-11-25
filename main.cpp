#ifndef LOGIC_ONLY

#include "sim_view.h"

#else

#include "simulation.h"

#endif

#include <cassert>
#include <iostream>
#include <string>

void test() {
    test_demographic_cycle();
    test_demographic_sim();
    test_env_grid_cell();
    test_environment();
    test_env_param();
    test_funder_data();
    test_funders();
    test_funders_success();
#ifndef LOGIC_ONLY
    test_grid_view();
#endif
    test_GRN();
    test_ind_param();
    test_individual();
    test_phenotype();
    test_meta_param();
    test_pop_param();
    test_population();
    test_simulation();
#ifndef LOGIC_ONLY
    test_sim_view();
#endif
    test_sim_param();
    test_utilities();
}


int main(int argc, char ** argv) //!OCLINT tests may be long
{

    const std::vector<std::string> args(argv, argv + argc);

#ifndef NDEBUG
    if (args.size() > 1 && args[1] == "--test")
    {
        test();
        // We've already tested, so the program is done
        return 0;
    }
#endif
    if (args.size() > 1 && args[1] == "--profile")
    {
        test_demographic_cycle();
        // We've already tested, so the program is done
        return 0;
    }

    int n_cycles = 500;
    int cycle_duration = 125;
    int seed = 0;
    int change_freq = 0;
    int n_random_conditions = 50;
    int pop_max = pow(10,4);
    double amplitude = 3;
    bool overwrite = false;
    int replay_cycle = 0;
    int seed_rand_cond;
    int rand_cond_n;

    check_for_cmd_param(args,
                        seed,
                        change_freq,
                        n_random_conditions,
                        replay_cycle,
                        amplitude,
                        overwrite,
                        seed_rand_cond,
                        rand_cond_n);

    meta_param m{n_cycles,
                 cycle_duration,
                 seed,
                 change_freq,
                 pop_max};

    ind_param ind{};

    pop_param p{1,
                100,
                0.1,
                0.025,
                0.1,
                0.01,
                10,
                0.0,
                ind
               };

    env_param e{200,
                0.1,
                10,
                0.1,
                0.1,
                0.1,
                0.01,
                0.01
               };


    if(args.size() > 1 && args[1] == "--visual")
    {
#ifndef LOGIC_ONLY
        run_visual_evo(e,m,p);
#endif
    }
    else if(args.size() > 1 && args[1] == "--replay")
    {
#ifndef LOGIC_ONLY
        replay_cycle_from_evo(change_freq,
                              seed,
                              replay_cycle);
#endif
    }
    else if(args.size() > 1 && args[1] == "--replay_rand_cond")
    {
#ifndef LOGIC_ONLY
        replay_rand_cond(change_freq,
                         seed,
                         n_random_conditions,
                         amplitude,
                         seed_rand_cond,
                         rand_cond_n);
#endif
    }
    else if(args.size() > 1 && args[1] == "--sim")
    {
        run_sim_evo(e,m,p,
                    change_freq,
                    seed,
                    overwrite);
    }
    else  if(args.size() > 1 && args[1] == "--rand")
    {
        run_sim_rand(amplitude,
                     change_freq,
                     n_random_conditions,
                     pop_max,
                     seed,
                     overwrite);
    }
    else  if(args.size() > 1 && args[1] == "--rand_best")
    {
        run_sim_best_rand(amplitude,
                          change_freq,
                          n_random_conditions,
                          pop_max,
                          seed,
                          overwrite);
    }
    else if(args.size() > 1 && args[1] == "--reac_norm")
    {
        ///!!!!!!!Here hard coded params!!!!!!
        double max_food = 20.0;
        auto max_energy = max_food;
        auto max_metabolite = max_energy;
        auto step = max_metabolite / 200.0;
        ///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        run_reac_norm_best(change_freq,
                           max_food,
                           max_energy,
                           max_metabolite,
                           step,
                           seed,
                           overwrite);
    }
    else if(args.size() > 1 && args[1] == "--create_rand_cond_vec")
    {
        create_vector_random_conditions(e,
                                        p.get_ind_param(),
                                        amplitude,
                                        n_random_conditions);

    }
    else
    {
        run_standard(e,m,p,
                     amplitude,
                     change_freq,
                     n_random_conditions,
                     pop_max,
                     seed);
    }

    return 0;
}

