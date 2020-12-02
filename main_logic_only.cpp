#include "simulation.h"
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
    test_GRN();
    test_ind_param();
    test_individual();
    test_phenotype();
    test_meta_param();
    test_pop_param();
    test_population();
    test_simulation();
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
        test();
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
    int collision_check_interval = 0;

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
                pop_max,
                collision_check_interval};

    ind_param i{};

    pop_param p{1,
                100,
                0.1,
                0.025,
                0.1,
                0.01,
                10,
                0.0,
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


    if(args.size() > 1 && args[1] == "--sim")
    {
        run_sim_evo(e, i, m, p,
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
    else if(args.size() > 1 && args[1] == "--rand_evo")
    {
        run_sim_evo_rand(amplitude,
                         change_freq,
                         n_random_conditions,
                         pop_max,
                         seed,
                         rand_cond_n,
                         overwrite);
    }
    else if(args.size() > 1 && args[1] == "--reac_norm")
    {
        ///!!!!!!!Here hard coded params!!!!!!
        double max_food = 20.0;
        auto max_energy = max_food;
        auto max_metabolite = max_energy;
        auto step = max_metabolite / 100.0;
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
                                        i,
                                        amplitude,
                                        n_random_conditions);

    }
    else
    {
        run_standard(e, i, m,p,
                     amplitude,
                     change_freq,
                     n_random_conditions,
                     pop_max,
                     seed);
    }

    return 0;
}

