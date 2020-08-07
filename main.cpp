#ifndef LOGIC_ONLY

#include "sim_view.h"

#else

#include "simulation.h"

#endif

#include <cassert>
#include <chrono>
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

    int seed = 0;
    int change_freq = 0;
    int n_random_conditions = 50;
    double amplitude = 3;

    check_for_cmd_param(args,
                        seed,
                        change_freq,
                        n_random_conditions,
                        amplitude);
    meta_param m{200,
                 125,
                 seed,
                         change_freq
                };

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
        sim_view v(sim_param{e, m, p});

        std::cout << "view: ";
        auto start = std::chrono::high_resolution_clock::now();
        auto rand_s = no_demographic_copy(load_sim_for_rand_cond(seed,change_freq));
        place_start_cells(rand_s.get_pop());
        v.run_random_conditions(rand_s, n_random_conditions, amplitude);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        std::cout << duration.count() << "s" << std::endl;
        std::cout << "n_cycles:" << v.get_sim().get_cycle() << std::endl;
#endif
    }
    else if(args.size() > 1 && args[1] == "--sim")
    {
        if(exists(create_sim_par_name(seed, change_freq)))
        {
            std::cout << "this simulation has already been run" << std::endl;
            return 0;
        }

        simulation s{sim_param{e, m, p}};

        auto start = std::chrono::high_resolution_clock::now();

        exec(s);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        std::cout << "simualtion :"<< duration.count() << "s" << std::endl;
    }
    else  if(args.size() > 1 && args[1] == "--rand")
    {
        auto rand_s = load_sim_for_rand_cond(seed,change_freq);

        if(exists(create_random_condition_name(rand_s,amplitude)))
        {
            std::cout << "The random conditions for this simulation"
                         " have already been tested!" << std::endl;
            return 0;
        }

        auto rand_start = std::chrono::high_resolution_clock::now();
        place_start_cells(rand_s.get_pop());

        save_demographic_sim(run_random_conditions(rand_s,
                                                   n_random_conditions,
                                                   amplitude),
                             create_random_condition_name(rand_s,
                                                          amplitude));

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - rand_start);
        std::cout<< "random condition test :" << duration.count() << "s" << std::endl;
    }
    else  if(args.size() > 1 && args[1] == "--rand_best")
    {
        auto rand_s = load_best_ind_for_rand_cond(seed,change_freq);

        if(exists(create_best_random_condition_name(rand_s,amplitude)))
        {
            std::cout << "The random conditions against the best for this simulation"
                         " have already been tested!" << std::endl;
            return 0;
        }

        auto rand_start = std::chrono::high_resolution_clock::now();
        place_start_cells(rand_s.get_pop());
        save_demographic_sim(run_random_conditions(rand_s,
                                                   n_random_conditions,
                                                   amplitude),
                             create_best_random_condition_name(rand_s,
                                                          amplitude));

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - rand_start);
        std::cout<< "random condition best test :" << duration.count() << "s" << std::endl;
    }
    else
    {
        simulation s{sim_param{e, m, p}};
        auto start = std::chrono::high_resolution_clock::now();

        exec(s);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        std::cout << "simualtion :"<< duration.count() << "s" << std::endl;

        auto rand_start = std::chrono::high_resolution_clock::now();

        auto rand_s = load_sim_for_rand_cond(seed,change_freq);
        run_random_conditions(rand_s, n_random_conditions, amplitude);

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration<float>(stop - rand_start);
        std::cout<< "random condition test :" << duration.count() << "s" << std::endl;

        duration = std::chrono::duration<float>(stop - start);
        std::cout<< "overall time :" << duration.count() << "s" << std::endl;
    }

    return 0;
}

