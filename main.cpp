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

    if (args.size() > 3
            && (args[1] == "--sim" || args[1] == "--rand")
            )
    {
        if(args[2][0] == 's' && std::isdigit(args[2][1]))
        {
            std::string s_seed;
            for(size_t i = 1; i != args[2].size(); i++)
            {
                s_seed += args[2][i];
            }
            seed = std::stoi(s_seed);
        }
        else
        {
            abort();
        }

        if(args[3][0] == 'f' && std::isdigit(args[3][1]))
        {
            std::string s_change_freq;
            for(size_t i = 1; i != args[3].size(); i++)
            {
                s_change_freq += args[3][i];
            }
            change_freq = std::stoi(s_change_freq);
        }
        else
        {
            abort();
        }
    }
    else if (args.size() > 2)
    {
        abort();
    }

    meta_param m{10,
                 50,
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

    int n_random_conditions = 5;
    double amplitude = 3;

    if (args.size() > 4
            &&  args[1] == "--rand"
            &&  args[4][0] == 'a'
            )
    {
        std::string s_amplitude;
        for(size_t i = 1; i != args[4].size(); i++)
        {
            s_amplitude += args[4][i];
        }
        amplitude = std::stod(s_amplitude);
    }

    if (args.size() > 5
            &&  args[1] == "--rand"
            &&  args[5][0] == 'n'
            )
    {
        std::string s_n_random_conditions;
        for(size_t i = 1; i != args[5].size(); i++)
        {
            s_n_random_conditions += args[5][i];
        }
        n_random_conditions = std::stoi(s_n_random_conditions);
    }

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
        simulation s{sim_param{e, m, p}};
        auto start = std::chrono::high_resolution_clock::now();

        exec(s);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        std::cout << "simualtion :"<< duration.count() << "s" << std::endl;
    }
    else  if(args.size() > 1 && args[1] == "--rand")
    {

        auto rand_start = std::chrono::high_resolution_clock::now();

        auto rand_s = load_sim_for_rand_cond(seed,change_freq);
        place_start_cells(rand_s.get_pop());
        run_random_conditions(rand_s,
                              n_random_conditions,
                              amplitude);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - rand_start);
        std::cout<< "random condition test :" << duration.count() << "s" << std::endl;
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

