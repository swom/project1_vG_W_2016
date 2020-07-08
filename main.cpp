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

#else
    // In release mode, all asserts are removed from the code
    assert(1 == 2);
#endif

    env_param e{200,
                0.1,
                10,
                0.1,
                10,
                0.05
               };

    meta_param m;

    if (args.size() > 1
            && args[1][0] == 's'
            && std::isdigit(args[1][1])
            )
    {
        std::string s_seed;
        for(size_t i = 1; i != args[1].size(); i++)
        {
            s_seed += args[1][i];
        }
        auto seed = std::stoi(s_seed);

        m = meta_param {500,
                125,
                seed,
                1};
    }
    else
    {
        m = meta_param {30,
                50,
                0,
                1};
    }

    ind_param i{};

    pop_param p{1,
                100,
                0.1,
                0.025,
                0.1,
                0.01,
                10,
                0.0,
                0.000001,
                i};

    if(args.size() > 1 && args[1] == "--visual")
    {
#ifndef LOGIC_ONLY
        sim_view v(sim_param{e, m, p});
        std::cout << "view: ";
        auto start = std::chrono::high_resolution_clock::now();
        v.exec();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        std::cout << duration.count() << "s" << std::endl;
        std::cout << "n_cycles:" << v.get_sim().get_cycle() << std::endl;
#endif
    }
    else
    {
        simulation s{sim_param{e, m, p}};
        std::cout << "logic: ";
        auto start = std::chrono::high_resolution_clock::now();
        exec(s);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        std::cout << duration.count() << "s" << std::endl;
    }

    return 0;
}

