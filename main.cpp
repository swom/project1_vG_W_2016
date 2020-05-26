#include "sim_view.h"

#include <cassert>
#include <chrono>
#include <iostream>
#include <string>


void test() {
    test_env_grid_cell();
    test_environment();
    test_env_param();
    test_grid_view();
    test_GRN();
    test_ind_param();
    test_individual();
    test_phenotype();
    test_meta_param();
    test_pop_param();
    test_population();
    test_simulation();
    test_sim_view();
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
    else
#else
    // In release mode, all asserts are removed from the code
    assert(1 == 2);
#endif
#ifndef LOGIC_ONLY
    if(args.size() > 1 && args[1] == "--visual")
    {
        env_param e{200, 0.01,10,0.1};

        meta_param m{50,
                     125};

        ind_param i{0.8,
                    10,
                    0.1,
                    0};

        pop_param p{1,
                    100,
                    0.1,
                    0.025,
                    0.1,
                    0.01,
                    10,
                    0.5,
                    0,
                    i};

        sim_view v(sim_param{e, m, p});
        std::cout << "view: ";
        auto start = std::chrono::high_resolution_clock::now();
        v.exec();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        std::cout << duration.count() << "s" << std::endl;

        simulation s{sim_param{e, m, p}};
        std::cout << "logic: ";
        start = std::chrono::high_resolution_clock::now();
        exec(s);
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration<float> (stop - start);
        std::cout << duration.count() << "s" << std::endl;

        //        int sim_time = 0;
        //        while (v.get_window().isOpen())
        //          {
        //            bool must_quit{v.process_events()};
        //            if (must_quit)
        //              return 0;
        //            if(sf::Keyboard::isKeyPressed(sf::Keyboard::A))
        //              {
        //                sim_time++;
        //                tick(v.get_sim());
        //              }
        //            if(sf::Keyboard::isKeyPressed(sf::Keyboard::R))
        //              {
        //               reset_sim(v.get_sim());
        //               v.prepare_pop();
        //              }
        //            v.show();
        //          }
    }
#endif
    {
        simulation s(sim_param{1,1,0.1,20,0.1,1});
        int time_now = 0;
        int time_limit = 200;
        while(time_now != time_limit)
        {
            tick(s);
            time_now++;
        }
    }

    return 0;

}

