#ifndef LOGIC_ONLY

#include "sim_view.h"

#else

#include "simulation.h"

#endif

#include <cassert>
#include <chrono>
#include <iostream>
#include <string>

int run_reac_norm_best(int change_freq,
                       double max_food,
                       double max_energy,
                       double max_metabolite,
                       double step,
                       int seed,
                       bool overwrite)
{
    auto name = create_reaction_norm_name(seed, change_freq);
    if(exists(name) && !overwrite)
    {
        std::cout<<"The reaction norm for the best individual"
                   "of this simulation has already been calculated";
        return 0;
    }

    auto funders_name = create_funders_success_name(seed, change_freq);
    funders_success funders_success;

    if(exists(funders_name))
        funders_success = load_funders_success(funders_name);
    else
        abort();

    auto best_ind_grn = find_last_gen_best_ind_grn(funders_success);

    auto reac_norm = calc_reaction_norm(best_ind_grn,
                                        max_food,
                                        max_energy,
                                        max_metabolite,
                                        step);

    save_reaction_norm(reac_norm, name);

    return 0;
}

int run_sim_best_rand(double amplitude,
                      int change_frequency,
                      int n_random_conditions,
                      int seed,
                      bool overwrite)
{
    auto rand_s = load_best_ind_for_rand_cond(seed,change_frequency);
    auto name = create_best_random_condition_name(rand_s,amplitude);

    if(exists(name) && !overwrite)
    {
        std::cout << "The random conditions against the best for this simulation"
                     " have already been tested!" << std::endl;
        return 0;
    }

    auto rand_start = std::chrono::high_resolution_clock::now();
    place_start_cells(rand_s.get_pop());
    run_random_conditions(rand_s,
                          n_random_conditions,
                          amplitude,
                          name);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - rand_start);
    std::cout<< "random condition best test :" << duration.count() << "s" << std::endl;
    return 0;
}

int run_sim_rand(double amplitude,
                 int change_frequency,
                 int n_random_conditions,
                 int seed,
                 bool overwrite)
{
    auto rand_s = load_sim_from_last_pop(seed,change_frequency);
    auto name = create_random_condition_name(rand_s,amplitude);

    if(exists(name) && !overwrite)
    {
        std::cout << "The random conditions for this simulation"
                     " have already been tested!" << std::endl;
        return 0;
    }

    auto rand_start = std::chrono::high_resolution_clock::now();
    place_start_cells(rand_s.get_pop());
    run_random_conditions(rand_s,
                          n_random_conditions,
                          amplitude,
                          name);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - rand_start);
    std::cout<< "random condition test :" << duration.count() << "s" << std::endl;
    return 0;
}

int run_sim_evo(const env_param& e,
                const meta_param& m,
                const pop_param& p,
                int change_frequency,
                int seed,
                bool overwrite)
{
    if(exists(create_sim_par_name(seed, change_frequency))
            && !overwrite)
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
    return 0;
}

int run_standard(const env_param& e,
                 const meta_param& m,
                 const pop_param& p,
                 double amplitude,
                 int change_frequency,
                 int n_random_conditions,
                 int seed)
{
    simulation s{sim_param{e, m, p}};
    auto start = std::chrono::high_resolution_clock::now();

    exec(s);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    std::cout << "simualtion :"<< duration.count() << "s" << std::endl;

    auto rand_start = std::chrono::high_resolution_clock::now();

    auto rand_s = load_sim_from_last_pop(seed,change_frequency);
    run_random_conditions(rand_s, n_random_conditions, amplitude, "standard_rand_run.csv");

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<float>(stop - rand_start);
    std::cout<< "random condition test :" << duration.count() << "s" << std::endl;

    duration = std::chrono::duration<float>(stop - start);
    std::cout<< "overall time :" << duration.count() << "s" << std::endl;
    return 0;
}

#ifndef LOGIC_ONLY
int run_visual_evo (const env_param& e,
                    const meta_param& m,
                    const pop_param& p,
                    double amplitude,
                    int change_frequency,
                    int n_random_conditions,
                    int seed)
{

    sim_view v(sim_param{e, m, p});
    v.exec();
    std::cout << "view: ";
    auto start = std::chrono::high_resolution_clock::now();
    auto rand_s = no_demographic_copy(load_sim_from_last_pop(seed,change_frequency));
    place_start_cells(rand_s.get_pop());
    v.run_random_conditions(rand_s, n_random_conditions, amplitude);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    std::cout << duration.count() << "s" << std::endl;
    std::cout << "n_cycles:" << v.get_sim().get_cycle() << std::endl;
    return 0;
}

int replay_cycle_from_evo (
                    int change_frequency,
                    int seed,
        int cycle)
{

    sim_view v(load_sim_parameters(create_sim_par_name(seed, change_frequency)));
    auto funders_success = load_funders_success(create_funders_success_name(seed, change_frequency));
    v.get_sim().set_funders_success(funders_success);
    v.get_sim().set_demo_sim(load_demographic_sim(create_sim_demo_name(seed,change_frequency)));
    pop_from_funders_gen(v.get_sim(),cycle);
    v.get_grid_view().prepare_grid(v.get_sim().get_env().get_grid());
    v.exec_cycle_visual(v.get_sim());
    std::cout << "view: ";
    return 0;
}
#endif


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
    test_utilities();
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
    bool overwrite = false;
    int replay_cycle = 0;

    check_for_cmd_param(args,
                        seed,
                        change_freq,
                        n_random_conditions,
                        replay_cycle,
                        amplitude,
                        overwrite);
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
        run_visual_evo(e,m,p,
                       amplitude,
                       change_freq,
                       n_random_conditions,
                       seed);
#endif
    }
    if(args.size() > 1 && args[1] == "--replay")
    {
#ifndef LOGIC_ONLY
        replay_cycle_from_evo(change_freq,
                              seed,
                              replay_cycle);
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
                     seed,
                     overwrite);
    }
    else  if(args.size() > 1 && args[1] == "--rand_best")
    {
        run_sim_best_rand(amplitude,
                          change_freq,
                          n_random_conditions,
                          seed,
                          overwrite);
    }
    else if(args.size() > 1 && args[1] == "--reac_norm")
    {
        ///Here hard coded params!
        double max_food = 20.0;
        auto max_energy = max_food;
        auto max_metabolite = max_energy;
        auto step = max_metabolite / 100.0;
        ///
        run_reac_norm_best(change_freq,
                           max_food,
                           max_energy,
                           max_metabolite,
                           step,
                           seed,
                           overwrite);
    }
    else
    {
        run_standard(e,m,p,
                     amplitude,
                     change_freq,
                     n_random_conditions,
                     seed);
    }

    return 0;
}

