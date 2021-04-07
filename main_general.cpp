#include "sim_view.h"
#include <cassert>
#include <iostream>
#include <string>



int main(int argc, char ** argv) //!OCLINT tests may be long
{

    const std::vector<std::string> args(argv, argv + argc);

    int n_cycles = 500;
    int cycle_duration = 125;
    int seed = 0;
    int change_freq = 0;
    int n_seq = 50;
    int pop_max = pow(10,4);
    double amplitude = 3;
    bool overwrite = false;
    int replay_cycle = 0;
    int conditions_per_seq;
    int seq_index;
    int collision_check_interval = 0;
    bool death = false;

    check_for_cmd_param(args,
                        seed,
                        change_freq,
                        n_seq,
                        replay_cycle,
                        amplitude,
                        overwrite,
                        conditions_per_seq,
                        seq_index,
                        death);

    meta_param m{n_cycles,
                cycle_duration,
                seed,
                change_freq,
                pop_max,
                collision_check_interval};

    ind_param i{0.8,
               10,
               0.1};

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

    if (args.size() > 1 && args[1] == "--profile")
    {
        meta_param meta{10,100};
        simulation s{sim_param{e, i, meta, p}};
        exec(s);
        return 0;
    }
    if(args.size() > 1 && args[1] == "--sim")
    {
        run_sim_evo(e, i, m, p,
                    overwrite);
    }
    else  if(args.size() > 1 && args[1] == "--visual")
    {
        run_visual_evo(e,i,m,p);
    }
    else  if(args.size() > 1 && args[1] == "--rand")
    {
        run_sim_rand(amplitude,
                     change_freq,
                     n_seq,
                     pop_max,
                     seed,
                     overwrite);
    }
    else  if(args.size() > 1 && args[1] == "--rand_best")
    {
        run_sim_best_rand(amplitude,
                          change_freq,
                          n_seq,
                          pop_max,
                          seed,
                          overwrite);
    }
    else if(args.size() > 1 && args[1] == "--rand_evo")
    {
        run_sim_evo_rand(amplitude,
                         change_freq,
                         n_seq,
                         conditions_per_seq,
                         pop_max,
                         seed,
                         seq_index,
                         overwrite,
                         death);
    }
    else if(args.size() > 1 && args[1] == "--rand_evo_visual")
    {
        run_sim_evo_rand_visual(amplitude,
                         change_freq,
                         n_seq,
                         conditions_per_seq,
                         pop_max,
                         seed,
                         seq_index,
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
    else if(args.size() > 1 && args[1] == "--replay_rand_cond")
    {
        replay_rand_cond_test(change_freq,
                              seed,
                              n_seq,
                              amplitude,
                              conditions_per_seq,
                              seq_index,
                              pop_max);
    }
    else if(args.size() > 1 && args[1] == "--replay_rand_cond_best")
    {
        replay_best_rand_cond(change_freq,
                              seed,
                              n_seq,
                              amplitude,
                              conditions_per_seq,
                              seq_index,
                              pop_max);
    }
    else if(args.size() > 1 && args[1] == "--replay_rand_cond_evo")
    {
        replay_rand_cond_evo(change_freq,
                             seed,
                             n_seq,
                             amplitude,
                             conditions_per_seq,
                             seq_index,
                             pop_max);
    }
    else if(args.size() > 1 && args[1] == "--create_rand_cond_vec")
    {
        create_vector_random_conditions(e,
                                        i,
                                        amplitude,
                                        n_seq);

    }
    else
    {
        run_standard(e, i, m,p,
                     amplitude,
                     change_freq,
                     n_seq,
                     pop_max,
                     seed);
    }

    return 0;
}

