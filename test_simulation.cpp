#include"test_simulation.h"

void test_activate_death()
{
    simulation s;

    //let's give a cycle some time
    s.get_meta_param().set_n_timesteps(50);
    //let's give individuals no energy, so they should die
    set_all_inds_en(s.get_pop(), 0);
    //let's give env no food
    set_all_cell_food(s.get_env(),0);
    assert(all_en_pop_equals(s.get_pop(),0));

    //With death not active no ind should starve
    assert(!death_is_active(s));
    auto n_ind_start = s.get_pop().get_pop_size();

    exec_cycle(s);
    s.reset_timesteps();
    s.tick_cycles();

    auto n_ind_end = s.get_pop().get_pop_size();

    assert(n_ind_end != 0);
    assert(n_ind_end == n_ind_start);
    assert(all_en_pop_NOT_equals(s.get_pop(),0));

    //With death active all inds should starve

    s.activate_death();

    //let's give individuals no energy, so they should die
    set_all_inds_en(s.get_pop(), 0);
    //let's give env no food
    set_all_cell_food(s.get_env(),0);

    exec_cycle(s);
    s.reset_timesteps();
    s.tick_cycles();

    n_ind_end = s.get_pop().get_pop_size();

    assert(n_ind_end < n_ind_start);
    assert(n_ind_end == 0);
}

void test_exec_change()
{
    simulation s;
    ///let' put more than one cycle in meta_param
    int n_cycles = 50;

    s.get_meta_param() =  meta_param{n_cycles};

    int n_rand_cond = 50;
    double amplitude = 3.0;
    int seed = 0;

    auto rand_cond = create_rand_conditions_unif(s.get_env().get_param(),
                                                 s.get_pop().get_ind(0).get_param(),
                                                 n_rand_cond,
                                                 amplitude,
                                                 seed);

    exec_change(s, rand_cond);

    auto demo_cycles = s.get_demo_sim().get_demo_cycles();

    for(size_t i = 0; i != demo_cycles.size() - 1 ; i++)
    {
        auto condition_index = i / (n_cycles / n_rand_cond);
        assert(demo_cycles[i].get_env_param() == rand_cond[condition_index].first);
        assert(demo_cycles[i].get_ind_param() == rand_cond[condition_index].second);
    }
}

void test_rand_cond_matrix()
{

    double amplitude = 3.0;
    int conditions_per_sequence = 3;
    int number_of_sequences = 2;


    auto rand_cond_matrix = create_rand_conditions_matrix(env_param{}, ind_param{},
                                                          number_of_sequences,
                                                          conditions_per_sequence,
                                                          amplitude);


    for(int i = 0; i != number_of_sequences; i++)
    {
        assert(rand_cond_matrix[i] == create_rand_conditions_unif(env_param{},ind_param{},conditions_per_sequence,amplitude, i));
    }

}

void test_rand_cond_matrix_extreme()
{

    double amplitude = 3.0;
    int conditions_per_sequence = 50;
    int number_of_sequences = 50;
    env_param e{};
    ind_param ip{};

    auto expected_max_diff_coeff = e.get_mean_diff_coeff() + (e.get_var_diff_coeff() * amplitude * 3);
    auto expected_min_diff_coeff = e.get_mean_diff_coeff() - (e.get_var_diff_coeff() * amplitude * 3);


    auto rand_cond_matrix = create_rand_conditions_matrix_extreme(e, ip,
                                                                  number_of_sequences,
                                                                  conditions_per_sequence,
                                                                  amplitude);


    for(int i = 0; i != number_of_sequences; i++)
    {
        assert(rand_cond_matrix[i] == create_rand_conditions_unif_extreme(env_param{},ind_param{},conditions_per_sequence,amplitude, i));
    }

}


///It is possible to create random conditions drawing from a uniform distribution
void test_create_rand_cond_unif()
{
    env_param e{};
    ind_param i{};
    int n_cond = 100;
    double amplitude = 3;
    int seed = 0;

    auto expected_mean_diff_coeff = e.get_mean_diff_coeff();
    auto expected_max_diff_coeff = e.get_mean_diff_coeff() + e.get_var_diff_coeff() * amplitude * 3;
    auto expected_min_diff_coeff = e.get_mean_diff_coeff() - (e.get_var_diff_coeff() * amplitude * 3);

    auto rand_cond = create_rand_conditions_unif(e, i, n_cond, amplitude, seed);
    auto max_diff_coeff = find_max_diff_coeff_rand_cond(rand_cond);
    auto min_diff_coeff = find_min_diff_coeff_rand_cond(rand_cond);
    auto mean_diff_coeff = mean_diff_coeff_rand_cond(rand_cond);

    assert(max_diff_coeff - expected_max_diff_coeff > -0.01 &&
           max_diff_coeff - expected_max_diff_coeff < 0.01);

    assert(min_diff_coeff - expected_min_diff_coeff > -0.01 &&
           min_diff_coeff - expected_min_diff_coeff < 0.01);

    assert(mean_diff_coeff - expected_mean_diff_coeff > -0.01 &&
           mean_diff_coeff - expected_mean_diff_coeff < 0.01);
}

void test_create_rand_cond_unif_extreme()
{
    env_param e{};
    ind_param i{};
    int n_cond = 200;
    double amplitude = 3;
    int seed = 0;

    auto expected_mean_diff_coeff = e.get_mean_diff_coeff();
    auto expected_max_diff_coeff = e.get_mean_diff_coeff() + e.get_var_diff_coeff() * amplitude * 3;
    auto expected_min_diff_coeff = e.get_mean_diff_coeff() - (e.get_var_diff_coeff() * amplitude * 3);

    auto rand_cond = create_rand_conditions_unif_extreme(e, i, n_cond, amplitude, seed);

    auto max_diff_coeff = find_max_diff_coeff_rand_cond(rand_cond);
    auto min_diff_coeff = find_min_diff_coeff_rand_cond(rand_cond);
    auto mean_diff_coeff = mean_diff_coeff_rand_cond(rand_cond);

    //Test that min max and mean are same as in create_rand_cond_unif()
    assert(max_diff_coeff - expected_max_diff_coeff > -0.01 && max_diff_coeff - expected_max_diff_coeff < 0.01);
    assert(min_diff_coeff - expected_min_diff_coeff > -0.01 && min_diff_coeff - expected_min_diff_coeff < 0.01);
    assert(mean_diff_coeff - expected_mean_diff_coeff > -0.01 && mean_diff_coeff - expected_mean_diff_coeff < 0.01);

    //Assert that all values are not in the inner range of the uniform distribution
    assert(std::all_of(rand_cond.begin(), rand_cond.end(),
                       [&expected_max_diff_coeff, &expected_min_diff_coeff, &amplitude]
                       (const std::pair<env_param, ind_param>& cond)
    {return !(cond.first.get_diff_coeff() > expected_min_diff_coeff + cond.first.get_var_diff_coeff() * amplitude &&
              cond.first.get_diff_coeff() < expected_max_diff_coeff - cond.first.get_var_diff_coeff() * amplitude);}));

}

void test_rand_cond_test_after_rand_evo()
{
    int grid_side = 10;
    int n_cycles = 10;
    int cycle_duration = 30;
    int seed = 1478;
    //Run a rand_evo sim
    sim_param s_p{env_param{grid_side},
              ind_param{},
                    meta_param{n_cycles,
                                cycle_duration,
                                seed},
                          pop_param{1,10}};
    simulation s{s_p};

    int n_seq = 1;
    int seq_index = n_seq - 1;
    int cond_per_seq = 2;
    double amplitude = 3;
    int pop_max = 100;
    auto prefix = "../"+create_rand_extreme_prefix(amplitude, cond_per_seq, seq_index);

    auto s_f = run_evo_random_conditions(s, n_seq, cond_per_seq, seq_index,
                                         pop_max, amplitude, prefix);

    auto original_seed = s.get_meta_param().get_seed();
    auto original_change = s.get_meta_param().get_change_freq();


    auto new_s = load_sim_no_pop(original_seed, original_change, prefix);
    assert( s_f == new_s.get_funders_success());

    int n_rand_cond = 2;

    int first_gen = 1;
    recreate_generation(new_s, first_gen);
    check_funder_equal_pop(new_s.get_pop(), s_f, first_gen);
    auto test_one = test_against_random_conditions(new_s, n_rand_cond, pop_max, amplitude, "test_1");

    int last_gen = n_cycles - 1;
    recreate_generation(new_s, last_gen);
    check_funder_equal_pop(new_s.get_pop(), s_f, last_gen);
    auto test_two = test_against_random_conditions(new_s, n_rand_cond, pop_max, amplitude, "test_2");

    assert(test_one != test_two);
}
void test_evolution_and_selection_for_spores()
{
    int n_cycles = 50;
    meta_param m{n_cycles};
    simulation s{sim_param{env_param{},
                           ind_param{},
                           m,
                           pop_param{}}
                };
    s.evolve_and_select_for_spores();
    assert(!s.is_selecting_only_spores());
    exec(s);
    assert(s.is_selecting_only_spores());
}
void test_simulation()//!OCLINT tests may be many
{
#ifndef NDEBUG

    ///It is possible in the exec function to activate the selection for spores after 50 cycles
    test_evolution_and_selection_for_spores();

    ///It is possible to test any generation of rand_evo against a set of random environment
    test_rand_cond_test_after_rand_evo();

    ///It is possible to create random conditions drawing from extremes of a uniform distribution
    /// The extremes are the intervals between 2 and 3 var * amplitude values
    test_create_rand_cond_unif_extreme();

    ///It is possible to create random matrix that samples from the extremes
    test_rand_cond_matrix_extreme();

    ///The function exec_change() takes random a random condition vector and
    /// changes the environment accordingly throughout the simulation
    test_exec_change();

    ///A matrix of random conditions can be created
    test_rand_cond_matrix();

    ///It is possible to create random conditions drawing from a uniform distribution
    test_create_rand_cond_unif();

    ///It is possible to activate the death mechanic in a simulation(in tick(simulatio&))
    test_activate_death();

    //A simulation can be initialized by getting a sim_parmater class as an argument
    {
        unsigned int pop_size = 1;
        unsigned int exp_new_pop_size = 1;
        double min_dist = 0.1;
        int grid_side = 1;
        double diff_coeff = 0.1;
        double init_food = 1.0;
        double mutation_prob = 0.01;
        double mutation_step = 0.1;
        double base_disp_prob = 0.01;
        double spore_advantage = 10.0;
        double repr_trsh = 0.1;
        double metab_degr_rate = 0.1;
        int n_cycles = 42;
        int cycle_duration = 5;

        sim_param s_p(pop_size,
                      exp_new_pop_size,
                      min_dist,
                      grid_side,
                      diff_coeff,
                      init_food,
                      mutation_prob,
                      mutation_step,
                      base_disp_prob,
                      spore_advantage,
                      repr_trsh,
                      metab_degr_rate,
                      n_cycles,
                      cycle_duration
                      );
        simulation s(s_p);
        //tests for pop_param
        assert(s.get_pop().get_v_ind().size() == pop_size);
        assert(s.get_pop().get_param().get_exp_new_pop_size() == exp_new_pop_size);
        assert(s.get_pop().get_param().get_min_dist() - min_dist < 0.0001 &&
               s.get_pop().get_param().get_min_dist() - min_dist > -0.0001);
        assert(s.get_pop().get_param().get_mu_p() - mutation_prob < 0.0001 &&
               s.get_pop().get_param().get_mu_p() - mutation_prob > -0.0001);
        assert(s.get_pop().get_param().get_mu_st() - mutation_step < 0.0001 &&
               s.get_pop().get_param().get_mu_st() - mutation_step > -0.0001);
        assert(s.get_pop().get_param().get_base_disp_prob() - base_disp_prob < 0.00001 &&
               s.get_pop().get_param().get_base_disp_prob() - base_disp_prob > -0.00001);
        assert(s.get_pop().get_param().get_spo_adv() - spore_advantage < 0.0001 &&
               s.get_pop().get_param().get_spo_adv() - spore_advantage > -0.0001);

        //tests for env param
        assert(s.get_env().get_param().get_diff_coeff() - diff_coeff < 0.00001 &&
               s.get_env().get_param().get_diff_coeff() - diff_coeff > -0.00001);
        assert(s.get_env().get_param().get_init_food() - init_food < 0.0001 &&
               s.get_env().get_param().get_init_food() - init_food > -0.0001);

        assert(s.get_env().get_param().get_degr_rate() - metab_degr_rate < 0.0001 &&
               s.get_env().get_param().get_degr_rate() - metab_degr_rate > -0.0001);
        //test for simulation meta param
        assert(s.get_meta_param().get_n_cycles() == n_cycles);
        assert(s.get_meta_param().get_cycle_duration() == cycle_duration);
    }


    //A simulation has an environment
    // The value -1234567890 is irrelevant: just get this to compile

    {
        simulation s;
        assert(s.get_env().get_param().get_grid_side() > -1234567890);
    }

    //A simulation can be initialized with a certain amount
    //of food in all the cells of its environment grid
    //1 by default
    {
        simulation s;
        for( auto& grid_cell : s.get_env().get_grid())
        {
            assert(grid_cell.get_food() - 1 < 0.000001
                   && grid_cell.get_food() - 1 > -0.0000001);
        }

        double starting_food = 3.14;
        simulation s1 (sim_param{1, 1, 1, 1, 0.1, starting_food});
        for( auto& grid_cell : s1.get_env().get_grid())
        {
            assert(grid_cell.get_food() - starting_food < 0.000001
                   && grid_cell.get_food() - starting_food > -0.0000001);
        }
    }

    //Individuals can take up energy from the environment
    {
        simulation s(sim_param{1,1,0.1,2});
        double total_food_init = std::accumulate
                (
                    s.get_env().get_grid().begin(),
                    s.get_env().get_grid().end(),0.0,
                    [](double sum,const env_grid_cell& c){return sum + c.get_food();}
        );

        double total_en_init = std::accumulate
                (
                    s.get_pop().get_v_ind().begin(),
                    s.get_pop().get_v_ind().end(),0.0,
                    [](double sum,const individual& i){return sum + i.get_energy();}
        );

        feeding(s);

        double total_food_after = std::accumulate
                (
                    s.get_env().get_grid().begin(),
                    s.get_env().get_grid().end(),0.0,
                    [](double sum, const env_grid_cell& c){return sum + c.get_food();}
        );

        double total_en_after =
                std::accumulate(
                    s.get_pop().get_v_ind().begin(),
                    s.get_pop().get_v_ind().end(),0.0,
                    [](double sum,const individual& i){
            return sum + i.get_energy();}
        );

        assert(total_food_init > total_food_after);
        assert(total_en_init < total_en_after);

        auto total_uptake =
                std::accumulate(s.get_pop().get_v_ind().begin(),
                                s.get_pop().get_v_ind().end(), 0.0,
                                [](double sum, const individual& i)
        {return sum + i.get_param().get_uptake_rate();});

        assert(total_en_after - (total_uptake + total_en_init) < 0.00000001 &&
               total_en_after - (total_uptake + total_en_init) > -0.00000001);
        assert(total_food_init - (total_food_after + total_uptake) < 0.000001 &&
               total_food_init - (total_food_after + total_uptake) > -0.000001);
    }

    //Individuals outside the grid do not feed
    {
        simulation s(sim_param{1,1,0.1,2});

        double total_food_init = std::accumulate
                (
                    s.get_env().get_grid().begin(),
                    s.get_env().get_grid().end(),0.0,
                    [](double sum,const env_grid_cell& c){return sum + c.get_food();}
        );

        double total_en_init = std::accumulate
                (
                    s.get_pop().get_v_ind().begin(),
                    s.get_pop().get_v_ind().end(),0.0,
                    [](double sum,const individual& i){return sum + i.get_energy();}
        );

        set_pos(s.get_pop().get_ind(0),std::pair<double,double>(-42,42));

        feeding(s);

        double total_food_after = std::accumulate
                (
                    s.get_env().get_grid().begin(),
                    s.get_env().get_grid().end(),0.0,
                    [](double sum, const env_grid_cell& c){return sum + c.get_food();}
        );

        double total_en_after = std::accumulate
                (
                    s.get_pop().get_v_ind().begin(),
                    s.get_pop().get_v_ind().end(),0.0,
                    [](double sum,const individual& i){
            return sum + i.get_energy();}
        );

        assert(total_food_init - total_food_after < 0.00001
               && total_food_init - total_food_after > -0.00001);
        assert(total_en_init - total_en_after < 0.000001
               && total_en_init - total_en_after > -0.000001);
    }

    //In one tick/timestep individuals take in input,
    //determine phenotype(based on previous timestep),
    // feed, than reproduce, than substances diffuse
    {
        simulation s(sim_param{1,1,0.1,3,1,1,0});
        //Set all the hid nodes and H2O and H2H weights to one so
        //that we are sure the phenotype will stay = active;
        for(auto& ind : s.get_pop().get_v_ind())
        {
            ind.get_grn().set_all_hid_nodes(1);
            ind.get_grn().set_all_out_nodes(1);
            ind.get_grn().set_all_H2O(1);
            ind.get_grn().set_all_H2H(1);
        }

        //The single individual in this population
        //after a tick should reproduce
        auto init_pop_size = s.get_pop().get_v_ind().size();
        auto ind_en = get_ind_tr_en( s.get_pop(), 0)
                + s.get_pop().get_ind(0).get_param().get_metabolic_rate() + 0.01
                - s.get_pop().get_ind(0).get_param().get_uptake_rate();
        s.get_pop().get_ind(0).set_energy(ind_en);

        //and the grid_cell where it is should recieve
        //food nutrients
        response(s);
        feeding(s);

        auto grid_index_ind =  find_grid_index(s.get_pop().get_ind(0),
                                               s.get_env().get_param().get_grid_side());

        double food_after_feeding = s.get_env().get_cell(grid_index_ind).get_food();

        metabolism_pop(s.get_pop());
        if(division(s.get_pop(), s.get_rng()))
        {
            manage_static_collisions(s.get_pop());
        }

        diffusion(s.get_env());
        double food_after_diffusion = s.get_env().get_cell(grid_index_ind).get_food();

        //The difference in food before and after the diffusion step is not 0
        assert(!(food_after_feeding - food_after_diffusion < 0.00001
                 && food_after_feeding - food_after_diffusion > -0.0001));
        assert(s.get_pop().get_v_ind().size() == 2 * init_pop_size);
        assert(!has_collision(s.get_pop() ));
    }

    //If nothing else happens, food should constantly decrease when cells are feeding
    {
        meta_param m{};
        pop_param p{2,
                    1,
                    0.1,
                    0.01,
                    0.1,
                    0.01,
                    10,
                    1
                   };
        ind_param indiviual{};
        env_param e{3,
                    0.1,
                    indiviual.get_treshold_energy() * 10
                   };
        simulation s (sim_param{e, indiviual, m,p});

        auto food_begin =
                std::accumulate(s.get_env().get_grid().begin(),
                                s.get_env().get_grid().end(),
                                0.0,
                                [](double sum, const env_grid_cell& c)
        {return sum + c.get_food();});

        //The simulation will last long enough  for the individuals to reproduce
        auto sim_time = s.get_pop().get_ind(0).get_param().get_treshold_energy() /
                s.get_pop().get_ind(0).get_param().get_uptake_rate() * 2;

        auto init_pop_size =  s.get_pop().get_v_ind().size();

        for( int i = 0; i != static_cast<int>(sim_time); i++)
        {

            auto food_before_feed =
                    std::accumulate(s.get_env().get_grid().begin(),
                                    s.get_env().get_grid().end(), 0.0,
                                    [](double sum, const env_grid_cell& c)
            {return sum + c.get_food();});

            feeding(s);

            auto food_after_feed =
                    std::accumulate(s.get_env().get_grid().begin(),
                                    s.get_env().get_grid().end(), 0.0,
                                    [](double sum, const env_grid_cell& c)
            {return sum + c.get_food();});

            double food_eaten = 0.0;
            for(const auto& ind : s.get_pop().get_v_ind())
            {
                auto grid_cell_ind = find_grid_index(ind, s.get_env().get_param().get_grid_side());
                if( grid_cell_ind != -100 && s.get_env().get_cell(grid_cell_ind).get_food() > 0)
                {
                    food_eaten += ind.get_param().get_uptake_rate();
                }
            }

            auto balance_uptake = food_before_feed - (food_after_feed + food_eaten);
            assert(balance_uptake < 0.0001 && balance_uptake > -0.0001);

            metabolism_pop(s.get_pop());
            division(s.get_pop(), s.get_rng());
            manage_static_collisions(s.get_pop());

            auto food_before_diff =
                    std::accumulate(s.get_env().get_grid().begin(),
                                    s.get_env().get_grid().end(), 0.0,
                                    [](double sum, const env_grid_cell& c)
            {return sum + c.get_food();});

            calc_diffusion_food(s.get_env());

            auto new_grid = s.get_env().get_grid();

            auto food_after_diff =
                    std::accumulate(new_grid.begin(),
                                    new_grid.end(), 0.0,
                                    [](double sum, const env_grid_cell& c)
            {return sum + c.get_food();});

            auto balance_diffusion = food_before_diff - food_after_diff;
            assert(balance_diffusion < 0.0001 && balance_diffusion > -0.0001);

            new_grid.swap(s.get_env().get_grid());

        }
        auto food_end =
                std::accumulate(s.get_env().get_grid().begin(),
                                s.get_env().get_grid().end(), 0.0,
                                [](double sum, const env_grid_cell& c)
        {return sum + c.get_food();});

        assert(food_end < food_begin);

        auto final_pop_size =  s.get_pop().get_v_ind().size();

        assert(init_pop_size < final_pop_size);
    }

    //A simulation is initiallized with a degradation rate
    {
        double degradation_rate = 0.314;
        simulation s(sim_param{0,0,0,0,1,0,0,0,0,0,0,degradation_rate});
        assert(s.get_env().get_param().get_degr_rate() - degradation_rate < 0.000001 &&
               s.get_env().get_param().get_degr_rate() - degradation_rate > -0.000001);
    }
    //Every time step the individuals produce new metabolite and metabolite degrades in grid_cells
    {
        double degradation_rate = 0.314;
        double init_metab = degradation_rate;
        simulation s(sim_param{1,0,0,1,1,0,0,0,0,0,0,degradation_rate});

        for(auto& grid_cell : s.get_env().get_grid())
        {
            grid_cell.set_metab(init_metab);
        }

        double tot_metab_before =
                std::accumulate(s.get_env().get_grid().begin(),
                                s.get_env().get_grid().end(), 0.0,
                                [](int sum, const env_grid_cell& g)
        {return sum + g.get_metab();});

        secretion_metabolite(s);
        degradation_metabolite(s.get_env());

        double tot_metab_after =
                std::accumulate(s.get_env().get_grid().begin(),
                                s.get_env().get_grid().end(), 0.0,
                                [](int sum, const env_grid_cell& g)
        {return sum + g.get_metab();});

        double tot_production =
                std::accumulate(s.get_pop().get_v_ind().begin(),
                                s.get_pop().get_v_ind().end(), 0.0,
                                [](int sum, const individual& i)
        {return sum + i.get_param().get_metab_secr_rate();});

        double tot_degradation =
                s.get_env().get_grid_size() * s.get_env().get_param().get_degr_rate();

        auto metab_balance =
                tot_metab_before - tot_degradation + tot_production - tot_metab_after;

        assert(metab_balance < 0.000001 && metab_balance > -0.000001);


    }

    //every timestep/tick collisions are handled
    {
        simulation s(sim_param{7,3});
        //The central individual in this population
        //after a tick should reproduce
        auto init_pop_size = s.get_pop().get_v_ind().size();
        s.get_pop().get_ind(1).set_energy(get_ind_tr_en(s.get_pop(), 1)
                                          + s.get_pop().get_ind(1).get_param().get_metabolic_rate()
                                          + 0.01
                                          - s.get_pop().get_ind(1).get_param().get_uptake_rate());
        feeding(s);
        metabolism_pop(s.get_pop());
        division(s.get_pop(), s.get_rng());
        manage_static_collisions(s.get_pop());

        assert(s.get_pop().get_v_ind().size() == init_pop_size + 1);
        assert(!has_collision(s.get_pop()));
    }

    //A simulation is initialized with a m_tick = 0;
    {
        simulation s;
        assert(s.get_timestep() == 0);
    }

    //After each tick the simulation updates its m_tick
    //by one
    {
        simulation s{};
        for(int i = 0; i != 3; i++)
        {
            assert(s.get_timestep() == i);
            tick(s);
        }
    }

    //If sparse collision resolution is used, collisions are only checked every n timesteps
    {
        //create a simulation with 100 individuals to ensure there will be collisions early on
        pop_param p{100};
        ind_param ind{};
        env_param e{100};
        meta_param m{};
        simulation s{sim_param{e, ind, m, p}};
        int n_ticks = 2;
        int n_total_ticks = 4;

        //create collisions
        while(!has_collision(s.get_pop()))
        {
            response(s);
            jordi_feeding(s);
            metabolism_pop(s.get_pop());
            secretion_metabolite(s);
            jordi_death(s.get_pop());
            division(s.get_pop(), s.get_rng());
        }

        //check that collisions are resolved only every n_ticks
        for(int i = 0; i != n_total_ticks; i++)
        {
            tick(s, n_ticks);

            //the timestep is updated in tick so you need to check
            //for timestep - 1
            if((s.get_timestep() - 1) % n_ticks == 0)
            {
                assert(!has_collision(s.get_pop()));
                //create collisions
                while(!has_collision(s.get_pop()))
                {
                    response(s);
                    jordi_feeding(s);
                    metabolism_pop(s.get_pop());
                    secretion_metabolite(s);
                    jordi_death(s.get_pop());
                    division(s.get_pop(), s.get_rng());
                }
            }
            else
            {
                assert(has_collision(s.get_pop()));
            }
        }


    }
    //A simulation runs a cycle for a certain amount of ticks stated in it s parameters
    {
        auto n_cycles = 1;
        auto cycle_duration = 1;
        meta_param m{ n_cycles, cycle_duration};
        sim_param sp{env_param{}, ind_param{}, m, pop_param{}};
        simulation s{sp};
        exec_cycle(s);
        assert(s.get_timestep() == cycle_duration);
    }

    //A simulation can be run for the amount of cycles stated in its parameters
    {
        auto n_cycles = 1;
        auto cycle_duration = 1;
        pop_param p{100};
        meta_param m{ n_cycles, cycle_duration};
        sim_param sp{env_param{}, ind_param{}, m, pop_param{}};
        simulation s{sp};
        exec(s);
        assert(s.get_cycle() == n_cycles);
    }

    //After the exec_cycle the time_step timer is reset
    {
        auto n_cycles = 1;
        auto cycle_duration = 1;
        pop_param p{100};
        meta_param m{ n_cycles, cycle_duration};
        sim_param sp{env_param{}, ind_param{}, m, pop_param{}};
        simulation s{sp};
        //As the parameters indicate only one cycle will be executed
        exec(s);
        assert(s.get_cycle() == n_cycles);
        assert(s.get_timestep() == 0);
    }

    //After each cycle a new population(max 100 individuals
    //is created from the previous one(see population tests)
    //and the environment is reset
    //(similar to dispersal() test)
    {
        auto n_cycles = 1;
        auto cycle_duration = 1;
        meta_param m{ n_cycles, cycle_duration};
        pop_param p{100};
        ind_param ind{};
        env_param e{3};
        sim_param sp{e, ind, m, p};
        simulation s{sp};
        auto init_env = s.get_env();
        //Run a few ticks to make sure env is different from original
        for(int i = 0; i != 1; i++)
        {
            tick(s);
        }
        assert(s.get_env() != init_env);
        //Reset tick timer to avoid crushing on assert  from prepare_funders()
        s.reset_timesteps();
        //As the parameters indicate only one cycle will be executed
        exec(s);
        assert(s.get_pop().get_pop_size() - s.get_pop().get_param().get_exp_new_pop_size() == 0);
        assert(s.get_env() == init_env);
    }

    //Spores do not feed
    {
        simulation s;
        auto init_food = s.get_env().get_cell(0).get_food();
        s.get_pop().get_ind(0).set_phen(phenotype::spore);
        feeding(s);
        assert(init_food - s.get_env().get_cell(0).get_food() < 0.000001
               && init_food - s.get_env().get_cell(0).get_food() > -0.000001);
    }


    //Sporulating individuals do not feed but they lose energy
    {
        simulation s;
        s.get_pop().get_ind(0).set_phen(phenotype::sporulating);
        set_ind_en(s.get_pop().get_ind(0),1);
        auto init_en_ind0 = get_ind_en(s.get_pop(), 0);
        auto init_food = s.get_env().get_cell(0).get_food();
        feeding(s);
        metabolism_pop(s.get_pop());
        assert(init_en_ind0 - get_ind_en(s.get_pop(), 0) -
               s.get_pop().get_ind(0).get_param().get_spor_metabolic_rate() < 0.000001
               &&
               init_en_ind0 - get_ind_en(s.get_pop(), 0) -
               s.get_pop().get_ind(0).get_param().get_spor_metabolic_rate() > -0.000001);

        assert(init_food - s.get_env().get_cell(0).get_food() < 0.000001
               && init_food - s.get_env().get_cell(0).get_food() > -0.000001);
    }

    //Max 100 ind, selected based on phenotype, are placed in a hex pattern,
    //in a new env after dispersal
    //Tests all of the above
    {
        unsigned int pop_size = 1000;
        unsigned int new_pop_size = 100;
        auto food = 42.1;
        auto metabolite = 42.1;
        simulation s(sim_param{pop_size,new_pop_size});
        for(auto& grid_cell : s.get_env().get_grid())
        {
            grid_cell.set_food(food);
            grid_cell.set_metab(metabolite);
        }
        environment ref_env = s.get_env();
        dispersal(s);
        //Max 100 ind
        assert(s.get_pop().get_v_ind().size() == new_pop_size);
        //Hex pattern
        auto n_hex_l = count_hex_layers(s.get_pop().get_pop_size());
        auto v_modulus = modulus_of_btw_ind_angles(s.get_pop(), M_PI/ (6 * n_hex_l));
        for(auto ind_modulus : v_modulus)
        {
            assert( ind_modulus < 0.0000000001 ||
                    (ind_modulus > M_PI / (6 * n_hex_l) - 0.1 &&
                     ind_modulus <= M_PI / (6 * n_hex_l) + 0.1));
        }
        assert(!has_collision(s.get_pop()));
        //Reset env
        assert( s.get_env() != ref_env );
        auto init_food = s.get_env().get_param().get_init_food();
        for(const auto& grid_cell : s.get_env().get_grid())
        {
            assert(grid_cell.get_food() - init_food <  0.000001 &&
                   grid_cell.get_food() - init_food >  -0.000001);
            assert(grid_cell.get_metab() <  0.000001 &&
                   grid_cell.get_metab() >  -0.000001);
        }
    }

    //All inds have the same size after dispersal
    {

        int time = 50;
        env_param e;
        ind_param ind;
        meta_param m;
        pop_param p;
        sim_param sp{e, ind, m, p};
        simulation s{sp};

        for( int i = 0; i != time; i++)
        {
            tick(s);
        }

        //        assert( s.get_pop().get_v_ind().end() !=
        //                std::adjacent_find( s.get_pop().get_v_ind().begin(),
        //                                    s.get_pop().get_v_ind().end(),
        //                                    [](const auto& lhs, const auto& rhs)
        //        {return lhs.get_radius() != rhs.get_radius();}));

        dispersal(s);

        assert( s.get_pop().get_v_ind().end() ==
                std::adjacent_find( s.get_pop().get_v_ind().begin(),
                                    s.get_pop().get_v_ind().end(),
                                    [](const auto& lhs, const auto& rhs)
        {return lhs.get_radius() != rhs.get_radius();}));
    }

    //All individuals in a simulation can respond to their environment and surrounding
    {
        double food_amount = 3.14;
        double metabolite_amount = 3.14;
        double energy_amount = 3.14;
        simulation s(sim_param{2,1,0.1,4});
        for(auto & grid_cell : s.get_env().get_grid())
        {
            grid_cell.set_food(food_amount);
            grid_cell.set_metab(metabolite_amount);
        }
        for(auto& ind : s.get_pop().get_v_ind())
        {
            ind.set_energy(energy_amount);
            //let's set all the weights of the network to 1
            //in this case we expect that the outputs will be one
            ind.get_grn().set_all_I2H(1);
            ind.get_grn().set_all_H2O(1);
        }
        //Since the networks react to input of timestpe t
        //at timestep t+1 I will run responds(i) 2 times
        //In this first response call, since the hidden nodes are 0
        //The output will be 0 and therefore the individuals phenotype
        //will be phenotype::sporulating
        response(s);
        for(auto& ind : s.get_pop().get_v_ind())
        {
            assert(is_active(ind));
        }

        //The values of food and metabolite in the grid_cell
        //as well as the energy in the individual
        //are changed so that all inputs are -1,
        //the network should therefore give outputs == 0/false
        //since the outputs node will recieve a signal
        //that is negative(below the treshold)
        //after responds(i) is called 2 more times
        for(auto& grid_cell : s.get_env().get_grid())
        {
            grid_cell.set_food(-1);
            grid_cell.set_metab(-1);
        }
        for(auto& ind : s.get_pop().get_v_ind())
        {
            ind.set_energy(-1);
        }
        //1st respose, the individuals respond to the initial parameter
        //expected output value == 1
        response(s);
        for(auto& ind : s.get_pop().get_v_ind())
        {
            assert(is_active(ind));
        }
        //2nd response, the individual responds to the changed parameter(all 0s)
        //expected output value == 0
        response(s);
        for(auto& ind : s.get_pop().get_v_ind())
        {
            assert(!is_active(ind));
            assert(is_sporulating(ind));
        }
    }

    //It is possible to extract the demographic state of a population
    {
        env_param e;
        ind_param i;
        meta_param m;
        pop_param p{0};
        sim_param sp{e,i,m,p};
        simulation s{sp};
        assert(s.get_pop().get_pop_size() == 0);

        int n_spores = 2;
        int n_sporulating = 3;
        int n_actives = 4;
        individual ind{ind_param{}};

        ind.set_phen(phenotype::spore);
        for(int i = 0; i != n_spores; i++)
        {
            s.get_pop().get_v_ind().push_back(ind);
        }

        ind.set_phen(phenotype::sporulating);
        for(int i = 0; i != n_sporulating; i++)
        {
            s.get_pop().get_v_ind().push_back(ind);
        }

        ind.set_phen(phenotype::active);
        for(int i = 0; i != n_actives; i++)
        {
            s.get_pop().get_v_ind().push_back(ind);
        }

        assert(s.get_pop().get_pop_size() == n_spores + n_sporulating + n_actives);

        demographic_cycle d_c = demographics(s, env_param{});

        assert(d_c.get_n_spores() == n_spores);
        assert(d_c.get_n_sporulating() == n_sporulating);
        assert(d_c.get_n_actives() == n_actives);
    }

    ///It is possible to set inds in a pop
    /// from a funders object so that pop has
    /// the same number of individuals
    /// as the ones in the funder object
    /// as well as with the same neural networks
    {
        int seed = 123;
        int change_freq = 4567;
        int funders_generation = 1;
        meta_param m{3,
                     2,
                     seed,
                             change_freq};

        env_param e;
        ind_param i;
        pop_param p;
        sim_param sp{e,i,m,p};
        simulation s{sp};

        exec(s);
        save_data(s);
        population pop{};

        auto funders_success = load_funders_success(create_funders_success_name(seed, change_freq));
        auto selected_funders = funders_success.get_v_funders()[funders_generation];

        auto demo_sim = load_demographic_sim(create_sim_demo_name(seed,change_freq));
        auto selected_conditions = demo_sim.get_demo_cycles()[funders_generation];

        pop.set_pop_inds(pop_from_funders(funders_success, demo_sim, funders_generation));

        for( int i = 0; i != pop.get_pop_size(); i++)
        {
            assert(pop.get_ind(i).get_param() == selected_conditions.get_ind_param());
            assert(pop.get_ind(i).get_grn() == selected_funders.get_v_funder_data()[i].get_grn());
        }

    }

    //It is possible to store the demographics
    //of the population contained in a simulation
    //at a certain point in time in a vector
    {
        simulation s;
        auto demo_sim_length = s.get_demo_sim().get_demo_cycles().size();
        store_demographics(s);
        auto demo_sim_length2 = s.get_demo_sim().get_demo_cycles().size();
        assert(demo_sim_length != demo_sim_length2);
        assert(demo_sim_length + 1 == demo_sim_length2);
    }

    //At the end of each cycle the demographics are stored
    {
        int n_cycles = 2;
        int cycle_duration = 1;
        meta_param m{n_cycles, cycle_duration};
        env_param e;
        pop_param p;
        ind_param i;
        sim_param s_p{e, i, m, p};

        simulation s{s_p};
        auto init_demo_cycle = s.get_demo_sim().get_demo_cycles().size();
        assert(init_demo_cycle == 0);
        exec_cycle(s);
        auto one_demo_cycle = s.get_demo_sim().get_demo_cycles().size();
        assert(one_demo_cycle == 1);
    }

    //At the end of the simulation the demographic is saved in a file
    {
        int n_cycles = 2;
        int cycle_duration = 1;
        meta_param m{n_cycles, cycle_duration};
        env_param e;
        pop_param p;
        ind_param i;
        sim_param s_p{e, i, m,p};

        simulation s{s_p};
        std::string expected_file_name = create_sim_demo_name(s);
        exec(s);
        save_data(s);

        assert(exists(expected_file_name));
        auto d_s = load_demographic_sim(expected_file_name);
        assert(d_s == s.get_demo_sim());

    }

    //A simulation is initialized with a funders_success  object
    {
        simulation s;

        assert(s.get_funders_success().get_v_funders().size() == 0);
    }

    //It is possible to store the ancestor_ID and GRN of the
    //individuals composing a population at the start of the cycle
    //updating the funders_success member of simulation
    {
        int n_cycles = 1;
        int cycle_duration = 50;
        meta_param m{n_cycles, cycle_duration};
        env_param e;
        pop_param p;
        ind_param i;
        sim_param s_p{e, i, m, p};

        simulation s{s_p};
        auto zero_cycle_funders = prepare_funders(s);
        exec_cycle(s);
        s.reset_timesteps();
        auto first_cycle_funders = prepare_funders(s);
        assert(zero_cycle_funders != first_cycle_funders);

    }

    //Every cycle the funders are stored
    {
        int n_cycles = 2;
        int cycle_duration = 50;
        meta_param m{n_cycles, cycle_duration};
        env_param e;
        pop_param p;
        ind_param i;
        sim_param s_p{e, i, m, p};

        simulation s{s_p};
        //the first cycle is funded by a single individual
        exec_cycle(s);
        s.reset_timesteps();
        auto first_cycle_funders_size = s.get_funders_success().get_v_funders().size();
        exec_cycle(s);
        s.reset_timesteps();
        auto second_cycle_funder_size = s.get_funders_success().get_v_funders().size();

        assert(first_cycle_funders_size != second_cycle_funder_size);
        assert(s.get_funders_success().get_v_funders().begin() !=
                s.get_funders_success().get_v_funders().begin() + 1);
    }

    //The success of each funder can be calculated
    //The success of a funder is equal to:
    // the fractions of individuals with its same ancestor_ID
    // over the total number of individuals of the population
    // considering the spore advantage in fitness
    {
        int n_cycles = 1;
        int cycle_duration = 1;
        meta_param m_p{n_cycles, cycle_duration};

        unsigned int pop_size = 3;
        pop_param p_p{pop_size};
        ind_param ind;
        sim_param s_p{env_param{}, ind, m_p, p_p};
        simulation s{s_p};

        add_new_funders(s);

        auto funders_with_success = calc_funders_success(s);

        double funder_total_success =
                std::accumulate(funders_with_success.get_v_funder_data().begin(),
                                funders_with_success.get_v_funder_data().end(), 0.0,
                                [](double sum, const funder_data& f)
        {return sum + f.get_success();});

        double total_fitness =
                std::accumulate(s.get_pop().get_v_ind().begin(),
                                s.get_pop().get_v_ind().end(),
                                0.0,
                                [&s](int sum, const individual& i)
        {return sum + (is_spore(i) ? s.get_pop().get_param().get_spo_adv() : 1);});

        total_fitness /= s.get_pop().get_pop_size();

        assert(funder_total_success - total_fitness < 0.0001
               && funder_total_success - total_fitness > -0.0001);
    }

    ///The success of each funder is calculated at the end of each cycle
    {
        int n_cycles = 1;
        int cycle_duration = 1;
        meta_param m_p{n_cycles, cycle_duration};

        unsigned int pop_size = 3;
        pop_param p_p{pop_size};
        ind_param i;
        sim_param s_p{env_param{}, i, m_p, p_p};
        simulation s{s_p};

        auto pre_cycle_funders = prepare_funders(s);
        exec_cycle(s);


        assert(s.get_funders_success().get_v_funders().back() !=
                pre_cycle_funders);

        const auto& funders =
                s.get_funders_success().get_v_funders().back().get_v_funder_data();

        for(const auto& funder : funders)
            assert(funder.get_success() - 1.0/pop_size < 0.00001
                   && funder.get_success() - 1.0/pop_size > -0.00001);
    }

    //At the end of the simulation the funders_success is saved in a file
    {
        int n_cycles = 5;
        int cycle_duration = 1;
        meta_param m{n_cycles, cycle_duration};
        env_param e;
        pop_param p;
        ind_param i;
        sim_param s_p{e, i, m, p};

        simulation s{s_p};

        exec(s);
        save_data(s);

        std::string expected_file_name = create_funders_success_name(s);
        assert(exists(expected_file_name));
        auto f_s = load_funders_success(expected_file_name);
        assert(f_s == s.get_funders_success());

    }

    //A simulation with a certain seed in metaparameters initializes the
    //Random number generator of population to that seed
    {
        int seed = 42;
        meta_param m{
            1,
            1,
            seed
        };
        env_param e;
        pop_param p;
        ind_param i;
        sim_param s_p{e, i, m, p};
        simulation s{s_p};
        auto ref_rng = std::minstd_rand();
        ref_rng.seed(seed);
        assert(s.get_pop().get_rng() == ref_rng);
    }

    //Every cycle multiple of the change_frequency metaparameter
    // environmental parameters change
    {

        env_param e{};
        pop_param p;
        ind_param i;
        int change_frequency = 1;
        meta_param m{change_frequency + 1,
                    1,
                    1,
                    change_frequency
                    };
        sim_param s_p{e, i, m, p};
        simulation s{s_p};

        exec(s);
        assert(e != s.get_env().get_param());
    }

    //Env param can be changed in a simulation
    //Based on env params
    {
        auto mean_diff_coeff = 0.1;
        auto mean_degr_coeff = 0.1;
        auto var_degr_rate = mean_degr_coeff / 3 - 0.01;
        auto var_diff_coeff = mean_diff_coeff / 3 - 0.01;
        env_param e_p{1, 0.1, 0.1, 0.1,
                      mean_diff_coeff,
                              mean_degr_coeff,
                              var_diff_coeff,
                              var_degr_rate};
        pop_param p;
        meta_param m;
        ind_param ind;
        sim_param s_p{e_p, ind, m, p};
        simulation s{s_p};
        int repeats = 100;

        for(int i  = 0 ; i != repeats; i++)
        {
            auto prev_env = s.get_env();
            double prev_deg_rate = prev_env.get_param().get_degr_rate();
            double prev_diff_coeff =  prev_env.get_param().get_diff_coeff();

            change_env(s);

            auto deg_rate = s.get_env().get_param().get_degr_rate();
            auto diff_coeff = s.get_env().get_param().get_diff_coeff();
            assert(prev_deg_rate != deg_rate);
            assert(prev_diff_coeff != diff_coeff);
        }
    }

    //Env_param and ind_param changes with a frequency dictated by the metaparameters
    //frequency_of_change value
    {
        auto frequency_of_change = 2;
        auto n_cycles = frequency_of_change +1;
        //environment will only change once

        meta_param m{n_cycles,
                    1,
                    1,
                    frequency_of_change};

        auto mean_diff_coeff = 0.1;
        auto mean_degr_coeff = 0.1;
        auto var_degr_rate = mean_degr_coeff / 3 - 0.01;
        auto var_diff_coeff = mean_diff_coeff / 3 - 0.01;
        env_param e{1, 0.1, 0.1, 0.1,
                    mean_diff_coeff,
                            mean_degr_coeff,
                            var_diff_coeff,
                            var_degr_rate};

        pop_param p;
        ind_param i;
        sim_param s_p{e, i, m, p};
        simulation s{s_p};
        auto prev_env = s.get_env();
        auto prev_pop = s.get_pop().get_v_ind();
        exec(s);
        assert(prev_env != s.get_env());
        assert(prev_pop != s.get_pop().get_v_ind());
    }

    //Each cycle env_param and ind_param will be stored in the demographic cycle
    //That will be use to store the particular parameters that can change throughout the simulation
    {
        simulation s;
        demographic_cycle d_c = demographics(s, s.get_env().get_param());
        assert(d_c.get_ind_param() == s.get_pop().get_v_ind().begin()->get_param());
        assert(d_c.get_env_param() == s.get_env().get_param());
        change_env(s);
        change_pop(s);
        demographic_cycle d_c1 = demographics(s, s.get_env().get_param());
        assert(d_c != d_c1);
    }


    //It is possible to save and load a random conditions in a vector
    {
        //Make a random conditions .csv file
        env_param env;
        ind_param ind;
        std::minstd_rand rng;
        int n_conditions = 2;
        int seed = 0;
        double amplitude = 1.1;
        std::string filename = create_name_vec_rand_cond(n_conditions, amplitude, seed);

        auto random_conditions = create_rand_conditions_unif(env, ind, n_conditions, amplitude, seed);
        save_vector_of_rand_cond(random_conditions, filename);
        auto random_conditions_loaded = load_random_conditions(filename);
        assert( random_conditions_loaded == random_conditions);

    }

    //It is possible to run a population against multiple random conditions
    //And store their demographics in a file
    {
        ///Run against the random conditions
        pop_param p{100};
        env_param e{100};
        ind_param i;
        int seed_ID = 1234567890;

        meta_param m{1,
                     1,
                     seed_ID};

        simulation s{sim_param{e, i, m, p}};

        double amplitude = 1;
        int repeats = 2;
        int pop_max = pow(10,4);

        auto rand_conditions_vector = create_rand_conditions_unif_extreme(
                    s.get_env().get_param(),
                    s.get_pop().get_v_ind().begin()->get_param(),
                    repeats,
                    amplitude);

        std::string filename = create_test_random_condition_name(s,amplitude);

        auto demo_sim = test_against_random_conditions(s, repeats, pop_max, amplitude, filename);

        assert(std::all_of(rand_conditions_vector.begin(), rand_conditions_vector.end(),
                            [&demo_sim](const std::pair<env_param, ind_param>& r) {
            return std::find_if(demo_sim.get_demo_cycles().begin(),
                             demo_sim.get_demo_cycles().end(),
                             [r](const demographic_cycle& d)
            {return r.first == d.get_env_param() && r.second == d.get_ind_param();}) !=  demo_sim.get_demo_cycles().end();
        }
        )
               );
//        assert(std::equal(rand_conditions_vector.begin(),rand_conditions_vector.end(),
//                          demo_sim.get_demo_cycles().begin(),
//                          [](const std::pair<env_param, ind_param>& r,
//                          const demographic_cycle& d)
//        {return r.first == d.get_env_param() && r.second == d.get_ind_param();})
//               );
        assert(exists(filename));
        auto demo_sim_load = load_demographic_sim(filename);
        assert(demo_sim == demo_sim_load);

    }
    //It is possible to change param of a simulation with other
    //taken from a random condition vector
    {

        //Create random conditions
        env_param env;
        ind_param ind;
        std::minstd_rand rng;
        std::vector<std::pair<env_param,ind_param>> rand_conditions;
        int repeats = 2;

        for(int i = 0; i != repeats; i++)
        {
            rand_conditions.push_back({change_env_param_norm(env, rng),
                                       change_ind_param_norm(ind, rng)});
        }

        simulation s;
        for(const auto& rand_condition : rand_conditions)
        {
            assert(s.get_env().get_param() != rand_condition.first
                    && s.get_pop().get_v_ind().begin()->get_param() != rand_condition.second);
            set_new_params(s,rand_condition.first, rand_condition.second);
            assert(s.get_env().get_param() == rand_condition.first
                   && s.get_pop().get_v_ind().begin()->get_param() == rand_condition.second);
        }
    }

    //A simulation can be instantiated given another simulation.
    //With the same population and environment,
    //but with empty data vectors for demographics
    {
        env_param e{5};
        pop_param p{2};
        ind_param i;
        meta_param m{2,1,5,2};
        simulation s{sim_param{e, i, m, p}};
        exec_cycle(s);
        assert(!s.get_demo_sim().get_demo_cycles().empty());
        simulation new_s = no_dem_and_fund_copy(s);
        assert(s.get_env() == new_s.get_env());
        assert(s.get_pop() == new_s.get_pop());
        assert(s.get_meta_param() == new_s.get_meta_param());
        assert(s.get_demo_sim() != new_s.get_demo_sim());
        assert(new_s.get_demo_sim().get_demo_cycles().empty());
    }

    /// It is possible to create and save a certain amount of random conditions
    /// giving the base env_ and ind_ param, and how much larger the range of variation
    /// is going to be proportionally.
    /// It is possible also to decide the seed
    /// of the random number generator for consistency across
    /// different runs.
    {
        int n_random_conditions = 2;
        int seed = 55;
        env_param e;
        ind_param i;
        auto amplitude = 1.0; //the range will stay the same
        auto rand_cond = create_rand_conditions_unif(e,
                                                     i,
                                                     n_random_conditions,
                                                     amplitude,
                                                     seed);

        auto rand_cond2 = create_rand_conditions_unif(e,
                                                      i,
                                                      n_random_conditions,
                                                      amplitude,
                                                      seed);

        assert(static_cast<unsigned int>(n_random_conditions) == rand_cond.size());
        assert(rand_cond == rand_cond2);
    }

    ///The final population is saved at the end of the run_sim_evo(s) function
    ///as funders
    {
        ind_param i;
        env_param e{5};
        meta_param m{1,
                     1};
        pop_param p{100,
                    100};
        simulation s{sim_param{e, i, m, p}};

        exec(s);
        save_data(s);

        std::string expected_filename = create_last_pop_name(s);
        assert(prepare_funders(s) == load_funders(expected_filename));
    }

    ///Simulation parameters are saved by the end of run_sim function
    {
        ind_param i;
        env_param e{5};
        meta_param m{1,
                     1};
        pop_param p{100,
                    100};
        sim_param s_p{e, i, m, p};
        simulation s{s_p};

        exec(s);
        save_data(s);

        std::string expected_filename = create_sim_par_name(s);
        assert(s_p == load_sim_parameters(expected_filename));
    }

    ///A simulation can be loaded given the seed and the change freq
    ///By loading the sim_params of a given simulation
    /// and instantiating the last population of that given simulation
    {
        ind_param i;
        env_param e{5};
        int seed = 23;
        int change_freq = 21;
        meta_param m{1,
                     1,
                     seed,
                             change_freq};
        pop_param p{100,
                    100};
        sim_param s_p{e, i, m, p};
        simulation s{s_p};

        exec(s);
        save_data(s);

        simulation s1 = load_sim(seed,change_freq);
        update_radius_pop(s1.get_pop());
        place_start_cells(s1.get_pop());

        assert(s.get_pop().get_v_ind() == s1.get_pop().get_v_ind());
        assert(s.get_pop().get_param() == s1.get_pop().get_param());
        assert(s.get_env().get_param() == s1.get_env().get_param());
        assert(s.get_meta_param() == s1.get_meta_param());
        assert(s.get_env() == s1.get_env());
        assert(s.get_funders_success() == s1.get_funders_success());
        assert(s.get_demo_sim() == s1.get_demo_sim());
    }

    ///A simulation last population of funders can be loaded given the seed and the change freq
    ///By loading the sim_params of a given simulation
    /// and instantiating the last population of that given simulation
    {
        ind_param i;
        env_param e{5};
        int seed = 23;
        int change_freq = 21;
        meta_param m{1,
                     50,
                     seed,
                             change_freq};
        pop_param p{1,
                    100};
        sim_param s_p{e, i, m,p};
        simulation s{s_p};

        exec(s);
        save_data(s);

        simulation s1 = load_sim_last_pop(seed,change_freq);
        place_start_cells(s.get_pop());
        assert(s.get_pop().get_v_ind() == s1.get_pop().get_v_ind());
        assert(s.get_pop().get_param() == s1.get_pop().get_param());
        assert(s.get_env().get_param() == s1.get_env().get_param());
        assert(s.get_meta_param() == s1.get_meta_param());
        assert(s.get_env() == s1.get_env());
    }

    /// It is possible to load a population of a number of individuals
    /// equal to the expexcted new population size specified in pop_param
    /// made only of the best final network of a population
    {
        ind_param i;
        pop_param p{1,100};
        int seed = 123;
        int change_freq = 10;
        meta_param m{2,50,seed,change_freq};
        env_param e{};
        simulation s{sim_param{e, i, m, p}};

        exec(s);
        save_data(s);

        auto filename = create_funders_success_name(s);
        auto fund_succ = load_funders_success(filename);
        assert(fund_succ == s.get_funders_success());
        auto best_grn = find_last_gen_best_ind_grn(fund_succ);
        auto rand_sim = load_best_ind_for_rand_cond(seed, change_freq);
        for(const auto& individual : rand_sim.get_pop().get_v_ind())
        {
            assert(individual.get_grn() == best_grn);
        }
    }

    ///It is possible to create a simulation where
    /// sim_parameters, funders_success and sim_demographics
    /// are loaded from another simulation
    {
        int seed = 123;
        int change_freq = 10;
        //The simualtion is already run in the file above
        auto funders_filename = create_funders_success_name(seed, change_freq);
        assert(exists(funders_filename));
        auto demo_sim_file_name = create_sim_demo_name(seed, change_freq);
        assert(exists(demo_sim_file_name));
        auto sim_param_filename = create_sim_par_name(seed, change_freq);
        assert(exists(sim_param_filename));


        auto fund_succ = load_funders_success(funders_filename);
        auto sim_dem = load_demographic_sim(demo_sim_file_name);
        auto sim_par = load_sim_parameters(sim_param_filename);

        simulation s = load_sim_no_pop(seed, change_freq);
        auto sp = sim_param{s.get_env().get_param(),
                s.get_pop().get_ind(0).get_param(),
                s.get_meta_param(),
                s.get_pop().get_param()};

        assert(s.get_funders_success() == fund_succ);
        assert(s.get_demo_sim() == sim_dem);
        assert(sp == sim_par );
    }

    ///It is possible to create a population
    /// with the same number of individuals
    /// as the ones in a funder object
    /// as well as with the same neural networks
    {
        int seed = 123;
        int change_freq = 10;
        int funders_generation = 1;

        auto funders_success = load_funders_success(create_funders_success_name(seed, change_freq));
        auto selected_funders = funders_success.get_v_funders()[funders_generation];

        auto demo_sim = load_demographic_sim(create_sim_demo_name(seed,change_freq));
        auto selected_conditions = demo_sim.get_demo_cycles()[funders_generation];

        auto s = load_sim_no_pop(seed, change_freq);

        reproduce_cycle(s, funders_generation);

        assert(s.get_env() == environment{selected_conditions.get_env_param()});

        for( int i = 0; i != s.get_pop().get_pop_size(); i++)
        {
            assert(s.get_pop().get_ind(i).get_grn() == selected_funders.get_v_funder_data()[i].get_grn());
        }

    }

    ///A cycle can be run until population reaches a certain number of individuals
    {
        int pop_max = 2;
        ind_param i;
        env_param e;
        pop_param p;
        meta_param m{1,100,1,1,
                     pop_max,
                             0};
        simulation s{sim_param{e, i, m, p}};
        exec_cycle(s);
        //The cycle will stop before it reaches the max
        //number of timesteps
        assert(s.get_timestep() < s.get_meta_param().get_cycle_duration());

    }

    /// When save_data is called
    /// the last and before last pop of funders is saved
    {
        auto seed = 4242;
        auto change_freq = 4848;
        env_param e;
        pop_param p;
        ind_param ind;
        meta_param m{1,50,
                     seed,
                             change_freq
                    };

        simulation s{sim_param{e, ind, m, p}};

        int n_cycles = 3;
        for(int i = 0; i != n_cycles; i++)
        {
            exec_cycle(s);
            s.reset_timesteps();
        }
        save_data(s);

        auto last_pop = load_sim_last_pop(seed, change_freq);
        auto before_last_pop = load_sim_before_last_pop(seed, change_freq);

        //Check last_pop
        assert( std::equal(last_pop.get_pop().get_v_ind().begin(),
                           last_pop.get_pop().get_v_ind().end(),
                           s.get_funders_success().get_v_funders().back().get_v_funder_data().begin(),
                           [](const individual& i, const funder_data& funder)
        {return i.get_grn() == funder.get_grn();}
        ));
        //Check before_last_pop
        assert( std::equal(before_last_pop.get_pop().get_v_ind().begin(),
                           before_last_pop.get_pop().get_v_ind().end(),
                           s.get_funders_success().get_v_funders().rbegin()[1].get_v_funder_data().begin(),
                [](const individual& i, const funder_data& funder)
        {return i.get_grn() == funder.get_grn();}
        ));

    }

    ///It is possible to add a prefix to the files
    /// where befoer last and last population are saved
    {
        env_param e;
        ind_param i;
        meta_param m{3,
                     20};
        pop_param p{2};
        sim_param sp{e,i,m,p};
        simulation s{sp};

        std::string prefix_3_cycles = "3_cycles";
        std::string prefix_6_cycles = "6_cycles";
        std::string expexcted_name_3_cycles_last_pop = prefix_3_cycles + "last_pop_s1change_0.csv";
        std::string expexcted_name_6_cycles_last_pop = prefix_6_cycles + "last_pop_s1change_0.csv";
        std::string expexcted_name_3_cycles_b_l_p = prefix_3_cycles + "before_last_pop_s1change_0.csv";
        std::string expexcted_name_6_cycles_b_l_p = prefix_6_cycles + "before_last_pop_s1change_0.csv";
        exec(s);
        save_before_last_pop(s,prefix_3_cycles);
        save_last_pop(s,prefix_3_cycles);
        assert(exists(expexcted_name_3_cycles_b_l_p));
        assert(exists(expexcted_name_3_cycles_last_pop));

        exec(s);
        save_before_last_pop(s,prefix_6_cycles);
        save_last_pop(s,prefix_6_cycles);
        assert(exists(expexcted_name_6_cycles_last_pop));
        assert(exists(expexcted_name_6_cycles_b_l_p));
    }
    ///The two last population of funders can be saved
    {
        env_param e;
        ind_param i;
        meta_param m{3,
                     20};
        pop_param p{2};
        sim_param sp{e,i,m,p};
        simulation s{sp};

        exec(s);
        auto last_2_pops = save_two_last_pops(s);

        auto last_pop = load_funders(create_last_pop_name(s));
        auto before_last_pop = load_funders(create_before_last_pop_name(s));

        std::vector<funders> load_last_2_pops{before_last_pop, last_pop};

        assert(last_2_pops == load_last_2_pops);

    }

#endif
}


