#include "simulation.h"
#include <cassert>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>

simulation::simulation(sim_param param):
    m_pop(param.get_pop_param()),
    m_e{param.get_env_param()},
    m_meta_param{param.get_meta_param()}
{
    m_rng.seed(m_meta_param.get_seed());
    m_pop.get_rng().seed(m_meta_param.get_seed());
    place_start_cells(m_pop);
}

void add_new_funders(simulation& s) noexcept
{
    s.get_funders_success().get_v_funders().push_back(prepare_funders(s));
}

void add_success_funders(simulation& s) noexcept
{
    s.get_funders_success().get_v_funders().back() = calc_funders_success(s);
}

funders calc_funders_success(const simulation& s)
{
    auto funders = s.get_funders_success().get_v_funders().back();

    double n_tot_spores =  std::count_if(
                s.get_pop().get_v_ind().begin(),
                s.get_pop().get_v_ind().end(),
                [](const individual i)
    {return  i.get_phen() == phenotype::spore;});

    double tot_fitness = n_tot_spores * s.get_pop().get_param().get_spo_adv()
            + s.get_pop().get_pop_size() - n_tot_spores;

    for(auto& funder : funders.get_v_funder_data())
    {
        assert(funder.get_success() == 0);

        double n_non_spore_descendants =
                std::count_if(
                    s.get_pop().get_v_ind().begin(),
                    s.get_pop().get_v_ind().end(),
                    [&funder](const individual i)
        {return (i.get_ancestor() == funder.get_ancestor_ID())
                    && i.get_phen() != phenotype::spore;});

        double n_spore_descendants =
                std::count_if(
                    s.get_pop().get_v_ind().begin(),
                    s.get_pop().get_v_ind().end(),
                    [&funder](const individual i)
        {return (i.get_ancestor() == funder.get_ancestor_ID())
                    && i.get_phen() == phenotype::spore;});

        auto fitness = n_non_spore_descendants +
                n_spore_descendants * s.get_pop().get_param().get_spo_adv();

        auto success = fitness / tot_fitness;

        funder.set_success(success);
    }
    return funders;
}

void change_conditions(simulation& s) noexcept
{
    change_env(s);
    change_pop(s);
}

void change_env(simulation& s) noexcept
{
    auto new_env_param = change_env_param_norm(s.get_env().get_param(),s.get_rng());
    s.get_env().set_param(new_env_param);
}

void change_pop( simulation& s)
{
    auto& p = s.get_pop();
    const auto new_ind_param = change_ind_param_norm(p.get_param().get_ind_param(), p.get_rng());

    ///Change ind_param object contained in pop_param
    p.get_param().get_ind_param() = new_ind_param;

    ///Change ind_params of all inds in pop
    p.get_v_ind() = change_inds(p,new_ind_param);
}

std::string create_best_random_condition_name(const simulation& s, double amplitude)
{
    return  std::string{
        "best_ind_random_cond_sim_demographic_s" +
        std::to_string(s.get_meta_param().get_seed()) +
                "_change_" +
                std::to_string(s.get_meta_param().get_change_freq())  +
                "_amplitude_"+
                std::to_string(amplitude)+
                ".csv"
    };
}

std::string create_best_random_condition_name(double amplitude, int change_freq, int seed)
{
    return  std::string{
        "best_ind_random_cond_sim_demographic_s" +
        std::to_string(seed) +
                "_change_" +
                std::to_string(change_freq)  +
                "_amplitude_"+
                std::to_string(amplitude)+
                ".csv"
    };
}

std::string create_funders_success_name(const simulation& s)
{
    return  std::string{
        "funders_success_s" +
        std::to_string(s.get_meta_param().get_seed()) +
                "change_" +
                std::to_string(s.get_meta_param().get_change_freq()) +
                ".csv"
    };
}

std::string create_last_pop_name(const simulation& s)
{
    return  std::string{
        "last_pop_s" +
        std::to_string(s.get_meta_param().get_seed()) +
                "change_" +
                std::to_string(s.get_meta_param().get_change_freq()) +
                ".csv"
    };
}

std::string create_last_pop_name(int seed, int change_freq)
{
    return  std::string{
        "last_pop_s" +
        std::to_string(seed) +
                "change_" +
                std::to_string(change_freq) +
                ".csv"
    };
}

std::string create_random_condition_name(const simulation& s, double amplitude)
{
    return  std::string{
        "random_cond_sim_demographic_s" +
        std::to_string(s.get_meta_param().get_seed()) +
                "_change_" +
                std::to_string(s.get_meta_param().get_change_freq())  +
                "_amplitude_"+
                std::to_string(amplitude)+
                ".csv"
    };
}

std::string create_random_condition_name(double amplitude, int change_freq, int seed)
{
    return  std::string{
        "random_cond_sim_demographic_s" +
        std::to_string(seed) +
                "_change_" +
                std::to_string(change_freq)  +
                "_amplitude_"+
                std::to_string(amplitude)+
                ".csv"
    };
}


std::string create_sim_demo_name(const simulation& s)
{
    return  std::string{
        "sim_demographic_s" +
        std::to_string(s.get_meta_param().get_seed()) +
                "change_" +
                std::to_string(s.get_meta_param().get_change_freq()) +
                ".csv"
    };
}

std::string create_sim_demo_name(int seed, int change_freq)
{
    return  std::string{
        "sim_demographic_s" +
        std::to_string(seed) +
                "change_" +
                std::to_string(change_freq) +
                ".csv"
    };
}

std::string create_sim_par_name(const simulation& s)
{
    return  std::string{
        "sim_par_s" +
        std::to_string(s.get_meta_param().get_seed()) +
                "change_" +
                std::to_string(s.get_meta_param().get_change_freq()) +
                ".csv"
    };
}

std::string create_sim_par_name(int seed, int change_freq)
{
    return  std::string{
        "sim_par_s" +
        std::to_string(seed) +
                "change_" +
                std::to_string(change_freq) +
                ".csv"
    };
}

std::vector<std::pair<env_param, ind_param>> create_rand_conditions_unif(const env_param& e,
                                                                         const ind_param& i,
                                                                         int n_rand_conditions,
                                                                         double amplitude,
                                                                         int seed)
{
    std::minstd_rand rng;
    rng.seed(seed);

    std::vector<std::pair<env_param, ind_param>> random_conditions;

    auto env = change_range_env_param(e, amplitude);
    auto ind = change_range_ind_param(i, amplitude);

    for(int r = 0; r != n_rand_conditions; r++)
    {
        random_conditions.push_back({change_env_param_unif(env, rng),
                                     change_ind_param_unif(ind, rng)});
    }

    return random_conditions;
}

void change_params(simulation& s, const env_param& e, const ind_param& i)
{
    s.get_env().set_param(e);

    ///Change ind_param object contained in pop_param
    s.get_pop().get_param().get_ind_param() = i;

    ///Change ind_params of all inds in pop
    s.get_pop().get_v_ind() = change_inds(s.get_pop(),i);

}

void dispersal(simulation &s)
{
    fund_new_pop(s.get_pop());
    reset_env(s.get_env());
}


void exec_cycle(simulation& s) noexcept
{

    add_new_funders(s);
    while(s.get_timestep() != s.get_meta_param().get_cycle_duration())
    {tick(s);}
    add_success_funders(s);
    store_demographics(s);
    dispersal(s);

}

void exec(simulation& s) noexcept
{
    while(s.get_cycle() != s.get_meta_param().get_n_cycles())
    {
        exec_cycle(s);
        if(s.get_cycle() != 0
                && s.get_meta_param().get_change_freq() != 0
                && s.get_cycle() % s.get_meta_param().get_change_freq() == 0
                )
        {
            change_env(s);
            change_pop(s);
        }
        s.reset_timesteps();
        s.tick_cycles();
    }

    save_data(s);
}


void feeding(simulation& s)
{
    for(auto& ind : s.get_pop().get_v_ind())
    {
        auto index_grid = find_grid_index(ind,s.get_env().get_param().get_grid_side());
        if(index_grid == -100 ||
                ind.get_phen() != phenotype::active)
        {continue;}
        feed(ind,s.get_env().get_cell(index_grid));
    }
}


void jordi_feeding(simulation& s)
{
    for(auto& ind : s.get_pop().get_v_ind())
    {
        auto index_grid = find_grid_index(ind,s.get_env().get_param().get_grid_side());
        if(index_grid == -100 ||
                ind.get_phen() != phenotype::active)
        {continue;}
        jordi_feed(ind,s.get_env().get_cell(index_grid));
    }
}

simulation load_sim_from_record(int seed, int change_freq)
{
    auto sim_par_name = create_sim_par_name(seed,change_freq);
    assert(exists(sim_par_name));

    simulation s{load_sim_parameters(sim_par_name)};

    auto funders_name = create_funders_success_name(seed,change_freq);
    assert(exists(funders_name));
    auto funders_success = load_funders_success(funders_name);

    auto demo_name = create_sim_demo_name(seed,change_freq);
    assert(exists(demo_name));
    auto demo_sim = load_demographic_sim(demo_name);

    s.set_funders_success(funders_success) ;
    s.set_demo_sim(demo_sim);

    return s;
}

std::vector<std::pair<env_param, ind_param>> load_random_conditions(const std::string& filename)
{
    std::vector<std::pair<env_param, ind_param>> random_conditions;
    env_param env;
    ind_param ind;
    std::string dummy;
    std::ifstream is{filename};
    while(is >> env
          && is >> dummy
          && is >> ind)
    {
        random_conditions.push_back({env,
                                     ind});
    }

    return random_conditions;
}

simulation load_sim_from_last_pop(int seed, int change_freq)
{
    auto filename = create_sim_par_name(seed,change_freq);

    assert(exists(filename));

    simulation s{load_sim_parameters(filename)};

    auto last_pop = load_funders(create_last_pop_name(seed,change_freq));
    s.get_pop().get_v_ind().resize(last_pop.get_v_funder_data().size());

    for(size_t i = 0; i != last_pop.get_v_funder_data().size(); i++)
    {
        s.get_pop().get_v_ind()[i].get_grn()
                = last_pop.get_v_funder_data()[i].get_grn();
    }

    return s;
}

simulation load_best_ind_for_rand_cond(int seed, int change_freq)
{
    simulation s{load_sim_parameters(create_sim_par_name(seed,change_freq))};

    auto best_ind_grn = find_best_ind_grn(load_funders_success(create_funders_success_name(seed, change_freq)));

    s.get_pop().get_v_ind().resize(s.get_pop().get_param().get_exp_new_pop_size());

    for( auto& individual : s.get_pop().get_v_ind())
    {
        individual.get_grn() = best_ind_grn;
    }

    return s;
}

simulation no_demographic_copy(const simulation& s)
{
    simulation new_s{sim_param{s.get_env().get_param(),
                    s.get_meta_param(),
                    s.get_pop().get_param()
                              }
                    };
    new_s.get_pop().get_v_ind() = s.get_pop().get_v_ind();
    return new_s;
}

funders prepare_funders(const simulation& s)
{
    assert(s.get_timestep() == 0);
    assert(s.get_pop().get_v_ind().size() <= 100);
    funders f;
    for(const auto& ind : s.get_pop().get_v_ind())
    {
        f.get_v_funder_data().push_back(funder_data{ind});
    }
    f.set_cycle(s.get_cycle());
    return f;
}

void reset_sim(simulation& s) noexcept
{
    reset_env(s.get_env());
    reset_pop(s.get_pop());
}

void response(simulation& s)
{
    for(auto& ind : s.get_pop().get_v_ind())
    {
        if(ind.get_phen() == phenotype::spore){continue;}
        auto index = find_grid_index(ind,s.get_env().get_param().get_grid_side());
        if(index == -100)//When individual is outside grid
        {
            responds(ind, env_grid_cell(0,0,0,0));
            continue;
        }
        responds(ind, s.get_env().get_cell(index));
    }
}

demographic_sim run_random_conditions(const simulation& s,
                                      int n_number_rand_cond,
                                      double amplitude,
                                      std::string name)
{
    auto random_conditions = create_rand_conditions_unif(
                s.get_env().get_param(),
                s.get_pop().get_param().get_ind_param(),
                n_number_rand_cond,
                amplitude,
                0);

    auto test_pop = s.get_pop().get_v_ind();

    simulation rand_s = no_demographic_copy(s);

    int counter = 0;

    for(const auto & condition : random_conditions)
    {
        assert(rand_s.get_pop().get_v_ind() == test_pop);
        assert(rand_s.get_env().get_grid() == s.get_env().get_grid());
        change_params(rand_s, condition.first, condition.second);
        auto start = std::chrono::high_resolution_clock::now();
        exec_cycle(rand_s);
        rand_s.tick_cycles();
        rand_s.reset_timesteps();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        std::cout<< "condition: " << duration.count() << "\n";
        rand_s.get_pop().get_v_ind() = test_pop;
        counter++;
        save_demographic_sim(rand_s.get_demo_sim(), name);
    }
    return rand_s.get_demo_sim();
}

void save_data(const simulation& s)
{
    std::string sim_param_name = create_sim_par_name(s);

    sim_param s_p{s.get_env().get_param(),
                s.get_meta_param(),
                s.get_pop().get_param()};

    save_sim_parameters(s_p, sim_param_name);

    std::string last_pop_name = create_last_pop_name(s);

    save_funders(prepare_funders(s), last_pop_name);

    std::string sim_demo_name = create_sim_demo_name(s);

    save_demographic_sim(s.get_demo_sim(), sim_demo_name);

    std::string funders_success_name = create_funders_success_name(s);

    save_funders_success(s.get_funders_success(), funders_success_name);
}

void secretion_metabolite(simulation& s)
{
    int index;
    for(const auto& ind : s.get_pop().get_v_ind())
    {
        index = find_grid_index(ind,s.get_env().get_param().get_grid_side());
        if(index == - 100)
        {
            continue;
        }
        secretes_metab(ind,s.get_env().get_cell(index));
    }
}

void store_demographics( simulation& s) noexcept
{
    s.set_demo_sim(update_demographics(s));
}

int tick(simulation& s)
{
    int time = 0;
    response(s);
    //feeding(s);
    jordi_feeding(s);
    metabolism_pop(s.get_pop());
    secretion_metabolite(s);
    //death(s.get_pop());
    jordi_death(s.get_pop());
    if(division(s.get_pop()))
    {
        time += manage_static_collisions(s.get_pop());
    }
    degradation_metabolite(s.get_env());
    diffusion(s.get_env());
    s.tick_timesteps();
    return time;
}

demographic_sim update_demographics(const simulation& s) noexcept
{
    auto d_s = s.get_demo_sim();
    d_s.get_demo_cycles().push_back(demographics(s.get_pop(), s.get_env().get_param()));
    return d_s;
}

void test_simulation()//!OCLINT tests may be many
{
#ifndef NDEBUG

//    //A simulation can be initialized by getting a sim_parmater class as an argument
//    {
//        unsigned int pop_size = 1;
//        unsigned int exp_new_pop_size = 1;
//        double min_dist = 0.1;
//        int grid_side = 1;
//        double diff_coeff = 0.1;
//        double init_food = 1.0;
//        double mutation_prob = 0.01;
//        double mutation_step = 0.1;
//        double base_disp_prob = 0.01;
//        double spore_advantage = 10.0;
//        double repr_trsh = 0.1;
//        double metab_degr_rate = 0.1;
//        int n_cycles = 42;
//        int cycle_duration = 5;

//        sim_param s_p(pop_size,
//                      exp_new_pop_size,
//                      min_dist,
//                      grid_side,
//                      diff_coeff,
//                      init_food,
//                      mutation_prob,
//                      mutation_step,
//                      base_disp_prob,
//                      spore_advantage,
//                      repr_trsh,
//                      metab_degr_rate,
//                      n_cycles,
//                      cycle_duration
//                      );
//        simulation s(s_p);
//        //tests for pop_param
//        assert(s.get_pop().get_v_ind().size() == pop_size);
//        assert(s.get_pop().get_param().get_exp_new_pop_size() == exp_new_pop_size);
//        assert(s.get_pop().get_param().get_min_dist() - min_dist < 0.0001 &&
//               s.get_pop().get_param().get_min_dist() - min_dist > -0.0001);
//        assert(s.get_pop().get_param().get_mu_p() - mutation_prob < 0.0001 &&
//               s.get_pop().get_param().get_mu_p() - mutation_prob > -0.0001);
//        assert(s.get_pop().get_param().get_mu_st() - mutation_step < 0.0001 &&
//               s.get_pop().get_param().get_mu_st() - mutation_step > -0.0001);
//        assert(s.get_pop().get_param().get_base_disp_prob() - base_disp_prob < 0.00001 &&
//               s.get_pop().get_param().get_base_disp_prob() - base_disp_prob > -0.00001);
//        assert(s.get_pop().get_param().get_spo_adv() - spore_advantage < 0.0001 &&
//               s.get_pop().get_param().get_spo_adv() - spore_advantage > -0.0001);

//        //tests for env param
//        assert(s.get_env().get_param().get_diff_coeff() - diff_coeff < 0.00001 &&
//               s.get_env().get_param().get_diff_coeff() - diff_coeff > -0.00001);
//        assert(s.get_env().get_param().get_init_food() - init_food < 0.0001 &&
//               s.get_env().get_param().get_init_food() - init_food > -0.0001);

//        assert(s.get_env().get_param().get_degr_rate() - metab_degr_rate < 0.0001 &&
//               s.get_env().get_param().get_degr_rate() - metab_degr_rate > -0.0001);
//        //test for simulation meta param
//        assert(s.get_meta_param().get_n_cycles() == n_cycles);
//        assert(s.get_meta_param().get_cycle_duration() == cycle_duration);
//    }


//    //A simulation has an environment
//    // The value -1234567890 is irrelevant: just get this to compile

//    {
//        simulation s;
//        assert(s.get_env().get_param().get_grid_side() > -1234567890);
//    }

//    //A simulation can be initialized with a certain amount
//    //of food in all the cells of its environment grid
//    //1 by default
//    {
//        simulation s;
//        for( auto& grid_cell : s.get_env().get_grid())
//        {
//            assert(grid_cell.get_food() - 1 < 0.000001
//                   && grid_cell.get_food() - 1 > -0.0000001);
//        }

//        double starting_food = 3.14;
//        simulation s1 (sim_param{1, 1, 1, 1, 0.1, starting_food});
//        for( auto& grid_cell : s1.get_env().get_grid())
//        {
//            assert(grid_cell.get_food() - starting_food < 0.000001
//                   && grid_cell.get_food() - starting_food > -0.0000001);
//        }
//    }

//    //Individuals can take up energy from the environment
//    {
//        simulation s(sim_param{1,1,0.1,2});
//        double total_food_init = std::accumulate
//                (
//                    s.get_env().get_grid().begin(),
//                    s.get_env().get_grid().end(),0.0,
//                    [](double sum,const env_grid_cell& c){return sum + c.get_food();}
//        );

//        double total_en_init = std::accumulate
//                (
//                    s.get_pop().get_v_ind().begin(),
//                    s.get_pop().get_v_ind().end(),0.0,
//                    [](double sum,const individual& i){return sum + i.get_energy();}
//        );

//        feeding(s);

//        double total_food_after = std::accumulate
//                (
//                    s.get_env().get_grid().begin(),
//                    s.get_env().get_grid().end(),0.0,
//                    [](double sum, const env_grid_cell& c){return sum + c.get_food();}
//        );

//        double total_en_after =
//                std::accumulate(
//                    s.get_pop().get_v_ind().begin(),
//                    s.get_pop().get_v_ind().end(),0.0,
//                    [](double sum,const individual& i){
//            return sum + i.get_energy();}
//        );

//        assert(total_food_init > total_food_after);
//        assert(total_en_init < total_en_after);

//        auto total_uptake =
//                std::accumulate(s.get_pop().get_v_ind().begin(),
//                                s.get_pop().get_v_ind().end(), 0.0,
//                                [](double sum, const individual& i)
//        {return sum + i.get_param().get_uptake_rate();});

//        assert(total_en_after - (total_uptake + total_en_init) < 0.00000001 &&
//               total_en_after - (total_uptake + total_en_init) > -0.00000001);
//        assert(total_food_init - (total_food_after + total_uptake) < 0.000001 &&
//               total_food_init - (total_food_after + total_uptake) > -0.000001);
//    }

//    //Individuals outside the grid do not feed
//    {
//        simulation s(sim_param{1,1,0.1,2});

//        double total_food_init = std::accumulate
//                (
//                    s.get_env().get_grid().begin(),
//                    s.get_env().get_grid().end(),0.0,
//                    [](double sum,const env_grid_cell& c){return sum + c.get_food();}
//        );

//        double total_en_init = std::accumulate
//                (
//                    s.get_pop().get_v_ind().begin(),
//                    s.get_pop().get_v_ind().end(),0.0,
//                    [](double sum,const individual& i){return sum + i.get_energy();}
//        );

//        set_pos(s.get_pop().get_ind(0),std::pair<double,double>(-42,42));

//        feeding(s);

//        double total_food_after = std::accumulate
//                (
//                    s.get_env().get_grid().begin(),
//                    s.get_env().get_grid().end(),0.0,
//                    [](double sum, const env_grid_cell& c){return sum + c.get_food();}
//        );

//        double total_en_after = std::accumulate
//                (
//                    s.get_pop().get_v_ind().begin(),
//                    s.get_pop().get_v_ind().end(),0.0,
//                    [](double sum,const individual& i){
//            return sum + i.get_energy();}
//        );

//        assert(total_food_init - total_food_after < 0.00001
//               && total_food_init - total_food_after > -0.00001);
//        assert(total_en_init - total_en_after < 0.000001
//               && total_en_init - total_en_after > -0.000001);
//    }

//    //In one tick/timestep individuals take in input,
//    //determine phenotype(based on previous timestep),
//    // feed, than reproduce, than substances diffuse
//    {
//        simulation s(sim_param{1,1,0.1,3,1,1,0});
//        //Set all the hid nodes and H2O and H2H weights to one so
//        //that we are sure the phenotype will stay = active;
//        for(auto& ind : s.get_pop().get_v_ind())
//        {
//            ind.get_grn().set_all_hid_nodes(1);
//            ind.get_grn().set_all_out_nodes(1);
//            ind.get_grn().set_all_H2O(1);
//            ind.get_grn().set_all_H2H(1);
//        }

//        //The single individual in this population
//        //after a tick should reproduce
//        auto init_pop_size = s.get_pop().get_v_ind().size();
//        auto ind_en = get_ind_tr_en( s.get_pop(), 0)
//                + s.get_pop().get_ind(0).get_param().get_metabolic_rate() + 0.01
//                - s.get_pop().get_ind(0).get_param().get_uptake_rate();
//        s.get_pop().get_ind(0).set_energy(ind_en);

//        //and the grid_cell where it is should recieve
//        //food nutrients
//        response(s);
//        feeding(s);

//        auto grid_index_ind =  find_grid_index(s.get_pop().get_ind(0),
//                                               s.get_env().get_param().get_grid_side());

//        double food_after_feeding = s.get_env().get_cell(grid_index_ind).get_food();

//        metabolism_pop(s.get_pop());
//        if(division(s.get_pop()))
//        {
//            manage_static_collisions(s.get_pop());
//        }

//        diffusion(s.get_env());
//        double food_after_diffusion = s.get_env().get_cell(grid_index_ind).get_food();

//        //The difference in food before and after the diffusion step is not 0
//        assert(!(food_after_feeding - food_after_diffusion < 0.00001
//                 && food_after_feeding - food_after_diffusion > -0.0001));
//        assert(s.get_pop().get_v_ind().size() == 2 * init_pop_size);
//        assert(!has_collision(s.get_pop() ));
//    }

//    //If nothing else happens, food should constantly decrease when cells are feeding
//    {
//        meta_param m{};
//        pop_param p{2,
//                    1,
//                    0.1,
//                    0.01,
//                    0.1,
//                    0.01,
//                    10,
//                    1
//                   };
//        env_param e{3,
//                    0.1,
//                    p.get_ind_param().get_treshold_energy() * 10
//                   };
//        simulation s (sim_param{e,m,p});

//        auto food_begin =
//                std::accumulate(s.get_env().get_grid().begin(),
//                                s.get_env().get_grid().end(),
//                                0.0,
//                                [](double sum, const env_grid_cell& c)
//        {return sum + c.get_food();});

//        //The simulation will last long enough  for the individuals to reproduce
//        auto sim_time = s.get_pop().get_ind(0).get_param().get_treshold_energy() /
//                s.get_pop().get_ind(0).get_param().get_uptake_rate() * 2;

//        auto init_pop_size =  s.get_pop().get_v_ind().size();

//        for( int i = 0; i != static_cast<int>(sim_time); i++)
//        {

//            auto food_before_feed =
//                    std::accumulate(s.get_env().get_grid().begin(),
//                                    s.get_env().get_grid().end(), 0.0,
//                                    [](double sum, const env_grid_cell& c)
//            {return sum + c.get_food();});

//            feeding(s);

//            auto food_after_feed =
//                    std::accumulate(s.get_env().get_grid().begin(),
//                                    s.get_env().get_grid().end(), 0.0,
//                                    [](double sum, const env_grid_cell& c)
//            {return sum + c.get_food();});

//            double food_eaten = 0.0;
//            for(const auto& ind : s.get_pop().get_v_ind())
//            {
//                auto grid_cell_ind = find_grid_index(ind, s.get_env().get_param().get_grid_side());
//                if( grid_cell_ind != -100 && s.get_env().get_cell(grid_cell_ind).get_food() > 0)
//                {
//                    food_eaten += ind.get_param().get_uptake_rate();
//                }
//            }

//            auto balance_uptake = food_before_feed - (food_after_feed + food_eaten);
//            assert(balance_uptake < 0.0001 && balance_uptake > -0.0001);

//            metabolism_pop(s.get_pop());
//            division(s.get_pop());
//            manage_static_collisions(s.get_pop());

//            auto food_before_diff =
//                    std::accumulate(s.get_env().get_grid().begin(),
//                                    s.get_env().get_grid().end(), 0.0,
//                                    [](double sum, const env_grid_cell& c)
//            {return sum + c.get_food();});

//            calc_diffusion_food(s.get_env());

//            auto new_grid = s.get_env().get_grid();

//            auto food_after_diff =
//                    std::accumulate(new_grid.begin(),
//                                    new_grid.end(), 0.0,
//                                    [](double sum, const env_grid_cell& c)
//            {return sum + c.get_food();});

//            auto balance_diffusion = food_before_diff - food_after_diff;
//            assert(balance_diffusion < 0.0001 && balance_diffusion > -0.0001);

//            new_grid.swap(s.get_env().get_grid());

//        }
//        auto food_end =
//                std::accumulate(s.get_env().get_grid().begin(),
//                                s.get_env().get_grid().end(), 0.0,
//                                [](double sum, const env_grid_cell& c)
//        {return sum + c.get_food();});

//        assert(food_end < food_begin);

//        auto final_pop_size =  s.get_pop().get_v_ind().size();

//        assert(init_pop_size < final_pop_size);
//    }

//    //A simulation is initiallized with a degradation rate
//    {
//        double degradation_rate = 0.314;
//        simulation s(sim_param{0,0,0,0,1,0,0,0,0,0,0,degradation_rate});
//        assert(s.get_env().get_param().get_degr_rate() - degradation_rate < 0.000001 &&
//               s.get_env().get_param().get_degr_rate() - degradation_rate > -0.000001);
//    }
//    //Every time step the individuals produce new metabolite and metabolite degrades in grid_cells
//    {
//        double degradation_rate = 0.314;
//        double init_metab = degradation_rate;
//        simulation s(sim_param{1,0,0,1,1,0,0,0,0,0,0,degradation_rate});

//        for(auto& grid_cell : s.get_env().get_grid())
//        {
//            grid_cell.set_metab(init_metab);
//        }

//        double tot_metab_before =
//                std::accumulate(s.get_env().get_grid().begin(),
//                                s.get_env().get_grid().end(), 0.0,
//                                [](int sum, const env_grid_cell& g)
//        {return sum + g.get_metab();});

//        secretion_metabolite(s);
//        degradation_metabolite(s.get_env());

//        double tot_metab_after =
//                std::accumulate(s.get_env().get_grid().begin(),
//                                s.get_env().get_grid().end(), 0.0,
//                                [](int sum, const env_grid_cell& g)
//        {return sum + g.get_metab();});

//        double tot_production =
//                std::accumulate(s.get_pop().get_v_ind().begin(),
//                                s.get_pop().get_v_ind().end(), 0.0,
//                                [](int sum, const individual& i)
//        {return sum + i.get_param().get_metab_secr_rate();});

//        double tot_degradation =
//                s.get_env().get_grid_size() * s.get_env().get_param().get_degr_rate();

//        auto metab_balance =
//                tot_metab_before - tot_degradation + tot_production - tot_metab_after;

//        assert(metab_balance < 0.000001 && metab_balance > -0.000001);


//    }

//    //every timestep/tick collisions are handled
//    {
//        simulation s(sim_param{7,3});
//        //The central individual in this population
//        //after a tick should reproduce
//        auto init_pop_size = s.get_pop().get_v_ind().size();
//        s.get_pop().get_ind(1).set_energy(get_ind_tr_en(s.get_pop(), 1)
//                                          + s.get_pop().get_ind(1).get_param().get_metabolic_rate()
//                                          + 0.01
//                                          - s.get_pop().get_ind(1).get_param().get_uptake_rate());
//        feeding(s);
//        metabolism_pop(s.get_pop());
//        division(s.get_pop());
//        manage_static_collisions(s.get_pop());

//        assert(s.get_pop().get_v_ind().size() == init_pop_size + 1);
//        assert(!has_collision(s.get_pop()));
//    }

//    //A simulation is initialized with a m_tick = 0;
//    {
//        simulation s;
//        assert(s.get_timestep() == 0);
//    }

//    //After each tick the simulation updates its m_tick
//    //by one
//    {
//        simulation s;
//        for(int i = 0; i != 3; i++)
//        {
//            assert(s.get_timestep() == i);
//            tick(s);
//        }
//    }

//    //A simulation runs a cycle for a certain amount of ticks stated in it s parameters
//    {
//        auto n_cycles = 1;
//        auto cycle_duration = 1;
//        meta_param m{ n_cycles, cycle_duration};
//        sim_param sp{ env_param(), m, pop_param()};
//        simulation s{sp};
//        exec_cycle(s);
//        assert(s.get_timestep() == cycle_duration);
//    }

//    //A simulation can be run for the amount of cycles stated in its parameters
//    {
//        auto n_cycles = 1;
//        auto cycle_duration = 1;
//        pop_param p{100};
//        meta_param m{ n_cycles, cycle_duration};
//        sim_param sp{ env_param(), m, p};
//        simulation s{sp};
//        exec(s);
//        assert(s.get_cycle() == n_cycles);
//    }

//    //After the exec_cycle the time_step timer is reset
//    {
//        auto n_cycles = 1;
//        auto cycle_duration = 1;
//        pop_param p{100};
//        meta_param m{ n_cycles, cycle_duration};
//        sim_param sp{ env_param(), m, p};
//        simulation s{sp};
//        //As the parameters indicate only one cycle will be executed
//        exec(s);
//        assert(s.get_cycle() == n_cycles);
//        assert(s.get_timestep() == 0);
//    }

//    //After each cycle a new population(max 100 individuals
//    //is created from the previous one(see population tests)
//    //and the environment is reset
//    //(similar to dispersal() test)
//    {
//        auto n_cycles = 1;
//        auto cycle_duration = 1;
//        meta_param m{ n_cycles, cycle_duration};
//        pop_param p{100};
//        env_param e{3};
//        sim_param sp{e, m, p};
//        simulation s{sp};
//        auto init_env = s.get_env();
//        //Run a few ticks to make sure env is different from original
//        for(int i = 0; i != 1; i++)
//        {
//            tick(s);
//        }
//        assert(s.get_env() != init_env);
//        //Reset tick timer to avoid crushing on assert  from prepare_funders()
//        s.reset_timesteps();
//        //As the parameters indicate only one cycle will be executed
//        exec(s);
//        assert(s.get_pop().get_pop_size() - s.get_pop().get_param().get_exp_new_pop_size() == 0);
//        assert(s.get_env() == init_env);
//    }

//    //Spores do not feed
//    {
//        simulation s;
//        auto init_food = s.get_env().get_cell(0).get_food();
//        s.get_pop().get_ind(0).set_phen(phenotype::spore);
//        feeding(s);
//        assert(init_food - s.get_env().get_cell(0).get_food() < 0.000001
//               && init_food - s.get_env().get_cell(0).get_food() > -0.000001);
//    }


//    //Sporulating individuals do not feed but they lose energy
//    {
//        simulation s;
//        s.get_pop().get_ind(0).set_phen(phenotype::sporulating);
//        set_ind_en(s.get_pop().get_ind(0),1);
//        auto init_en_ind0 = get_ind_en(s.get_pop(), 0);
//        auto init_food = s.get_env().get_cell(0).get_food();
//        feeding(s);
//        metabolism_pop(s.get_pop());
//        assert(init_en_ind0 - get_ind_en(s.get_pop(), 0) -
//               s.get_pop().get_ind(0).get_param().get_spor_metabolic_rate() < 0.000001
//               &&
//               init_en_ind0 - get_ind_en(s.get_pop(), 0) -
//               s.get_pop().get_ind(0).get_param().get_spor_metabolic_rate() > -0.000001);

//        assert(init_food - s.get_env().get_cell(0).get_food() < 0.000001
//               && init_food - s.get_env().get_cell(0).get_food() > -0.000001);
//    }

//    //Max 100 ind, selected based on phenotype, are placed in a hex pattern,
//    //in a new env after dispersal
//    //Tests all of the above
//    {
//        unsigned int pop_size = 1000;
//        unsigned int new_pop_size = 100;
//        auto food = 42.1;
//        auto metabolite = 42.1;
//        simulation s(sim_param{pop_size,new_pop_size});
//        for(auto& grid_cell : s.get_env().get_grid())
//        {
//            grid_cell.set_food(food);
//            grid_cell.set_metab(metabolite);
//        }
//        environment ref_env = s.get_env();
//        dispersal(s);
//        //Max 100 ind
//        assert(s.get_pop().get_v_ind().size() == new_pop_size);
//        //Hex pattern
//        auto n_hex_l = count_hex_layers(s.get_pop().get_pop_size());
//        auto v_modulus = modulus_of_btw_ind_angles(s.get_pop(), M_PI/ (6 * n_hex_l));
//        for(auto ind_modulus : v_modulus)
//        {
//            assert( ind_modulus < 0.0000000001 ||
//                    (ind_modulus > M_PI / (6 * n_hex_l) - 0.1 &&
//                     ind_modulus <= M_PI / (6 * n_hex_l) + 0.1));
//        }
//        assert(!has_collision(s.get_pop()));
//        //Reset env
//        assert( s.get_env() != ref_env );
//        auto init_food = s.get_env().get_param().get_init_food();
//        for(const auto& grid_cell : s.get_env().get_grid())
//        {
//            assert(grid_cell.get_food() - init_food <  0.000001 &&
//                   grid_cell.get_food() - init_food >  -0.000001);
//            assert(grid_cell.get_metab() <  0.000001 &&
//                   grid_cell.get_metab() >  -0.000001);
//        }
//    }

//    //All individuals in a simulation can respond to their environment and surrounding
//    {
//        double food_amount = 3.14;
//        double metabolite_amount = 3.14;
//        double energy_amount = 3.14;
//        simulation s(sim_param{2,1,0.1,4});
//        for(auto & grid_cell : s.get_env().get_grid())
//        {
//            grid_cell.set_food(food_amount);
//            grid_cell.set_metab(metabolite_amount);
//        }
//        for(auto& ind : s.get_pop().get_v_ind())
//        {
//            ind.set_energy(energy_amount);
//            //let's set all the weights of the network to 1
//            //in this case we expect that the outputs will be one
//            ind.get_grn().set_all_I2H(1);
//            ind.get_grn().set_all_H2O(1);
//        }
//        //Since the networks react to input of timestpe t
//        //at timestep t+1 I will run responds(i) 2 times
//        //In this first response call, since the hidden nodes are 0
//        //The output will be 0 and therefore the individuals phenotype
//        //will be phenotype::sporulating
//        response(s);
//        for(auto& ind : s.get_pop().get_v_ind())
//        {
//            assert(is_active(ind));
//        }

//        //The values of food and metabolite in the grid_cell
//        //as well as the energy in the individual
//        //are changed so that all inputs are -1,
//        //the network should therefore give outputs == 0/false
//        //since the outputs node will recieve a signal
//        //that is negative(below the treshold)
//        //after responds(i) is called 2 more times
//        for(auto& grid_cell : s.get_env().get_grid())
//        {
//            grid_cell.set_food(-1);
//            grid_cell.set_metab(-1);
//        }
//        for(auto& ind : s.get_pop().get_v_ind())
//        {
//            ind.set_energy(-1);
//        }
//        //1st respose, the individuals respond to the initial parameter
//        //expected output value == 1
//        response(s);
//        for(auto& ind : s.get_pop().get_v_ind())
//        {
//            assert(is_active(ind));
//        }
//        //2nd response, the individual responds to the changed parameter(all 0s)
//        //expected output value == 0
//        response(s);
//        for(auto& ind : s.get_pop().get_v_ind())
//        {
//            assert(!is_active(ind));
//            assert(is_sporulating(ind));
//        }
//    }
//    //It is possible to store the demographics
//    //of the population contained in a simulation
//    //at a certain point in time in a vector
//    {
//        simulation s;
//        auto demo_sim_length = s.get_demo_sim().get_demo_cycles().size();
//        store_demographics(s);
//        auto demo_sim_length2 = s.get_demo_sim().get_demo_cycles().size();
//        assert(demo_sim_length != demo_sim_length2);
//        assert(demo_sim_length + 1 == demo_sim_length2);
//    }

//    //At the end of each cycle the demographics are stored
//    {
//        int n_cycles = 2;
//        int cycle_duration = 1;
//        meta_param m{n_cycles, cycle_duration};
//        env_param e;
//        pop_param p;
//        sim_param s_p{e,m,p};

//        simulation s{s_p};
//        auto init_demo_cycle = s.get_demo_sim().get_demo_cycles().size();
//        assert(init_demo_cycle == 0);
//        exec_cycle(s);
//        auto one_demo_cycle = s.get_demo_sim().get_demo_cycles().size();
//        assert(one_demo_cycle == 1);
//    }

//    //At the end of the simulation the demographic is saved in a file
//    {
//        int n_cycles = 2;
//        int cycle_duration = 1;
//        meta_param m{n_cycles, cycle_duration};
//        env_param e;
//        pop_param p;
//        sim_param s_p{e,m,p};

//        simulation s{s_p};
//        std::string expected_file_name = create_sim_demo_name(s);
//        exec(s);

//        assert(exists(expected_file_name));
//        auto d_s = load_demographic_sim(expected_file_name);
//        assert(d_s == s.get_demo_sim());

//    }

//    //A simulation is initialized with a funders_success  object
//    {
//        simulation s;
//        assert(s.get_funders_success().get_v_funders().size() >= 0u);
//    }

//    //It is possible to store the ancestor_ID and GRN of the
//    //individuals composing a population at the start of the cycle
//    //updating the funders_success member of simulation
//    {
//        int n_cycles = 1;
//        int cycle_duration = 50;
//        meta_param m{n_cycles, cycle_duration};
//        env_param e;
//        pop_param p;
//        sim_param s_p{e,m,p};

//        simulation s{s_p};
//        auto zero_cycle_funders = prepare_funders(s);
//        exec_cycle(s);
//        s.reset_timesteps();
//        auto first_cycle_funders = prepare_funders(s);
//        assert(zero_cycle_funders != first_cycle_funders);

//    }

//    //Every cycle the funders are stored
//    {
//        int n_cycles = 2;
//        int cycle_duration = 50;
//        meta_param m{n_cycles, cycle_duration};
//        env_param e;
//        pop_param p;
//        sim_param s_p{e,m,p};

//        simulation s{s_p};
//        //the first cycle is funded by a single individual
//        exec_cycle(s);
//        s.reset_timesteps();
//        auto first_cycle_funders_size = s.get_funders_success().get_v_funders().size();
//        exec_cycle(s);
//        s.reset_timesteps();
//        auto second_cycle_funder_size = s.get_funders_success().get_v_funders().size();

//        assert(first_cycle_funders_size != second_cycle_funder_size);
//        assert(s.get_funders_success().get_v_funders().begin() !=
//                s.get_funders_success().get_v_funders().begin() + 1);
//    }

//    //The success of each funder can be calculated
//    //The success of a funder is equal to:
//    // the fractions of individuals with its same ancestor_ID
//    // over the total number of individuals of the population
//    // considering the spore advantage in fitness
//    {
//        int n_cycles = 1;
//        int cycle_duration = 1;
//        meta_param m_p{n_cycles, cycle_duration};

//        unsigned int pop_size = 3;
//        pop_param p_p{pop_size};
//        sim_param s_p{env_param{}, m_p, p_p};
//        simulation s{s_p};

//        add_new_funders(s);

//        auto funders_with_success = calc_funders_success(s);

//        double funder_total_success =
//                std::accumulate(funders_with_success.get_v_funder_data().begin(),
//                                funders_with_success.get_v_funder_data().end(), 0.0,
//                                [](double sum, const funder_data& f)
//        {return sum + f.get_success();});

//        double total_fitness =
//                std::accumulate(s.get_pop().get_v_ind().begin(),
//                                s.get_pop().get_v_ind().end(),
//                                0.0,
//                                [&s](int sum, const individual& i)
//        {return sum + (is_spore(i) ? s.get_pop().get_param().get_spo_adv() : 1);});

//        total_fitness /= s.get_pop().get_pop_size();

//        assert(funder_total_success - total_fitness < 0.0001
//               && funder_total_success - total_fitness > -0.0001);
//    }

//    ///The success of each funder is calculated at the end of each cycle
//    {
//        int n_cycles = 1;
//        int cycle_duration = 1;
//        meta_param m_p{n_cycles, cycle_duration};

//        unsigned int pop_size = 3;
//        pop_param p_p{pop_size};
//        sim_param s_p{env_param{}, m_p, p_p};
//        simulation s{s_p};

//        auto pre_cycle_funders = prepare_funders(s);
//        exec_cycle(s);


//        assert(s.get_funders_success().get_v_funders().back() !=
//                pre_cycle_funders);

//        const auto& funders =
//                s.get_funders_success().get_v_funders().back().get_v_funder_data();

//        for(const auto& funder : funders)
//            assert(funder.get_success() - 1.0/pop_size < 0.00001
//                   && funder.get_success() - 1.0/pop_size > -0.00001);
//    }

//    //At the end of the simulation the funders_success is saved in a file
//    {
//        int n_cycles = 5;
//        int cycle_duration = 1;
//        meta_param m{n_cycles, cycle_duration};
//        env_param e;
//        pop_param p;
//        sim_param s_p{e,m,p};

//        simulation s{s_p};
//        exec(s);

//        std::string expected_file_name = create_funders_success_name(s);
//        assert(exists(expected_file_name));
//        auto f_s = load_funders_success(expected_file_name);
//        assert(f_s == s.get_funders_success());

//    }

//    //A simulation with a certain seed in metaparameters initializes the
//    //Random number generator of population to that seed
//    {
//        int seed = 42;
//        meta_param m{
//            1,
//            1,
//            seed
//        };
//        env_param e;
//        pop_param p;
//        sim_param s_p{e, m, p};
//        simulation s{s_p};
//        auto ref_rng = std::minstd_rand();
//        ref_rng.seed(seed);
//        assert(s.get_pop().get_rng() == ref_rng);
//    }

//    //Every cycle multiple of the change_frequency metaparameter
//    // environmental parameters change
//    {

//        env_param e{};
//        pop_param p;
//        int change_frequency = 1;
//        meta_param m{change_frequency + 1,
//                    1,
//                    1,
//                    change_frequency
//                    };
//        sim_param s_p{e, m, p};
//        simulation s{s_p};

//        exec(s);
//        assert(e != s.get_env().get_param());
//    }

//    //Env param can be changed in a simulation
//    //Based on env params
//    {
//        auto mean_diff_coeff = 0.1;
//        auto mean_degr_coeff = 0.1;
//        auto var_degr_rate = mean_degr_coeff / 3 - 0.01;
//        auto var_diff_coeff = mean_diff_coeff / 3 - 0.01;
//        env_param e_p{1, 0.1, 0.1, 0.1,
//                      mean_diff_coeff,
//                              mean_degr_coeff,
//                              var_diff_coeff,
//                              var_degr_rate};
//        pop_param p;
//        meta_param m;
//        sim_param s_p{e_p, m, p};
//        simulation s{s_p};
//        int repeats = 100;

//        for(int i  = 0 ; i != repeats; i++)
//        {
//            auto prev_env = s.get_env();
//            double prev_deg_rate = prev_env.get_param().get_degr_rate();
//            double prev_diff_coeff =  prev_env.get_param().get_diff_coeff();

//            change_env(s);

//            auto deg_rate = s.get_env().get_param().get_degr_rate();
//            auto diff_coeff = s.get_env().get_param().get_diff_coeff();
//            assert(prev_deg_rate != deg_rate);
//            assert(prev_diff_coeff != diff_coeff);
//        }
//    }

//    //Env_param and ind_param changes with a frequency dictated by the metaparameters
//    //frequency_of_change value
//    {
//        auto frequency_of_change = 2;
//        auto n_cycles = frequency_of_change +1;
//        //environment will only change once

//        meta_param m{n_cycles,
//                    1,
//                    1,
//                    frequency_of_change};

//        auto mean_diff_coeff = 0.1;
//        auto mean_degr_coeff = 0.1;
//        auto var_degr_rate = mean_degr_coeff / 3 - 0.01;
//        auto var_diff_coeff = mean_diff_coeff / 3 - 0.01;
//        env_param e{1, 0.1, 0.1, 0.1,
//                    mean_diff_coeff,
//                            mean_degr_coeff,
//                            var_diff_coeff,
//                            var_degr_rate};

//        pop_param p{};
//        sim_param s_p{e,m,p};
//        simulation s{s_p};
//        auto prev_env = s.get_env();
//        auto prev_pop = s.get_pop().get_v_ind();
//        exec(s);
//        assert(prev_env != s.get_env());
//        assert(prev_pop != s.get_pop().get_v_ind());
//    }

//    //Each cycle env_param and ind_param will be stored in the demographic cycle
//    //That will be use to store the particular parameters that can change throughout the simulation
//    {
//        simulation s;
//        demographic_cycle d_c = demographics(s.get_pop(), s.get_env().get_param());
//        assert(d_c.get_ind_param() == s.get_pop().get_param().get_ind_param());
//        assert(d_c.get_env_param() == s.get_env().get_param());
//        change_env(s);
//        change_pop(s);
//        demographic_cycle d_c1 = demographics(s.get_pop(), s.get_env().get_param());
//        assert(d_c != d_c1);
//    }


//    //It is possible to load the random conditions in a vector
//    {
//        //Make a random conditions .csv file
//        env_param env;
//        ind_param ind;
//        std::minstd_rand rng;
//        std::string filename{"test_random_conditions.csv"};
//        std::ofstream os{filename};

//        std::vector<std::pair<env_param,ind_param>> rand_conditions_saved;

//        int repeats = 2;

//        for(int i = 0; i != repeats; i++)
//        {
//            rand_conditions_saved.push_back({change_env_param_norm(env, rng),
//                                             change_ind_param_norm(ind, rng)});
//            os << rand_conditions_saved.back().first << " , "
//               << rand_conditions_saved.back().second << std::endl;
//        }



//        auto random_conditions_loaded = load_random_conditions(filename);
//        assert( random_conditions_loaded == rand_conditions_saved);

//    }

//    //It is possible to run a population against multiple random conditions
//    //And store their demographics in a file
//    {
//        ///Run against the random conditions
//        pop_param p{100};
//        env_param e{100};
//        int seed_ID = 1234567890;

//        meta_param m{1,
//                     1,
//                     seed_ID};

//        simulation s{sim_param{e,m,p}};

//        double amplitude = 1;
//        int repeats = 2;

//        auto rand_conditions_vector = create_rand_conditions_unif(
//                    s.get_env().get_param(),
//                    s.get_pop().get_param().get_ind_param(),
//                    repeats,
//                    amplitude,
//                    0);

//        std::string filename = create_random_condition_name(s,amplitude);

//        auto demo_sim = run_random_conditions(s, repeats, amplitude, filename);
//        assert(std::equal(rand_conditions_vector.begin(),rand_conditions_vector.end(),
//                          demo_sim.get_demo_cycles().begin(),
//                          [](const std::pair<env_param, ind_param>& r,
//                          const demographic_cycle& d)
//        {return r.first == d.get_env_param() && r.second == d.get_ind_param();})
//               );
//        assert(exists(filename));
//        auto demo_sim_load = load_demographic_sim(filename);
//        assert(demo_sim == demo_sim_load);

//    }
//    //It is possible to change param of a simulation with other
//    //taken from a random condition vector
//    {

//        //Create random conditions
//        env_param env;
//        ind_param ind;
//        std::minstd_rand rng;
//        std::vector<std::pair<env_param,ind_param>> rand_conditions;
//        int repeats = 2;

//        for(int i = 0; i != repeats; i++)
//        {
//            rand_conditions.push_back({change_env_param_norm(env, rng),
//                                       change_ind_param_norm(ind, rng)});
//        }

//        simulation s;
//        for(const auto& rand_condition : rand_conditions)
//        {
//            assert(s.get_env().get_param() != rand_condition.first
//                    && s.get_pop().get_param().get_ind_param() != rand_condition.second);
//            change_params(s,rand_condition.first, rand_condition.second);
//            assert(s.get_env().get_param() == rand_condition.first
//                   && s.get_pop().get_param().get_ind_param() == rand_condition.second);
//        }
//    }

//    //A simulation can be instantiated given another simulation.
//    //With the same population and environment,
//    //but with empty data vectors for demographics
//    {
//        env_param e{5};
//        pop_param p{2};
//        meta_param m{2,1,5,2};
//        simulation s{sim_param{e,m,p}};
//        exec_cycle(s);
//        assert(!s.get_demo_sim().get_demo_cycles().empty());
//        simulation new_s = no_demographic_copy(s);
//        assert(s.get_env() == new_s.get_env());
//        assert(s.get_pop() == new_s.get_pop());
//        assert(s.get_meta_param() == new_s.get_meta_param());
//        assert(s.get_demo_sim() != new_s.get_demo_sim());
//        assert(new_s.get_demo_sim().get_demo_cycles().empty());
//    }

//    /// It is possible to create and save a certain amount of random conditions
//    /// giving the base env_ and ind_ param, and how much larger the range of variation
//    /// is going to be proportionally.
//    /// It is possible also to decide the seed
//    /// of the random number generator for consistency across
//    /// different runs.
//    {
//        int n_random_conditions = 2;
//        int seed = 55;
//        env_param e;
//        ind_param i;
//        auto amplitude = 1.0; //the range will stay the same
//        auto rand_cond = create_rand_conditions_unif(e,
//                                                     i,
//                                                     n_random_conditions,
//                                                     amplitude,
//                                                     seed);

//        auto rand_cond2 = create_rand_conditions_unif(e,
//                                                      i,
//                                                      n_random_conditions,
//                                                      amplitude,
//                                                      seed);

//        assert(static_cast<unsigned int>(n_random_conditions) == rand_cond.size());
//        assert(rand_cond == rand_cond2);
//    }

//    ///The final population is saved at the end of the exec(s) function
//    ///as funders
//    {
//        env_param e{5};
//        meta_param m{1,
//                     1};
//        pop_param p{100,
//                    100};
//        simulation s{sim_param{e,m,p}};
//        exec(s);
//        std::string expected_filename = create_last_pop_name(s);
//        assert(prepare_funders(s) == load_funders(expected_filename));
//    }

//    ///Simulation parameters are saved by the end of exec(s) function
//    {
//        env_param e{5};
//        meta_param m{1,
//                     1};
//        pop_param p{100,
//                    100};
//        sim_param s_p{e,m,p};
//        simulation s{s_p};
//        exec(s);
//        std::string expected_filename = create_sim_par_name(s);
//        assert(s_p == load_sim_parameters(expected_filename));
//    }

    ///A simulation can be loaded given the seed and the change freq
    ///By loading the sim_params of a given simulation
    /// and instantiating the last population of that given simulation
    {
        env_param e{5};
        int seed = 23;
        int change_freq = 21;
        meta_param m{1,
                     1,
                     seed,
                             change_freq};
        pop_param p{100,
                    100};
        sim_param s_p{e,m,p};
        simulation s{s_p};
        exec(s);
        simulation s1 = load_sim_from_last_pop(seed,change_freq);
        assert(s.get_pop().get_v_ind() == s1.get_pop().get_v_ind());
        assert(s.get_pop().get_param() == s1.get_pop().get_param());
        assert(s.get_env().get_param() == s1.get_env().get_param());
        assert(s.get_meta_param() == s1.get_meta_param());
        assert(s.get_env() == s.get_env());
    }

    /// It is possible to load a population of a number of individuals
    /// equal to the expexcted new population size specified in pop_param
    /// made only of the best final network of a population
    {
        int seed = 123;
        int change_freq = 10;
        pop_param p{1,100};
        meta_param m{10,125,seed,change_freq};
        env_param e{};
        simulation s{sim_param{e, m, p}};
        exec(s);
        auto filename = create_funders_success_name(s);
        auto fund_succ = load_funders_success(filename);
        assert(fund_succ == s.get_funders_success());
        auto best_grn = find_best_ind_grn(fund_succ);
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

        simulation s = load_sim_from_record(seed, change_freq);
        auto sp = sim_param{s.get_env().get_param(),
                                 s.get_meta_param(),
                                 s.get_pop().get_param()};

        assert(s.get_funders_success() == fund_succ);
        assert(s.get_demo_sim() == sim_dem);
        assert(sp == sim_par );

    }
#endif
}
















