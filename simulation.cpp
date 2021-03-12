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
    m_pop(param.get_pop_param(), param.get_ind_param()),
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
    s.get_env().set_new_env_param(new_env_param);
}

void change_pop( simulation& s)
{
    auto& p = s.get_pop();
    const auto new_ind_param = change_ind_param_norm(p.get_v_ind().begin()->get_param(), p.get_rng());

    ///Change ind_params of all inds in pop
    p.get_v_ind() = set_new_ind_par(p.get_v_ind(),new_ind_param);
}


int continue_evo(int seed, int change_freq)
{
    auto s = load_sim_last_pop(seed, change_freq);
    auto start = std::chrono::high_resolution_clock::now();

    exec(s);
    save_data(s);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    std::cout << "simualtion :"<< duration.count() << "s" << std::endl;
    return 0;
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

std::string create_funders_success_name(const simulation& s,
                                        std::string prefix,
                                        std::string suffix)
{
    return  std::string{
        prefix +
                "funders_success_s" +
                std::to_string(s.get_meta_param().get_seed()) +
                "change_" +
                std::to_string(s.get_meta_param().get_change_freq()) +
                suffix +
                ".csv"
    };
}

std::string create_last_pop_name(const simulation& s, const std::string &prefix)
{
    return prefix + std::string{
        "last_pop_s" +
        std::to_string(s.get_meta_param().get_seed()) +
                "change_" +
                std::to_string(s.get_meta_param().get_change_freq()) +
                ".csv"
    };
}

std::string create_before_last_pop_name(const simulation& s, const std::string &prefix)
{
    return prefix + std::string{
        "before_last_pop_s" +
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

std::string create_name_vec_rand_cond(int n_of_conditions, double amplitude, int seed)
{
    return  std::string{
        "v_rand_cond_n" +
        std::to_string(n_of_conditions) +
                "a_" +
                std::to_string(amplitude) +
                "s_" +
                std::to_string(seed) +
                ".csv"
    };
}

std::string create_test_random_condition_name(const simulation& s, double amplitude)
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

std::string create_test_random_condition_name(double amplitude, int change_freq, int seed)
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

std::string create_sim_demo_name(const simulation& s, std::string prefix, std::string suffix)
{
    return  std::string{
        prefix +
                "sim_demographic_s" +
                std::to_string(s.get_meta_param().get_seed()) +
                "change_" +
                std::to_string(s.get_meta_param().get_change_freq()) +
                suffix +
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

std::vector<std::pair<env_param, ind_param>> create_rand_conditions_unif_extreme(const env_param& e,
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
         random_conditions.push_back({change_env_param_unif_extreme(env, rng),
                                      change_ind_param_unif_extreme(ind, rng)});
     }

     return random_conditions;
 }

std::vector<std::pair<env_param, ind_param>> create_vector_random_conditions(const env_param& e,
                                                                             const ind_param& i,
                                                                             double amplitude,
                                                                             int n_conditions,
                                                                             int seed)
{
    auto name = create_name_vec_rand_cond(n_conditions,amplitude,seed);
    auto rand_cond_vector = create_rand_conditions_unif(e,i,n_conditions,amplitude,seed);
    save_vector_of_rand_cond(rand_cond_vector, name);
    return rand_cond_vector;
}

std::vector<std::vector<std::pair<env_param, ind_param>>> create_rand_conditions_matrix(const env_param& ep,
                                                                                        const ind_param& ip,
                                                                                        int number_of_sequences,
                                                                                        int conditions_per_sequence,
                                                                                        double amplitude)
{
    std::vector<std::vector<std::pair<env_param, ind_param>>> condition_matrix;

    for(int i = 0; i != number_of_sequences; i++)
    {
        auto condition_sequence = create_rand_conditions_unif(ep, ip, conditions_per_sequence, amplitude, i);
        condition_matrix.push_back(condition_sequence);
    }

    return condition_matrix;
}

std::vector<std::vector<std::pair<env_param, ind_param>>> create_rand_conditions_matrix_extreme(const env_param& ep,
                                                                                         const ind_param& ip,
                                                                                         int number_of_sequences,
                                                                                         int conditions_per_sequence,
                                                                                         double amplitude)
 {
     std::vector<std::vector<std::pair<env_param, ind_param>>> condition_matrix;

     for(int i = 0; i != number_of_sequences; i++)
     {
         auto condition_sequence = create_rand_conditions_unif_extreme(ep, ip, conditions_per_sequence, amplitude, i);
         condition_matrix.push_back(condition_sequence);
     }

     return condition_matrix;
 }

demographic_cycle demographics(const simulation &s, const env_param &e) noexcept
{
    const auto& p = s.get_pop();
    return demographic_cycle{
        count_actives(p),
                count_spores(p),
                count_sporulating(p),
                s.get_timestep(),
                e,
                p.get_v_ind().begin()->get_param()};
}

void dispersal(simulation &s)
{
    fund_new_pop(s.get_pop());
    reset_env(s.get_env());
}


void exec_cycle(simulation& s) noexcept
{
    auto rand_start = std::chrono::high_resolution_clock::now();

    add_new_funders(s);
    while(s.get_timestep() != s.get_meta_param().get_cycle_duration() &&
          s.get_pop().get_pop_size() < s.get_meta_param().get_pop_max())
    {
        tick(s, s.get_meta_param().get_collision_check_interval());
    }
    add_success_funders(s);
    store_demographics(s);
    dispersal(s);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - rand_start);
    std::cout<< "random_cond_evo cycle n " << s.get_cycle() <<
                ":" << duration.count() << "s" << std::endl;
}

void exec(simulation& s) noexcept
{
    while(s.get_cycle() != s.get_meta_param().get_n_cycles())
    {
        auto start = std::chrono::high_resolution_clock::now();

        exec_cycle(s);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        std::cout<< "cycle n " << s.get_cycle() << ":" << std::endl <<
                    "time: " << duration.count() << "s" << std::endl <<
                    "n_individuals: " << s.get_pop().get_pop_size() << std::endl <<
                    std::endl;

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
}

void exec_change(simulation& s,
                 const std::vector<std::pair<env_param, ind_param>>& rand_conditions)
{
    try {
        if(rand_conditions.size() > s.get_meta_param().get_n_cycles())
        {
            throw std::string{"More random conditions than cycles!"};
        }
        else if(s.get_meta_param().get_n_cycles() % rand_conditions.size())
        {
            throw std::string{"conditions are not a dividend of the number of cycles!"};
        }
    }
    catch (std::string& e) {
        std::cout << e << std::endl;
#ifdef NDEBUG
        abort();
#endif
    }

    size_t condition_index = 0;
    auto condition_duration = s.get_meta_param().get_n_cycles() / rand_conditions.size();

    while(s.get_cycle() != s.get_meta_param().get_n_cycles())
    {
        auto start = std::chrono::high_resolution_clock::now();
        auto modulus =  s.get_cycle() % condition_duration;

        if( modulus == 0 && condition_index < rand_conditions.size())
        {
            reproduce_rand_cond(s, rand_conditions, condition_index);
            condition_index++;
        }

        exec_cycle(s);

        s.reset_timesteps();
        s.tick_cycles();

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        std::cout<< "cycle n " << s.get_cycle() << ":" << std::endl <<
                    "time: " << duration.count() << "s" << std::endl <<
                    "n_individuals: " << s.get_pop().get_pop_size() << std::endl << std::endl;
    }
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

double find_max_diff_coeff_rand_cond(const std::vector<std::pair<env_param,ind_param>>& rand_cond)
{
    auto max = std::max_element(rand_cond.begin(), rand_cond.end(),
                     [](const std::pair<env_param,ind_param>& lhs, const std::pair<env_param, ind_param>& rhs)
    {return lhs.first.get_diff_coeff() < rhs.first.get_diff_coeff();});

    return max->first.get_diff_coeff();

}

double find_min_diff_coeff_rand_cond(const std::vector<std::pair<env_param,ind_param>>& rand_cond)
{
    auto min = std::max_element(rand_cond.begin(), rand_cond.end(),
                     [](const std::pair<env_param,ind_param>& lhs, const std::pair<env_param, ind_param>& rhs)
    {return lhs.first.get_diff_coeff() > rhs.first.get_diff_coeff();});

    return min->first.get_diff_coeff();
}

double mean_diff_coeff_rand_cond(const std::vector<std::pair<env_param,ind_param>>& rand_cond)
{
    auto sum = std::accumulate(rand_cond.begin(), rand_cond.end(), 0.0,
                               [] (const double& s, const std::pair<env_param,ind_param>& lhs)
    {return lhs.first.get_diff_coeff() + s;});

    return sum / rand_cond.size();
}

simulation load_sim_no_pop(int seed, int change_freq)
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
    if(!exists(filename))
    {
        std::cerr << "random condition vector filername does not exist";
        abort();
    }

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

simulation load_sim(int seed, int change_freq)
{
    auto sim_par_filename = create_sim_par_name(seed,change_freq);
    assert(exists(sim_par_filename));
    simulation s{load_sim_parameters(sim_par_filename)};

    auto sim_demographic_filename = create_sim_demo_name(seed, change_freq);
    assert(exists(sim_demographic_filename));
    s.set_demo_sim(load_demographic_sim(sim_demographic_filename));

    auto funders_success_filename = create_funders_success_name(seed, change_freq);
    assert(exists(funders_success_filename));
    s.get_funders_success() = load_funders_success(funders_success_filename);

    auto last_pop = s.get_funders_success().get_v_funders().back();
    s.get_pop().get_v_ind().resize(last_pop.get_v_funder_data().size());

    for(size_t i = 0; i != last_pop.get_v_funder_data().size(); i++)
    {
        s.get_pop().get_v_ind()[i].get_grn()
                = last_pop.get_v_funder_data()[i].get_grn();
    }

    update_radius_pop(s.get_pop());
    place_start_cells(s.get_pop());

    return s;
}

simulation load_sim_last_pop(int seed, int change_freq)
{
    auto sim_par_filename = create_sim_par_name(seed,change_freq);
    assert(exists(sim_par_filename));
    simulation s{load_sim_parameters(sim_par_filename)};

    auto last_pop_name = create_last_pop_name(seed, change_freq);
    auto last_pop = load_funders(last_pop_name);
    s.get_pop().get_v_ind().resize(last_pop.get_v_funder_data().size());

    for(size_t i = 0; i != last_pop.get_v_funder_data().size(); i++)
    {
        s.get_pop().get_v_ind()[i].get_grn()
                = last_pop.get_v_funder_data()[i].get_grn();
    }

    update_radius_pop(s.get_pop());
    place_start_cells(s.get_pop());

    //Change internal state of rng member of simulation
    //to avoid pseudo rng to replicate same exact mutations
    //in different kinds runs.
    //I do it based on the weight of the last connection of
    //last ind in the pop, so that runs of the same type
    //replicate consistently.
    auto scramble_seed = s.get_pop().get_v_ind().back().get_grn().get_H2O().back().back();
    auto scramble_rng = std::minstd_rand(scramble_seed);
    auto scramble_n = std::uniform_int_distribution(50,100)(scramble_rng);
    for(int i = 0; i != scramble_n; i++)
    {
        s.get_rng()();
    }

    return s;
}

simulation load_sim_before_last_pop(int seed, int change_freq)
{
    auto sim_par_filename = create_sim_par_name(seed,change_freq);
    assert(exists(sim_par_filename));
    simulation s{load_sim_parameters(sim_par_filename)};

    auto last_pop_name = create_last_pop_name(seed, change_freq);
    auto last_pop = load_funders(last_pop_name);
    s.get_pop().get_v_ind().resize(last_pop.get_v_funder_data().size());

    for(size_t i = 0; i != last_pop.get_v_funder_data().size(); i++)
    {
        s.get_pop().get_v_ind()[i].get_grn()
                = last_pop.get_v_funder_data()[i].get_grn();
    }

    update_radius_pop(s.get_pop());
    place_start_cells(s.get_pop());

    //Change internal state of rng member of simulation
    //to avoid pseudo rng to replicate same exact mutations
    //in different kinds runs.
    //I do it based on the weight of the last connection of
    //last ind in the pop, so that runs of the same type
    //replicate consistently.
    auto scramble_seed = s.get_pop().get_v_ind().back().get_grn().get_H2O().back().back();
    auto scramble_rng = std::minstd_rand(scramble_seed);
    auto scramble_n = std::uniform_int_distribution(50,100)(scramble_rng);
    for(int i = 0; i != scramble_n; i++)
    {
        s.get_rng()();
    }

    return s;
}

simulation load_best_ind_for_rand_cond(int seed, int change_freq)
{
    auto name = create_sim_par_name(seed,change_freq);

    simulation s{load_sim_parameters(name)};

    auto best_ind_grn = find_last_gen_best_ind_grn(load_funders_success(create_funders_success_name(seed, change_freq)));

    s.get_pop().get_v_ind().resize(s.get_pop().get_param().get_exp_new_pop_size());

    for( auto& individual : s.get_pop().get_v_ind())
    {
        individual.get_grn() = best_ind_grn;
    }
    update_radius_pop(s.get_pop());

    place_start_cells(s.get_pop());

    return s;
}

simulation no_dem_and_fund_copy(const simulation& s)
{
    simulation new_s{sim_param{s.get_env().get_param(),
                    s.get_pop().get_v_ind().begin()->get_param(),
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

void reproduce_cycle_env(simulation&s, int cycle)
{
    assert(static_cast<size_t>(cycle) < s.get_demo_sim().get_demo_cycles().size());
    environment env_generation{s.get_demo_sim().get_demo_cycles()[cycle].get_env_param()};
    s.get_env() = env_generation;
}

void reproduce_cycle(simulation&s, int cycle)
{
    reproduce_cycle_env(s, cycle);
    s.get_pop().set_pop_inds(pop_from_funders(s.get_funders_success(), s.get_demo_sim(), cycle));

    update_radius_pop(s.get_pop());
    place_start_cells(s.get_pop());

}

void reproduce_rand_cond(simulation&s, const std::vector<std::pair<env_param, ind_param>>& rand_cond, int n_rand_cond)
{
    if(n_rand_cond < static_cast<int>(rand_cond.size()))
    {
        set_new_params(s, rand_cond[n_rand_cond].first, rand_cond[n_rand_cond].second);
    }
    else
    {
        std::cerr << "Random condition n: " << n_rand_cond << "does not exist";
        abort();
    }
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

        auto index = find_grid_index(ind,s.get_env().get_param().get_grid_side());

        if(ind.get_phen() == phenotype::spore)
        {
        }
        else if(index == -100)//When individual is outside grid
        {
            responds(ind, env_grid_cell(0,0,0,0));
        }
        else
        {
            responds(ind, s.get_env().get_cell(index));
        }
    }
}

demographic_sim run_test_random_conditions(const simulation& s,
                                           int n_number_rand_cond,
                                           int pop_max,
                                           double amplitude,
                                           std::string name)
{
    auto random_conditions = create_rand_conditions_unif(
                s.get_env().get_param(),
                s.get_pop().get_v_ind().begin()->get_param(),
                n_number_rand_cond,
                amplitude,
                0);

    auto test_pop = s.get_pop().get_v_ind();

    simulation rand_s = no_dem_and_fund_copy(s);
    rand_s.get_meta_param().get_pop_max() = pop_max;

    int counter = 0;

    for(const auto & condition : random_conditions)
    {
        assert(rand_s.get_pop().get_v_ind() == test_pop);
        assert(rand_s.get_env().get_grid() == s.get_env().get_grid());
        set_new_params(rand_s, condition.first, condition.second);

        auto start = std::chrono::high_resolution_clock::now();

        exec_cycle(rand_s);
        rand_s.tick_cycles();
        rand_s.reset_timesteps();

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration<float>(stop - start);
        std::cout<< "condition n"<< counter <<": " << duration.count() << std::endl;
        auto tot_inds = rand_s.get_demo_sim().get_demo_cycles().back().get_n_actives() +
                rand_s.get_demo_sim().get_demo_cycles().back().get_n_spores() +
                rand_s.get_demo_sim().get_demo_cycles().back().get_n_sporulating();
        std::cout<< "n individuals:" << tot_inds << std::endl << std::endl;

        rand_s.get_pop().get_v_ind() = test_pop;
        counter++;
        //resave every time all currently obtained results from test over all currently tested conditions
        save_demographic_sim(rand_s.get_demo_sim(), name);
    }
    return rand_s.get_demo_sim();
}

demographic_sim run_evo_random_conditions(const simulation& s,
                                          int number_of_sequences,
                                          int cond_per_seq,
                                          int seq_index,
                                          int pop_max,
                                          double amplitude,
                                          std::string prefix)
{
    auto random_conditions = create_rand_conditions_matrix_extreme(
                s.get_env().get_param(),
                ind_param{},
                number_of_sequences,
                cond_per_seq,
                amplitude);

    simulation rand_s = no_dem_and_fund_copy(s);

    ///For now do not allow the env and ind param
    /// to change, no matter the freq of change.
    rand_s.get_meta_param().get_change_freq() = 0;

    rand_s.get_meta_param().get_pop_max() = pop_max;

    auto start = std::chrono::high_resolution_clock::now();

    exec_change(rand_s, random_conditions[seq_index]);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    std::cout<< "duration : " << duration.count() << "s" << std::endl <<
                "n_individuals: " << s.get_pop().get_pop_size() << std::endl <<
                std::endl;

    std::cout << "saving demographics" << std::endl;
    save_demographic_sim(rand_s.get_demo_sim(), create_sim_demo_name(s, prefix));
    std::cout << "saving funders" << std::endl;
    save_funders_success(rand_s.get_funders_success(), create_funders_success_name(rand_s, prefix));

    return rand_s.get_demo_sim();
}

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
                      int pop_max,
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
    run_test_random_conditions(rand_s,
                               n_random_conditions,
                               pop_max,
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
                 int pop_max,
                 int seed,
                 bool overwrite)
{
    auto rand_s = load_sim_last_pop(seed,change_frequency);
    auto name = create_test_random_condition_name(rand_s,amplitude);

    if(exists(name) && !overwrite)
    {
        std::cout << "The random conditions for this simulation"
                     " have already been tested!" << std::endl;
        return 0;
    }

    auto rand_start = std::chrono::high_resolution_clock::now();
    run_test_random_conditions(rand_s,
                               n_random_conditions,
                               pop_max,
                               amplitude,
                               name);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - rand_start);
    std::cout<< "random condition test :" << duration.count() << "s" << std::endl;
    return 0;
}

int run_sim_evo_rand(double amplitude,
                     int change_frequency,
                     int num_of_sequences,
                     int conditions_per_seq,
                     int pop_max,
                     int seed,
                     int seq_index,
                     bool overwrite)
{

    auto s = load_sim_last_pop(seed,change_frequency);

    auto prefix = "rand_evo_extreme_a" + std::to_string(amplitude) +
            "seq_" + std::to_string(seq_index) +
            "cond_per_seq" + std::to_string(conditions_per_seq);

    if(exists(prefix + create_sim_demo_name(s)) && !overwrite)
    {
        std::cout << "The evolution in random conditions for this simulation"
                     " have already been tested!" << std::endl;
        return 0;
    }

    auto start = std::chrono::high_resolution_clock::now();
    run_evo_random_conditions(s,
                              num_of_sequences,
                              conditions_per_seq,
                              seq_index,
                              pop_max,
                              amplitude,
                              prefix);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    std::cout<< "random condition evo :" << duration.count() << "s" << std::endl;
    return 0;
}

int run_sim_evo(const env_param& e,
                const ind_param& i,
                const meta_param& m,
                const pop_param& p,
                bool overwrite)
{
    if(exists(create_sim_par_name(m.get_seed(),
                                  m.get_change_freq()))
            && !overwrite)
    {
        std::cout << "this simulation has already been run" << std::endl;
        return 0;
    }

    simulation s{sim_param{e, i, m, p}};

    auto start = std::chrono::high_resolution_clock::now();

    exec(s);
    save_data(s);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    std::cout << "simualtion :"<< duration.count() << "s" << std::endl;
    return 0;
}

int run_standard(const env_param& e,
                 const ind_param& i,
                 const meta_param& m,
                 const pop_param& p,
                 double amplitude,
                 int change_frequency,
                 int n_random_conditions,
                 int pop_max,
                 int seed)
{
    simulation s{sim_param{e, i, m, p}};
    auto start = std::chrono::high_resolution_clock::now();

    exec(s);
    save_data(s);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration<float>(stop - start);
    std::cout << "simualtion :"<< duration.count() << "s" << std::endl;

    auto rand_start = std::chrono::high_resolution_clock::now();

    auto rand_s = load_sim_last_pop(seed,change_frequency);
    run_test_random_conditions(rand_s, n_random_conditions, pop_max, amplitude, "standard_rand_run.csv");

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration<float>(stop - rand_start);
    std::cout<< "random condition test :" << duration.count() << "s" << std::endl;

    duration = std::chrono::duration<float>(stop - start);
    std::cout<< "overall time :" << duration.count() << "s" << std::endl;
    return 0;
}

void save_vector_of_rand_cond(const std::vector<std::pair<env_param, ind_param>>& rand_cond_v,
                              const std::string& filename)
{
    std::ofstream os{filename};
    for(const auto& condition : rand_cond_v)
    {
        os << condition.first << " , "
           << condition.second << std::endl;
    }
}

void save_data(const simulation& s)
{
    std::string sim_param_name = create_sim_par_name(s);

    sim_param s_p{s.get_env().get_param(),
                s.get_pop().get_ind(0).get_param(),
                s.get_meta_param(),
                s.get_pop().get_param()};

    save_sim_parameters(s_p, sim_param_name);

    save_two_last_pops(s);

    std::string sim_demo_name = create_sim_demo_name(s);

    save_demographic_sim(s.get_demo_sim(), sim_demo_name);

    std::string funders_success_name = create_funders_success_name(s);

    save_funders_success(s.get_funders_success(), funders_success_name);
}


funders save_before_last_pop(const simulation& s, std::string prefix)
{
    std::string before_last_pop_name = create_before_last_pop_name(s, prefix);
    save_funders(s.get_funders_success().get_v_funders().back(), before_last_pop_name);
    return s.get_funders_success().get_v_funders().back();
}

funders save_last_pop(const simulation& s, const std::string& prefix)
{
    std::string last_pop_name =  create_last_pop_name(s, prefix);
    auto real_last_new_funders = prepare_funders(s);
    save_funders(real_last_new_funders, last_pop_name);
    return real_last_new_funders;
}

std::vector<funders> save_two_last_pops(const simulation& s, const std::string& prefix)
{
    auto before_last = save_before_last_pop(s, prefix);
    auto last = save_last_pop(s, prefix);
    return std::vector<funders>{before_last, last};
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

void set_new_params(simulation& s, const env_param& e, const ind_param& i)
{
    s.get_env().set_new_env_param(e);
    s.get_pop().set_pop_inds(set_new_ind_par(s.get_pop().get_v_ind(),i));
}

void store_demographics( simulation& s) noexcept
{
    s.set_demo_sim(update_demographics(s));
}

int tick(simulation& s, int n_ticks)
{
    int time = 0;
    response(s);
    //feeding(s);
    jordi_feeding(s);
    metabolism_pop(s.get_pop());
    secretion_metabolite(s);
    //death(s.get_pop());
    //jordi_death(s.get_pop());
    update_radius_pop(s.get_pop());
    auto division_happens = division(s.get_pop());
    if( n_ticks > 0 && s.get_timestep() % n_ticks == 0)
    {
        manage_static_collisions(s.get_pop());
    }
    else if(n_ticks == 0 && division_happens)
    {
        manage_static_collisions(s.get_pop());
    }
    degradation_metabolite(s.get_env());
    diffusion(s.get_env());
    s.tick_timesteps();
    return time;
}

demographic_sim update_demographics(const simulation& s) noexcept
{
    auto d_s = s.get_demo_sim();
    d_s.get_demo_cycles().push_back(demographics(s, s.get_env().get_param()));
    return d_s;
}















