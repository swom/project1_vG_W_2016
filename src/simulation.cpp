#define _USE_MATH_DEFINES
#include <cmath>
#include <cassert>
#include "simulation.h"
#include <cassert>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <sys/stat.h>
//#include <unistd.h>
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

    const auto n_tot_spores =  std::count_if(
                s.get_pop().get_v_ind().begin(),
                s.get_pop().get_v_ind().end(),
                [](const individual i)
    {return  i.get_phen() == phenotype::spore;});

    double tot_fitness = n_tot_spores * s.get_pop().get_param().get_spo_adv()
            + s.get_pop().get_pop_size() - n_tot_spores;

    for(auto& funder : funders.get_v_funder_data())
    {
        assert(funder.get_success() == 0);

        const auto n_non_spore_descendants = 
                std::count_if(
                    s.get_pop().get_v_ind().begin(),
                    s.get_pop().get_v_ind().end(),
                    [&funder](const individual i)
        {return (i.get_ancestor() == funder.get_ancestor_ID())
                    && i.get_phen() != phenotype::spore; });

        const auto n_spore_descendants =
                std::count_if(
                    s.get_pop().get_v_ind().begin(),
                    s.get_pop().get_v_ind().end(),
                    [&funder](const individual i)
        {return (i.get_ancestor() == funder.get_ancestor_ID())
                    && i.get_phen() == phenotype::spore;});

        auto fitness = static_cast<double>(n_non_spore_descendants +
                n_spore_descendants * s.get_pop().get_param().get_spo_adv());

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

std::string create_funder_success_name(const simulation& s)
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
    std::cerr << s.get_meta_param().get_n_cycles() << ' ';
    while(s.get_cycle() != s.get_meta_param().get_n_cycles())
    {
        std::cerr << s.get_cycle() << ' ';
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
    std::cerr << '\n';
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

simulation load_sim_for_rand_cond(int seed, int change_freq)
{
    simulation s{load_sim_parameters(create_sim_par_name(seed,change_freq))};

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

    auto best_ind_grn = find_best_ind_grn(load_funders_success(create_funder_success_name(seed, change_freq)));

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
                                      double amplitude)
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

    auto name = create_sim_demo_name(rand_s);

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

    std::string funders_success_name = create_funder_success_name(s);

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

















