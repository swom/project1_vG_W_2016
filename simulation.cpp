#include "simulation.h"
#include <cassert>
#include <numeric>
#include <algorithm>
#include <cmath>

simulation::simulation(sim_param param):
    m_pop(param.get_pop_param()),
    m_e{param.get_env_param()}

{

}

void dispersal(simulation &s)
{
    select_new_pop(s.pop());
    fund_pop(s.pop());
    place_start_cells(s.pop());
    reset_env(s.get_env());
}


void exec(simulation& s, int n_tick) noexcept
{
    while(s.get_tick() != n_tick){tick(s);}
}

void feeding(simulation& s)
{
    for(auto& ind : s.pop().get_v_ind())
    {
        auto index_grid = find_grid_index(ind,s.get_env().get_param().get_grid_side());
        if(index_grid == -100 ||
                ind.get_phen() != phenotype::active)
        {continue;}
        feed(ind,s.get_env().get_cell(index_grid));
    }
}

void response(simulation& s)
{
    for(auto& ind : s.pop().get_v_ind())
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

void secretion_metabolite(simulation& s)
{
    int index;
    for(const auto& ind : s.pop().get_v_ind())
    {
        index = find_grid_index(ind,s.get_env().get_param().get_grid_side());
        if(index == - 100)
        {
            continue;
        }
        secretes_metab(ind,s.get_env().get_cell(index));
    }
}


int tick(simulation& s)
{
    int time = 0;
    response(s);
    feeding(s);
    metabolism_pop(s.pop());
    secretion_metabolite(s);
    death(s.pop());
    if(division(s.pop()))
    {
        time += manage_static_collisions(s.pop());
    }
    degradation_metabolite(s.get_env());
    diffusion(s.get_env());
    s.update_sim_timer();
    return time;
}


void test_simulation()//!OCLINT tests may be many
{
#ifndef NDEBUG

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
        double metab_degr_rate = 0.01;

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
                      metab_degr_rate);
        simulation s(s_p);
        //tests for pop_param
        assert(s.pop().get_v_ind().size() == pop_size);
        assert(s.pop().get_param().get_exp_new_pop_size() == exp_new_pop_size);
        assert(s.pop().get_param().get_min_dist() - min_dist < 0.0001 &&
               s.pop().get_param().get_min_dist() - min_dist > -0.0001);
        assert(s.pop().get_param().get_mu_p() - mutation_prob < 0.0001 &&
               s.pop().get_param().get_mu_p() - mutation_prob > -0.0001);
        assert(s.pop().get_param().get_mu_st() - mutation_step < 0.0001 &&
               s.pop().get_param().get_mu_st() - mutation_step > -0.0001);
        assert(s.pop().get_param().get_base_disp_prob() - base_disp_prob < 0.00001 &&
               s.pop().get_param().get_base_disp_prob() - base_disp_prob > -0.00001);
        assert(s.pop().get_param().get_spo_adv() - spore_advantage < 0.0001 &&
               s.pop().get_param().get_spo_adv() - spore_advantage > -0.0001);
        assert(s.pop().get_param().get_repr_p() - repr_trsh < 0.0001 &&
               s.pop().get_param().get_repr_p() - repr_trsh > -0.0001);
        //tests for env param
        assert(s.get_env().get_param().get_diff_coeff() - diff_coeff < 0.00001 &&
               s.get_env().get_param().get_diff_coeff() - diff_coeff > -0.00001);
        assert(s.get_env().get_param().get_init_food() - init_food < 0.0001 &&
               s.get_env().get_param().get_init_food() - init_food > -0.0001);

        assert(s.get_env().get_param().get_metab_degr_rate() - metab_degr_rate < 0.0001 &&
               s.get_env().get_param().get_metab_degr_rate() - metab_degr_rate > -0.0001);
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
        s = simulation(sim_param{1, 1, 1, 1, 0.1, starting_food});
        for( auto& grid_cell : s.get_env().get_grid())
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
                    s.pop().get_v_ind().begin(),
                    s.pop().get_v_ind().end(),0.0,
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
                    s.pop().get_v_ind().begin(),
                    s.pop().get_v_ind().end(),0.0,
                    [](double sum,const individual& i){
            return sum + i.get_energy();}
        );

        assert(total_food_init > total_food_after);
        assert(total_en_init < total_en_after);

        auto total_uptake =
                std::accumulate(s.pop().get_v_ind().begin(),
                                s.pop().get_v_ind().end(), 0.0,
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
                    s.pop().get_v_ind().begin(),
                    s.pop().get_v_ind().end(),0.0,
                    [](double sum,const individual& i){return sum + i.get_energy();}
        );

        set_pos(s.pop().get_ind(0),std::pair<double,double>(-42,42));

        feeding(s);

        double total_food_after = std::accumulate
                (
                    s.get_env().get_grid().begin(),
                    s.get_env().get_grid().end(),0.0,
                    [](double sum, const env_grid_cell& c){return sum + c.get_food();}
        );

        double total_en_after = std::accumulate
                (
                    s.pop().get_v_ind().begin(),
                    s.pop().get_v_ind().end(),0.0,
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
        for(auto& ind : s.pop().get_v_ind())
        {
            ind.get_grn().set_all_hid_nodes(1);
            ind.get_grn().set_all_out_nodes(1);
            ind.get_grn().set_all_H2O(1);
            ind.get_grn().set_all_H2H(1);
        }

        //The single individual in this population
        //after a tick should reproduce
        auto init_pop_size = s.pop().get_v_ind().size();
        auto ind_en = get_ind_tr_en( s.pop(), 0)
                + s.pop().get_ind(0).get_param().get_metabolic_rate() + 0.01
                - s.pop().get_ind(0).get_param().get_uptake_rate();
        s.pop().get_ind(0).set_energy(ind_en);

        //and the grid_cell where it is should recieve
        //food nutrients
        response(s);
        feeding(s);

        auto grid_index_ind =  find_grid_index(s.pop().get_ind(0),
                                               s.get_env().get_param().get_grid_side());

        double food_after_feeding = s.get_env().get_cell(grid_index_ind).get_food();

        metabolism_pop(s.pop());
        death(s.pop());
        if(division(s.pop()))
        {
            manage_static_collisions(s.pop());
        }

        diffusion(s.get_env());
        double food_after_diffusion = s.get_env().get_cell(grid_index_ind).get_food();

        //The difference in food before and after the diffusion step is not 0
        assert(!(food_after_feeding - food_after_diffusion < 0.00001
                 && food_after_feeding - food_after_diffusion > -0.0001));
        assert(s.pop().get_v_ind().size() == 2 * init_pop_size);
        assert(!has_collision(s.pop() ));
    }

    //If nothing else happens, food should constantly decrease when cells are feeding
    {
        simulation s (sim_param{2, 1, 0, 3, 0.1, 5});
        auto food_begin =
                std::accumulate(s.get_env().get_grid().begin(),
                                s.get_env().get_grid().end(), 0.0,
                                [](double sum, const env_grid_cell& c)
        {return sum + c.get_food();});

        //The simulation will last long enough  for the individuals to reproduce
        auto sim_time = s.pop().get_ind(0).get_param().get_treshold_energy() /
                s.pop().get_ind(0).get_param().get_uptake_rate() + 5;

        auto init_pop_size =  s.pop().get_v_ind().size();

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
            for(const auto& ind : s.pop().get_v_ind())
            {
                auto grid_cell_ind = find_grid_index(ind, s.get_env().get_param().get_grid_side());
                if( grid_cell_ind != -100 && s.get_env().get_cell(grid_cell_ind).get_food() > 0)
                {
                    food_eaten += ind.get_param().get_uptake_rate();
                }
            }

            auto balance_uptake = food_before_feed - (food_after_feed + food_eaten);
            assert(balance_uptake < 0.0001 && balance_uptake > -0.0001);

            metabolism_pop(s.pop());
            death(s.pop());
            division(s.pop());
            manage_static_collisions(s.pop());

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

        auto final_pop_size =  s.pop().get_v_ind().size();

        assert(init_pop_size < final_pop_size);
    }

    //A simulation is initiallized with a degradation rate
    {
        double degradation_rate = 3.14;
        simulation s(sim_param{0,0,0,0,0,0,0,0,0,0,0,degradation_rate});
        assert(s.get_env().get_param().get_metab_degr_rate() - degradation_rate < 0.000001 &&
               s.get_env().get_param().get_metab_degr_rate() - degradation_rate > -0.000001);
    }
    //Every time step the individuals produce new metabolite and metabolite degrades in grid_cells
    {
        double degradation_rate = 3.14;
        double init_metab = degradation_rate;
        simulation s(sim_param{1,0,0,1,0,0,0,0,0,0,0,degradation_rate});

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
                std::accumulate(s.pop().get_v_ind().begin(),
                                s.pop().get_v_ind().end(), 0.0,
                                [](int sum, const individual& i)
        {return sum + i.get_param().get_metab_secr_rate();});

        double tot_degradation =
                s.get_env().get_grid_size() * s.get_env().get_param().get_metab_degr_rate();

        auto metab_balance =
                tot_metab_before - tot_degradation + tot_production - tot_metab_after;

        assert(metab_balance < 0.000001 && metab_balance > -0.000001);


    }

    //every timestep/tick collisions are handled
    {
        simulation s(sim_param{7,3});
        //The central individual in this population
        //after a tick should reproduce
        auto init_pop_size = s.pop().get_v_ind().size();
        s.pop().get_ind(1).set_energy(get_ind_tr_en(s.pop(), 1)
                                      + s.pop().get_ind(1).get_param().get_metabolic_rate()
                                      + 0.01
                                      - s.pop().get_ind(1).get_param().get_uptake_rate());
        feeding(s);
        metabolism_pop(s.pop());
        division(s.pop());
        manage_static_collisions(s.pop());

        assert(s.pop().get_v_ind().size() == init_pop_size + 1);
        assert(!has_collision(s.pop()));
    }

    //A simulation is initialized with a m_tick = 0;
    {
        simulation s;
        assert(s.get_tick() == 0);
    }

    //After each tick the simulation updates its m_tick
    //by one
    {
        simulation s;
        for(int i = 0; i != 3; i++)
        {
            assert(s.get_tick() == i);
            tick(s);
        }
    }

    //A simulation can be run for a certain amount of ticks
    {
        simulation s;
        auto n_ticks = 100;
        exec(s, n_ticks);
        assert(s.get_tick() - n_ticks == 0);
    }


    //Spores do not feed
    {
        simulation s;
        auto init_food = s.get_env().get_cell(0).get_food();
        s.pop().get_ind(0).set_phen(phenotype::spore);
        feeding(s);
        assert(init_food - s.get_env().get_cell(0).get_food() < 0.000001
               && init_food - s.get_env().get_cell(0).get_food() > -0.000001);
    }


    //Sporulating individuals do not feed but they lose energy
    {
        simulation s;
        s.pop().get_ind(0).set_phen(phenotype::sporulating);
        set_ind_en(s.pop().get_ind(0),1);
        auto init_en_ind0 = get_ind_en(s.pop(), 0);
        auto init_food = s.get_env().get_cell(0).get_food();
        feeding(s);
        metabolism_pop(s.pop());
        assert(init_en_ind0 - get_ind_en(s.pop(), 0) -
               s.pop().get_ind(0).get_param().get_metabolic_rate() < 0.000001
               &&
               init_en_ind0 - get_ind_en(s.pop(), 0) -
               s.pop().get_ind(0).get_param().get_metabolic_rate() > -0.000001);

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
        assert(s.pop().get_v_ind().size() == new_pop_size);
        //Hex pattern
        auto n_hex_l = count_hex_layers(s.pop().get_pop_size());
        auto v_modulus = modulus_of_btw_ind_angles(s.pop(), M_PI/ (6 * n_hex_l));
        for(auto ind_modulus : v_modulus)
        {
            assert( ind_modulus < 0.0000000001 ||
                    (ind_modulus > M_PI / (6 * n_hex_l) - 0.1 &&
                     ind_modulus <= M_PI / (6 * n_hex_l) + 0.1));
        }
        assert(!has_collision(s.pop()));
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
        for(auto& ind : s.pop().get_v_ind())
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
        for(auto& ind : s.pop().get_v_ind())
        {
            assert(is_sporulating(ind));
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
        for(auto& ind : s.pop().get_v_ind())
        {
            ind.set_energy(-1);
        }
        //1st respose, the individuals respond to the initial parameter
        //expected output value == 1
        response(s);
        for(auto& ind : s.pop().get_v_ind())
        {
            assert(is_active(ind));
        }
        //2nd response, the individual responds to the changed parameter(all 0s)
        //expected output value == 0
        response(s);
        for(auto& ind : s.pop().get_v_ind())
        {
            assert(!is_active(ind));
            assert(is_sporulating(ind));
        }
    }
#endif
}
















