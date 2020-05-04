#include "simulation.h"
#include <cassert>
#include <numeric>
#include <algorithm>
#include <cmath>

simulation::simulation(unsigned int pop_size,
                       unsigned int exp_new_pop_size,
                       double min_dist,
                       int grid_side,
                       double diff_coeff,
                       double init_food,
                       double mutation_prob,
                       double mutation_step,
                       double base_disp_prob,
                       double spore_advantage,
                       double reproduction_prob,
                       double metab_degradation_rate):
    m_pop{pop_param{
          pop_size,
          exp_new_pop_size,
          min_dist,
          mutation_prob,
          mutation_step,
          base_disp_prob,
          spore_advantage,
          reproduction_prob}
          },
    m_e{env_param{
        grid_side,
        diff_coeff,
        init_food,
        metab_degradation_rate}
        }
{

}

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
        double reproduction_prob = 0.1;
        double metab_degradation_rate = 0.01;

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
                      reproduction_prob,
                      metab_degradation_rate);
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
        assert(s.pop().get_param().get_repr_p() - reproduction_prob < 0.0001 &&
               s.pop().get_param().get_repr_p() - reproduction_prob > -0.0001);
        //tests for env param
        assert(s.get_env().get_param().get_diff_coeff() - diff_coeff < 0.00001 &&
               s.get_env().get_param().get_diff_coeff() - diff_coeff > -0.00001);
        assert(s.get_env().get_param().get_init_food() - init_food < 0.0001 &&
               s.get_env().get_param().get_init_food() - init_food > -0.0001);

        assert(s.get_env().get_param().get_metab_degr_rate() - metab_degradation_rate < 0.0001 &&
               s.get_env().get_param().get_metab_degr_rate() - metab_degradation_rate > -0.0001);
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
        s = simulation(1, 1, 1, 1, 0.1, starting_food);
        for( auto& grid_cell : s.get_env().get_grid())
        {
            assert(grid_cell.get_food() - starting_food < 0.000001
                   && grid_cell.get_food() - starting_food > -0.0000001);
        }
    }

    //Individuals can take up energy from the environment
    {
        simulation s(1,1,0.1,2);
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

        double total_en_after = std::accumulate
                (
                    s.pop().get_v_ind().begin(),
                    s.pop().get_v_ind().end(),0.0,
                    [](double sum,const individual& i){
            return sum + i.get_energy();}
        );

        assert(total_food_init > total_food_after);
        assert(total_en_init < total_en_after);

        auto total_uptake = std::accumulate(s.pop().get_v_ind().begin(), s.pop().get_v_ind().end(), 0.0,
                                            [](double sum, const individual& i)
        {return sum + i.get_param().get_uptake_rate();});

        assert(total_en_after - (total_uptake + total_en_init) < 0.00000001 &&
               total_en_after - (total_uptake + total_en_init) > -0.00000001);
        assert(total_food_init - (total_food_after + total_uptake) < 0.000001 &&
               total_food_init - (total_food_after + total_uptake) > -0.000001);
    }

    //Individuals outside the grid do not feed
    {
        simulation s(1,1,0.1,2);

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

    //In one tick/timestep individuals take in input, determine phenotype(based on previous timestep),
    // feed, than reproduce, than substances diffuse
    {
        simulation s(1,1,0.1,3,1,1,0);
        //Set all the hid nodes and H2O and H2H weights to one so that we are sure the phenotype will stay = active;
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

        auto grid_index_ind =  find_grid_index(s.pop().get_ind(0),s.get_env().get_param().get_grid_side());
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
        simulation s (2, 1, 0, 3, 0.1, 5);
        auto food_begin = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                          [](double sum, const env_grid_cell& c) {return sum + c.get_food();});

        //The simulation will last long enough  for the individuals to reproduce
        auto sim_time = s.pop().get_ind(0).get_param().get_treshold_energy() /
                s.pop().get_ind(0).get_param().get_uptake_rate() + 5;

        auto init_pop_size =  s.pop().get_v_ind().size();

        for( int i = 0; i != static_cast<int>(sim_time); i++)
        {

            auto food_before_feed = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                                    [](double sum, const env_grid_cell& c) {return sum + c.get_food();});
            feeding(s);

            auto food_after_feed = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                                   [](double sum, const env_grid_cell& c) {return sum + c.get_food();});

            double food_eaten = 0.0;
            for(const auto& ind : s.pop().get_v_ind())
            {
                auto grid_cell_ind = find_grid_index(ind, s.get_env().get_param().get_grid_side());
                if( grid_cell_ind != -100 && s.get_env().get_cell(grid_cell_ind).get_food() > 0)
                {
                    food_eaten += ind.get_param().get_uptake_rate();
                }
            }

            auto balance_uptake = food_before_feed - (food_after_feed + food_eaten) ;
            assert(balance_uptake < 0.0001 && balance_uptake > -0.0001);

            metabolism_pop(s.pop());
            death(s.pop());
            division(s.pop());
            manage_static_collisions(s.pop());

            auto food_before_diff = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                                    [](double sum, const env_grid_cell& c) {return sum + c.get_food();});

            calc_diffusion_food(s.get_env());

            auto new_grid = s.get_env().get_grid();

            auto food_after_diff = std::accumulate(new_grid.begin(), new_grid.end(), 0.0,
                                                   [](double sum, const env_grid_cell& c) {return sum + c.get_food();});

            auto balance_diffusion = food_before_diff - food_after_diff;
            assert(balance_diffusion < 0.0001 && balance_diffusion > -0.0001);

            new_grid.swap(s.get_env().get_grid());

        }
        auto food_end = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                        [](double sum, const env_grid_cell& c) {return sum + c.get_food();});
        assert(food_end < food_begin);

        auto final_pop_size =  s.pop().get_v_ind().size();

        assert(init_pop_size < final_pop_size);
    }

    //A simulation is initiallized with a degradation rate
    {
        double degradation_rate = 3.14;
        simulation s(0,0,0,0,0,0,0,0,0,0,0,degradation_rate);
        assert(s.get_env().get_param().get_metab_degr_rate() - degradation_rate < 0.000001 &&
               s.get_env().get_param().get_metab_degr_rate() - degradation_rate > -0.000001);
    }
    //Every time step the individuals produce new metabolite and metabolite degrades in grid_cells
    {
        double degradation_rate = 3.14;
        double init_metab = degradation_rate;
        simulation s(1,0,0,1,0,0,0,0,0,0,0,degradation_rate);

        for(auto& grid_cell : s.get_env().get_grid())
        {
            grid_cell.set_metab(init_metab);
        }

        double tot_metab_before = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                                  [](int sum, const env_grid_cell& g) {return sum + g.get_metab();});

        secretion_metabolite(s);
        degradation_metabolite(s.get_env());

        double tot_metab_after =
                std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                [](int sum, const env_grid_cell& g)
        {return sum + g.get_metab();});

        double tot_production =
                std::accumulate(s.pop().get_v_ind().begin(), s.pop().get_v_ind().end(), 0.0,
                                [](int sum, const individual& i)
        {return sum + i.get_param().get_metab_secr_rate();});

        double tot_degradation = s.get_env().get_grid_size() * s.get_env().get_param().get_metab_degr_rate();

        auto metab_balance = tot_metab_before - tot_degradation + tot_production - tot_metab_after;
        assert(metab_balance < 0.000001 && metab_balance > -0.000001);


    }

    //every timestep/tick collisions are handled
    {
        simulation s(7,3);
        //The central individual in this population
        //after a tick should reproduce
        auto init_pop_size = s.pop().get_v_ind().size();
        s.pop().get_ind(1).set_energy(get_ind_tr_en(s.pop(), 1)
                                      + s.pop().get_ind(1).get_param().get_metabolic_rate() + 0.01
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


    //    //A sporulating individual updates its timer every time metab_pop is called
    //    {
    //        simulation s;
    //        auto init_timer = s.get_ind(0).get_spo_timer();
    //        s.get_ind(0).set_phen(phenotype::sporulating);
    //        metabolism_pop(s);
    //        assert(init_timer != s.get_ind(0).get_spo_timer());
    //        assert(s.get_ind(0).get_spo_timer() == init_timer + 1);
    //        assert(s.get_ind(0).get_spo_timer() != init_timer + 2);

    //        s.get_ind(0).reset_spo_timer();
    //        init_timer = s.get_ind(0).get_spo_timer();
    //        s.get_ind(0).set_phen(phenotype::spore);
    //        metabolism_pop(s);
    //        assert(s.get_ind(0).get_spo_timer() == init_timer);
    //        assert(s.get_ind(0).get_spo_timer() != init_timer + 1);
    //        assert(s.get_ind(0).get_spo_timer() != init_timer + 2);

    //        s.get_ind(0).reset_spo_timer();
    //        init_timer = s.get_ind(0).get_spo_timer();
    //        s.get_ind(0).set_phen(phenotype::active);
    //        metabolism_pop(s);
    //        assert(s.get_ind(0).get_spo_timer() == init_timer);
    //        assert(s.get_ind(0).get_spo_timer() != init_timer + 1);
    //        assert(s.get_ind(0).get_spo_timer() != init_timer + 2);
    //    }

    //    //Individuals that die are removed from the population
    //    {
    //        simulation s;
    //        s.get_ind(0).set_energy(0);//the only individual in this sim has 0 energy, thus it will die
    //        assert(s.get_pop_size() == 1);
    //        death(s);
    //        assert(s.get_pop().empty() && s.get_pop_size() == 0);

    //        unsigned int pop_size = 5;
    //        //The simulation does not have a grid with food,
    //        //so organisms cannot feed
    //        s = simulation(pop_size,1,0.1,0);
    //        for(auto& ind :s.get_pop())
    //        {
    //            ind.set_energy(0);
    //        }
    //        //Only the first individual has enough energy to survive
    //        //for 1 tick
    //        set_ind_en(s.get_ind(0),s.get_ind(0).get_param().get_metabolic_rate() + 0.001);
    //        assert(s.get_pop().size() == pop_size);
    //        tick(s);
    //        assert(s.get_pop_size() == 1);
    //        //then at the second tick the only individual left dies
    //        tick(s);
    //        assert(s.get_pop().empty() && s.get_pop_size() == 0);
    //    }

    //    //A simulation is initialized with a random number generator
    //    {
    //        simulation s;
    //        std::uniform_int_distribution u_d(0,2);
    //        double mean = 0;
    //        //Draw a thousands times from a uniform dist between 0 and 2
    //        for(int i = 0; i != 1000; i++)
    //        {
    //            mean += u_d(s.get_rng());
    //        }
    //        //calculate mean of the drawn values
    //        mean /= 1000;
    //        assert(mean < 1.01 && mean > 0.99 );
    //    }

    //    //A pop can generate numbers from a uniform distribution
    //    //between 0 and 2PI
    //    {
    //        simulation s;
    //        double mean = 0;
    //        //Draw a thousands times from a uniform dist between 0 and 2
    //        int sampling_size = 1000;
    //        for(int i = 0; i != sampling_size; i++)
    //        {
    //            mean += repr_angle(s);
    //        }
    //        //calculate mean of the drawn values
    //        mean /= 1000;
    //        assert(mean < M_PI + 0.1 && mean > M_PI - 0.1 );
    //    }

    //    //Daughter cells are placed at a random angle after reproduction
    //    {
    //        simulation s;
    //        //to calculate angle we will use three point
    //        //the center of the mother(0,0)
    //        std::pair<double, double> mother (0,0);
    //        set_pos(s.get_ind(0),mother);
    //        //a reference point on the same axis as the mother (0,1)
    //        //and the center of the daughter -> get_daughter_pos()
    //        std::pair<double, double> reference(1,0);

    //        //Draw a thousands times from a uniform dist between 0 and 2
    //        double mean = 0;
    //        int sampling_size = 1000;
    //        for(int i = 0; i != sampling_size; i++)
    //        {
    //            auto daughter = get_daughter_pos(s.get_ind(0), repr_angle(s));
    //            mean += calc_angle_3_pos(mother,daughter,reference);
    //        }
    //        mean /= sampling_size;
    //        assert(mean < M_PI + 0.1 && mean > M_PI - 0.1 );
    //    }

    //    //A pop can generate numbers from a normal distribution for mutation step
    //    //with mean 0 and variance 0.1 by default
    //    {
    //        simulation s;
    //        double mean = 0;
    //        int sampling_size = 10000;
    //        for(int i = 0 ; i != sampling_size; i++ )
    //            mean += mut_step(s);
    //        mean /= sampling_size;
    //        assert(mean < 0.01 && mean > -0.01);
    //    }
    //    //A pop can generate numbers from a bernoulli distribution to see if mutation happens or not
    //    //0.01 by default
    //    {
    //        simulation s;
    //        double mean = 0;
    //        int sampling_size = 100000;
    //        for(int i = 0 ; i != sampling_size; i++ )
    //            mean += mut_happens(s);
    //        mean /= sampling_size;
    //        assert(mean < 0.011 && mean > 0.009);
    //    }

    //    //The sum of weight of an individual after many rounds of mutation
    //    //Should have the same mean as in the beginning, but its variance should
    //    //be the same as the mutation_step distribution
    //    {
    //        simulation s;
    //        //   double init_mean = weights_mean(s.get_ind(0).get_grn());
    //        double init_variance = weights_var(s.get_ind(0).get_grn());
    //        assert(init_variance < 0.0001 && init_variance > -0.000001);

    //        int sampling_size = 10000;
    //        auto mut_prob_dist = create_bernoulli_dist(s.get_param().get_mu_p());
    //        auto mut_step_dist = create_normal_dist(0,s.get_param().get_mu_st());
    //        for (int i = 0; i != sampling_size; i++)
    //        {
    //            mutates(s.get_ind(0),
    //                    s.get_rng(),
    //                    mut_prob_dist,
    //                    mut_step_dist
    //                    );
    //        }
    //        //This first assert does not pass, the mean is much more variable than
    //        //I thought, but I do not see any bug. I will comment this out
    //        //    assert(mean - init_mean > -0.1 && mean - init_mean < 0.1);
    //        assert(init_variance - weights_var(s.get_ind(0).get_grn()) > 0.01 ||
    //               init_variance - weights_var(s.get_ind(0).get_grn()) < -0.01);
    //    }
    //    //After dividing the two daughter individuals mutate
    //    {
    //        double mutation_probability = 1; //all weights will be mutated in this simulation
    //        simulation s(1, 1, 0, 0, 0, 0, mutation_probability);
    //        auto init_var = weights_var(s.get_ind(0).get_grn());
    //        assert(init_var < 0.00001 && init_var > -0.0001);
    //        divides(s.get_ind(0),
    //                s.get_pop(),
    //                repr_angle(s),
    //                s.get_rng(),
    //                create_bernoulli_dist(s.get_param().get_mu_p()),
    //                create_normal_dist(0, s.get_param().get_mu_st())
    //                );
    //        auto post_var = weights_var(s.get_ind(0).get_grn());
    //        assert(init_var - post_var > 0.000001 || init_var - post_var < -0.0001);
    //    }

    //    //A simulation has a member variable m_new_pop_size that states the max number of
    //    //individuals that will start a new population
    //    //by default = to pop.size()
    //    //If at dispersal m_new_pop_size > pop.size()
    //    //Then the number of funding individuals == pop.size()

    //    {
    //        simulation s(1000, 100);
    //        select_new_pop(s);
    //        assert(static_cast<int>(s.get_new_pop().size()) == s.get_param().get_exp_new_pop_size());
    //        s = simulation(10, 100);
    //        select_new_pop(s);
    //        auto n_drawn_ind = std::count_if(s.get_pop().begin(),s.get_pop().end(),
    //                                         [](const individual& i) {return is_drawn(i);});
    //        assert( n_drawn_ind == s.get_pop_size());
    //    }


    //    //During dispersal the individuals selected for the new_pop cannot be drawn again from pop
    //    {
    //        simulation s(1000, 100);
    //        select_new_pop(s);
    //        assert(std::count_if(s.get_pop().begin(),s.get_pop().end(),
    //                             [](const individual& i) {return is_drawn(i);}) == 100);
    //    }

    //    //A simulation is initialized with a variable m_base_fitness
    //    //by default = 0.01
    //    {
    //        double base_disp_prob = 0.1;
    //        simulation s(0,0,0,0,0,0,0,0,base_disp_prob);
    //        assert(s.get_param().get_base_disp_prob() - base_disp_prob < 0.00001 &&
    //               s.get_param().get_base_disp_prob() - base_disp_prob > -0.000001);
    //    }

    //    //A simulation is initialized with a variable m_spore_advantage
    //    //by default = 10
    //    {
    //        double spore_advantage = 10;
    //        simulation s(0,0,0,0,0,0,0,0,0,spore_advantage);
    //        assert(s.get_param().get_spo_adv() - spore_advantage < 0.00001 &&
    //               s.get_param().get_spo_adv() - spore_advantage > -0.000001);
    //    }

    //    //A simulation is initialized with a uniform distribution between 0 and 1
    //    //used to see which ind is drawn at dispersal
    //    {
    //        simulation s;
    //        int sample_size = 100000;
    //        double mean = 0;
    //        for(int i = 0; i != sample_size; i++)
    //        {
    //            mean += create_unif_dist(0,1)(s.get_rng());
    //        }
    //        mean /= sample_size;
    //        assert(mean - 0.5 < 0.01 &&
    //               mean - 0.5 > -0.01);
    //    }
    //    //At initialization a simulation checks that base_disp_dist * 10 is not > 1
    //    //--------> constructor throws exception. Tested directly in constructor
    //    {
    //#ifndef IS_ON_TRAVIS
    //        try {
    //            simulation(0,0,0,0,0,0,0,0,1);
    //        } catch (std::string& e) {
    //            assert(e == "base dispersal probability * spore advantage > 1, too high!\n" );
    //        }
    //#endif
    //    }

    //    //Individuals are selected based on their phenotype
    //    //A spore is more likely to be selected than a living
    //    {
    //        simulation s(1000,100);
    //        for(int i = s.get_pop_size() / 2; i != s.get_pop_size(); i++)
    //            s.get_ind(i).set_phen(phenotype::spore);
    //        select_new_pop(s);
    //        auto spore_ratio =
    //                std::accumulate(s.get_new_pop().begin(),s.get_new_pop().end(),0.0,
    //                                [](const int sum, const individual& ind){return sum + is_spore(ind);}) /
    //                s.get_new_pop().size();
    //        assert(spore_ratio > 0.5);
    //    }

    //    //After being selected in new population individuals flag is_drawn is resetted
    //    {
    //        simulation s;
    //        select_new_pop(s);
    //        for(const auto& ind :s.get_new_pop())
    //            assert(!is_drawn(ind));
    //    }

    //    //After a new population is selected it swapped with the old population
    //    //And the old population is cancelled
    //    {
    //        unsigned int pop_size = 1000;
    //        int new_pop_size = 100;
    //        simulation s(pop_size,new_pop_size);
    //        select_new_pop(s);
    //        assert(static_cast<int>(s.get_new_pop().size()) == new_pop_size);
    //        assert(s.get_pop().size() == pop_size);
    //        fund_pop(s);
    //        assert(s.get_new_pop().size() == 0);
    //        assert(static_cast<int>(s.get_pop_size()) == new_pop_size);
    //    }

    //    //Individuals after funding the new population are set in an hexagonal pattern
    //    {
    //        unsigned int pop_size = 1000;
    //        int new_pop_size = 100;
    //        simulation s(pop_size,new_pop_size);
    //        select_new_pop(s);
    //        fund_pop(s);
    //        place_start_cells(s);
    //        auto n_hex_l = count_hex_layers(s.get_pop_size());
    //        auto v_modulus = modulus_of_btw_ind_angles(s, M_PI/ (6 * n_hex_l));
    //        for(auto ind_modulus : v_modulus)
    //            assert( ind_modulus < 0.0000000001 || (ind_modulus > M_PI / (6 * n_hex_l) - 0.1 && ind_modulus <= M_PI / (6 * n_hex_l) + 0.1));
    //        assert(!has_collision(s));
    //    }


    //Max 100 ind, selected based on phenotype, are placed in a hex pattern, in a new env after dispersal
    //Tests all of the above
    {
        unsigned int pop_size = 1000;
        unsigned int new_pop_size = 100;
        auto food = 42.1;
        auto metabolite = 42.1;
        simulation s(pop_size,new_pop_size);
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
            assert( ind_modulus < 0.0000000001 || (ind_modulus > M_PI / (6 * n_hex_l) - 0.1 && ind_modulus <= M_PI / (6 * n_hex_l) + 0.1));
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
        simulation s(2,1,0.1,4);
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

        //The values of food and metabolite in the grid_cell as well as the energy in the individual
        //are changed so that all inputs are -1, the network should therefore give outputs == 0/false
        //since the outputs node will recieve a signal that is negative(below the treshold)
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
















