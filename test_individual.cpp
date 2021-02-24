#include"tests.h"

void test_individual()//!OCLINT tests may be many
{
#ifndef NDEBUG
    //An individual should be initialized with the defined starting size
    {
        double starting_size = 14.0;
        individual i(ind_param{starting_size});
        assert(i.get_radius() - starting_size < 0.0000001);
    }

    //An individual should be initialized at a certain position
    {
        double x = 100;
        double y = 100;
        individual i(ind_param{}, x, y);
        assert(i.get_x() - x < 0.0000001);
        assert(i.get_y() - y < 0.0000001);
    }

    //An individual is initialized with an uptake value
    //0.1 by default
    {
        individual i{ind_param{}};
        assert(i.get_param().get_uptake_rate() - 0.1 < 0.000001);

        double uptake_rate = 0.3;
        individual i2(ind_param{0,0,uptake_rate});
        assert(i2.get_param().get_uptake_rate() - uptake_rate < 0.000001);
    }

    //An individual is initialized with a m_metab_rate
    //0.01 by default
    {
        individual i{ind_param{}};
        assert(i.get_param().get_metabolic_rate() - 0.01 < 0.000001
               && i.get_param().get_metabolic_rate() - 0.01 > -0.000001);
        double metabolic_rate = 2;
        individual i2(ind_param{0.8,10,0.1,0.1,0.002,metabolic_rate});
        assert(i2.get_param().get_metabolic_rate() - metabolic_rate < 0.000001
               && i2.get_param().get_metabolic_rate() - metabolic_rate > -0.000001);

    }
    //Individuals should be initialized with a internal energy value
    //by default 0.1
    {
        double internal_energy = 3.14;
        individual i{ind_param{},0,0,internal_energy};
        assert(i.get_energy() - internal_energy < 0.000001 &&
               i.get_energy() - internal_energy > -0.000001);
    }

    //An individual is initialized with a treshold level of energy
    {
        double treshold_energy = 3;
        individual i(ind_param{0, treshold_energy});
        assert(i.get_param().get_treshold_energy() - treshold_energy < 0.00000001);
    }

    // an individual's energy can be set
    {
        individual i{ind_param{}};
        double lhs = i.get_energy();
        double new_energy = 3 + i.get_energy();
        i.set_energy(new_energy);
        double rhs = i.get_energy();
        assert(std::abs(rhs - lhs)>0.0000001);
    }

    //an individual energy can be changed
    {
        individual i{ind_param{}};
        double init_en = i.get_energy();
        double en_change = 3.14;
        i.change_en(en_change);
        assert(i.get_energy() - en_change - init_en < 0.000001);
    }

    //Energy after reproduction should be half of the excess of energy
    {
        auto treshold = 2.0;
        individual i{ind_param{0,
                               treshold},
                     0,
                     0,
                     treshold * 2};//this individual should have energy in excess = treshold
        //after division
        double excess_energy = i.get_energy() - i.get_param().get_treshold_energy();
        auto split_energy = i.split_excess_energy();

        assert(split_energy - excess_energy/2 < 0.000000001 &&
               split_energy - excess_energy/2 > -0.000000001);

    }
    // An individual position can be extracted as a pair of doubles x,y
    {
        auto x = 0.0;
        auto y = 0.0;
        individual i(ind_param{}, x, y);
        assert(get_pos(i).first - i.get_x() < 0.00001);
        assert(get_pos(i).second - i.get_y() < 0.00001);
    }

    //The distance between two individuals can be calculated
    {
        std::vector<individual> pop(2, individual(ind_param{}, 0, 0));
        //place individuals along x axis at 1 distance from each other
        for (auto i = 0; i != static_cast<int>(pop.size()); i++)
        {
            pop[static_cast<unsigned int>(i)].set_x(i);
        }
        for (unsigned int i = 0; i != pop.size(); i++)
        {
            for (unsigned int j = 0; j != pop.size(); j++)
            {
                auto distance_i_j = distance(pop[i],pop[j]);
                if(i == j)
                {
                    assert(distance_i_j < 0.000001 && distance_i_j > -0.00001);
                }
                else if(i != j)
                {
                    assert(distance_i_j - 1 < 0.000001 && distance_i_j - 1 > -0.00001);
                }
            }
        }
    }

    //The position of the second daughter cell is just outside the mother cell
    {
        double wiggle_room = 0.00001;
        individual mom(ind_param{}, 0, 0);
        auto daughter_pos = calculate_daughter_pos(mom,0);
        individual daughter2(ind_param{}, daughter_pos.first, daughter_pos.second);
        assert(!are_colliding(mom,daughter2, wiggle_room));
    }

    //The overlap of two individuals can be calculated
    {
        double wiggle_room = 0.00001;
        //Two individuals overlap by half their radius
        individual lhs(ind_param{}, 0, 0);
        individual rhs(ind_param{}, 0,lhs.get_radius());
        assert(half_overlap(lhs, rhs, wiggle_room) - (lhs.get_radius() / 2 + wiggle_room) < 0.000001 &&
               half_overlap(lhs, rhs, wiggle_room) - (lhs.get_radius() / 2 + wiggle_room) > -0.000001 );
    }

    //Given  2 overlapping individuals lhs and rhs:
    //the displacement components that the lhs need to move away
    //by half the overlap between lhs and rhs
    //are the opposite of the ones needed by rhs to do the same
    {
        double wiggle_room = 0.00001;
        //2 overlapping individuals
        double radius = 0.5;
        individual lhs(ind_param{radius}, 0, 0);
        individual rhs(ind_param{radius}, 0, radius);
        auto lhs_disp = get_displacement(lhs, rhs, wiggle_room);
        auto rhs_disp = get_displacement(rhs, lhs, wiggle_room);
        auto x_component = lhs_disp.first + rhs_disp.first;
        auto y_component = lhs_disp.second + rhs_disp.second;
        assert(x_component < 0.00001 && x_component > -0.00001);
        assert(y_component < 0.00001 && y_component > -0.00001);
    }

    //Two cells can be displaced so not to overlap anymore
    {
        double wiggle_room = 0.00001;

        individual lhs(ind_param{}, 0.0001,0);
        individual rhs(ind_param{}, 0, 0);
        add_displacement(lhs, rhs, wiggle_room);
        add_displacement(rhs, lhs, wiggle_room);

        lhs.displace();
        rhs.displace();
        auto dist = distance(lhs,rhs);
        auto radii_sum = lhs.get_radius() + rhs.get_radius();
        auto space_btw_circles = dist - radii_sum;
        assert( space_btw_circles < 0.00001
                && space_btw_circles > -0.00001);
    }

    //The grid cell where an individual should be
    //can be deduced from individual's coordinates
    {
        individual i(ind_param{},0,0);
        int grid_side = 1;
        assert(find_grid_index(i, grid_side) == 0);
        grid_side = 3;
        assert(find_grid_index(i, grid_side) == 4);
        auto pos = std::pair<double,double>(-0.2, 1.2);
        set_pos(i,pos);
        assert(find_grid_index(i, grid_side) == 7);
    }

    //If an individual is past the grid side
    //find_grid_index() will return -100
    {
        individual i(0,0);
        int grid_side = 3;
        auto pos = std::pair<double,double>(1.7, 0);
        set_pos(i,pos);
        assert( find_grid_index(i, grid_side) == -100);

        pos = std::pair<double,double>(-1.7, 0);
        set_pos(i,pos);
        assert(find_grid_index(i, grid_side) == -100);
    }

    //If an individual is above or below the grid limit
    //find_grid_index() will return -100
    {
        individual i(0,0);
        int grid_side = 3;
        auto pos = std::pair<double,double>(0, -1.7);
        set_pos(i,pos);
        assert(find_grid_index(i, grid_side) == -100);

        pos = std::pair<double,double>(0, 1.7);
        set_pos(i,pos);
        assert(find_grid_index(i, grid_side) == -100);
    }

    //An individual can increase its energy by eating food
    //It also depletes the food
    {
        individual i;
        double init_en = i.get_energy();
        double init_food = 0.1;
        env_grid_cell c(0,init_food);
        feed(i,c);
        assert(i.get_energy() > init_en);
        assert(init_food - i.get_param().get_uptake_rate() < 0.00001);
        assert((i.get_energy() + i.get_param().get_uptake_rate() + c.get_food()) -
               (init_en + init_food) > 0.000001);

    }
    //if there is less food than a individual
    //can normally eat, it eats what is available
    {
        individual i;
        double init_en = i.get_energy();
        double init_food = 0.01;
        env_grid_cell c(0,init_food);
        assert(i.get_param().get_uptake_rate() > c.get_food());

        feed(i, c);
        assert(c.get_food() < 0.000001 && c.get_food() > -0.000001);
        assert(i.get_energy() - init_food - init_en < 0.000001);
    }


    //Active individuals lose energy due to their metabolism
    {
        individual i(0,0,0,1);
        auto init_en = i.get_energy();
        active_metabolism(i);
        auto en_after = i.get_energy();
        assert(init_en > en_after);
    }

    //Sporulating individuals lose energy due to their metabolism
    {
        individual i(0,0,0,1);
        auto init_en = i.get_energy();
        i.set_phen(phenotype::sporulating);
        sporulating_metabolism(i);
        auto en_after = i.get_energy();
        assert(init_en > en_after);
    }

    //An active individual's energy cannot go below 0
    {
        individual i{ind_param{},0,0,0};
        assert(i.get_energy() - i.get_param().get_metabolic_rate() < 0);
        active_metabolism(i);
        assert(i.get_energy() < 0.000001
               && i.get_energy() > -0.0000001);
    }

    //A sporulating individual's energy cannot go below 0
    {
        individual i(ind_param{},0,0,0);
        i.set_phen(phenotype::sporulating);
        assert(i.get_energy() - i.get_param().get_spor_metabolic_rate() < 0);
        sporulating_metabolism(i);
        assert(i.get_energy() < 0.000001
               && i.get_energy() > -0.0000001);
    }

    //An individual can release metabolite in a grid_cell
    {
        individual i;
        env_grid_cell c;
        auto init_metab = c.get_metab();
        secretes_metab(i,c);
        assert(c.get_metab() > init_metab);
    }


    //An individual is initialized with a individual_type(living, sporulating, spore)
    //By default it is initialized with  individual_type::living
    {
        individual i;
        assert(to_str(i.get_phen()) == "living");
        individual i2{ind_param{},0,0,0,phenotype::spore};
        assert(to_str(i2.get_phen()) == "spore");
        individual i3(ind_param{},0,0,0,phenotype::sporulating);
        assert(to_str(i3.get_phen()) == "sporulating");
    }

    //An individuals type can be changed after initialization
    {
        individual i;
        assert(to_str(i.get_phen()) == "living");
        i.set_phen(phenotype::spore);
        assert(to_str(i.get_phen()) == "spore");
        assert(to_str(i.get_phen()) != "living");
    }

    //Individuals are initialized with a sporulating timer
    //By default is 0
    {
        individual i;
        assert(i.get_spo_timer() == 0);
    }

    //The timer can be increased by one
    {
        individual i;
        auto timer_init = i.get_spo_timer();
        assert( timer_init == 0);
        i.tick_spo_timer();
        assert(i.get_spo_timer() == timer_init + 1);
    }

    //The sporulation timer can be reset to 0
    {
        individual ind;
        auto timer_value = 42;
        for (int i = 0; i != timer_value; i++){ind.tick_spo_timer(); }
        assert(ind.get_spo_timer() == timer_value);
        ind.reset_spo_timer();
        assert(ind.get_spo_timer() == 0);
    }

    //An individual is initialized with a sporulation time
    //by default 5, the value will be constant
    {
        individual i;
        assert(i.get_param().get_transformation_time() == 5);
        auto transformation_time = 42;
        individual i2(ind_param{0.1, 0.1, 0.1, 0.1,
                                0.01, 0.01, 0.5, 0.5,
                                0.07,0.5,0.5,0.07,
                                transformation_time},
                      0,0,0,
                      phenotype::sporulating);
        assert(i2.get_param().get_transformation_time() == transformation_time);
    }

    //When the timer_for_sporulation reaches the transformation time
    //The sporulating individual will turn into a spore
    {
        individual i;
        i.set_phen(phenotype::sporulating);
        i.set_energy(i.get_param().get_metabolic_rate() * i.get_param().get_transformation_time());
        for(int j = 0; j != i.get_param().get_transformation_time(); j++)
        {
            assert(i.get_phen() == phenotype::sporulating);
            sporulation(i);
        }
        assert(i.get_phen() == phenotype::spore);
    }

    //Sporulating individuals that revert back to living
    //get their timer reset
    {
        individual i;
        i.set_phen(phenotype::sporulating);
        int time = 42;
        for(int j = 0; j != time; j++)
        {
            i.tick_spo_timer();
        }
        assert(i.get_spo_timer() == time);
        reverts(i);
        assert(i.get_spo_timer() == 0 && i.get_phen() == phenotype::active);
    }

    //It is possible to determine if an individual is sporulating
    {
        individual i;
        assert(!is_sporulating(i));
        i.set_phen(phenotype::sporulating);
        assert(is_sporulating(i));
    }

    //It is possible to determine if an individual is living
    {
        individual i;
        assert(is_active(i));
        i.set_phen(phenotype::sporulating);
        assert(!is_active(i));
    }

    //It is possible to determine if an individual is a spore
    {
        individual i;
        assert(!is_spore(i));
        i.set_phen(phenotype::spore);
        assert(is_spore(i));
    }

    //An individual can start sporulating
    {
        individual i;
        assert(i.get_phen() != phenotype::sporulating);
        assert(i.get_phen() == phenotype::active);
        starts_sporulation(i);
        assert(i.get_phen() == phenotype::sporulating);
        assert(i.get_phen() != phenotype::active);

    }
    //An individual with 0 energy is signaled to be destroyed
    {
        individual i{0,0,0,0};
        assert(i.get_energy() < 0.00001 && i.get_energy() > -0.00001);
        assert(is_dead(i));
    }

    //After feeding and metabolism if energy is 0 individual is destroyed
    {
        individual i{0,0,0,0};
        assert(is_dead(i));//individual is initialized with 0 energy
        //Let's create an env grid cell that will feed the individual
        //the exact quantity it will lose with metabolism
        env_grid_cell c(0,i.get_param().get_metabolic_rate());
        feed(i,c);
        assert(!is_dead(i));
        active_metabolism(i);
        assert(is_dead(i));
    }

    //Spores do not die even if their energy is 0
    {
        individual i{0,0,0,0};
        assert(i.get_phen() != phenotype::spore);
        assert(is_dead(i));
        i.set_phen(phenotype::spore);
        assert(!is_dead(i));
    }

    //The fitness of an individual is determined by its phenotype
    //Spores dispersal probability == Living dispersal probability * spore_advantage
    {
        individual i;
        double base_disp_prob = 0.1;
        double spore_advantage = 10;
        assert(get_fitness(i, base_disp_prob, spore_advantage) - base_disp_prob < 0.0001 &&
               get_fitness(i, base_disp_prob, spore_advantage) - base_disp_prob > -0.0001);
        individual i2;
        i2.set_phen(phenotype::spore);
        assert(get_fitness(i2, base_disp_prob, spore_advantage) -
               get_fitness(i,base_disp_prob, spore_advantage) > 0.0001 ||
               get_fitness(i2, base_disp_prob, spore_advantage) -
               get_fitness(i,base_disp_prob, spore_advantage)  < -0.0001);
    }
    //An individual is initialized with a GRN
    //The get_input_nodes()..... etc part is irrelevant
    //Just get the function to run
    {
        individual i;
        assert(i.get_grn().get_input_nodes().size() == 3);
    }

    //An individual can sense food, metabolite and energy amounts
    //and use them as inputs to its GRN
    {
        individual i{0,0,0,0};
        env_grid_cell g;
        sense(i, g);
        for(const auto& input : i.get_grn().get_input_nodes())
            assert(input < 0.0001 && input > -0.0001);
    }

    //An individual can elaborate internal and external signals
    {
        double food_amount = 3.14;
        double metabolite_amount = 3.14;
        double energy_amount = 3.14;
        individual i(0,0,0,energy_amount);
        env_grid_cell c(metabolite_amount, food_amount);
        //let's set all the weights of the network to 1
        //in this case we expect that the outputs will be one
        i.get_grn().set_all_I2H(1);
        i.get_grn().set_all_H2O(1);
        //Since the network reacts to input of timestpe t
        //at timestep t+1 I will run responds(i) 2 times
        //so we can read the output
        responds(i, c);
        responds(i, c);
        for(int j = 0; j != static_cast<int>(i.get_grn().get_output_nodes().size()); j++)
        {
            assert(i.get_grn().get_output_nodes()[static_cast<unsigned int>(j)] == 1);
        }
    }
    //An individual determines its type according to its GRN outputs
    //if output[0] == false then it is sporulating
    //if output[0] == true then it is living
    {
        individual i;
        assert(is_active(i));
        i.get_grn().set_all_out_nodes(false);
        determine_phen(i);
        assert(!is_active(i));
        assert(is_sporulating(i));

        i.get_grn().set_all_out_nodes(true);
        determine_phen(i);
        assert(is_active(i));
        assert(!is_sporulating(i));
    }
    // An individual can read cues and respond by switching phenotype
    {
        double food_amount = 3.14;
        double metabolite_amount = 3.14;
        double energy_amount = 3.14;
        individual i(0,0,0,energy_amount);
        env_grid_cell c(metabolite_amount, food_amount);
        //let's set all the weights of the network to 1
        //in this case we expect that the outputs will be one
        i.get_grn().set_all_I2H(1);
        i.get_grn().set_all_H2O(1);
        //Since the network reacts to input of timestpe t
        //at timestep t+1 I will run responds(i) 2 times
        //In this first response call, since the hidden nodes are 1
        //The output will be 1 and therefore the individual phenotype
        //will be phenotype::sporulating
        responds(i, c);
        assert(is_active(i));
        //The values of food and metabolite in the grid_cell as well as the energy in the individual
        //are changed so that all inputs are -1, the network should therefore give outputs == 0/false
        //since the outputs node will recieve a signal that is negative(below the treshold)
        //after responds(i) is called 2 more times
        c.set_food(-1);
        c.set_metab(-1);
        i.set_energy(-1);
        //1st responds(i), the individual responds to the initial parameter
        //expected output value == 1
        responds(i, c);
        assert(is_active(i));
        //2nd responds(i), the individual responds to the changed parameter(all 0s)
        //expected output value == 0
        responds(i, c);
        assert(!is_active(i));
        assert(is_sporulating(i));
    }

    //An individual is initialized with a m_is_drawn member
    //0 by default
    {
        individual i;
        assert(!is_drawn(i));
    }

    //And individual can be drawn to fund the new population
    {
        individual i;
        assert(!is_drawn(i));
        draw(i);
        assert(is_drawn(i));
    }

    //A drawn individual can be reset to be drwan again
    {
        individual i;
        draw(i);
        assert(is_drawn(i));
        draw_flag_reset(i);
        assert(!is_drawn(i));
    }

    //An individual has a parameter member ind_param
    //The value -123456789 is irrelevant, used only to use assert
    {
        individual i;
        assert(i.get_param().get_metabolic_rate() > -123456789);
    }

    //An individual has a vector that stores the IDs of its ancestors
    //It is possible to retrieve the line of funders from which an individual derives
    {
        individual i;
        assert(i.get_ancestor().size() >= 0);
    }


    //It is possible to check if two individuals have the same ancestor
    {
        individual i1;
        auto i2 = i1;
        assert(have_same_ancestor(i1, i2));
    }

#endif
}
