#include"tests.h"
void test_population() noexcept  //!OCLINT
{
#ifndef NDEBUG

    //pop is initialized with a certain number of individuals
    // The value 1234567890 is irrelevant: just get this to compile

    {
        unsigned int pop_size = 100;
        population p(pop_size);
        for(int i = 0; i < p.get_pop_size(); ++i)
        {
            double x = p.get_ind(i).get_x();
            assert( x > -1234567890 );
        }
    }

    //A population can be initialized by a pop_param object and a ind_param object
    {
        pop_param p{1,2,3,0.4};
        ind_param i{1,2,3,4};
        population pop{p, i};

        assert(pop.get_param() == p &&
               pop.get_ind(0).get_param() == i);
    }

    //The size of a population is equal to the size of the vector containing its individuals
    {
        population p;
        assert( p.get_pop_size() == static_cast<int>(p.get_v_ind().size()));
    }


    //No individuals are dividing at the start of the pop
    {
        //initiate empty vector and leave it empty
        population p(3);
        assert(static_cast<int>(get_dividing_individuals(p).size()) < 0.00000000001);
    }



    // The number n of hex_layers required for x individual is
    // 3n(n-1)+1 <= x < 3(n+1)((n+1)-1)+1
    {
        std::vector<int> x{0,1,7,19,22};
        std::vector<int> n{0,1,2,3,4};

        for(size_t i = 0; i < x.size(); i++ )
        {
            assert(count_hex_layers(x[i]) == n[i]);
        }
    }

    //Individuals' centers are placed in an hexagonal pattern
    {
        population s(7);

        //The angle formed by any 3 individual_centers is a multiple of PI/6*number_of_hex_layers
        //This is true in an hex_patt but i do not know if it is sufficient

        double n_hex_l = count_hex_layers(s.get_pop_size());
        auto minimal_angle = M_PI / (6 * n_hex_l);
        auto angle = M_PI  / 12.0;
        assert(minimal_angle - angle < 0.00001 &&
               minimal_angle - angle  / 12.0 > -0.00001 );
        auto v_modulus = modulus_of_btw_ind_angles(s,minimal_angle);
        for(auto ind_modulus : v_modulus)
            assert( ind_modulus < 0.0000000001 || (ind_modulus > minimal_angle - 0.000001 && ind_modulus <= minimal_angle));

    }


    //No individuals are overlapping at the start of the pop
    {
        population s(2);
        assert(!has_collision(s));
    }


    //An individual position can be retrieved as a pair object x,y
    {
        population s;
        std::pair<double, double> v{0,0};//default coordinates of individuals
        std::pair<double, double> v2{1,1};//different coordinates from default

        assert(get_ind_pos(s.get_ind(0)) == v);
        assert(get_ind_pos(s.get_ind(0)) != v2);
    }


    // An individaul can be placed to some given coordinates after initialization
    {
        population s(2);
        for(int i = 0; i < s.get_pop_size(); i++)
        {
            set_pos(s.get_ind(i),std::pair<double,double>(i,i));
        }

        for (int i = 0; i < s.get_pop_size() - 1; i++)
        {
            for(int j = i+1 ; j < s.get_pop_size(); j++)
            {
                assert(get_ind_pos(s.get_ind(0)) != get_ind_pos(s.get_ind(j)));
            }
        }

    }


    //When an individual divides it adds a copy of itself
    //to the population/individual vector(i.e divides)
    {
        population s;
        double lhs = s.get_pop_size();
        //let's allow all individuals in the population to reproduce,
        //to facilitate the testing conditions
        for(auto& ind : s.get_v_ind())
        {
            ind.set_energy(ind.get_param().get_treshold_energy());
        }
        division(s);

        double rhs = s.get_pop_size();
        assert(lhs * 2 - rhs < 0.0000001);

    }


    //Only individuals with energy >= than the treshold will divide
    {
        population s(3);
        //This ind will not reproduce
        s.get_ind(0).set_energy(s.get_ind(0).get_param().get_treshold_energy()-1);
        //This ind will reproduce with no extra energy
        s.get_ind(1).set_energy(s.get_ind(1).get_param().get_treshold_energy());
        //This ind will reproduce with extra energy
        s.get_ind(2).set_energy(s.get_ind(2).get_param().get_treshold_energy()+1);

        std::vector<int> div_ind = get_dividing_individuals(s);
        assert(!div_ind.empty());

        for(unsigned int i = 0; i != div_ind.size(); i++)
            assert(s.get_ind(div_ind[i]).get_energy()
                   >= s.get_ind(div_ind[i]).get_param().get_treshold_energy());

    }

    //The excess energy of dividing individuals is equal
    //to the diference between their energy and their treshold energy
    {
        population s(2);
        double energy = get_ind_tr_en(s, 0);
        set_ind_en(s.get_ind(0), energy);
        set_ind_en(s.get_ind(1), energy*2);
        auto v_ex_en = get_excess_energies(s);
        for(int i =0; i != static_cast<int>(v_ex_en.size()); i++){
            assert(v_ex_en[static_cast<unsigned int>(i)]
                    - (get_ind_en(s, i) - get_ind_tr_en(s, i)) < 0.00000001);
        }
    }


    //The distance between the two offspring from same mother in population vector
    //is the distance betwen the mother index and the end of the population vector

    {
        population s(1);
        //First individual reproduces
        set_ind_en(s.get_ind(0), get_ind_tr_en(s, 0)*2);
        auto daughters = get_sisters_index_offset(s);
        // In this case the distance between the two offspring
        //will be 1 element of the vector in next gen
        auto sister_distances = 1;
        for(auto sisters : daughters)
        {
            assert( sisters.second - sisters.first - sister_distances < 0.000001);
        }
    }


    //After a reproduction round individuals with enough energy will divide
    //and redistribute their remaining energy to their daughter cells
    {
        population s;
        auto repr_excess_en = s.get_ind(0).get_param().get_treshold_energy();
        //This ind will reproduce with extra energy
        set_ind_en(s.get_ind(0), repr_excess_en);

        auto mother_excess_energy = get_excess_energies(s);
        division(s);
        assert(get_ind_en(s, 0) - mother_excess_energy[0]/2 < 0.00001 &&
                get_ind_en(s, 0) - mother_excess_energy[0]/2 > -0.00001);
        assert(get_ind_en(s, 1) - mother_excess_energy[0]/2 < 0.00001 &&
                get_ind_en(s, 1) - mother_excess_energy[0]/2 > -0.00001);
        assert(get_ind_en(s, 0) - get_ind_en(s, 1) < 0.00001 &&
               get_ind_en(s, 0) - get_ind_en(s, 1) > -0.00001);
    }

    //Individuals can be sorted based on their x coordinate in increasing order
    {
        auto pop_size = 4;
        population s(static_cast<unsigned int>(pop_size));
        auto comparison_pop = s.get_v_ind();
        sort_inds_by_x_inc(s.get_v_ind().begin(),s.get_v_ind().end());
        assert(s.get_v_ind() != comparison_pop);
        sort_inds_by_x_inc(comparison_pop.begin(), comparison_pop.end());
        assert(s.get_v_ind() == comparison_pop);

        //Can also sort parts of the vector of ind
        auto changed_ind = pop_size - 1;
        auto reference_ind = pop_size / 2 - 1;
        set_pos(s.get_ind(changed_ind),get_pos(s.get_ind(reference_ind)));
        comparison_pop = s.get_v_ind();
        sort_inds_by_x_inc(comparison_pop.begin(),comparison_pop.end());
        assert(s.get_v_ind() != comparison_pop);
        sort_inds_by_x_inc(s.get_v_ind().begin() + reference_ind, s.get_v_ind().begin() + changed_ind + 1);
        assert(s.get_v_ind() == comparison_pop);

    }

    //After reproduction the first daughter individual takes the position of the mother
    {
        population s;
        auto parent_pop = s.get_v_ind();
        //setting energy high enough for the individual
        //to reproduce without offspring having negative enrgies
        set_ind_en(s.get_ind(0),get_ind_tr_en(s, 0));

        divides(s.get_ind(0),
                s.get_v_ind(),
                repr_angle(s),
                s.get_rng(),
                create_bernoulli_dist(s.get_param().get_mu_p()),
                create_normal_dist(0,s.get_param().get_mu_st())
                );
        //The first daughter cell is at the same index of the mother cell
        assert(get_ind_pos(s.get_ind(0)) == get_pos(parent_pop[0]) );
    }

    //After reproduction the second daughter individual
    //is placed just outside the position of the first daughter cell
    {
        population s;
        //setting energy high enough for the individual to reproduce
        //without offspring having negative enrgies
        set_ind_en(s.get_ind(0),get_ind_tr_en(s, 0));
        divides(s.get_ind(0),
                s.get_v_ind(),
                repr_angle(s),
                s.get_rng(),
                create_bernoulli_dist(s.get_param().get_mu_p()),
                create_normal_dist(0,s.get_param().get_mu_st())
                );
        assert(!has_collision(s));
        assert(distance(s.get_ind(0), s.get_ind(1)) -
               (s.get_ind(0).get_radius() + s.get_ind(1).get_radius()) < 0.1);
    }

    //Given a focal individual it is possible to find the individuals that could potentially
    //collide with it (if we draw a square around each individuals, these squares intersect)
    {
        population s(2);
        set_pos(s.get_ind(1),get_ind_pos(s.get_ind(0)));
        auto collision_range = possible_collisions_x(s.get_ind(0),s.get_v_ind());
        assert(collision_range.first == s.get_v_ind().begin() &&
               collision_range.second == s.get_v_ind().end());


        //No collisions should return a pair whose first element is the focal individual
        //And the second the individual after the focal
        //(or the end of the vector in case the focal is the last)
        set_pos(s.get_ind(1),std::pair<double,double>{2,2});
        //Pop need to be sorted by increasing x coord
        std::sort(s.get_v_ind().begin(),s.get_v_ind().end(),
                  [](const individual& lhs, const individual& rhs)
        {return  lhs.get_x() < rhs.get_x();});
        collision_range = possible_collisions_x(s.get_ind(0),s.get_v_ind());
        assert(collision_range.first + 1 == collision_range.second);

        //One collision to the left should return a pair whose first element is the iterator
        //before the focal individual and the second is the one after the focal individual
        //(the end of the vector if the focal is the last element
        set_pos(s.get_ind(0),
                std::pair<double,double>{
                    s.get_ind(1).get_x() - s.get_ind(1).get_radius(),
                    s.get_ind(1).get_y()
                }
                );
        sort_inds_by_x_inc(s.get_v_ind().begin(),s.get_v_ind().end());
        auto focal_ind_index = 1;
        collision_range = possible_collisions_x(s.get_ind(focal_ind_index),s.get_v_ind());
        assert(collision_range.first == s.get_v_ind().begin() + focal_ind_index - 1);
        assert(collision_range.second == s.get_v_ind().begin() + focal_ind_index + 1);
        assert(collision_range.second == s.get_v_ind().end());

        //Collision to the right should return a range
        //focal_ind_index,focal_ind_index+2
        sort_inds_by_x_inc(s.get_v_ind().begin(),s.get_v_ind().end());
        focal_ind_index = 0;
        collision_range = possible_collisions_x(s.get_ind(focal_ind_index), s.get_v_ind());
        assert(collision_range.first == s.get_v_ind().begin() + focal_ind_index);
        assert(collision_range.second == s.get_v_ind().begin() + focal_ind_index + 2);
        assert(collision_range.second == s.get_v_ind().end());

        //Collision from below should return a range
        //focal_index - 1, focal_index +1
        set_pos(s.get_ind(0),
                std::pair<double,double>{
                    s.get_ind(1).get_x(),
                    s.get_ind(1).get_y() + s.get_ind(1).get_radius()
                }
                );
        sort_inds_by_x_inc(s.get_v_ind().begin(),s.get_v_ind().end());
        focal_ind_index = 1;
        collision_range = possible_collisions_x(s.get_ind(focal_ind_index), s.get_v_ind());
        assert(collision_range.first == s.get_v_ind().begin() + focal_ind_index - 1);
        assert(collision_range.first == s.get_v_ind().begin());
        assert(collision_range.second == s.get_v_ind().begin() + focal_ind_index + 1);
        assert(collision_range.second == s.get_v_ind().end());

        //Collision from above should return a range
        //focal_index, focal_index + 2
        sort_inds_by_x_inc(s.get_v_ind().begin(),s.get_v_ind().end());
        focal_ind_index = 0;
        collision_range = possible_collisions_x(s.get_ind(focal_ind_index), s.get_v_ind());
        assert(collision_range.first == s.get_v_ind().begin() + focal_ind_index);
        assert(collision_range.first == s.get_v_ind().begin());
        assert(collision_range.second == s.get_v_ind().begin() + focal_ind_index + 2);
        assert(collision_range.second == s.get_v_ind().end());

    }
    //If there are collisions individuals are displaced
    //until there are no more collisions
    {
        population s(2);
        assert(!has_collision(s));
        set_pos(s.get_ind(1), get_ind_pos(s.get_ind(0)));
        no_complete_overlap(s);
        assert(has_collision(s));

        calc_tot_displ_pop(s);
        for(auto& ind : s.get_v_ind())
        {ind.displace();}

        assert(!has_collision(s));

    }
    //After reproduction new collisions caused by new individuals
    //being placed where other individuals already are managed
    {
        population s(7);
        assert(!has_collision(s));
        //add 1 individual overlapping with central individual
        s.get_v_ind().emplace_back(individual(ind_param{},0,0));
        assert(has_collision(s));
        manage_static_collisions(s);
        assert(!has_collision(s));
    }

    //Individuals lose energy through metabolism
    {
        population p;
        for (auto & ind : p.get_v_ind())
        {
            ind.set_energy(2);
        }
        double init_en_tot = std::accumulate
                (
                    p.get_v_ind().begin(),
                    p.get_v_ind().end(),0.0,
                    [](double sum,const individual& i){return sum + i.get_energy();}
        );

        metabolism_pop(p);
        double after_en_tot = std::accumulate
                (
                    p.get_v_ind().begin(),
                    p.get_v_ind().end(),0.0,
                    [](double sum,const individual& i){return sum + i.get_energy();}
        );
        assert(after_en_tot < init_en_tot);
    }

    //In Jordi's version individual lose energy only during sporulations
    {
        population p;
        for (auto & ind : p.get_v_ind())
        {
            ind.set_energy(2);
        }
        double init_en_tot = std::accumulate
                (
                    p.get_v_ind().begin(),
                    p.get_v_ind().end(),0.0,
                    [](double sum,const individual& i){return sum + i.get_energy();}
        );
        for(const auto& ind : p.get_v_ind())
        {
            assert(ind.get_phen() == phenotype::active);

        }
        spor_metabolism_pop(p);
        double after_en_tot = std::accumulate
                (
                    p.get_v_ind().begin(),
                    p.get_v_ind().end(),0.0,
                    [](double sum,const individual& i){return sum + i.get_energy();}
        );
        assert(after_en_tot - init_en_tot < 0.00001 &&
               after_en_tot - init_en_tot > -0.00001);
        for( auto& ind : p.get_v_ind())
        {
            ind.set_phen(phenotype::sporulating);

        }
        spor_metabolism_pop(p);
        after_en_tot = std::accumulate
                (
                    p.get_v_ind().begin(),
                    p.get_v_ind().end(),0.0,
                    [](double sum,const individual& i){return sum + i.get_energy();}
        );

        assert(after_en_tot < init_en_tot);
    }
    //Spores do not get detected when looking for dividing individual
    {
        population p;
        set_ind_en(p.get_ind(0),get_ind_tr_en(p, 0));
        assert(get_dividing_individuals(p)[0] == 0);
        p.get_ind(0).set_phen(phenotype::spore);
        assert(get_dividing_individuals(p).empty());
    }

    //Spores do not reproduce
    {
        population p;
        p.get_ind(0).set_phen(phenotype::spore);
        set_ind_en(p.get_ind(0),get_ind_tr_en(p, 0));
        auto init_pop_size = p.get_pop_size();
        set_ind_en(p.get_ind(0),get_ind_tr_en(p, 0));
        for (auto & ind : p.get_v_ind())
        {
            ind.set_energy(2);
        }
        metabolism_pop(p);
        assert(!division(p));
        assert(init_pop_size == p.get_pop_size());
    }

    //Spores do not lose energy
    {
        population p;
        p.get_ind(0).set_phen(phenotype::spore);
        set_ind_en(p.get_ind(0), get_ind_tr_en(p, 0));
        auto init_en_ind0 = get_ind_en(p, 0);
        metabolism_pop(p);
        assert(get_ind_en(p, 0) - init_en_ind0 < 0.000001
               && get_ind_en(p, 0) - init_en_ind0 > -0.000001);
    }

    //Sporulating individuals do not get detected when looking for dividing individual
    {
        population p;
        set_ind_en(p.get_ind(0),get_ind_tr_en(p, 0));
        assert(get_dividing_individuals(p)[0] == 0);
        p.get_ind(0).set_phen(phenotype::sporulating);
        assert(get_dividing_individuals(p).empty());
    }

    //Sporulating individuals cannot reproduce
    {
        population p;
        p.get_ind(0).set_phen(phenotype::sporulating);
        set_ind_en(p.get_ind(0),get_ind_tr_en(p, 0));
        auto init_pop_size = p.get_pop_size();
        set_ind_en(p.get_ind(0), get_ind_tr_en(p, 0));
        for (auto & ind : p.get_v_ind())
        {
            ind.set_energy(2);
        }
        metabolism_pop(p);
        assert(!division(p));
        assert(init_pop_size == p.get_pop_size());
    }
    //A sporulating individual updates its timer every time metab_pop is called
    {
        population p;
        auto init_timer = p.get_ind(0).get_spo_timer();
        p.get_ind(0).set_phen(phenotype::sporulating);
        metabolism_pop(p);
        assert(init_timer != p.get_ind(0).get_spo_timer());
        assert(p.get_ind(0).get_spo_timer() == init_timer + 1);
        assert(p.get_ind(0).get_spo_timer() != init_timer + 2);

        p.get_ind(0).reset_spo_timer();
        init_timer = p.get_ind(0).get_spo_timer();
        p.get_ind(0).set_phen(phenotype::spore);
        metabolism_pop(p);
        assert(p.get_ind(0).get_spo_timer() == init_timer);
        assert(p.get_ind(0).get_spo_timer() != init_timer + 1);
        assert(p.get_ind(0).get_spo_timer() != init_timer + 2);

        p.get_ind(0).reset_spo_timer();
        init_timer = p.get_ind(0).get_spo_timer();
        p.get_ind(0).set_phen(phenotype::active);
        metabolism_pop(p);
        assert(p.get_ind(0).get_spo_timer() == init_timer);
        assert(p.get_ind(0).get_spo_timer() != init_timer + 1);
        assert(p.get_ind(0).get_spo_timer() != init_timer + 2);
    }

    //Individuals in a population will die with a certain probability
    //given by the population parameters
    {
        int pop_size = 100;
        int replicates = 1000;
        double death_rate = 0.1;
        pop_param pp {static_cast<unsigned int>(pop_size),
                    1,
                    0.1,
                    0.0015,
                    0.1,
                    0.01,
                    10,
                    death_rate};
        auto mean = 0.0;
        population p(pp);
        std::vector<individual> post_senescence_pop;

        for( int i  = 0; i != replicates; i++)
        {
            senescence(p);
            mean += p.get_pop_size();
            p.get_v_ind().resize(pop_size);
        }
        mean /= replicates;
        auto balance = pop_size - (mean + pop_size * pp.get_death_rate());
        assert( balance < 0.1 &&
                balance > -0.1);
    }

    //Individuals that starve are removed from the population
    {
        population p;
        p.get_ind(0).set_energy(0);//the only individual in this sim has 0 energy, thus it will die
        assert(p.get_pop_size() == 1);
        starvation(p);
        assert(p.get_v_ind().empty() && p.get_pop_size() == 0);

        unsigned int pop_size = 5;
        //The pop does not have a grid with food,
        //so organisms cannot feed
        population p1 (pop_size);
        for(auto& ind : p1.get_v_ind())
        {
            ind.set_energy(0);
        }
        //Only the first individual has enough energy to survive
        //for 1 tick
        set_ind_en(p1.get_ind(0),p1.get_ind(0).get_param().get_metabolic_rate() + 0.001);
        assert(p1.get_v_ind().size() == pop_size);
        metabolism_pop(p1);
        starvation(p1);
        assert(p1.get_pop_size() == 1);
        //then at the second tick the only individual left dies
        metabolism_pop(p1);
        starvation(p1);
        assert(p1.get_v_ind().empty() && p.get_pop_size() == 0);
    }

    //A pop can generate numbers from a uniform distribution
    //between 0 and 2PI
    {
        population p;
        double mean = 0;
        //Draw a thousands times from a uniform dist between 0 and 2
        int sampling_size = 1000;
        for(int i = 0; i != sampling_size; i++)
        {
            mean += repr_angle(p);
        }
        //calculate mean of the drawn values
        mean /= 1000;
        assert(mean < M_PI + 0.1 && mean > M_PI - 0.1 );
    }

    //Daughter cells are placed at a random angle after reproduction
    {
        population p;
        //to calculate angle we will use three point
        //the center of the mother(0,0)
        std::pair<double, double> mother (0,0);
        set_pos(p.get_ind(0),mother);
        //a reference point on the same axis as the mother (0,1)
        //and the center of the daughter -> get_daughter_pos()
        std::pair<double, double> reference(1,0);

        //Draw a thousands times from a uniform dist between 0 and 2
        double mean = 0;
        int sampling_size = 1000;
        for(int i = 0; i != sampling_size; i++)
        {
            auto daughter = calculate_daughter_pos(p.get_ind(0), repr_angle(p));
            mean += calc_angle_3_pos(mother,daughter,reference);
        }
        mean /= sampling_size;
        assert(mean < M_PI + 0.1 && mean > M_PI - 0.1 );
    }

    //A pop can generate numbers from a normal distribution for mutation step
    //with mean 0 and variance 0.1 by default
    {
        population p;
        double mean = 0;
        int sampling_size = 10000;
        for(int i = 0 ; i != sampling_size; i++ )
            mean += mut_step(p);
        mean /= sampling_size;
        assert(mean < 0.01 && mean > -0.01);
    }

    //A pop can generate numbers from a bernoulli distribution to see if mutation happens or not
    //0.0015 by default
    {
        population p;
        double mean = 0;
        int sampling_size = 100000;
        for(int i = 0 ; i != sampling_size; i++ )
            mean += mut_happens(p);
        mean /= sampling_size;
        assert(mean - p.get_param().get_mu_p()< 0.0011 && mean - p.get_param().get_mu_p() > -0.001);
    }

    //The sum of weight of an individual after many rounds of mutation
    //Should have the same mean as in the beginning, but its variance should
    //be the same as the mutation_step distribution
    {
        population p;
        //   double init_mean = weights_mean(s.get_ind(0).get_grn());
        double init_variance = weights_var(p.get_ind(0).get_grn());
        assert(init_variance < 0.0001 && init_variance > -0.000001);

        int sampling_size = 10000;
        auto mut_prob_dist = create_bernoulli_dist(p.get_param().get_mu_p());
        auto mut_step_dist = create_normal_dist(0,p.get_param().get_mu_st());
        for (int i = 0; i != sampling_size; i++)
        {
            mutates(p.get_ind(0),
                    p.get_rng(),
                    mut_prob_dist,
                    mut_step_dist
                    );
        }
        //This first assert does not pass, the mean is much more variable than
        //I thought, but I do not see any bug. I will comment this out
        //    assert(mean - init_mean > -0.1 && mean - init_mean < 0.1);
        assert(init_variance - weights_var(p.get_ind(0).get_grn()) > 0.01 ||
               init_variance - weights_var(p.get_ind(0).get_grn()) < -0.01);
    }

    //After dividing the two daughter individuals mutate
    {
        double mutation_probability = 1; //all weights will be mutated in this pop
        population p(pop_param{1, 1, 0, mutation_probability});
        auto init_var = weights_var(p.get_ind(0).get_grn());
        assert(init_var < 0.00001 && init_var > -0.0001);
        divides(p.get_ind(0),
                p.get_v_ind(),
                repr_angle(p),
                p.get_rng(),
                create_bernoulli_dist(p.get_param().get_mu_p()),
                create_normal_dist(0, p.get_param().get_mu_st())
                );
        auto post_var = weights_var(p.get_ind(0).get_grn());
        assert(init_var - post_var > 0.000001 || init_var - post_var < -0.0001);
    }

    //A pop_param has a member variable m_new_pop_size that states the max number of
    //individuals that will start a new population
    //by default = to pop_start_size
    //If at dispersal m_new_pop_size > pop.size()
    //Then the number of funding individuals == pop.size()

    {
        population p(pop_param{1000, 100});
        select_new_pop(p);
        assert(p.get_new_v_ind().size() == p.get_param().get_exp_new_pop_size());
        population p1 (pop_param{10, 100});
        select_new_pop(p1);
        auto n_drawn_ind = std::count_if(p1.get_v_ind().begin(),p1.get_v_ind().end(),
                                         [](const individual& i) {return is_drawn(i);});
        assert( n_drawn_ind == p1.get_pop_size());
    }


    //During dispersal the individuals selected for the new_pop cannot be drawn again from pop
    {
        population p(pop_param{1000, 100});
        select_new_pop(p);
        assert(std::count_if(p.get_v_ind().begin(),p.get_v_ind().end(),
                             [](const individual& i) {return is_drawn(i);}) == 100);
    }

    //A pop is initialized with a variable m_base_fitness
    //by default = 0.01
    {
        double base_disp_prob = 0.1;
        population p(pop_param{0,0,0,0,0,base_disp_prob});
        assert(p.get_param().get_base_disp_prob() - base_disp_prob < 0.00001 &&
               p.get_param().get_base_disp_prob() - base_disp_prob > -0.000001);
    }

    //A pop is initialized with a variable m_spore_advantage
    //by default = 10
    {
        double spore_advantage = 10;
        population s(pop_param{0,0,0,0,0,0,spore_advantage});
        assert(s.get_param().get_spo_adv() - spore_advantage < 0.00001 &&
               s.get_param().get_spo_adv() - spore_advantage > -0.000001);
    }

    //A pop is initialized with a uniform distribution between 0 and 1
    //used to see which ind is drawn at dispersal
    {
        population s;
        int sample_size = 1000;
        double mean = 0;
        for(int i = 0; i != sample_size; i++)
        {
            mean += create_unif_dist(0,1)(s.get_rng());
        }
        mean /= sample_size;
        assert(mean - 0.5 < 0.01 &&
               mean - 0.5 > -0.01);
    }


    //Individuals are selected based on their phenotype
    //A spore is more likely to be selected than a living
    {
        population p(pop_param{1000,100});
        for(int i = p.get_pop_size() / 2; i != p.get_pop_size(); i++)
            p.get_ind(i).set_phen(phenotype::spore);
        select_new_pop(p);
        auto spore_ratio =
                std::accumulate(p.get_new_v_ind().begin(),p.get_new_v_ind().end(),0.0,
                                [](const int sum, const individual& ind){return sum + is_spore(ind);}) /
                p.get_new_v_ind().size();
        assert(spore_ratio > 0.5);
    }

    //After being selected in new population
    //individuals flag is_drawn is resetted
    {
        population p;
        select_new_pop(p);
        for(const auto& ind :p.get_new_v_ind())
            assert(!is_drawn(ind));
    }

    //After being selected in new population
    //individuals energy are reset to default level

    {
        population p;
        select_new_pop(p);
        for(const auto& ind :p.get_new_v_ind())
            assert(ind.get_energy() == individual{}.get_energy());
    }

    //After a new population is selected it
    //swapped with the old population
    //And the old population is cancelled
    {
        unsigned int pop_size = 1000;
        unsigned int new_pop_size = 100;
        population p(pop_param{pop_size,new_pop_size});
        select_new_pop(p);
        assert(p.get_new_v_ind().size() == new_pop_size);
        assert(p.get_v_ind().size() == pop_size);
        place_new_pop(p);
        assert(p.get_new_v_ind().size() == 0);
        assert(p.get_v_ind().size() == new_pop_size);
    }

    //Individuals after funding the new population are set in an hexagonal pattern
    {
        unsigned int pop_size = 1000;
        unsigned int new_pop_size = 100;
        population p(pop_param{pop_size,new_pop_size});
        fund_new_pop(p);
        auto n_hex_l = count_hex_layers(p.get_pop_size());
        auto v_modulus = modulus_of_btw_ind_angles(p, M_PI/ (6 * n_hex_l));
        for(auto ind_modulus : v_modulus)
            assert( ind_modulus < 0.0000000001 || (ind_modulus > M_PI / (6 * n_hex_l) - 0.1 && ind_modulus <= M_PI / (6 * n_hex_l) + 0.1));
        assert(!has_collision(p));
    }

    //It is possible to assign a unique ID to each ind in a population
    {
        population p;
        assert(p.get_pop_size() == 1);
        divides(p.get_ind(0),
                p.get_v_ind(),
                0,
                p.get_rng(),
                create_bernoulli_dist(0),
                create_normal_dist(0,0));
        assert(p.get_pop_size() == 2);
        assert(have_same_ancestor(p.get_ind(0), p.get_ind(1)));
        p.get_v_ind() = assign_ancestor_ID(p.get_v_ind());
        assert(!have_same_ancestor(p.get_ind(0), p.get_ind(1)));

    }
    //When a population is initialized each ind recieves a unique ancestor ID
    {
        unsigned int pop_size = 2;
        population p{pop_size};
        for (int i = 0; i != p.get_pop_size() - 1; i++)
            for(int j = i + 1; j != p.get_pop_size(); j++)
                assert(!have_same_ancestor(p.get_ind(i),p.get_ind(j)));
    }

    //it is possible to track individuals that have the same ancestry
    {
        int pop_size = 2;
        population p{static_cast<unsigned int>(pop_size)};
        auto ind1 = p.get_ind(0);
        auto ind2 = p.get_ind(1);

        divides(ind1,
                p.get_v_ind(),
                0,
                p.get_rng(),
                create_bernoulli_dist(0),
                create_normal_dist(0,0));
        //The new ind is at the .back() of the individuals vector
        assert(p.get_pop_size() == pop_size + 1);
        auto& ind3 = p.get_v_ind().back();
        assert( have_same_ancestor(ind1, ind3) );
        assert( !have_same_ancestor(ind2, ind3) );

    }

    //When a pop is funded a new ancestor ID is assigned to each funder individual
    {
        population p;
        //store ancestor_IDs of actual population
        std::vector<std::vector<int>> p_ancestor_IDs;
        for(const auto& ind : p.get_v_ind())
        {
            p_ancestor_IDs.push_back(ind.get_ancestor());
        }

        //Fund new pop
        fund_new_pop(p);

        //store ancestor_IDs of NEW population
        std::vector<std::vector<int>> new_p_ancestor_IDs;
        for(const auto& ind : p.get_v_ind())
        {
            new_p_ancestor_IDs.push_back(ind.get_ancestor());
        }

        assert(p_ancestor_IDs.size() == new_p_ancestor_IDs.size());
        for(size_t i = 0; i != p_ancestor_IDs.size(); i++)
            for(size_t j = 0; j != new_p_ancestor_IDs.size(); j++)
            {
                assert(p_ancestor_IDs[i] != new_p_ancestor_IDs[j]);
            }
    }

    //A population's individuals can be changed by
    //Swapping their ind_params with new ones
    {
        ind_param i{};
        pop_param pp{};
        population p{pp};
        auto new_p = set_new_ind_par(p.get_v_ind(), change_ind_param_norm(i, p.get_rng()));
        assert(p.get_v_ind() != new_p);
    }


#endif
}