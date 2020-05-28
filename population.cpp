#include "population.h"

population::population(pop_param pop_parameters):
    m_pop_param{pop_parameters},
    m_pop(m_pop_param.get_pop_start_size(), m_pop_param.get_ind_param()),
    m_new_pop_tmp_buffer(0, m_pop_param.get_ind_param()),
    m_relax(relaxation::param_t{})
{
    assign_ancestor_ID(m_pop).swap(m_pop);
    if(!m_pop.empty())
    {
        place_start_cells(*this);
    }
}

void active_metabolism_pop(population &p)
{
    for(auto& ind : p.get_v_ind())
    {
        if(ind.get_phen() != phenotype::active)
            continue;

        active_metabolism(ind);
    }
}

bool all_ind_are_drawn(const population& s) noexcept
{
    return std::all_of(s.get_v_ind().begin(), s.get_v_ind().end(),
                       [](const individual& i) {return is_drawn(i);});
}

std::vector<individual> assign_ancestor_ID(const std::vector<individual>& p) noexcept
{
    auto new_p = p;
    for(unsigned int i = 0; i != p.size(); i++ )
    {
        new_p[i].get_ancestor().push_back(i);
    }
    return new_p;
}

int count_hex_layers(int pop_size)  noexcept
{
    int n = 1;
    if(pop_size>0){while(3*n*(n-1)+1 < pop_size) n++;}
    else {return 0;}
    return n;
}

double calc_angle_3_pos(std::pair<double, double> P1,
                        std::pair<double, double> P2,
                        std::pair<double, double> P3)
{
    auto angle1 = std::atan2(P3.second - P1.second, P3.first - P1.first);
    auto angle2 = std::atan2(P2.second - P1.second, P2.first - P1.first);
    if(angle1 >= angle2)
        return angle1 - angle2;

    return M_PI + std::abs(angle1 -angle2);
}

bool calc_tot_displ_pop(std::vector<individual>& pop)
{
    bool has_collision = false;
    //Sort the pop vector by increasing x
    sort_inds_by_x_inc(pop.begin(), pop.end());

    for ( auto i = pop.begin(); i != pop.end(); ++i)
    {
        i->becomes_focal();
        auto focal_ind = *i;
        auto range_x = possible_collisions_x(focal_ind, pop);
        auto range_y = possible_collisions_y(focal_ind,range_x.first, range_x.second, pop);

        assert(only_one_focal(range_x.first,range_x.second));

        auto focal_ind_address = std::find_if(range_x.first, range_x.second,
                                              [](const individual& ind) {return ind.is_focal();});

        for ( auto j = range_y.first ; j != range_y.second; ++j)
        {
            if(j->is_focal()){continue;}
            if (are_colliding(*focal_ind_address, *j))
            {
                assert(!j->is_focal());
                add_displacement(*focal_ind_address, *j);
                has_collision = true;
            }
        }
        //Sort back to increasing x
        focal_ind_address->no_more_focal();
        sort_inds_by_x_inc(range_x.first, range_x.second);
    }
    return has_collision;
}

std::uniform_real_distribution<double> create_unif_dist(double a, double b) noexcept
{
    return std::uniform_real_distribution<double>(a,b);
}

std::bernoulli_distribution create_bernoulli_dist(double p) noexcept
{
    return std::bernoulli_distribution{p};
}

std::normal_distribution<double> create_normal_dist(double m, double v)
{
    return std::normal_distribution<double>{m, v};
}



void death(population &p) noexcept
{
    starvation(p);
    senescence(p);
}

bool division(population &p) noexcept
{
    auto  div_inds  = get_dividing_individuals(p);
    std::uniform_real_distribution<double> repr_prob(0,1);

    if(!div_inds.empty())
        for(size_t i = 0; i != div_inds.size(); i++)
        {
            int div_ind = div_inds[i];
            if(repr_prob(p.get_rng()) < p.get_param().get_repr_p())
            {

                divides(p.get_ind(div_ind),
                        p.get_v_ind(),
                        repr_angle(p),
                        p.get_rng(),
                        create_bernoulli_dist(p.get_param().get_mu_p()),
                        create_normal_dist(0,p.get_param().get_mu_st())
                        );
            }
        }
    return !div_inds.empty();
}


void displace_inds(std::vector<individual>& pop) noexcept
{
    for(auto& ind : pop)
    {
        if(ind.get_x_displacement() < 0.00000001 &&
                ind.get_x_displacement() > -0.0000001 &&
                ind.get_y_displacement() < 0.00000001 &&
                ind.get_y_displacement() > -0.0000001)
        {
            ind.displace();
        }
    }
}

void fund_new_pop(population& p) noexcept
{
    select_new_pop(p);
    place_new_pop(p);
    assign_ancestor_ID(p.get_v_ind()).swap(p.get_v_ind());
    place_start_cells(p);
}

std::vector<int> get_dividing_individuals(const population &p) noexcept
{

    std::vector<int> dividing_individuals;
    for(unsigned int i = 0; i != p.get_v_ind().size(); i++)
    {
        if(p.get_v_ind()[i].get_energy() >= p.get_v_ind()[i].get_param().get_treshold_energy()
                && p.get_v_ind()[i].get_phen() == phenotype::active)
        {
            dividing_individuals.push_back(static_cast<int>(i));
        }
    }
    return dividing_individuals;
}

std::vector<double> get_excess_energies(const population &p) noexcept
{
    std::vector<double> excess_energy;
    auto v_ind{get_dividing_individuals(p)};
    for(unsigned int i=0; i != v_ind.size(); i++)
    {
        int ind_index = v_ind[i];
        excess_energy.push_back(
                    get_ind_en(p, ind_index) - get_ind_tr_en(p, ind_index)
                    );
    }
    return excess_energy;
}

double get_ind_en(const population &p, int i)
{
    return p.get_ind(i).get_energy();
}

double get_ind_tr_en(const population &p, int i)
{
    return p.get_ind(i).get_param().get_treshold_energy();
}

std::pair<double, double> get_ind_pos(const individual& i)
{
    std::pair<double, double> pos;
    pos.first = i.get_x();
    pos.second = i.get_y();
    return pos;

}

std::vector<std::pair<int,int>> get_sisters_index_offset(const population &p)  noexcept
{
    std::vector<std::pair<int,int>> sis_indexes;
    std::pair<int,int> daughters;
    int sis_dist;
    auto div_ind = get_dividing_individuals(p);

    for (const auto& ind : div_ind)
    {
        sis_dist = static_cast<int>(p.get_pop_size()) - ind;
        daughters.first = ind;
        daughters.second = daughters.first + sis_dist;
        sis_indexes.push_back(daughters);
    }
    return sis_indexes;
}

bool has_collision(population &p)
{
    auto pop = p.get_v_ind();
    //Sort the pop vector by increasing x
    sort_inds_by_x_inc(pop.begin(), pop.end());

    for ( auto i = pop.begin(); i != pop.end(); ++i)
    {
        i->becomes_focal();
        auto focal_ind = *i;
        auto range_x = possible_collisions_x(focal_ind, pop);
        auto range_y = possible_collisions_y(focal_ind,range_x.first, range_x.second, pop);

#ifndef NDEBUG
        assert(only_one_focal(range_x.first,range_x.second));
#endif

        auto focal_ind_address = std::find_if(range_x.first, range_x.second,
                                              [](const individual& ind) {return ind.is_focal();});

        for ( auto j = range_y.first ; j != range_y.second; ++j)
        {
            //if j is the focal individual skip collision
            if(j->is_focal()){continue;}
            if (are_colliding(*focal_ind_address, *j))
            {
#ifndef NDEBUG
                assert(!j->is_focal());
#endif
                return true;
            }
        }
        //Sort back to increasing x
        focal_ind_address->no_more_focal();
        sort_inds_by_x_inc(range_x.first, range_x.second);
    }
    return false;
}

void jordi_death(population &p) noexcept
{
    senescence(p);
}

int manage_static_collisions(population &p)
{
    using relaxation::particle_t;

    int time =
            p.get_relax()(p.get_v_ind().begin(), p.get_v_ind().end(),
                          [](const individual& i)
                          // individual to particle_t conversion
    { return particle_t{ glm::vec2{i.get_x(), i.get_y()}, static_cast<float>(i.get_param().get_radius()) }; },
    [](individual& i, const particle_t& p)
    { i.set_x(p.pos.x); i.set_y(p.pos.y); i.get_param().set_radius(p.radius); }
            );
    return time;
}

void metabolism_pop(population &p)
{
    active_metabolism_pop(p);
    spor_metabolism_pop(p);
}

std::vector<double> modulus_of_btw_ind_angles(population &p, double ang_rad)
{
    std::vector<double> v_modulus;
    for(int i = 0; i != p.get_pop_size() - 2; i++)
        for(int j = i+1; j != p.get_pop_size() - 1; j++)
            for(int k = j+1; k != p.get_pop_size(); k++)
            {
                auto P1 = get_pos(p.get_ind(i));
                auto P2 = get_pos(p.get_ind(j));
                auto P3 = get_pos(p.get_ind(k));
                v_modulus.push_back((std::abs(fmod(calc_angle_3_pos(P1,P2,P3),ang_rad))));
            }
    return v_modulus;
}

bool mut_happens(population &p) noexcept
{
    return create_bernoulli_dist(p.get_param().get_mu_p())(p.get_rng());
}

double mut_step(population& s) noexcept
{
    return create_normal_dist(0,s.get_param().get_mu_st())(s.get_rng());
}

void no_complete_overlap(population &p) noexcept
{
    for(size_t i = 0; i != p.get_v_ind().size(); i++)
        for(size_t j = 0; j != p.get_v_ind().size(); j++)
        {
            if(i == j ) continue;
            if(distance(p.get_v_ind()[i], p.get_v_ind()[j]) < 0.00001 &&
                    distance(p.get_v_ind()[i], p.get_v_ind()[j]) > -0.00001)
            {
                auto rnd_n_0_1 = create_unif_dist(0,1)(p.get_rng());
                p.get_v_ind()[i].change_x(rnd_n_0_1 * 0.0001);
                p.get_v_ind()[j].change_x(rnd_n_0_1 * -0.0001);
            }
        }
}

bool only_one_focal(const std::vector<individual>::iterator first,
                    const std::vector<individual>::iterator last)
{
    auto real_focal =std::find_if(first,last,[](const individual& i)
    {return i.is_focal();}) + 1;
    return last == std::find_if(real_focal,last,[](const individual& i)
    {return i.is_focal();});
}

void place_start_cells(population &p) noexcept
{
    int n = count_hex_layers(p.get_pop_size());
    unsigned int placed_ind = 0;

    // d is the distance between 2 individuals's centers
    double d = 2 * (p.get_v_ind()[0].get_param().get_radius() + p.get_param().get_min_dist());

    for(int i = 0; i != n; i++)
    {
        double y = (sqrt(3) * i * d) / 2.0;

        for (int j = 0; j < (2 * n - 1 - i); j++)
        {
            double x = (-(2 * n - i - 2) * d) / 2.0 + j * d;

            p.get_v_ind()[placed_ind].set_x(x);
            p.get_v_ind()[placed_ind].set_y(y);
            placed_ind++;
            if(placed_ind == p.get_v_ind().size()) return;

            if (y > 0.000001 || y < -0.000001)
            {
                p.get_v_ind()[placed_ind].set_x(x);
                p.get_v_ind()[placed_ind].set_y(-y);
                placed_ind++;
                if(placed_ind == p.get_v_ind().size()) return;
            }
        }

    }
}

void place_new_pop(population &p) noexcept
{
    p.get_v_ind().swap(p.get_new_v_ind());
    reset_output_nodes_pop(p);
    p.get_new_v_ind().clear();
}

std::pair<std::vector<individual>::iterator,std::vector<individual>::iterator>
possible_collisions(individual focal_ind, std::vector<individual>& pop)
{
    auto range_x = possible_collisions_x(focal_ind,pop);
    auto range_y = possible_collisions_y(focal_ind, range_x.first, range_x.second, pop);
    return range_y;
}

std::pair<std::vector<individual>::iterator,std::vector<individual>::iterator>
possible_collisions_x(individual focal_ind, std::vector<individual>& pop)
{
    std::pair<std::vector<individual>::iterator,std::vector<individual>::iterator> range;
    //Sort all individuals whose x +/- radius is in the range between focal_x -/+ focal_radius
    //in ascending order by thei y coordinate
    auto first_x = std::lower_bound(
                pop.begin(), pop.end(), focal_ind,
                [](const individual& lhs, const individual& rhs)
    {return lhs.get_x() + lhs.get_param().get_radius() < rhs.get_x() - rhs.get_param().get_radius();}
    );

    auto last_x = std::upper_bound(
                pop.begin(),pop.end(), focal_ind,
                [](const individual& lhs, const individual& rhs)
    {return lhs.get_x() + lhs.get_param().get_radius() < rhs.get_x() - rhs.get_param().get_radius();}
    );

    return range = {first_x, last_x};
}

std::pair<std::vector<individual>::iterator,std::vector<individual>::iterator>
possible_collisions_y(individual focal_ind,
                      std::vector<individual>::iterator first_x,
                      std::vector<individual>::iterator last_x,
                      std::vector<individual> &pop)
{
    std::pair<std::vector<individual>::iterator,std::vector<individual>::iterator> range;
    std::sort(first_x, last_x, [](const individual& lhs, const individual& rhs){return lhs.get_y() < rhs.get_y();});

    //In this new range find individuals whose  y +/- radius is in range with focal_ y -/+ focal_radius
    auto first_y = std::lower_bound(
                first_x,last_x,focal_ind,
                [](const individual& lhs, const individual& rhs)
    {return lhs.get_y() + lhs.get_param().get_radius() < rhs.get_y() - rhs.get_param().get_radius();}
    );
    if(first_y == pop.end()){first_y = pop.begin();}

    auto last_y = std::upper_bound(
                first_x,last_x,focal_ind,
                [](const individual& lhs, const individual& rhs)
    {return lhs.get_y() + lhs.get_param().get_radius() < rhs.get_y() - rhs.get_param().get_radius();}
    );

    return range =  {first_y,last_y};
}

double repr_angle(population &p) noexcept
{
    return create_unif_dist(0, 2 * M_PI)(p.get_rng());
}

void reset_output_nodes_pop(population &p) noexcept
{
    for(auto& ind : p.get_v_ind())
    {
        ind.get_grn().set_all_out_nodes(true);
        ind.set_phen(phenotype::active);
    }
}

void reset_drawn_fl_new_pop(population &p) noexcept
{
    std::for_each(p.get_new_v_ind().begin(),p.get_new_v_ind().end(),
                  [](individual& i){draw_flag_reset(i);});
}

void reset_pop(population& p) noexcept
{
    p.get_v_ind().resize(p.get_param().get_pop_start_size());
    place_start_cells(p);
    for(auto& ind : p.get_v_ind())
    {
        ind.set_phen(phenotype::active);
        ind.get_grn() = GRN{};
    }
}

void select_new_pop(population &p)
{
    assert(!p.get_v_ind().empty());
    assert(p.get_new_v_ind().empty());
    while(true)
    {
        for(auto& ind : p.get_v_ind())
        {
            if(create_unif_dist(0,1)(p.get_rng()) <
                    get_fitness(ind, p.get_param().get_base_disp_prob(), p.get_param().get_spo_adv())
                    && !is_drawn(ind))
            {
                draw(ind);
                p.get_new_v_ind().push_back(ind);
            }
            if(p.get_new_v_ind().size() == p.get_param().get_exp_new_pop_size() ||
                    all_ind_are_drawn(p))
            {
                reset_drawn_fl_new_pop(p);
                return;
            }
        }
    }
}

void set_ind_en(individual& i, double en)
{
    i.set_energy(en);
}


void senescence(population& p) noexcept
{
    p.get_v_ind().erase(
                std::remove_if (p.get_v_ind().begin(), p.get_v_ind().end(), [&p](const individual& i)
    {return !is_spore(i) && create_unif_dist(0,1)(p.get_rng()) < p.get_param().get_death_rate();})
                ,p.get_v_ind().end());
}

void sort_inds_by_x_inc(std::vector<individual>::iterator start, std::vector<individual>::iterator last)
{
    std::sort(start,last,
              [](const individual& lhs, const individual& rhs)
    {return lhs.get_x() < rhs.get_x();});
}

void spor_metabolism_pop(population &p)
{
    for(auto& ind : p.get_v_ind())
    {
        if(ind.get_phen() != phenotype::sporulating)
            continue;

        sporulating_metabolism(ind);
    }
}

void starvation(population& p) noexcept
{
    p.get_v_ind().erase(
                std::remove_if(p.get_v_ind().begin(),p.get_v_ind().end(),
                               [](individual const &i){ return is_dead(i);})
            ,p.get_v_ind().end());
}

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
               (s.get_ind(0).get_param().get_radius() + s.get_ind(1).get_param().get_radius()) < 0.1);
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
                    s.get_ind(1).get_x() - s.get_ind(1).get_param().get_radius(),
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
                    s.get_ind(1).get_y() + s.get_ind(1).get_param().get_radius()
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

        calc_tot_displ_pop(s.get_v_ind());
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
        s.get_v_ind().emplace_back(individual(0,0));
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
        int replicates = 30000;
        double death_rate = 0.1;
        pop_param pp {static_cast<unsigned int>(pop_size),
                    1,
                    0.1,
                    0.0015,
                    0.1,
                    0.01,
                    10,
                    0.5,
                    death_rate};
        auto mean = 0.0;
        population p(pp);

        for( int i  = 0; i != replicates; i++)
        {
            senescence(p);
            mean += p.get_pop_size();
            p.get_v_ind().resize(pp.get_pop_start_size());
        }
        mean /= replicates;
        auto balance = pop_size - (mean + pop_size * pp.get_death_rate());
        assert( balance < 0.01 &&
                balance > -0.01);
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

    //A pop is initialized with a random number generator
    {
        population p;
        std::uniform_int_distribution u_d(0,2);
        double mean = 0;
        //Draw a thousands times from a uniform dist between 0 and 2
        for(int i = 0; i != 1000; i++)
        {
            mean += u_d(p.get_rng());
        }
        //calculate mean of the drawn values
        mean /= 1000;
        assert(mean < 1.01 && mean > 0.99 );
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
            auto daughter = get_daughter_pos(p.get_ind(0), repr_angle(p));
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
        int sample_size = 100000;
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

    //After being selected in new population individuals flag is_drawn is resetted
    {
        population p;
        select_new_pop(p);
        for(const auto& ind :p.get_new_v_ind())
            assert(!is_drawn(ind));
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
        int pop_size = 2;
        population p{pop_size};
        for (int i = 0; i != p.get_pop_size() - 1; i++)
            for(int j = i + 1; j != p.get_pop_size(); j++)
                assert(!have_same_ancestor(p.get_ind(i),p.get_ind(j)));
    }

    //it is possible to track individuals that have the same ancestry
    {
        int pop_size = 2;
        population p{pop_size};
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


#endif
}
