#include "population.h"
#include "utilities.h"

population::population(const pop_param& pop_parameters, const ind_param& ind_parameters):
    m_pop_param{pop_parameters},
    m_pop(m_pop_param.get_pop_start_size(), ind_parameters),
    m_new_pop_tmp_buffer(0, ind_parameters),
    m_relax(relaxation::param_t{})
{
    assign_ancestor_ID(m_pop).swap(m_pop);
    if(!m_pop.empty())
    {
        place_start_cells(*this);
    }
}

bool operator== (const population& lhs, const population& rhs)
{
    return lhs.get_param() == rhs.get_param()
            && lhs.get_v_ind() == rhs.get_v_ind();
}

bool operator!=  (const population& lhs, const population& rhs)
{
    return !(lhs == rhs);
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

bool all_en_pop_equals(const population& p, double en)
{
    return std::all_of(p.get_v_ind().begin(), p.get_v_ind().end(),
                       [&en](const individual& i)
    {return i.get_energy() - en > -0.0001 && i.get_energy() - en < 0.00001;});
}


bool all_en_pop_NOT_equals(const population& p, double en)
{
    return !all_en_pop_equals(p,en);
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

int count_actives(const population& pop)
{
    const auto& p = pop.get_v_ind();
    return  std::count_if(p.begin(),p.end(),
                          [](const individual& ind){ return ind.get_phen() == phenotype::active;});
}

int count_spores(const population& pop)
{
    const auto& p = pop.get_v_ind();
    return  std::count_if(p.begin(),p.end(),
                          [](const individual& ind){ return ind.get_phen() == phenotype::spore;});
}

int count_sporulating(const population& pop)
{
    const auto& p = pop.get_v_ind();
    return  std::count_if(p.begin(),p.end(),
                          [](const individual& ind){ return ind.get_phen() == phenotype::sporulating;});
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

bool calc_tot_displ_pop(population& population)
{
    std::vector<individual>& pop = population.get_v_ind();
    bool has_collision = false;
    auto wiggle_room = population.get_param().get_wiggle_room();
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
            if (are_colliding(*focal_ind_address, *j, wiggle_room))
            {
                assert(!j->is_focal());
                add_displacement(*focal_ind_address, *j, wiggle_room);
                has_collision = true;
            }
        }
        //Sort back to increasing x
        focal_ind_address->no_more_focal();
        sort_inds_by_x_inc(range_x.first, range_x.second);
    }
    return has_collision;
}

std::vector<individual> set_new_ind_par(std::vector<individual> new_p_v, const ind_param& new_ind_params)
{
    for(auto& individual : new_p_v)
    {
        individual.set_param(new_ind_params);
    }
    return new_p_v;
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
            if(repr_prob(p.get_rng()) < p.get_ind(div_ind).get_param().get_repr_prob())
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
    reset_output_nodes_pop(p);
    update_radius_pop(p);
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
    auto wiggle_room =  p.get_param().get_wiggle_room();

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
            if (are_colliding(*focal_ind_address, *j, wiggle_room))
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
    { return particle_t{ glm::vec2{i.get_x(), i.get_y()}, static_cast<float>(i.get_radius()) }; },
    [](individual& i, const particle_t& part)
    { i.set_x(part.pos.x); i.set_y(part.pos.y); i.get_radius() = part.radius; }
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
    if(p.get_v_ind().empty())
        return;

    int n = count_hex_layers(p.get_pop_size());
    unsigned int placed_ind = 0;

    // d is the distance between 2 individuals's centers
    double d = 2 * (p.get_v_ind()[0].get_radius() + p.get_param().get_min_dist());

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
    {return lhs.get_x() + lhs.get_radius() < rhs.get_x() - rhs.get_radius();}
    );

    auto last_x = std::upper_bound(
                pop.begin(),pop.end(), focal_ind,
                [](const individual& lhs, const individual& rhs)
    {return lhs.get_x() + lhs.get_radius() < rhs.get_x() - rhs.get_radius();}
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
    {return lhs.get_y() + lhs.get_radius() < rhs.get_y() - rhs.get_radius();}
    );
    if(first_y == pop.end()){first_y = pop.begin();}

    auto last_y = std::upper_bound(
                first_x,last_x,focal_ind,
                [](const individual& lhs, const individual& rhs)
    {return lhs.get_y() + lhs.get_radius() < rhs.get_y() - rhs.get_radius();}
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
        ind.reset_spo_timer();
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
    if(p.get_v_ind().empty())
    {
        return;
    }
    assert(p.get_new_v_ind().empty());
    std::shuffle(p.get_v_ind().begin(),p.get_v_ind().end(),p.get_rng());
    //The energy level at which an individual will be initialized
    auto default_energy = individual{}.get_energy();

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
                p.get_new_v_ind().back().set_energy(default_energy);
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

std::vector<individual> pop_from_funders(const funders_success &f_s,
                                         const demographic_sim& d_s,
                                         int generation)
{
    auto funders = f_s.get_v_funders()[generation];

    individual base_ind{d_s.get_demo_cycles()[generation].get_ind_param()} ;
    std::vector<individual> pop;

    for(const auto& funder : funders.get_v_funder_data())
    {
        base_ind.get_grn() = funder.get_grn();
        pop.push_back(base_ind);
    }

    return pop;
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

bool someone_starved(const population& p)
{
    return std::any_of(p.get_v_ind().begin(), p.get_v_ind().end(), [](const individual& i)
    {return is_dead(i);});
}


void starvation( population& p) noexcept
{

    p.get_v_ind().erase(
                std::remove_if(
                p.get_v_ind().begin(), p.get_v_ind().end(), [](const individual &i){return is_dead(i);}
            )
            ,p.get_v_ind().end());
}

void update_radius_pop(population& p)
{
    for(auto& ind : p.get_v_ind())
    {
        update_radius(ind);
    }
}

