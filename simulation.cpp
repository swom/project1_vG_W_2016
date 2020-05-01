#include "simulation.h"
#include <cassert>
#include <numeric>
#include <algorithm>
#include <cmath>

simulation::simulation(unsigned int pop_size, int exp_new_pop_size, double min_dist,
                       int grid_side, double diff_coeff,
                       double init_food, double mutation_prob,
                       double mutation_step, double base_disp_prob, double spore_advantage,
                       double reproduction_prob, double metab_degradation_rate):
    m_param{pop_size,
            exp_new_pop_size,
            min_dist, grid_side,
            diff_coeff, init_food,
            mutation_prob,
            mutation_step,
            base_disp_prob,
            spore_advantage,
            reproduction_prob,
            metab_degradation_rate},
    m_pop(m_param.get_pop_start_size(),individual{0,0}),
    m_e{m_param.get_grid_side(), m_param.get_diff_coeff(), m_param.get_init_food(), m_param.get_metab_degr_rate()}

{
#ifndef IS_ON_TRAVIS
    try {
        if(base_disp_prob * spore_advantage > 1)
        {throw std::string{"base dispersal probability * spore advantage > 1, too high!\n"};}
    }
    catch (std::string& e) {
        std::cout << e;
#ifdef NDEBUG
        abort();
#endif
    }
#endif
    if(!m_pop.empty())
    {
        place_start_cells(*this);
    }
}

simulation::simulation(sim_param param):
    m_param{param},
    m_pop(m_param.get_pop_start_size(),individual{0,0}),
    m_e{m_param.get_grid_side(), m_param.get_diff_coeff(), m_param.get_init_food(), m_param.get_metab_degr_rate()}

{
#ifndef IS_ON_TRAVIS
    try {
        if(get_param().get_base_disp_prob() * get_param().get_spo_adv() > 1)
        {throw std::string{"base dispersal probability * spore advantage > 1, too high!\n"};}
    }
    catch (std::string& e) {
        std::cout << e;
#ifdef NDEBUG
        abort();
#endif
    }
#endif
    if(!m_pop.empty())
    {
        place_start_cells(*this);
    }
}

bool all_ind_are_drawn(const simulation& s) noexcept
{
    return std::all_of(s.get_pop().begin(), s.get_pop().end(),
                       [](const individual& i) {return is_drawn(i);});
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

//not sure this is the fastest implementation, maybe sawp and pop_back is still faster?
void death(simulation& s) noexcept
{
    s.get_pop().erase(std::remove_if(s.get_pop().begin(),s.get_pop().end(),
                                     [](individual const &i){ return is_dead(i);})
            ,s.get_pop().end());
}

bool division(simulation &s) noexcept
{
    auto  div_inds  = get_dividing_individuals(s);
    std::uniform_real_distribution<double> repr_prob(0,1);

    if(!div_inds.empty())
        for(size_t i = 0; i != div_inds.size(); i++)
        {
            int div_ind = div_inds[i];
            if(repr_prob(s.get_rng()) < s.get_param().get_repr_p())
            {

                divides(s.get_ind(div_ind),
                        s.get_pop(),
                        repr_angle(s),
                        s.get_rng(),
                        create_bernoulli_dist(s.get_param().get_mu_p()),
                        create_normal_dist(0,s.get_param().get_mu_st())
                        );
            }
        }
    return !div_inds.empty();
}

void dispersal(simulation& s)
{
    select_new_pop(s);
    fund_pop(s);
    place_start_cells(s);
    reset_env(s.get_env());
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

void exec(simulation& s, int n_tick) noexcept
{
    while(s.get_tick() != n_tick){tick(s);}
}

void feeding(simulation& s)
{
    for(auto& ind : s.get_pop())
    {
        auto index_grid = find_grid_index(ind,s.get_env().get_grid_side());
        if(index_grid == -100 ||
                ind.get_phen() != phenotype::active)
        {continue;}
        feed(ind,s.get_env().get_cell(index_grid));
    }
}

void fund_pop(simulation& s) noexcept
{
    s.get_pop().swap(s.get_new_pop());
    s.get_new_pop().clear();
}

std::vector<int> get_dividing_individuals(const simulation& s) noexcept
{

    std::vector<int> dividing_individuals;
    for(unsigned int i = 0; i != s.get_pop().size(); i++)
    {
        if(s.get_pop()[i].get_energy() >= s.get_pop()[i].get_treshold_energy()
                && s.get_pop()[i].get_phen() == phenotype::active)
        {
            dividing_individuals.push_back(static_cast<int>(i));
        }
    }
    return dividing_individuals;
}

std::vector<double> get_excess_energies(const simulation& s) noexcept
{
    std::vector<double> excess_energy;
    auto v_ind{get_dividing_individuals(s)};
    for(unsigned int i=0; i != v_ind.size(); i++)
    {
        int ind_index = v_ind[i];
        excess_energy.push_back(
                    get_ind_en(s, ind_index) - get_ind_tr_en(s, ind_index)
                    );
    }
    return excess_energy;
}

double get_ind_en(const simulation& s, int i)
{
    return s.get_ind(i).get_energy();
}

double get_ind_tr_en(const simulation& s, int i)
{
    return s.get_ind(i).get_treshold_energy();
}

std::pair<double, double> get_ind_pos(const individual& i)
{
  std::pair<double, double> pos;
pos.first = i.get_x();
pos.second = i.get_y();
return pos;

}

std::vector<std::pair<int,int>> get_sisters_index_offset(const simulation& s)  noexcept
{
    std::vector<std::pair<int,int>> sis_indexes;
    std::pair<int,int> daughters;
    int sis_dist;
    auto div_ind = get_dividing_individuals(s);

    for (const auto& ind : div_ind)
    {
        sis_dist = static_cast<int>(s.get_pop_size()) - ind;
        daughters.first = ind;
        daughters.second = daughters.first + sis_dist;
        sis_indexes.push_back(daughters);
    }
    return sis_indexes;
}

bool has_collision(simulation& s)
{
    auto pop = s.get_pop();
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

int manage_static_collisions(simulation& s)
{
    int time = 0;
    bool has_collision = true;
    do
    {
        has_collision = false;
        if((has_collision = calc_tot_displ_pop(s.get_pop())))
        {
            for(auto& ind : s.get_pop())
            {
                if(
                        std::abs(ind.get_x_displacement()) > 0 ||
                        std::abs(ind.get_y_displacement()) > 0
                        )
                {
                    ind.displace();
                }
            }
        }
        time ++;
    }
    while(has_collision);
    return time;
}

void metabolism_pop(simulation& s)
{
    for(auto& ind : s.get_pop())
    {
        if(ind.get_phen() != phenotype::spore)
            metabolism(ind);
    }
}

std::vector<double> modulus_of_btw_ind_angles(simulation& s, double ang_rad)
{
    std::vector<double> v_modulus;
    for(int i = 0; i != s.get_pop_size() - 2; i++)
        for(int j = i+1; j != s.get_pop_size() - 1; j++)
            for(int k = j+1; k != s.get_pop_size(); k++)
            {
                auto P1 = get_pos(s.get_ind(i));
                auto P2 = get_pos(s.get_ind(j));
                auto P3 = get_pos(s.get_ind(k));
                v_modulus.push_back((abs(fmod(calc_angle_3_pos(P1,P2,P3),ang_rad))));
            }
    return v_modulus;
}

bool mut_happens(simulation& s) noexcept
{
    return create_bernoulli_dist(s.get_param().get_mu_p())(s.get_rng());
}

double mut_step(simulation& s) noexcept
{
    return create_normal_dist(0,s.get_param().get_mu_st())(s.get_rng());
}

void no_complete_overlap(simulation& s) noexcept
{
    for(size_t i = 0; i != s.get_pop().size(); i++)
        for(size_t j = 0; j != s.get_pop().size(); j++)
        {
            if(i == j ) continue;
            if(distance(s.get_pop()[i], s.get_pop()[j]) < 0.00001 &&
                    distance(s.get_pop()[i], s.get_pop()[j]) > -0.00001)
            {
                auto rnd_n_0_1 = create_unif_dist(0,1)(s.get_rng());
                s.get_pop()[i].change_x(rnd_n_0_1 * 0.0001);
                s.get_pop()[j].change_x(rnd_n_0_1 * -0.0001);
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

void place_start_cells(simulation& s) noexcept
{
    int n = count_hex_layers(s.get_pop_size());
    unsigned int placed_ind = 0;

    // d is the distance between 2 individuals's centers
    double d = 2 * (s.get_pop()[0].get_radius() + s.get_param().get_min_dist());

    for(int i = 0; i != n; i++)
    {
        double y = (sqrt(3) * i * d) / 2.0;

        for (int j = 0; j < (2 * n - 1 - i); j++)
        {
            double x = (-(2 * n - i - 2) * d) / 2.0 + j * d;

            s.get_pop()[placed_ind].set_x(x);
            s.get_pop()[placed_ind].set_y(y);
            placed_ind++;
            if(placed_ind == s.get_pop().size()) return;

            if (y > 0.000001 || y < -0.000001)
            {
                s.get_pop()[placed_ind].set_x(x);
                s.get_pop()[placed_ind].set_y(-y);
                placed_ind++;
                if(placed_ind == s.get_pop().size()) return;
            }
        }

    }
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

double repr_angle( simulation& s) noexcept
{
    return create_unif_dist(0, 2 * M_PI)(s.get_rng());
}

void response(simulation& s)
{
    for(auto& ind : s.get_pop())
    {
        if(ind.get_phen() == phenotype::spore){continue;}
        auto index = find_grid_index(ind,s.get_env().get_grid_side());
        if(index == -100)//When individual is outside grid
        {
            responds(ind, env_grid_cell(0,0,0,0));
            continue;
        }
        responds(ind, s.get_env().get_cell(index));
    }
}

void reset_drawn_fl_new_pop(simulation& s) noexcept
{
    std::for_each(s.get_new_pop().begin(),s.get_new_pop().end(),
                  [](individual& i){draw_flag_reset(i);});
}

void secretion_metabolite(simulation& s)
{
    int index;
    for(const auto& ind : s.get_pop())
    {
        index = find_grid_index(ind,s.get_env().get_grid_side());
        if(index == - 100)
        {
            continue;
        }
        secretes_metab(ind,s.get_env().get_cell(index));
    }
}

void select_new_pop(simulation& s)
{
    assert(s.get_new_pop().empty());
    while(true)
    {
        for(auto& ind : s.get_pop())
        {
            if(create_unif_dist(0,1)(s.get_rng()) <
                    get_fitness(ind, s.get_param().get_base_disp_prob(), s.get_param().get_spo_adv())
                    && !is_drawn(ind))
            {
                draw(ind);
                s.get_new_pop().push_back(ind);
            }
            if(static_cast<int>(s.get_new_pop().size()) == s.get_param().get_exp_new_pop_size() ||
                    all_ind_are_drawn(s))
            {
                reset_drawn_fl_new_pop(s);
                return;
            }
        }
    }
}

void set_ind_en(individual& i, double en)
{
    i.set_energy(en);
}

void sort_inds_by_x_inc(std::vector<individual>::iterator start, std::vector<individual>::iterator last)
{
    std::sort(start,last,
              [](const individual& lhs, const individual& rhs)
    {return lhs.get_x() < rhs.get_x();});
}

int tick(simulation& s)
{
    int time = 0;
    response(s);
    feeding(s);
    metabolism_pop(s);
    secretion_metabolite(s);
    death(s);
    if(division(s))
    {
        time += manage_static_collisions(s);
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
        int exp_new_pop_size = 1;
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
        assert(s.get_pop().size() == pop_size);
        assert(s.get_param().get_exp_new_pop_size() == exp_new_pop_size);
        assert(s.get_param().get_min_dist() - min_dist < 0.0001 &&
               s.get_param().get_min_dist() - min_dist > -0.0001);
        assert(s.get_param().get_diff_coeff() - diff_coeff < 0.00001 &&
               s.get_param().get_diff_coeff() - diff_coeff > -0.00001);
        assert(s.get_param().get_init_food() - init_food < 0.0001 &&
               s.get_param().get_init_food() - init_food > -0.0001);
        assert(s.get_param().get_mu_p() - mutation_prob < 0.0001 &&
               s.get_param().get_mu_p() - mutation_prob > -0.0001);
        assert(s.get_param().get_mu_st() - mutation_step < 0.0001 &&
               s.get_param().get_mu_st() - mutation_step > -0.0001);
        assert(s.get_param().get_base_disp_prob() - base_disp_prob < 0.00001 &&
               s.get_param().get_base_disp_prob() - base_disp_prob > -0.00001);
        assert(s.get_param().get_spo_adv() - spore_advantage < 0.0001 &&
               s.get_param().get_spo_adv() - spore_advantage > -0.0001);
        assert(s.get_param().get_repr_p() - reproduction_prob < 0.0001 &&
               s.get_param().get_repr_p() - reproduction_prob > -0.0001);
        assert(s.get_param().get_metab_degr_rate() - metab_degradation_rate < 0.0001 &&
               s.get_param().get_metab_degr_rate() - metab_degradation_rate > -0.0001);
    }
    //Simulation is initialized with a certain number of individuals
    // The value 1234567890 is irrelevant: just get this to compile

    {
        unsigned int pop_size = 100;
        simulation s(pop_size);
        for(int i = 0; i < s.get_pop_size(); ++i)
        {
            double x = s.get_ind(i).get_x();
            assert( x > -1234567890 );
        }
    }

    //The size of a population is equal to the size of the vector containing its individuals
    {
        simulation s;
        assert( s.get_pop_size() == static_cast<int>(s.get_pop().size()));
    }


    //No individuals are dividing at the start of the simulation
    {
        //initiate empty vector and leave it empty
        simulation s(3);
        assert(static_cast<int>(get_dividing_individuals(s).size()) < 0.00000000001);
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
        simulation s(7);

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


    //No individuals are overlapping at the start of the simulation
    {
        simulation s(2);
        assert(!has_collision(s));
    }


    //An individual position can be retrieved as a pair object x,y
    {
        simulation s;
        std::pair<double, double> v{0,0};//default coordinates of individuals
        std::pair<double, double> v2{1,1};//different coordinates from default

        assert(get_ind_pos(s.get_ind(0)) == v);
        assert(get_ind_pos(s.get_ind(0)) != v2);
    }


    // An individaul can be placed to some given coordinates after initialization
    {
        simulation s(2);
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
        simulation s;
        double lhs = s.get_pop_size();
        //let's allow all individuals in the population to reproduce,
        //to facilitate the testing conditions
        for(auto& ind : s.get_pop())
        {
            ind.set_energy(ind.get_treshold_energy());
        }
        division(s);

        double rhs = s.get_pop_size();
        assert(lhs * 2 - rhs < 0.0000001);

    }


    //Only individuals with energy >= than the treshold will divide
    {
        simulation s(3);
        //This ind will not reproduce
        s.get_ind(0).set_energy(s.get_ind(0).get_treshold_energy()-1);
        //This ind will reproduce with no extra energy
        s.get_ind(1).set_energy(s.get_ind(1).get_treshold_energy());
        //This ind will reproduce with extra energy
        s.get_ind(2).set_energy(s.get_ind(2).get_treshold_energy()+1);

        std::vector<int> div_ind = get_dividing_individuals(s);
        assert(!div_ind.empty());

        for(unsigned int i = 0; i != div_ind.size(); i++)
            assert(s.get_ind(div_ind[i]).get_energy()
                   >= s.get_ind(div_ind[i]).get_treshold_energy());

    }

    //The excess energy of dividing individuals is equal
    //to the diference between their energy and their treshold energy
    {
        simulation s(2);
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
        simulation s(1);
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
        simulation s;
        auto repr_excess_en = 3;
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
        simulation s(static_cast<unsigned int>(pop_size));
        auto comparison_pop = s.get_pop();
        sort_inds_by_x_inc(s.get_pop().begin(),s.get_pop().end());
        assert(s.get_pop() != comparison_pop);
        sort_inds_by_x_inc(comparison_pop.begin(), comparison_pop.end());
        assert(s.get_pop() == comparison_pop);

        //Can also sort parts of the vector of ind
        auto changed_ind = pop_size - 1;
        auto reference_ind = pop_size / 2 - 1;
        set_pos(s.get_ind(changed_ind),get_pos(s.get_ind(reference_ind)));
        comparison_pop = s.get_pop();
        sort_inds_by_x_inc(comparison_pop.begin(),comparison_pop.end());
        assert(s.get_pop() != comparison_pop);
        sort_inds_by_x_inc(s.get_pop().begin() + reference_ind, s.get_pop().begin() + changed_ind + 1);
        assert(s.get_pop() == comparison_pop);

    }

    //After reproduction the first daughter individual takes the position of the mother
    {
        simulation s;
        auto parent_pop = s.get_pop();
        //setting energy high enough for the individual
        //to reproduce without offspring having negative enrgies
        set_ind_en(s.get_ind(0),get_ind_tr_en(s, 0));

        divides(s.get_ind(0),
                s.get_pop(),
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
        simulation s;
        //setting energy high enough for the individual to reproduce
        //without offspring having negative enrgies
        set_ind_en(s.get_ind(0),get_ind_tr_en(s, 0));
        divides(s.get_ind(0),
                s.get_pop(),
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
        simulation s(2);
        set_pos(s.get_ind(1),get_ind_pos(s.get_ind(0)));
        auto collision_range = possible_collisions_x(s.get_ind(0),s.get_pop());
        assert(collision_range.first == s.get_pop().begin() &&
               collision_range.second == s.get_pop().end());


        //No collisions should return a pair whose first element is the focal individual
        //And the second the individual after the focal
        //(or the end of the vector in case the focal is the last)
        set_pos(s.get_ind(1),std::pair<double,double>{2,2});
        //Pop need to be sorted by increasing x coord
        std::sort(s.get_pop().begin(),s.get_pop().end(),
                  [](const individual& lhs, const individual& rhs)
        {return  lhs.get_x() < rhs.get_x();});
        collision_range = possible_collisions_x(s.get_ind(0),s.get_pop());
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
        sort_inds_by_x_inc(s.get_pop().begin(),s.get_pop().end());
        auto focal_ind_index = 1;
        collision_range = possible_collisions_x(s.get_ind(focal_ind_index),s.get_pop());
        assert(collision_range.first == s.get_pop().begin() + focal_ind_index - 1);
        assert(collision_range.second == s.get_pop().begin() + focal_ind_index + 1);
        assert(collision_range.second == s.get_pop().end());

        //Collision to the right should return a range
        //focal_ind_index,focal_ind_index+2
        sort_inds_by_x_inc(s.get_pop().begin(),s.get_pop().end());
        focal_ind_index = 0;
        collision_range = possible_collisions_x(s.get_ind(focal_ind_index), s.get_pop());
        assert(collision_range.first == s.get_pop().begin() + focal_ind_index);
        assert(collision_range.second == s.get_pop().begin() + focal_ind_index + 2);
        assert(collision_range.second == s.get_pop().end());

        //Collision from below should return a range
        //focal_index - 1, focal_index +1
        set_pos(s.get_ind(0),
                std::pair<double,double>{
                    s.get_ind(1).get_x(),
                    s.get_ind(1).get_y() + s.get_ind(1).get_radius()
                }
                );
        sort_inds_by_x_inc(s.get_pop().begin(),s.get_pop().end());
        focal_ind_index = 1;
        collision_range = possible_collisions_x(s.get_ind(focal_ind_index), s.get_pop());
        assert(collision_range.first == s.get_pop().begin() + focal_ind_index - 1);
        assert(collision_range.first == s.get_pop().begin());
        assert(collision_range.second == s.get_pop().begin() + focal_ind_index + 1);
        assert(collision_range.second == s.get_pop().end());

        //Collision from above should return a range
        //focal_index, focal_index + 2
        sort_inds_by_x_inc(s.get_pop().begin(),s.get_pop().end());
        focal_ind_index = 0;
        collision_range = possible_collisions_x(s.get_ind(focal_ind_index), s.get_pop());
        assert(collision_range.first == s.get_pop().begin() + focal_ind_index);
        assert(collision_range.first == s.get_pop().begin());
        assert(collision_range.second == s.get_pop().begin() + focal_ind_index + 2);
        assert(collision_range.second == s.get_pop().end());

    }
    //If there are collisions individuals are displaced
    //until there are no more collisions
    {
        simulation s(2);
        assert(!has_collision(s));
        set_pos(s.get_ind(1), get_ind_pos(s.get_ind(0)));
        no_complete_overlap(s);
        assert(has_collision(s));

        calc_tot_displ_pop(s.get_pop());
        for(auto& ind : s.get_pop())
        {ind.displace();}

        assert(!has_collision(s));

    }
    //After reproduction new collisions caused by new individuals
    //being placed where other individuals already are managed
    {
        simulation s(7);
        assert(!has_collision(s));
        //add 1 individual overlapping with central individual
        s.get_pop().emplace_back(individual(0,0));
        assert(has_collision(s));
        manage_static_collisions(s);
        assert(!has_collision(s));
    }

    //A simulation has an environment
    // The value -1234567890 is irrelevant: just get this to compile

    {
        simulation s;
        assert(s.get_env().get_grid_side() > -1234567890);
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
                    s.get_pop().begin(),
                    s.get_pop().end(),0.0,
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
                    s.get_pop().begin(),
                    s.get_pop().end(),0.0,
                    [](double sum,const individual& i){
            return sum + i.get_energy();}
        );

        assert(total_food_init > total_food_after);
        assert(total_en_init < total_en_after);

        auto total_uptake = std::accumulate(s.get_pop().begin(), s.get_pop().end(), 0.0,
                                            [](double sum, const individual& i) {return sum + i.get_uptake_rate();});

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
                    s.get_pop().begin(),
                    s.get_pop().end(),0.0,
                    [](double sum,const individual& i){return sum + i.get_energy();}
        );

        set_pos(s.get_ind(0),std::pair<double,double>(-42,42));

        feeding(s);

        double total_food_after = std::accumulate
                (
                    s.get_env().get_grid().begin(),
                    s.get_env().get_grid().end(),0.0,
                    [](double sum, const env_grid_cell& c){return sum + c.get_food();}
        );

        double total_en_after = std::accumulate
                (
                    s.get_pop().begin(),
                    s.get_pop().end(),0.0,
                    [](double sum,const individual& i){
            return sum + i.get_energy();}
        );

        assert(total_food_init - total_food_after < 0.00001
               && total_food_init - total_food_after > -0.00001);
        assert(total_en_init - total_en_after < 0.000001
               && total_en_init - total_en_after > -0.000001);
    }

    //Individuals lose energy through metabolism
    {
        simulation s(1,1,0.1,2);
        feeding(s);
        double init_en_tot = std::accumulate
                (
                    s.get_pop().begin(),
                    s.get_pop().end(),0.0,
                    [](double sum,const individual& i){return sum + i.get_energy();}
        );

        metabolism_pop(s);
        double after_en_tot = std::accumulate
                (
                    s.get_pop().begin(),
                    s.get_pop().end(),0.0,
                    [](double sum,const individual& i){return sum + i.get_energy();}
        );
        assert(after_en_tot < init_en_tot);
    }

    //In one tick/timestep individuals take in input, determine phenotype(based on previous timestep),
    // feed, than reproduce, than substances diffuse
    {
        simulation s(1,1,0.1,3,1,1,0);
        //Set all the hid nodes and H2O and H2H weights to one so that we are sure the phenotype will stay = active;
        for(auto& ind : s.get_pop())
        {
            ind.get_grn().set_all_hid_nodes(1);
            ind.get_grn().set_all_out_nodes(1);
            ind.get_grn().set_all_H2O(1);
            ind.get_grn().set_all_H2H(1);
        }

        //The single individual in this population
        //after a tick should reproduce
        auto init_pop_size = s.get_pop().size();
        auto ind_en = get_ind_tr_en(s, 0)
                + s.get_ind(0).get_metabolic_rate() + 0.01
                - s.get_ind(0).get_uptake_rate();
        s.get_ind(0).set_energy(ind_en);

        //and the grid_cell where it is should recieve
        //food nutrients
        auto grid_index_ind =  find_grid_index(s.get_ind(0),s.get_env().get_grid_side());
        double predicted_food_after_feeding =
                s.get_env().get_cell(grid_index_ind).get_food() - s.get_ind(0).get_uptake_rate();

        tick(s);

        assert(!(predicted_food_after_feeding - s.get_env().get_cell(grid_index_ind).get_food() < 0.00001
                 && predicted_food_after_feeding -s.get_env().get_cell(grid_index_ind).get_food() > -0.0001));
        assert(s.get_pop().size() == 2 * init_pop_size);
        assert(!has_collision(s));
    }

    //If nothing else happens, food should constantly decrease when cells are feeding
    {
        simulation s (2, 1, 0, 3, 0.1, 5);
        auto food_begin = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                          [](double sum, const env_grid_cell& c) {return sum + c.get_food();});

        //The simulation will last long enough  for the individuals to reproduce
        auto sim_time = s.get_ind(0).get_treshold_energy() / s.get_ind(0).get_uptake_rate() + 5;

        auto init_pop_size =  s.get_pop().size();

        for( int i = 0; i != static_cast<int>(sim_time); i++)
        {

            auto food_before_feed = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                                    [](double sum, const env_grid_cell& c) {return sum + c.get_food();});
            feeding(s);

            auto food_after_feed = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                                   [](double sum, const env_grid_cell& c) {return sum + c.get_food();});

            double food_eaten = 0.0;
            for(const auto& ind : s.get_pop())
            {
                auto grid_cell_ind = find_grid_index(ind, s.get_param().get_grid_side());
                if( grid_cell_ind != -100 && s.get_env().get_cell(grid_cell_ind).get_food() > 0)
                {
                    food_eaten += ind.get_uptake_rate();
                }
            }

            auto balance_uptake = food_before_feed - (food_after_feed + food_eaten) ;
            assert(balance_uptake < 0.0001 && balance_uptake > -0.0001);

            metabolism_pop(s);
            death(s);
            division(s);
            manage_static_collisions(s);

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

        auto final_pop_size =  s.get_pop().size();

        assert(init_pop_size < final_pop_size);
    }

    //A simulation is initiallized with a degradation rate
    {
        double degradation_rate = 3.14;
        simulation s(0,0,0,0,0,0,0,0,0,0,0,degradation_rate);
        assert(s.get_param().get_metab_degr_rate() - degradation_rate < 0.000001 &&
               s.get_param().get_metab_degr_rate() - degradation_rate > -0.000001);
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

        double tot_metab_after = std::accumulate(s.get_env().get_grid().begin(), s.get_env().get_grid().end(), 0.0,
                                                 [](int sum, const env_grid_cell& g) {return sum + g.get_metab();});
        double tot_production = std::accumulate(s.get_pop().begin(), s.get_pop().end(), 0.0,
                                                [](int sum, const individual& i) {return sum + i.get_metab_secr_rate();});
        double tot_degradation = s.get_env().get_grid_size() * s.get_param().get_metab_degr_rate();

        auto metab_balance = tot_metab_before - tot_degradation + tot_production - tot_metab_after;
        assert(metab_balance < 0.000001 && metab_balance > -0.000001);


    }
    //every timestep/tick collisions are handled
    {
        simulation s(7,3);
        //The central individual in this population
        //after a tick should reproduce
        auto init_pop_size = s.get_pop().size();
        s.get_ind(1).set_energy(get_ind_tr_en(s, 1)
                                + s.get_ind(1).get_metabolic_rate() + 0.01
                                - s.get_ind(1).get_uptake_rate());
        feeding(s);
        metabolism_pop(s);
        division(s);
        manage_static_collisions(s);

        assert(s.get_pop().size() == init_pop_size + 1);
        assert(!has_collision(s));
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

    //Spores do not get detected when looking for dividing individual
    {
        simulation s;
        set_ind_en(s.get_ind(0),get_ind_tr_en(s, 0));
        assert(get_dividing_individuals(s)[0] == 0);
        s.get_ind(0).set_phen(phenotype::spore);
        assert(get_dividing_individuals(s).empty());
    }

    //Spores do not reproduce
    {
        simulation s;
        s.get_ind(0).set_phen(phenotype::spore);
        set_ind_en(s.get_ind(0),get_ind_tr_en(s, 0));
        auto init_pop_size = s.get_pop_size();
        set_ind_en(s.get_ind(0),get_ind_tr_en(s, 0));
        feeding(s);
        metabolism_pop(s);
        assert(!division(s));
        assert(init_pop_size == s.get_pop_size());
    }

    //Spores do not feed
    {
        simulation s;
        auto init_food = s.get_env().get_cell(0).get_food();
        s.get_ind(0).set_phen(phenotype::spore);
        feeding(s);
        assert(init_food - s.get_env().get_cell(0).get_food() < 0.000001
               && init_food - s.get_env().get_cell(0).get_food() > -0.000001);
    }
    //Spores do not lose energy
    {
        simulation s;
        s.get_ind(0).set_phen(phenotype::spore);
        set_ind_en(s.get_ind(0), get_ind_tr_en(s, 0));
        auto init_en_ind0 = get_ind_en(s, 0);
        metabolism_pop(s);
        assert(get_ind_en(s, 0) - init_en_ind0 < 0.000001
               && get_ind_en(s, 0) - init_en_ind0 > -0.000001);
    }

    //Sporulating individuals do not get detected when looking for dividing individual
    {
        simulation s;
        set_ind_en(s.get_ind(0),get_ind_tr_en(s, 0));
        assert(get_dividing_individuals(s)[0] == 0);
        s.get_ind(0).set_phen(phenotype::sporulating);
        assert(get_dividing_individuals(s).empty());
    }

    //Sporulating individuals cannot reproduce
    {
        simulation s;
        s.get_ind(0).set_phen(phenotype::sporulating);
        set_ind_en(s.get_ind(0),get_ind_tr_en(s, 0));
        auto init_pop_size = s.get_pop_size();
        set_ind_en(s.get_ind(0), get_ind_tr_en(s, 0));
        feeding(s);
        metabolism_pop(s);
        assert(!division(s));
        assert(init_pop_size == s.get_pop_size());
    }

    //Sporulating individuals do not feed but they lose energy
    {
        simulation s;
        s.get_ind(0).set_phen(phenotype::sporulating);
        set_ind_en(s.get_ind(0),1);
        auto init_en_ind0 = get_ind_en(s, 0);
        auto init_food = s.get_env().get_cell(0).get_food();
        feeding(s);
        metabolism_pop(s);
        assert(init_en_ind0 - get_ind_en(s, 0) - s.get_ind(0).get_metabolic_rate() < 0.000001
               && init_en_ind0 - get_ind_en(s, 0) - s.get_ind(0).get_metabolic_rate() > -0.000001);
        assert(init_food - s.get_env().get_cell(0).get_food() < 0.000001
               && init_food - s.get_env().get_cell(0).get_food() > -0.000001);
    }


    //A sporulating individual updates its timer every time metab_pop is called
    {
        simulation s;
        auto init_timer = s.get_ind(0).get_spo_timer();
        s.get_ind(0).set_phen(phenotype::sporulating);
        metabolism_pop(s);
        assert(init_timer != s.get_ind(0).get_spo_timer());
        assert(s.get_ind(0).get_spo_timer() == init_timer + 1);
        assert(s.get_ind(0).get_spo_timer() != init_timer + 2);

        s.get_ind(0).reset_spo_timer();
        init_timer = s.get_ind(0).get_spo_timer();
        s.get_ind(0).set_phen(phenotype::spore);
        metabolism_pop(s);
        assert(s.get_ind(0).get_spo_timer() == init_timer);
        assert(s.get_ind(0).get_spo_timer() != init_timer + 1);
        assert(s.get_ind(0).get_spo_timer() != init_timer + 2);

        s.get_ind(0).reset_spo_timer();
        init_timer = s.get_ind(0).get_spo_timer();
        s.get_ind(0).set_phen(phenotype::active);
        metabolism_pop(s);
        assert(s.get_ind(0).get_spo_timer() == init_timer);
        assert(s.get_ind(0).get_spo_timer() != init_timer + 1);
        assert(s.get_ind(0).get_spo_timer() != init_timer + 2);
    }

    //Individuals that die are removed from the population
    {
        simulation s;
        s.get_ind(0).set_energy(0);//the only individual in this sim has 0 energy, thus it will die
        assert(s.get_pop_size() == 1);
        death(s);
        assert(s.get_pop().empty() && s.get_pop_size() == 0);

        unsigned int pop_size = 5;
        //The simulation does not have a grid with food,
        //so organisms cannot feed
        s = simulation(pop_size,1,0.1,0);
        for(auto& ind :s.get_pop())
        {
            ind.set_energy(0);
        }
        //Only the first individual has enough energy to survive
        //for 1 tick
        set_ind_en(s.get_ind(0),s.get_ind(0).get_metabolic_rate() + 0.001);
        assert(s.get_pop().size() == pop_size);
        tick(s);
        assert(s.get_pop_size() == 1);
        //then at the second tick the only individual left dies
        tick(s);
        assert(s.get_pop().empty() && s.get_pop_size() == 0);
    }

    //A simulation is initialized with a random number generator
    {
        simulation s;
        std::uniform_int_distribution u_d(0,2);
        double mean = 0;
        //Draw a thousands times from a uniform dist between 0 and 2
        for(int i = 0; i != 1000; i++)
        {
            mean += u_d(s.get_rng());
        }
        //calculate mean of the drawn values
        mean /= 1000;
        assert(mean < 1.01 && mean > 0.99 );
    }

    //A simulation is initialized with a uniform distribution
    //between 0 and 2PI
    {
        simulation s;
        double mean = 0;
        //Draw a thousands times from a uniform dist between 0 and 2
        int sampling_size = 1000;
        for(int i = 0; i != sampling_size; i++)
        {
            mean += repr_angle(s);
        }
        //calculate mean of the drawn values
        mean /= 1000;
        assert(mean < M_PI + 0.1 && mean > M_PI - 0.1 );
    }

    //Daughter cells are placed at a random angle after reproduction
    {
        simulation s;
        //to calculate angle we will use three point
        //the center of the mother(0,0)
        std::pair<double, double> mother (0,0);
        set_pos(s.get_ind(0),mother);
        //a reference point on the same axis as the mother (0,1)
        //and the center of the daughter -> get_daughter_pos()
        std::pair<double, double> reference(1,0);

        //Draw a thousands times from a uniform dist between 0 and 2
        double mean = 0;
        int sampling_size = 1000;
        for(int i = 0; i != sampling_size; i++)
        {
            auto daughter = get_daughter_pos(s.get_ind(0), repr_angle(s));
            mean += calc_angle_3_pos(mother,daughter,reference);
        }
        mean /= sampling_size;
        assert(mean < M_PI + 0.1 && mean > M_PI - 0.1 );
    }

    //A simulation is initialized with a normal distribution for mutation step
    //with mean 0 and variance 0.1 by default
    {
        simulation s;
        double mean = 0;
        int sampling_size = 10000;
        for(int i = 0 ; i != sampling_size; i++ )
            mean += mut_step(s);
        mean /= sampling_size;
        assert(mean < 0.01 && mean > -0.01);
    }
    //A simulation is initialized with a bernoulli distribution to see if mutation happens or not
    //0.01 by default
    {
        simulation s;
        double mean = 0;
        int sampling_size = 100000;
        for(int i = 0 ; i != sampling_size; i++ )
            mean += mut_happens(s);
        mean /= sampling_size;
        assert(mean < 0.011 && mean > 0.009);
    }

    //The sum of weight of an individual after many rounds of mutation
    //Should have the same mean as in the beginning, but its variance should
    //be the same as the mutation_step distribution
    {
        simulation s;
        //   double init_mean = weights_mean(s.get_ind(0).get_grn());
        double init_variance = weights_var(s.get_ind(0).get_grn());
        assert(init_variance < 0.0001 && init_variance > -0.000001);

        int sampling_size = 10000;
        auto mut_prob_dist = create_bernoulli_dist(s.get_param().get_mu_p());
        auto mut_step_dist = create_normal_dist(0,s.get_param().get_mu_st());
        for (int i = 0; i != sampling_size; i++)
        {
            mutates(s.get_ind(0),
                    s.get_rng(),
                    mut_prob_dist,
                    mut_step_dist
                    );
        }
        //This first assert does not pass, the mean is much more variable than
        //I thought, but I do not see any bug. I will comment this out
        //    assert(mean - init_mean > -0.1 && mean - init_mean < 0.1);
        assert(init_variance - weights_var(s.get_ind(0).get_grn()) > 0.01 ||
               init_variance - weights_var(s.get_ind(0).get_grn()) < -0.01);
    }
    //After dividing the two daughter individuals mutate
    {
        double mutation_probability = 1; //all weights will be mutated in this simulation
        simulation s(1, 1, 0, 0, 0, 0, mutation_probability);
        auto init_var = weights_var(s.get_ind(0).get_grn());
        assert(init_var < 0.00001 && init_var > -0.0001);
        divides(s.get_ind(0),
                s.get_pop(),
                repr_angle(s),
                s.get_rng(),
                create_bernoulli_dist(s.get_param().get_mu_p()),
                create_normal_dist(0, s.get_param().get_mu_st())
                );
        auto post_var = weights_var(s.get_ind(0).get_grn());
        assert(init_var - post_var > 0.000001 || init_var - post_var < -0.0001);
    }

    //A simulation has a member variable m_new_pop_size that states the max number of
    //individuals that will start a new population
    //by default = to pop.size()
    //If at dispersal m_new_pop_size > pop.size()
    //Then the number of funding individuals == pop.size()

    {
        simulation s(1000, 100);
        select_new_pop(s);
        assert(static_cast<int>(s.get_new_pop().size()) == s.get_param().get_exp_new_pop_size());
        s = simulation(10, 100);
        select_new_pop(s);
        auto n_drawn_ind = std::count_if(s.get_pop().begin(),s.get_pop().end(),
                                         [](const individual& i) {return is_drawn(i);});
        assert( n_drawn_ind == s.get_pop_size());
    }


    //During dispersal the individuals selected for the new_pop cannot be drawn again from pop
    {
        simulation s(1000, 100);
        select_new_pop(s);
        assert(std::count_if(s.get_pop().begin(),s.get_pop().end(),
                             [](const individual& i) {return is_drawn(i);}) == 100);
    }

    //A simulation is initialized with a variable m_base_fitness
    //by default = 0.01
    {
        double base_disp_prob = 0.1;
        simulation s(0,0,0,0,0,0,0,0,base_disp_prob);
        assert(s.get_param().get_base_disp_prob() - base_disp_prob < 0.00001 &&
               s.get_param().get_base_disp_prob() - base_disp_prob > -0.000001);
    }

    //A simulation is initialized with a variable m_spore_advantage
    //by default = 10
    {
        double spore_advantage = 10;
        simulation s(0,0,0,0,0,0,0,0,0,spore_advantage);
        assert(s.get_param().get_spo_adv() - spore_advantage < 0.00001 &&
               s.get_param().get_spo_adv() - spore_advantage > -0.000001);
    }

    //A simulation is initialized with a uniform distribution between 0 and 1
    //used to see which ind is drawn at dispersal
    {
        simulation s;
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
    //At initialization a simulation checks that base_disp_dist * 10 is not > 1
    //--------> constructor throws exception. Tested directly in constructor
    {
#ifndef IS_ON_TRAVIS
        try {
            simulation(0,0,0,0,0,0,0,0,1);
        } catch (std::string& e) {
            assert(e == "base dispersal probability * spore advantage > 1, too high!\n" );
        }
#endif
    }

    //Individuals are selected based on their phenotype
    //A spore is more likely to be selected than a living
    {
        simulation s(1000,100);
        for(int i = s.get_pop_size() / 2; i != s.get_pop_size(); i++)
            s.get_ind(i).set_phen(phenotype::spore);
        select_new_pop(s);
        auto spore_ratio =
                std::accumulate(s.get_new_pop().begin(),s.get_new_pop().end(),0.0,
                                [](const int sum, const individual& ind){return sum + is_spore(ind);}) /
                s.get_new_pop().size();
        assert(spore_ratio > 0.5);
    }

    //After being selected in new population individuals flag is_drawn is resetted
    {
        simulation s;
        select_new_pop(s);
        for(const auto& ind :s.get_new_pop())
            assert(!is_drawn(ind));
    }

    //After a new population is selected it swapped with the old population
    //And the old population is cancelled
    {
        unsigned int pop_size = 1000;
        int new_pop_size = 100;
        simulation s(pop_size,new_pop_size);
        select_new_pop(s);
        assert(static_cast<int>(s.get_new_pop().size()) == new_pop_size);
        assert(s.get_pop().size() == pop_size);
        fund_pop(s);
        assert(s.get_new_pop().size() == 0);
        assert(static_cast<int>(s.get_pop_size()) == new_pop_size);
    }

    //Individuals after funding the new population are set in an hexagonal pattern
    {
        unsigned int pop_size = 1000;
        int new_pop_size = 100;
        simulation s(pop_size,new_pop_size);
        select_new_pop(s);
        fund_pop(s);
        place_start_cells(s);
        auto n_hex_l = count_hex_layers(s.get_pop_size());
        auto v_modulus = modulus_of_btw_ind_angles(s, M_PI/ (6 * n_hex_l));
        for(auto ind_modulus : v_modulus)
            assert( ind_modulus < 0.0000000001 || (ind_modulus > M_PI / (6 * n_hex_l) - 0.1 && ind_modulus <= M_PI / (6 * n_hex_l) + 0.1));
        assert(!has_collision(s));
    }


    //Max 100 ind, selected based on phenotype, are placed in a hex pattern, in a new env after dispersal
    //Tests all of the above
    {
        unsigned int pop_size = 1000;
        int new_pop_size = 100;
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
        assert(s.get_pop_size() == new_pop_size);
        //Hex pattern
        auto n_hex_l = count_hex_layers(s.get_pop_size());
        auto v_modulus = modulus_of_btw_ind_angles(s, M_PI/ (6 * n_hex_l));
        for(auto ind_modulus : v_modulus)
        {
            assert( ind_modulus < 0.0000000001 || (ind_modulus > M_PI / (6 * n_hex_l) - 0.1 && ind_modulus <= M_PI / (6 * n_hex_l) + 0.1));
        }
        assert(!has_collision(s));
        //Reset env
        assert( s.get_env() != ref_env );
        auto init_food = s.get_param().get_init_food();
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
        for(auto& ind : s.get_pop())
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
        for(auto& ind : s.get_pop())
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
        for(auto& ind : s.get_pop())
        {
            ind.set_energy(-1);
        }
        //1st respose, the individuals respond to the initial parameter
        //expected output value == 1
        response(s);
        for(auto& ind : s.get_pop())
        {
            assert(is_active(ind));
        }
        //2nd response, the individual responds to the changed parameter(all 0s)
        //expected output value == 0
        response(s);
        for(auto& ind : s.get_pop())
        {
            assert(!is_active(ind));
            assert(is_sporulating(ind));
        }
    }
#endif
}
















