#include "individual.h"
#include <cassert>
#include <cmath>

individual::individual(ind_param ind_parameters,
                       double x_pos,
                       double y_pos,
                       double energy,
                       phenotype phenotype,
                       int sporulation_timer):
    m_ind_param{ind_parameters},
    m_x{x_pos},
    m_y{y_pos},
    m_energy{energy},
    m_phenotype{phenotype},
    m_radius{m_ind_param.get_base_radius()},
    m_sporulation_timer{sporulation_timer}
{

}

bool operator==(const individual& lhs, const individual& rhs)
{
    auto grn = lhs.get_grn() == rhs.get_grn();
    auto par = lhs.get_param() == rhs.get_param();
    auto phen = lhs.get_phen() == rhs.get_phen();
    auto x = lhs.get_x() - rhs.get_x() < 0.00001 &&
            lhs.get_x() - rhs.get_x() > -0.00001;
    auto y = lhs.get_y() -  rhs.get_y() < 0.00001 &&
            lhs.get_y() - rhs.get_y() > -0.00001;
    return grn && par && phen && x && y;
}

bool operator!=(const individual& lhs, const individual& rhs)
{
    return !(lhs == rhs);
}

void active_metabolism(individual& i) noexcept
{
    if(i.get_energy() >= i.get_param().get_metabolic_rate())
    {i.change_en(-i.get_param().get_metabolic_rate());}
    else
    {i.set_energy(0);}

    if(i.get_phen() == phenotype::sporulating )
    {
        sporulation(i);
    }

}

bool are_colliding(individual &lhs, individual &rhs, double wiggle_room) noexcept
{
    auto radii_sum = lhs.get_radius() + rhs.get_radius();
    const double sqrd_actual_distance = squared_distance(lhs,rhs) ;
    //Almost never called but to make sure we do not divide by 0
    //In  get_displacement()
    if(sqrd_actual_distance < wiggle_room)
    {
        std::minstd_rand rng;
        auto rnd_n_0_1 =static_cast<double>(std::uniform_real_distribution(0.f,1.f)(rng));
        lhs.change_x(rnd_n_0_1 * 0.1);
        rhs.change_x(rnd_n_0_1 * -0.1);
    }
    return sqrd_actual_distance + 0.000001 <
            radii_sum * radii_sum - wiggle_room;
}

void add_displacement(individual& lhs, individual& rhs, double wiggle_room) noexcept
{
    auto displacement = get_displacement(lhs,rhs, wiggle_room);
    lhs.x_displacement(displacement.first);
    lhs.y_displacement(displacement.second);
}

double distance(const individual& lhs, const individual& rhs) noexcept
{
    return sqrt((lhs.get_x() - rhs.get_x()) * (lhs.get_x() - rhs.get_x())
                + (lhs.get_y() - rhs.get_y()) * (lhs.get_y() - rhs.get_y()));
}

void determine_phen(individual& i) noexcept
{
    if(will_sporulate(i) && !is_sporulating(i))
    {
        starts_sporulation(i);
        assert(i.get_spo_timer() == 0);
    }
    else if(!will_sporulate(i) && is_sporulating(i))
    {
        reverts(i);
        assert(i.get_spo_timer() == 0);
    }
}



void divides(individual& i,
             std::vector<individual>& pop,
             double repr_angle,
             std::minstd_rand& rng,
             std::bernoulli_distribution mu_p,
             std::normal_distribution<double> mu_st)
{
    assert( i.get_phen() == phenotype::active);
    double offs_init_en = i.split_excess_energy();
    i.set_energy(offs_init_en);
    individual daughter = i;

    update_radius(i);
    update_radius(daughter);

    mutates(i, rng, mu_p, mu_st);
    mutates(daughter, rng, mu_p, mu_st);
    set_pos(daughter,calculate_daughter_pos(daughter, repr_angle));
    pop.push_back(daughter);

}

void draw(individual& i)
{
    assert(!is_drawn(i));
    i.set_drawn_flag(true);
}

void draw_flag_reset(individual& i)
{
    assert(is_drawn(i));
    i.set_drawn_flag(false);
}

void feed(individual& i, env_grid_cell& c) noexcept
{
    if(c.get_food() > i.get_param().get_uptake_rate())
    {
        c.set_food_change(- i.get_param().get_uptake_rate());
        c.increment_food();
        i.change_en(i.get_param().get_uptake_rate());
    }
    else
    {
        i.change_en(c.get_food());
        c.set_food(0);
    }
}

void jordi_feed(individual& i, env_grid_cell& c) noexcept
{

    auto food_moved = i.get_param().get_uptake_rate() * c.get_food();
    c.set_food_change(-food_moved);
    c.increment_food();
    i.change_en(food_moved);

}


int find_grid_index(const individual& i, double grid_side)
{
    auto x_offset = i.get_x() + grid_side/2;
    auto y_offset = i.get_y() + grid_side/2;
    //cast correctly negative numbers so to keep them outside the grid
    int x_index_offset = x_offset < 0 ?
                static_cast<int>(x_offset - 1) : static_cast<int>(x_offset);
    int y_index_offset = y_offset < 0 ?
                static_cast<int>(y_offset - 1) : static_cast<int>(y_offset);

    if(x_offset >= grid_side || x_offset < 0
            || y_offset >= grid_side || y_offset < 0)
    {return  -100;}
    return x_index_offset + y_index_offset * static_cast<int>(grid_side);
}

std::pair<double,double> calculate_daughter_pos(individual& i, double rnd_angle) noexcept
{
    std::pair<double, double> pos;
    pos.first = i.get_x() + cos(rnd_angle) * 2 * i.get_radius();
    pos.second += i.get_y() + sin(rnd_angle) * 2 * i.get_radius();
    return pos;
}

double get_fitness(const individual& i, double disp_prob, double spore_advantage) noexcept
{
    if(i.get_phen() == phenotype::spore)
        return disp_prob * spore_advantage;
    return disp_prob;
}

std::pair<double, double> get_pos(individual& i)  noexcept
{
    std::pair<double, double> pos;
    pos.first = i.get_x();
    pos.second = i.get_y();
    return pos;
}

std::pair<double, double> get_displacement(const individual& lhs, const individual& rhs, double wiggle_room) noexcept
{
    std::pair<double, double> displ_x_y;
    auto dist = distance(lhs, rhs);

    auto half_overlapping = half_overlap(lhs,rhs, wiggle_room);
    displ_x_y.first = half_overlapping * (lhs.get_x() - rhs.get_x()) / dist;
    displ_x_y.second = half_overlapping * (lhs.get_y() - rhs.get_y()) / dist;
    return  displ_x_y;
}

bool is_dead(individual const&  i) noexcept
{ if(i.get_phen() == phenotype::spore)
    {return false;}
    return i.get_energy() <= 0;
}

bool is_drawn(const individual& i) noexcept
{
    return i.get_drawn_flag();
}

bool is_active(const individual& i) noexcept
{
    return i.get_phen() == phenotype::active;
}

bool is_sporulating(const individual& i) noexcept
{
    return i.get_phen() == phenotype::sporulating;
}

bool is_spore(const individual& i) noexcept
{
    return i.get_phen() == phenotype::spore;
}

double half_overlap(const individual& lhs, const individual& rhs , double wiggle_room) noexcept
{
    auto dist = distance(lhs,rhs);
    auto sum_of_radiuses = lhs.get_radius() + rhs.get_radius();

    //if circles are one into the other
    if (dist < std::max(lhs.get_radius(), rhs.get_radius()))
    {
        return (std::max(lhs.get_radius(),rhs.get_radius()) -
                dist +
                std::min(lhs.get_radius(),rhs.get_radius())) / 2;
    }

    //if circles partially overlap or
    //their dist is equal to sum of radiuses
    //(in which case they should not pass through this function)
    return wiggle_room + (sum_of_radiuses - dist)  / 2;
}

bool have_same_ancestor(const individual& lhs, const individual& rhs) noexcept
{
    return lhs.get_ancestor() == rhs.get_ancestor();
}

void mutates(individual& i, std::minstd_rand& rng,
             std::bernoulli_distribution& mu_p,
             std::normal_distribution<double> mu_st) noexcept
{
    mutation(i.get_grn(), rng, mu_p, mu_st);
}

void responds(individual& i, const env_grid_cell& c)
{
    sense(i,c);
    jordi_response_mech(i.get_grn());
    determine_phen(i);
}

void reverts(individual& i) noexcept
{
    assert(i.get_phen() == phenotype::sporulating);
    i.set_phen(phenotype::active);
    i.reset_spo_timer();
}

void secretes_metab(const individual& i, env_grid_cell& c)
{
    assert(c.get_metab_change() < 0.000001 && c.get_metab_change() > -0.00001);
    c.set_metab_change(i.get_param().get_metab_secr_rate());
    c.increment_metabolite();
    c.set_metab_change(0);
}

void sense(individual& i, const env_grid_cell& c)
{
    auto& inputs = i.get_grn().get_input_nodes();
    inputs[0] = c.get_food();
    inputs[1] = c.get_metab();
    inputs[2] = i.get_energy();
}


void set_pos(individual& i, std::pair<double, double> pos)  noexcept
{
    i.set_x(pos.first);
    i.set_y(pos.second);
}

void sporulation(individual& i) noexcept
{

    assert(i.get_phen() == phenotype::sporulating);
    i.tick_spo_timer();
    assert(i.get_spo_timer() <= i.get_param().get_transformation_time());
    if(i.get_spo_timer() == i.get_param().get_transformation_time())
    {
        i.set_phen(phenotype::spore);
        i.reset_spo_timer();
    }

}

void sporulating_metabolism(individual& i) noexcept
{

    if(i.get_energy() >= i.get_param().get_spor_metabolic_rate())
    {i.change_en(-i.get_param().get_spor_metabolic_rate());}
    else
    {i.set_energy(0);}
    sporulation(i);

}

double squared_distance(const individual& lhs, const individual& rhs) noexcept
{
    return (lhs.get_x() - rhs.get_x()) * (lhs.get_x() - rhs.get_x())
            + (lhs.get_y() - rhs.get_y()) * (lhs.get_y() - rhs.get_y());
}

void starts_sporulation(individual& i)
{
    assert(is_active(i));
    assert(i.get_spo_timer() == 0);
    i.set_phen(phenotype::sporulating);
}

void update_radius(individual& i)
{
    i.get_radius() = i.get_param().get_base_radius()/* +
            i.get_param().get_uptake_mean() * i.get_energy()*/;
}

bool will_sporulate(individual& i) noexcept
{
    return i.get_grn().get_output_spo() == 0;
}


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
