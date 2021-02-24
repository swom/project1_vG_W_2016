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

