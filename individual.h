

 #ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include "env_grid_cell.h"
#include "grn.h"
#include "ind_param.h"
#include "phenotype.h"

#include <vector>
#include <utility>
#include <algorithm>
#include <random>
#include <cassert>

class individual
{
public:

    individual(double x_pos = 0, double y_pos = 0, double radius  = 0.5,
               double energy = 0.1, double treshold_energy = 3,
               double uptake_rate = 0.1, double metabolic_rate = 0.01,
               phenotype phenotype = phenotype::active,
               int sporulation_timer = 0, int transformation_time = 5,
               double wiggle_room = 0.01, double metab_secretion_rate = 0.1,
               double spor_metab_rate = 0.5);

    individual(ind_param ind_param,
               double x_pos = 0,
               double y_pos = 0,
               double energy = 0.1,
               phenotype phenotype = phenotype::active,
               int sporulation_timer = 0);


    ///Turns the flag (signalling that this is the focal individual during collision check
    void becomes_focal() noexcept {m_is_focal = true;}

    ///Changes x of an individual
    void change_x(double x) noexcept {m_x += x;}

    ///Changes y of an individual
    void change_y(double y) noexcept {m_y += y;}

    ///Displaces an individual based on its m_x_displacement and m_y_displacement values
    void displace() noexcept { change_x(m_x_displacement); change_y(m_y_displacement); reset_displacement();}

    ///Gets the flag that signal if the individual has been drawn to be member of the new pop
    bool get_drawn_flag() const noexcept {return m_is_drawn;}

    ///Changes the energy level of an individual
    void change_en(double en_change) noexcept {m_energy += en_change;}

    ///gets energy of individual
    double get_energy() const noexcept {return m_energy;}

    ///Gets the const ref to the m_grn
    const GRN& get_grn() const noexcept {return m_grn;}

    ///Gets the ref to the m_grn
    GRN& get_grn() noexcept {return m_grn;}

    ///Gets ref to the parametersof the individual
    ind_param& get_param() noexcept {return m_ind_param;}

    ///Gets const ref to the parametersof the individual
    const ind_param& get_param() const noexcept {return m_ind_param;}

    ///Gets the type of the individual
    phenotype get_phen() const noexcept {return m_phenotype;}

    ///gets x position of the individual
    double get_x() const noexcept {return m_x;}

    ///Gets the amount of displacement that will be applied on the x coordinate
    double get_x_displacement() const noexcept {return m_x_displacement;}

    ///Gets the amount of displacement that will be applied on the x coordinate
    double get_y_displacement() const noexcept {return m_y_displacement;}

    ///gets y position of the individual
    double get_y() const noexcept {return m_y;}

    ///Gets the sporulation timer
    int get_spo_timer() const noexcept {return m_sporulation_timer;}

    ///Returns the flag signalling if this is the focal individual during the collision check
    bool is_focal() const noexcept {return m_is_focal;}

    ///Resets the flag signalling if this is the focal individual in collision check to flase
    void no_more_focal() noexcept {
#ifndef NDEBUG
        assert(m_is_focal);
#endif
        m_is_focal = false;}

    ///resets the values of m_x_displacement and m_y_displacement to 0;
    void reset_displacement() {m_x_displacement = 0; m_y_displacement = 0;}

    ///Resets the sporulation timer
    void reset_spo_timer() noexcept {m_sporulation_timer = 0;}

    ///Sets the m_is_drawn flag to a bool value
    void set_drawn_flag(bool b) noexcept {m_is_drawn = b;}

    ///sets the energy of an individual
    void set_energy(double new_energy) {m_energy = new_energy;}

    ///Sets the type of an individual
    void set_phen(phenotype type) {m_phenotype = type;}

    ///Sets the x of an individual
    void set_x(double x) noexcept {m_x = x;}

    ///Sets the y of an individual
    void set_y(double y) noexcept {m_y = y;}

    ///Ticks the sporulation timer by one
    void tick_spo_timer() noexcept {m_sporulation_timer++;}

    ///Splits the excess energy not required for division in two9to be then
    ///(to be then assigned to the two daughter cells by
    /// simulation::reproduce/cells_divide)
    double split_excess_energy() const noexcept {return (m_energy - m_ind_param.get_treshold_energy())/2;}

    ///Changes x of an individual
    void x_displacement(double x_displacement) noexcept {m_x_displacement += x_displacement;}

    ///Changes y of an individual
    void y_displacement(double y_displacement) noexcept {m_y_displacement += y_displacement;}

private:

    ///The parameters of the individual
    ind_param m_ind_param;

    ///X coord of an individual
    double m_x;

    ///Y coord of an individual
    double m_y;

    ///Change that will be applued to the X coord
    double m_x_displacement = 0;

    ///Change that will be applued to the Y coord
    double m_y_displacement = 0;

    ///Level of internal energy gained from feeding
    double m_energy;

    ///The gene regulatory network of an individual,
    ///driving the decision of which phenotype to assume
    GRN m_grn;

    ///Flag that signals if an individual has already been drawn to be part
    /// of the new population
    bool m_is_drawn = false;

    ///Flag signalling if this is the focal individual
    ///  during loop of collision detection
    bool m_is_focal = false;

    ///The phenotype of the individual
    phenotype m_phenotype;

    ///The amount of timesteps an individual has been sporulating
    int m_sporulation_timer;
};

///Returns true if two individuals are in the same position
bool operator==(const individual& lhs, const individual& rhs);

///Returns true if two individuals are not in the same position
bool operator!=(const individual& lhs, const individual& rhs);

///Checks if two individuals are colliding
bool are_colliding(individual &lhs, individual &rhs) noexcept;

///Calculates how much individuals need to be displaced to not overlap
std::pair<double, double> get_displacement(const individual &lhs, const individual &rhs) noexcept;

///Adds a x and y displacement component to lhs and rhs displacement so that lhs and rhs do not overlap
void add_displacement(individual& lhs, individual& rhs) noexcept;

/// Finds the distance between two individuals
double distance(const individual& lhs, const individual& rhs) noexcept;

///Find the squared distance between 2 inds, faster than distance
double squared_distance(const individual& lhs, const individual& rhs) noexcept;

///Divides an individual at a given index
void divides(individual &i, std::vector<individual> &pop, double repr_angle,
             std::minstd_rand& rng, std::bernoulli_distribution mu_p,
             std::normal_distribution<double> mu_st);

///Sets the type of the individual accordingly to its GRN's outputs
void determine_phen(individual& i) noexcept;

///Signal that an individual has been drawn to be a funder of the new population
void draw(individual& i);

///Resets the flag of a drawn individual to 0
void draw_flag_reset(individual& i);

///Finds the grid cell where an individual is on,
///given the side of the environment grid
///the grid_cell at grid[0] is at position(x,y) = (-grid_side/2, -grid_side/2)
int find_grid_index(const individual &i, double grid_side) ;

///An individual increases its energy depleting food
void feed(individual& i, env_grid_cell& food) noexcept;

///returns the position of the second daughter cell just outside the mother
std::pair<double, double> get_daughter_pos(individual& i, double rnd_angle) noexcept;

///Returns the fitness of an individual based on its phenotype
///And the base_fitness declared in simulation
double get_fitness(const individual& i, double base_disp_prob, double spore_advantage) noexcept;

///Gets the x,y coordinates as a pair
std::pair<double, double> get_pos(individual& i) noexcept;

///Signals if an individual has to be destroyed
bool is_dead(const individual &i) noexcept;

///Checks if an individual has been drawn for becoming part of the new population
bool is_drawn(const individual& i) noexcept;

///Checks if an individual is active
bool is_active(const individual& i) noexcept;

///Checks if an individual is sporulating
bool is_sporulating(const individual& i) noexcept;

///Checks if an individual is a spore
bool is_spore(const individual& i) noexcept;

///Individuals lose energry due to metabolism
void active_metabolism(individual& i) noexcept;

///Mutates the GRN of an individual
/// ///For now requires to get distribution and rng from somewhere else
/// (simulation)
void mutates(individual& i, std::minstd_rand& rng,
             std::bernoulli_distribution& mu_p,
             std::normal_distribution<double> mu_st) noexcept;

///Finds the overlap between two individuals
double half_overlap(const individual& lhs, const individual& rhs) noexcept;

///Reverts a sporulating individual back to living (and resets the timer)
void reverts(individual& i) noexcept;

///An individual senses cues, internal and on its grid_cell
///  and responds to them
void responds(individual& i, const env_grid_cell& c);

///Adds metabolite to the cell where the individual's center is on
void secretes_metab(const individual& i, env_grid_cell& c);

///Takes food, metabolite from the grid cell and energy from individual
/// and sets them as inputs of the grn
void sense(individual& i, const env_grid_cell& c);

///Sets x,y given a pair of doubles(x,y)
void set_pos(individual& i, std::pair<double, double> pos)  noexcept;

///Changes an ind_type of an individual from living to sporulating
void starts_sporulation(individual& i);

///Applies metabolism of sporulating individuals
void sporulating_metabolism(individual& i) noexcept;

///Checks and updates the timer of sporulating individuals
/// and changes them into spores if they have sporulated long enough
void sporulation(individual& i) noexcept;

///Determines if the output of the network will make the individual
/// become of ind_type::sporulating
bool will_sporulate(individual& i) noexcept;

void test_individual();

#endif // INDIVIDUAL_H
