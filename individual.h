#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include "env_grid_cell.h"
#include "phenotype.h"
#include "grn.h"
#include <vector>
#include <utility>
#include <algorithm>
#include <random>

class individual
{
public:

  individual(double x_pos = 0, double y_pos = 0, double radius  = 0.5,
             double energy = 0, double treshold_energy = 2,
             double uptake_rate = 0.1, double metabolic_rate = 0.01,
             phenotype phenotype = phenotype::active,
             int sporulation_timer = 0, int transformation_time = 5);

  ///Changes x of an individual
  void change_x(double x) noexcept {m_x += x;}

  ///Changes y of an individual
  void change_y(double y) noexcept {m_y += y;}

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

  ///gets the metabolic rate of an individual
  double get_metab_rate() const noexcept {return m_metab_rate;}

  ///gets radius of individual
  double get_radius() const noexcept {return m_radius;}

  ///Gets the type of the individual
  phenotype get_phen() const noexcept {return m_individual_type;}

  ///Gets the time it takes to transform into a spore
  int get_transformation_time() const noexcept {return m_transformation_time;}

  ///gets the food uptake_rate of an individual
  double get_uptake_rate() const noexcept {return m_uptake_rate;}

  ///gets x position of the individual
  double get_x() const noexcept {return m_x;}

  ///gets y position of the individual
  double get_y() const noexcept {return m_y;}

  ///gets the treshold energy at which an individual reproduces
  double get_treshold_energy() const noexcept {return m_treshold_energy;}

  ///Gets the sporulation timer
  int get_spo_timer() const noexcept {return m_sporulation_timer;}

  ///Resets the sporulation timer
  void reset_spo_timer() noexcept {m_sporulation_timer = 0;}

  ///Sets the m_is_drawn flag to a bool value
  void set_drawn_flag(bool b) noexcept {m_is_drawn = b;}

  ///sets the energy of an individual
  void set_energy(double new_energy) {m_energy = new_energy;}

  ///Sets the type of an individual
  void set_phen(phenotype type) {m_individual_type = type;}

  ///Sets the x of an individual
  void set_x(double x) noexcept {m_x = x;}

  ///Sets the y of an individual
  void set_y(double y) noexcept {m_y = y;}

  ///Ticks the sporulation timer by one
  void tick_spo_timer() noexcept {m_sporulation_timer++;}

  ///Splits the excess energy not required for division in two9to be then
  ///(to be then assigned to the two daughter cells by
  /// simulation::reproduce/cells_divide)
  double split_excess_energy() const noexcept {return (m_energy - m_treshold_energy)/2;}

private:

  double m_x;
  double m_y;
  double m_radius;
  double m_energy;
  double m_treshold_energy;
  double m_uptake_rate;
  double m_metab_rate;
  phenotype m_individual_type;
  int m_sporulation_timer;
  int m_transformation_time; //number of time steps the individual needs
                             //to go through to sporulate, if it is alive at
                             //m_transformation time + 1 it will become a spore
  GRN m_grn;
  bool m_is_drawn = false;
};


///Checks if two individuals are colliding
bool are_colliding(const individual &lhs, const individual &rhs) noexcept;

///Calculates how much individuals need to be displaced to not overlap
std::pair<double, double> get_displacement(const individual &lhs, const individual &rhs) noexcept;

///Dislpaces two individuals so that they do not overlap
void displace(individual& lhs, individual& rhs) noexcept;

/// Finds the distance between two individuals
double distance(const individual& lhs, const individual& rhs) noexcept;

///Divides an individual at a given index
void divides(individual &i, std::vector<individual> &pop, double repr_angle,
             std::minstd_rand& rng, std::bernoulli_distribution& mu_p,
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
int find_grid_index( individual& i, double grid_side) ;

///An individual increases its energy depleting food
void feed(individual& i, env_grid_cell& food) noexcept;

///returns the position of the second daughter cell just outside the mother
const std::pair<double,double> get_daughter_pos(individual& i, double rnd_angle) noexcept;

///Returns the fitness of an individual based on its phenotype
///And the base_fitness declared in simulation
double get_fitness(const individual& i, double base_disp_prob, double spore_advantage) noexcept;

///Gets the x,y coordinates as a pair
const std::pair<double,double> get_pos(individual& i) noexcept;

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
void metabolism(individual& i) noexcept;

///Mutates the GRN of an individual
/// ///For now requires to get distribution and rng from somewhere else
/// (simulation)
void mutates(individual& i, std::minstd_rand& rng,
            std::bernoulli_distribution& mu_p,
            std::normal_distribution<double> mu_st) noexcept;

///Finds the overlap between two individuals
double overlap(const individual& lhs, const individual& rhs) noexcept;

///Reverts a sporulating individual back to living (and resets the timer)
void reverts(individual& i) noexcept;

///An individual senses cues, internal and on its grid_cell
///  and responds to them
void responds(individual& i, const env_grid_cell& c);

///Takes food, metabolite from the grid cell and energy from individual
/// and sets them as inputs of the grn
void sense(individual& i, const env_grid_cell& c);

///Sets x,y given a pair of doubles(x,y)
void set_pos(individual& i, std::pair<double, double> pos)  noexcept;

///Changes an ind_type of an individual from living to sporulating
void starts_sporulation(individual& i);

///Checks and updates the timer of sporulating individuals
/// and changes them into spores if they have sporulated long enough
void sporulation(individual& i) noexcept;

///Determines if the output of the network will make the individual
/// become of ind_type::sporulating
bool will_sporulate(individual& i) noexcept;

void test_individual();

#endif // INDIVIDUAL_H
