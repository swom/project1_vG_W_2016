#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include <random>
#include "individual.h"
#include "environment.h"

class simulation
{
public:
  simulation(int pop_size = 1, int exp_new_pop_size = 1,double min_dist = 0.1, int grid_side = 1,
             double diff_coeff = 0.1, double init_food = 1, double mutation_prob = 0.01, double mutation_step = 0.1,
             double base_disp_prob = 0.01, double spore_advantage = 10);

  ///Returns the value of the variable m_base_fitness that indicates
  /// the basal fitness/dispersal probability of an individual
  double get_base_disp_prob() const noexcept {return m_base_disp_prob;}

  ///Returns the diffusion coefficent of the environment that will be built in this simulation
  double get_diff_coeff() const noexcept {return m_diff_coeff;}

  ///Returns the reference to the uniform distribution m_disp_dist
  std::uniform_real_distribution<double> get_unif_dist() const noexcept {return m_disp_dist;}

  ///Gets the environment of a simulation
  const environment& get_env() const noexcept {return m_e;}

  ///Gets the reference to environment of a simulation
  environment& get_env() noexcept {return m_e;}

  ///Gets the dimension of the side of the env_grid that is built
  /// in this simulation
  int get_grid_side() const noexcept {return m_grid_side;}

  ///Gets an inidividual at a certain index in the vector
  const individual get_ind(int i) const
  {
    return m_pop[static_cast<unsigned int>(i)];
  }

  ///Gets an inidividual at a certain index in the vector and allows to change it
  individual& get_ind(int i)
  {
    return m_pop[static_cast<unsigned int>(i)];
  }

  ///Gets an individual's energy
  double get_ind_en(int i) const
  {
    return m_pop[static_cast<unsigned int>(i)].get_energy();
  }

  ///Gets an individual's treshold energy
  double get_ind_tr_en(int i) const
  {
    return m_pop[static_cast<unsigned int>(i)].get_treshold_energy();
  }

  ///Gets the position of an individual as a vector x,y
  std::pair<double, double> get_ind_pos(int i);

  ///Gets the initial food that will be provided in each grid_cell of the environment
  double get_init_food() const noexcept {return m_init_food;}

  ///Get minimum distance between individuals at the start of the simulation
  double get_min_dist() const noexcept {return m_min_init_dist_btw_cells;}

  ///Gets the reference to the mutation probability distribution
  std::bernoulli_distribution& get_mu_p() noexcept {return m_mutation_prob;}

  ///Gets the reference to the mutation step distribution
  std::normal_distribution<double>& get_mu_st() noexcept {return m_mutation_step;}

  ///Gets the vector containing all the individuals of the population
  const std::vector<individual>& get_new_pop() const noexcept {return m_new_pop;}

  ///Gets the vector containing all the individuals of the population
  std::vector<individual>& get_new_pop() noexcept {return m_new_pop;}

  ///Gets the size of the population
  int get_new_pop_size() const noexcept {return static_cast<int>(m_new_pop.size());}

  ///Gets the expected size of the new population
  int get_exp_new_pop_size() const noexcept {return m_exp_new_pop_size;}

  ///Gets the size of the population
  int get_pop_size() const noexcept {return static_cast<int>(m_pop.size());}

  ///Gets the vector containing all the individuals of the population
  const std::vector<individual>& get_pop() const noexcept {return m_pop;}

  ///Gets the vector containing all the individuals of the population
  std::vector<individual>& get_pop() noexcept {return m_pop;}

  ///Gets the reference to m_rng
  std::minstd_rand& get_rng() noexcept {return  m_rng;}

  ///Gets reference to m_reproduction_angle
  std::uniform_real_distribution<double>& get_repr_angle() noexcept {return m_reproduction_angle;}

  ///Returns the variable m_spore_advantage that indicates
  ///how many times a spore is more likely to get dispersed
  ///than the other phenotypes
  double get_spo_adv() const noexcept {return m_spore_advantage;}

  ///Gets the number of ticks in the simulation
  const int& get_tick() const noexcept {return m_sim_timer;}

  ///Gives back the mutation size based on initialization parameters
  double mut_step() noexcept {return m_mutation_step(m_rng);}

  ///Calculates if a mutation happens or not
  bool mut_happens() noexcept {return m_mutation_prob(m_rng);}

  ///Returns a random repr angle
  double repr_angle() noexcept {return m_reproduction_angle(m_rng);}

  ///Places an individual of index i at position x,y
  void set_ind_pos(individual& i, double x, double y);

  ///Places an individual of index i at position x,y
  void set_ind_pos(individual& i, std::pair<double, double> pos);

  ///Sets and individual's energy
  void set_ind_en(int i, double en)
  {
    m_pop[static_cast<unsigned int>(i)].set_energy(en);
  }

  ///Updates tick counter by one
  void update_sim_timer() noexcept {m_sim_timer++;}


private:
  ///Made public for testing reasons
  //  ///Returns a random angle
  //  double rnd_angle() noexcept { return m_reproduction_angle(m_rng);}

  std::vector<individual> m_pop;
  int m_exp_new_pop_size;
  std::vector<individual> m_new_pop;
  double m_min_init_dist_btw_cells;
  int m_grid_side; //possibly useless, think of taking out
  double m_diff_coeff; //possibly useless, think of taking out
  double m_init_food; //possibly useless, think of taking out
  environment m_e;
  int m_sim_timer = 0;
  std::minstd_rand m_rng;
  std::uniform_real_distribution<double> m_reproduction_angle;
  std::bernoulli_distribution m_mutation_prob;
  std::normal_distribution<double> m_mutation_step;
  std::uniform_real_distribution<double> m_disp_dist;
  double m_base_disp_prob;
  double m_spore_advantage;
};

///Checks if the entire population has been already drawn for funding the new population
bool all_ind_are_drawn(const simulation& s) noexcept;

///Counts the number of hexagonal layers necessary to place all individuals in hex pattern
int count_hex_layers(int pop_size)  noexcept ;

///Calculate the distance between two cells given their positions
double calculate_distance(std::pair<double, double> lhs, std::pair<double, double> rhs) noexcept;

///Calculates angle between 3 positions(first position is the vertex of the angle)
///returns an angle in radiants betweem 0 and 2PI
double calc_angle_3_pos(std::pair<double, double> P1,
                        std::pair<double, double> P2,
                        std::pair<double, double> P3);

///Calculates the total displacement of each individual in the population if there are collisions
bool calc_tot_displ_pop(std::vector<individual>& pop);

///Calculates the modulus of the angles between all the individuals of a population
///and a given angle in radiants
std::vector<double> modulus_of_btw_ind_angles(simulation& s, double ang_rad);

///Removes dead inidviduals from population vector
void death(simulation& s) noexcept;

///The individuals in the vector are copied at the end of the pop vector
bool division(simulation& s) noexcept;

///Selects a new population and places it in a new environment
void dispersal(simulation& s);

///Displaces the individuals after their total displacement
///has been calculated with calc_tot_disp_pop()
void displace_inds(std::vector<individual>& pop) noexcept;

/// All the individuals feed
void feeding(simulation& s);

///The new population becomes the actual population, the old population is cancelled
void fund_pop(simulation& s) noexcept;

///Returns the iterators to the range of individuals in a population
/// that could possibly be colliding with a focal individual along the x axis
std::pair<std::vector<individual>::iterator,std::vector<individual>::iterator>
possible_collisions_x(individual focal_ind, std::vector<individual> &pop);

///Returns the iterators to the range of individuals that could
/// collide with a focal individual on the x_axis and as well on the y axis
/// !!! It changes the order of the elemnts in pop vector
std::pair<std::vector<individual>::iterator,std::vector<individual>::iterator>
possible_collisions_y(individual focal_ind,
                      std::vector<individual>::iterator first_x,
                      std::vector<individual>::iterator last_x,
                      std::vector<individual> &pop);

///finds the indexes of the individuals ready to divide
std::vector<int> get_dividing_individuals(const simulation& s) noexcept;

///Gets the excess energy from a referenced vector of index of individuals in population vector
std::vector<double> get_excess_energies(const simulation& s) noexcept;

///Gets distance in vector elements between two duaghters of same mother cell
std::vector<std::pair<int, int> > get_sisters_index_offset(const simulation& s) noexcept;

///Checks if any two cells are colliding, return an empty vector in this case
/// or a vecto of the indexes of the first two colliding cells
bool has_collision(simulation &s);

/// Displaces colliding cells so that they do not collide anymore
int manage_static_collisions(simulation& s);

///Moves individuals that are perfectly on top of each other a little bit
///to allow correct displacement
void no_complete_overlap(simulation& s) noexcept;

///Check that in a certain range there is only one individual
bool only_one_focal(const std::vector<individual>::iterator first,
                    const std::vector<individual>::iterator last);

///All inidviduals lose energy due to metabolism
void metabolism_pop(simulation& s);

/// Places cells in an hexagonal grid
/// https://stackoverflow.com/questions/14280831/algorithm-to-generate-2d-magic-hexagon-lattice
void place_start_cells(simulation& s) noexcept;

///Resets the drawn flag for all individuals in the new_population vector
void reset_drawn_fl_new_pop(simulation& s) noexcept;

///Resets the environment in the simulation to be identical to the previous one at the start
void reset_env(simulation& s);

///Draws a 100 individual to fund the new population and puts them in m_new_pop
void select_new_pop(simulation& s);

///Runs all the necessary actions for a timestep to happen
int tick(simulation& s);

///Runs a simulation or a given amount of time
void exec(simulation& s, int n_tick) noexcept;

///Sorts individuals in a given range of a vector by increasing x coordinate
void sort_inds_by_x_inc(std::vector<individual>::iterator start, std::vector<individual>::iterator end);

void test_simulation();

#endif // SIMULATION_H
