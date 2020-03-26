#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include <random>
#include "individual.h"
#include "environment.h"

class simulation
{
public:
  simulation(int pop_size = 1, int grid_side = 1, double min_dist = 0.1, double starting_food = 1,
             double mutation_prob = 0.01, double mutation_step = 0.1);

  ///finds the indexes of the individuals ready to divide
  std::vector<int> get_dividing_individuals() const noexcept;

  ///Gets the environment of a simulation
  const environment& get_env() const noexcept {return m_e;}

  ///Gets the reference to environment of a simulation
  environment& get_env() noexcept {return m_e;}

  ///Gets an inidividual at a certain index in the vector
  const individual get_ind(int i) const
  {
    return m_population[static_cast<unsigned int>(i)];
  }

  ///Gets an inidividual at a certain index in the vector and allows to change it
  individual& get_ind(int i)
  {
    return m_population[static_cast<unsigned int>(i)];
  }

  ///Gets the excess energy from a referenced vector of index of individuals in population vector
  std::vector<double> get_excess_energies(std::vector<int> &v_ind) const noexcept;

  ///Gets the excess energy for a vector of index of individuals in population vector
  std::vector<double> get_excess_energies(std::vector<int> v_ind) const noexcept;

  ///Gets an individual's energy
  double get_ind_en(int i) const
  {
    return m_population[static_cast<unsigned int>(i)].get_energy();
  }

  ///Gets an individual's treshold energy
  double get_ind_tr_en(int i) const
  {
    return m_population[static_cast<unsigned int>(i)].get_treshold_energy();
  }

  ///Gets the position of an individual as a vector x,y
  const std::pair<double,double> get_ind_pos(int i);

  ///Get minimum distance between individuals at the start of the simulation
  double get_min_dist() const noexcept {return m_min_init_dist_btw_cells;}

  ///Gets the reference to the mutation probability distribution
  std::bernoulli_distribution& get_mu_p() noexcept {return m_mutation_prob;}

  ///Gets the reference to the mutation step distribution
  std::normal_distribution<double>& get_mu_st() noexcept {return m_mutation_step;}

  ///Gets the size of the population
  int get_pop_size() const noexcept {return static_cast<int>(m_population.size());}

  ///Gets the vector containing all the individuals of the population
  const std::vector<individual>& get_pop() const noexcept {return m_population;}

  ///Gets the vector containing all the individuals of the population
  std::vector<individual>& get_pop() noexcept {return m_population;}

  ///Gets the reference to m_rng
  std::minstd_rand& get_rng() noexcept {return  m_rng;}

  ///Gets reference to m_reproduction_angle
  std::uniform_real_distribution<double>& get_repr_angle() noexcept {return m_reproduction_angle;}

  ///Gets distance in vector elements between two duaghters of same mother cell
  std::vector<std::pair<int, int> > get_sisters_index_offset() const noexcept;

  ///Gets the number of ticks in the simulation
  const int& get_tick() const noexcept {return m_sim_timer;}

  ///Divides an individual at a given index
  void divide_ind(individual &i);

  ///The individuals in the vector are copied at the end of the pop vector
  // To be chaged to take vec of references as an argument
  void do_division(std::vector<int> dividing_individuals);

  ///Individuals with enough energy divide and redistribute energy
  void do_reprduction() noexcept;

  ///Gives back the mutation size based on initialization parameters
  double mut_step() noexcept {return m_mutation_step(m_rng);}

  ///Calculates if a mutation happens or not
  bool mut_happens() noexcept {return m_mutation_prob(m_rng);}

  /// Places cells in an hexagonal grid
  /// https://stackoverflow.com/questions/14280831/algorithm-to-generate-2d-magic-hexagon-lattice
  void place_start_cells() noexcept;

  ///Returns a random repr angle
  double rnd_repr_angle() noexcept {return m_reproduction_angle(m_rng);}

  ///Places an individual of index i at position x,y
  void set_ind_pos(individual& i, double x, double y);

  ///Places an individual of index i at position x,y
  void set_ind_pos(individual& i, const std::pair<double, double> pos);

  ///Sets and individual's energy
  void set_ind_en(int i, double en)
  {
    m_population[static_cast<unsigned int>(i)].set_energy(en);
  }

  ///Updates tick counter by one
  void update_sim_timer() noexcept {m_sim_timer++;}


private:
///Made public for testing reasons
//  ///Returns a random angle
//  double rnd_angle() noexcept { return m_reproduction_angle(m_rng);}

  std::vector<individual> m_population;
  double m_min_init_dist_btw_cells;
  environment m_e;
  int m_sim_timer = 0;
  std::minstd_rand m_rng;
  std::uniform_real_distribution<double> m_reproduction_angle;
  std::bernoulli_distribution m_mutation_prob;
  std::normal_distribution<double> m_mutation_step;
};

///Counts the number of hexagonal layers necessary to place all individuals in hex pattern
int count_hex_layers(int pop_size)  noexcept ;

///Calculate the distance between two cells given their positions
double calculate_distance(std::pair<double, double> lhs, std::pair<double, double> rhs) noexcept;

///Calculates angle between 3 positions(first position is the vertex of the angle)
///returns an angle in radiants betweem 0 and 2PI
double calc_angle_3_pos(std::pair<double, double> P1,
                       std::pair<double, double> P2,
                       std::pair<double, double> P3);

///Calculates the modulus of the angles between all the individuals of a population
///and a given angle in radiants
std::vector<double> modulus_of_btw_ind_angles(simulation& s, double ang_rad);

///Removes dead inidviduals from population vector
void death(simulation& s) noexcept;

/// All the individuals feed
void feeding(simulation& s);

///Checks if any two cells are colliding
bool has_collision(const simulation& s);

/// Displaces colliding cells so that they do not collide anymore
void manage_static_collisions(simulation& s);

///All inidviduals lose energy due to metabolism
void metabolism_pop(simulation& s);

///Runs all the necessary actions for a timestep to happen
void tick(simulation& s);

///Runs a simulation or a given amount of time
void exec(simulation& s, int n_tick) noexcept;

void test_simulation();

#endif // SIMULATION_H
