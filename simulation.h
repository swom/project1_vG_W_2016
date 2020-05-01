#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include <random>
#include "environment.h"
#include "individual.h"
#include "parameters.h"

class simulation
{
public:
    simulation(unsigned int pop_size = 1, int exp_new_pop_size = 1, double min_dist = 0.1, int grid_side = 1,
               double diff_coeff = 0.1, double init_food = 1.0, double mutation_prob = 0.01, double mutation_step = 0.1,
               double base_disp_prob = 0.01, double spore_advantage = 10.0, double reproduction_prob = 0.1,
               double metab_degradation_rate = 0.01);

    simulation(sim_param param, ind_param inds_param);

    ///Gets the environment of a simulation
    const environment& get_env() const noexcept {return m_e;}

    ///Gets the reference to environment of a simulation
    environment& get_env() noexcept {return m_e;}

    ///Gets an inidividual at a certain index in the vector
    const individual& get_ind(int i) const
    {
        assert(i >= 0);
        assert(static_cast<unsigned int>(i)< m_pop.size());
        return m_pop[static_cast<unsigned int>(i)];
    }

    ///Gets an inidividual at a certain index in the vector and allows to change it
    individual& get_ind(int i)
    {
        assert(i >= 0);
        assert(static_cast<unsigned int>(i)< m_pop.size());
        return m_pop[static_cast<unsigned int>(i)];
    }

    ///Gets the vector containing all the individuals of the population
    const std::vector<individual>& get_new_pop() const noexcept {return m_new_pop_tmp_buffer;}

    ///Gets the vector containing all the individuals of the population
    std::vector<individual>& get_new_pop() noexcept {return m_new_pop_tmp_buffer;}

     ///Gets const ref to parameters
    const sim_param& get_param() const noexcept {return m_param;}

    ///Gets the vector containing all the individuals of the population
    const std::vector<individual>& get_pop() const noexcept {return m_pop;}

    ///Gets the vector containing all the individuals of the population
    std::vector<individual>& get_pop() noexcept {return m_pop;}

    ///Gets the size of the population
    int get_pop_size() const noexcept {return static_cast<int>(m_pop.size());}

    ///Gets the reference to m_rng
    auto& get_rng() noexcept {return  m_rng;}

    ///Gets the number of ticks in the simulation
    int get_tick() const noexcept {return m_sim_timer;}

    ///Updates tick counter by one
    void update_sim_timer() noexcept {++m_sim_timer;}

private:

    ///The parameters of the simulation
    sim_param m_param;

    ///The vector of individuals representing a population
    std::vector<individual> m_pop;

    ///The buffer vector on which the new population will be copied and that
    /// will be then swapped with the current population
    /// Used to reduce memory allocation time
    std::vector<individual> m_new_pop_tmp_buffer;

    ///The environment class containing the grid where substances(food, metabolite) are
    environment m_e;

    ///The timer that keeps track of how many timesteps we are in the simulation
    int m_sim_timer = 0;

    ///The random number generator of simulation(used for everything)
    std::minstd_rand m_rng;

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

///Creates a uniform distribution
std::uniform_real_distribution<double> create_unif_dist(double a, double b) noexcept;

///Creates a bernoulli distribution
std::bernoulli_distribution create_bernoulli_dist(double p) noexcept;

///Creates a normal distribution
std::normal_distribution<double> create_normal_dist(double m, double v);

///Removes dead inidviduals from population vector
void death(simulation& s) noexcept;

//Pop
///The individuals in the vector are copied at the end of the pop vector
bool division(simulation& s) noexcept;

//Pop
///Selects a new population and places it in a new environment
void dispersal(simulation& s);

//Pop
///Displaces the individuals after their total displacement
///has been calculated with calc_tot_disp_pop()
void displace_inds(std::vector<individual>& pop) noexcept;

///Runs a simulation or a given amount of time
void exec(simulation& s, int n_tick) noexcept;

/// All the individuals feed
void feeding(simulation& s);

//Pop
///The new population becomes the actual population, the old population is cancelled
void fund_pop(simulation& s) noexcept;

//Pop
///finds the indexes of the individuals ready to divide
std::vector<int> get_dividing_individuals(const simulation& s) noexcept;

//Pop
///Gets the excess energy from a referenced vector of index of individuals in population vector
std::vector<double> get_excess_energies(const simulation& s) noexcept;

//Pop
///Gets an individual's energy
double get_ind_en(const simulation& s, int i);


//Pop
///Gets an individual's treshold energy
double get_ind_tr_en(const simulation& s, int i);

//Pop
///Gets the position of an individual as a vector x,y
std::pair<double, double> get_ind_pos(int i);

//Pop
///Gets distance in vector elements between two duaghters of same mother cell
std::vector<std::pair<int, int> > get_sisters_index_offset(const simulation& s) noexcept;

//Pop
///Checks if any two cells are colliding, return an empty vector in this case
/// or a vecto of the indexes of the first two colliding cells
bool has_collision(simulation &s);

//Pop
/// Displaces colliding cells so that they do not collide anymore
int manage_static_collisions(simulation& s);

//Pop
///All inidviduals lose energy due to metabolism
void metabolism_pop(simulation& s);

//Pop
///Calculates the modulus of the angles between all the individuals of a population
///and a given angle in radiants
std::vector<double> modulus_of_btw_ind_angles(simulation& s, double ang_rad);

//Pop
///Checks if a mutation happens given the mutation probability of a simulation
bool mut_happens(simulation& s) noexcept;

//Pop
///Gives back the mutation size based on initialization parameters
double mut_step() noexcept;

//Pop
///Moves individuals that are perfectly on top of each other a little bit
///to allow correct displacement
void no_complete_overlap(simulation& s) noexcept;

//Pop
///Check that in a certain range there is only one individual
bool only_one_focal(const std::vector<individual>::iterator first,
                    const std::vector<individual>::iterator last);

//Pop
/// Places cells in an hexagonal grid
/// https://stackoverflow.com/questions/14280831/algorithm-to-generate-2d-magic-hexagon-lattice
void place_start_cells(simulation& s) noexcept;

//Pop
///Returns the iterators to the range of individuals in a population
/// that could possibly be colliding with a focal individual along the x axis
std::pair<std::vector<individual>::iterator,std::vector<individual>::iterator>
possible_collisions_x(individual focal_ind, std::vector<individual> &pop);

//Pop
///Returns the iterators to the range of individuals that could
/// collide with a focal individual on the x_axis and as well on the y axis
/// !!! It changes the order of the elemnts in pop vector
std::pair<std::vector<individual>::iterator,std::vector<individual>::iterator>
possible_collisions_y(individual focal_ind,
                      std::vector<individual>::iterator first_x,
                      std::vector<individual>::iterator last_x,
                      std::vector<individual> &pop);

///Returns a random angle between 0 and and 2PI (in radians)
double repr_angle(simulation &s) noexcept;

//Pop
///Resets the drawn flag for all individuals in the new_population vector
void reset_drawn_fl_new_pop(simulation& s) noexcept;

///Resets the environment in the simulation to be identical to the previous one at the start
void reset_env(simulation& s);

///All individuals secrete metabolite
void secretion_metabolite(simulation& s);

//Pop
///Draws a 100 individual to fund the new population and puts them in m_new_pop
void select_new_pop(simulation& s);

//Pop
///Places an individual of index i at position x,y
void set_ind_pos(individual& i, double x, double y);

//Pop
///Places an individual of index i at position x,y
void set_ind_pos(individual& i, std::pair<double, double> pos);

//Pop
///Sets and individual's energy
void set_ind_en(individual& i, double en);

//Pop
///Sorts individuals in a given range of a vector by increasing x coordinate
void sort_inds_by_x_inc(std::vector<individual>::iterator start, std::vector<individual>::iterator end);

///Runs all the necessary actions for a timestep to happen
int tick(simulation& s);

void test_simulation();

#endif // SIMULATION_H
