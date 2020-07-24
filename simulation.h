#ifndef SIMULATION_H
#define SIMULATION_H
#include "demographic_sim.h"
#include "environment.h"
#include "funders_success.h"
#include "individual.h"
#include "population.h"
#include "sim_parameters.h"

#include <vector>
#include <random>

class simulation
{
public:

    simulation(sim_param param = sim_param());

    ///Returns const ref to the data structure containing the
    /// demographics of the population at different points in time
    /// of the simulation
    const demographic_sim& get_demo_sim() const noexcept {return m_demo_sim;}

    ///Gets the number of cycle counters
    int get_cycle() const noexcept {return m_executed_cycles;}

    ///Gets the environment of a simulation
    const environment& get_env() const noexcept {return m_e;}

    ///Gets the reference to environment of a simulation
    environment& get_env() noexcept {return m_e;}

    ///Gets const ref to the funders_success object
    const funders_success& get_funders_success() const noexcept {return m_funder_success;}

    ///Gets ref to the funders_success object
    funders_success& get_funders_success() noexcept {return m_funder_success;}

    ///Gets the const reference to meta parameters
    const meta_param get_meta_param() const noexcept {return m_meta_param;}

    ///Gets the vector containing all the individuals of the population
    const population& get_pop() const noexcept {return m_pop;}

    ///Gets the vector containing all the individuals of the population
     population& get_pop() noexcept {return m_pop;}

     ///Returns ref to rng
     std::minstd_rand& get_rng() noexcept  {return m_rng;}

    ///Gets the number of ticks in the simulation
    int get_timestep() const noexcept {return m_sim_timesteps;}

    ///Resets the timesteps to 0
    void reset_timesteps() noexcept {m_sim_timesteps = 0;}

    ///Sets the m_demo_sim to a new demographic_sim object
    void set_demo_sim( const demographic_sim& d_s) noexcept {m_demo_sim = d_s;}

    ///Sets the metaparameters of the simualtion given a metaparameter object
    void set_meta_param(const meta_param& m) noexcept {m_meta_param = m;}

    ///Ticks the counter of cycles up by oneù
    void tick_cycles() noexcept {++m_executed_cycles;}

    ///Ticks the number of timesteps counter up by one
    void tick_timesteps() noexcept {++m_sim_timesteps;}

private:

    ///Data structure that contains data regarding the demographics
    /// of a population at different points in time of the simulation
    demographic_sim m_demo_sim;

    ///Counts the number of cycles the simulation has gone through
    int m_executed_cycles = 0;

    ///Data structures storing the funders of various cycles:
    /// their ancestor_ID
    /// their GRNs
    /// their success(how much of their lineage constitutes the final population)
    funders_success m_funder_success;

    ///The vector of individuals representing a population
    population m_pop;

    ///The environment class containing the grid where substances(food, metabolite) are
    environment m_e;

    ///The metaparameters of a simulation: how many cycles, how long each cycle etc.
    meta_param m_meta_param;

    ///The timer that keeps track of how many timesteps we are in the simulation
    int m_sim_timesteps = 0;

    ///Random number generator
    std::minstd_rand m_rng;

};

///Adds a new vector of funders to the funders_success vector
void add_new_funders(simulation& s) noexcept;

///Calculates the success of each funder at the end of a cycle
funders calc_funders_success(const simulation& s);

///Changes the conditions for the individuals in the simulation
/// Both environemntal and individual
void change_conditions(simulation& s) noexcept;

///Changes the environemtnal parameters of a simulation based on its metaparameters
/// drawing from a random distribution
void change_env(simulation& s) noexcept;

///Changes the population(only the individual) parameters of a simulation based on its metaparameters
void change_pop( simulation& s);

///Changes the parameters of the simulation given two objects of the env_ and ind_
/// parameter classes
void change_params(simulation& s, const env_param& e, const ind_param& i);

///Creates a name for the file where the run for random conditions is saved
std::string create_random_condition_name(const simulation& s, double amplitude);

///Creates a vector of random conditions given parameter classes and how much wider
/// can be the oscillation of the new random conditions
std::vector<std::pair<env_param, ind_param>> create_rand_conditions_unif(const env_param& e,
                                                                         const ind_param& i,
                                                                         int n_rand_conditions,
                                                                         double amplitude,
                                                                         int seed);
///Selects a new population and places it in a new environment
void dispersal(simulation& s);

///Runs a simulation for a given amount of cycles
/// stated by the simulation meta parameters
void exec(simulation& s) noexcept;

///Runs a cycle for a given amount of timesteps
/// stated by the simulation meta parameters
void exec_cycle(simulation& s) noexcept;

///Checks if a file with a given name exists(not tested=
///Taken from:
///https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
inline bool exists (const std::string& name);

/// All the individuals feed from environment
void feeding(simulation& s);

///Individuals feed a proportion of total food in their cell not a fixed rate
void jordi_feeding(simulation& s);

////Loads sims from a sim_param and last_pop funders objects
/// saved with a given seed and change freq parameters
simulation load_sim_for_rand_cond(int seed, int change_freq);

/// Makes a copy of a simulation but erases
/// the funder_success and demographic sim vectors
simulation no_demographic_copy(const simulation& s);

///Stores ancestor_ID and GRN of funders of a cycle in funders_success
funders prepare_funders(const simulation& s);

///Resets a simulation to its initial conditions
void reset_sim(simulation& s) noexcept;

///Individuals read input from environment and determine their own phenotype
void response(simulation& s);

///Runs a poplation from a simulation against a series of random conditions
demographic_sim run_random_conditions(const simulation &s,
                           int n_number_rand_cond, double amplitude);

///Save all necessary data at the end of the simulation
void save_data(const simulation& s);

///All individuals secrete metabolite into environment
void secretion_metabolite(simulation& s);

///Changes the demographic cycle object with a mroe recent one
void store_demographics(simulation &s) noexcept;

 ///Runs all the necessary actions for a timestep to happen
int tick(simulation& s);

///Stores the demographic of the population in the simulation
/// at that point in time
demographic_sim update_demographics(const simulation& s) noexcept;

void test_simulation();

#endif // SIMULATION_H
