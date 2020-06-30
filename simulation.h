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

    ///Gets the number of ticks in the simulation
    int get_timestep() const noexcept {return m_sim_timesteps;}

    ///Resets the timesteps to 0
    void reset_timesteps() noexcept {m_sim_timesteps = 0;}

    ///Sets the m_demo_sim to a new demographic_sim object
    void set_demo_sim( const demographic_sim& d_s) noexcept {m_demo_sim = d_s;}

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

};

///Adds a new vector of funders to the funders_success vector
void add_new_funders(simulation& s) noexcept;

///Calculates the success of each funder at the end of a cycle
funders calc_funders_success(const simulation& s);

///Changes the environemtnal parameters of a simulation based on its metaparameters
void change_env(simulation& s) noexcept;

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

///Stores ancestor_ID and GRN of funders of a cycle in funders_success
funders prepare_funders(const simulation& s);

///Resets a simulation to its initial conditions
void reset_sim(simulation& s) noexcept;

///Individuals read input from environment and determine their own phenotype
void response(simulation& s);

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
