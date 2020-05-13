#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include <random>
#include "environment.h"
#include "individual.h"
#include "population.h"
#include "sim_parameters.h"

class simulation
{
public:

    simulation(sim_param param = sim_param());

    ///Gets the number of cycle counters
    int get_cycle() const noexcept {return m_executed_cycles;}

    ///Gets the environment of a simulation
    const environment& get_env() const noexcept {return m_e;}

    ///Gets the reference to environment of a simulation
    environment& get_env() noexcept {return m_e;}

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

    ///Ticks the counter of cycles up by oneù
    void tick_cycles() noexcept {++m_executed_cycles;}

    ///Ticks the number of timesteps counter up by one
    void tick_timesteps() noexcept {++m_sim_timesteps;}

private:

    ///Counts the number of cycles the simulation has gone through
    int m_executed_cycles = 0;

    ///The vector of individuals representing a population
    population m_pop;

    ///The environment class containing the grid where substances(food, metabolite) are
    environment m_e;

    ///The metaparameters of a simulation: how many cycles, how long each cycle etc.
    meta_param m_meta_param;

    ///The timer that keeps track of how many timesteps we are in the simulation
    int m_sim_timesteps = 0;

    ///The random number generator of simulation(used for everything)
    std::minstd_rand m_rng;

};


///Selects a new population and places it in a new environment
void dispersal(simulation& s);

///Runs a simulation for a given amount of cycles
/// stated by the simulation meta parameters
void exec(simulation& s) noexcept;

///Runs a cycle for a given amount of timesteps
/// stated by the simulation meta parameters
void exec_cycle(simulation& s) noexcept;

///Overload that uses the simulation meta parameters
void exec_cycle(simulation& s) noexcept;

/// All the individuals feed from environment
void feeding(simulation& s);

///All individuals secrete metabolite into environment
void secretion_metabolite(simulation& s);

///Resets a simulation to its initial conditions
void reset_sim(simulation& s) noexcept;

///Individuals read input from environment and determine their own phenotype
void response(simulation& s);

 ///Runs all the necessary actions for a timestep to happen
int tick(simulation& s);

void test_simulation();

#endif // SIMULATION_H
