#ifndef SIMULATION_H
#define SIMULATION_H
#include <vector>
#include <random>
#include "environment.h"
#include "individual.h"
#include "sim_parameters.h"
#include "population.h"

class simulation
{
public:

    simulation(sim_param param = sim_param());

    ///Gets the environment of a simulation
    const environment& get_env() const noexcept {return m_e;}

    ///Gets the reference to environment of a simulation
    environment& get_env() noexcept {return m_e;}

    ///Gets the vector containing all the individuals of the population
    const population& pop() const noexcept {return m_pop;}

    ///Gets the vector containing all the individuals of the population
     population& pop() noexcept {return m_pop;}

    ///Gets the number of ticks in the simulation
    int get_tick() const noexcept {return m_sim_timer;}

    ///Updates tick counter by one
    void update_sim_timer() noexcept {++m_sim_timer;}

private:

    ///The vector of individuals representing a population
    population m_pop;

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


///Selects a new population and places it in a new environment
void dispersal(simulation& s);

///Runs a simulation or a given amount of time
void exec(simulation& s, int n_tick) noexcept;

/// All the individuals feed from environment
void feeding(simulation& s);

///All individuals secrete metabolite into environment
void secretion_metabolite(simulation& s);

///Individuals read input from environment and determine their own phenotype
void response(simulation& s);

 ///Runs all the necessary actions for a timestep to happen
int tick(simulation& s);

void test_simulation();

#endif // SIMULATION_H
