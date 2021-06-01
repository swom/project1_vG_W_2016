#ifndef SIMULATION_H
#define SIMULATION_H

#include "env_changer.h"
#include "environment.h"
#include "individual.h"
#include "population.h"
#include "sim_parameters.h"
#include "utilities.h"
#include <vector>
#include <random>

class simulation
{
public:

    simulation(sim_param param = sim_param());

    ///Makes death active
    void activate_death() noexcept {m_death_is_active = true;}

    ///Checks that the selection for spores flag is active and enough time has passed (50 cycles)
    bool can_activate_selection_for_spores() noexcept {return m_select_for_spores && (m_executed_cycles > 50);}

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
    /// !!!ATTENTION!!! This is not super safe
    /// But is necessary to add funders to the vector
    /// of funder success and then update their success
    /// during that generation
    funders_success& get_funders_success() noexcept {return m_funder_success;}

    ///Gets the const reference to meta parameters
    const meta_param& get_meta_param() const noexcept {return m_meta_param;}

    ///Gets the reference to meta parameters
    meta_param& get_meta_param() noexcept {return m_meta_param;}

    ///Returns the sim_param
    sim_param get_sim_param() const noexcept {return  sim_param{m_e.get_param(),
                    m_pop.get_ind(0).get_param(),
                    m_meta_param,
                    m_pop.get_param()};}

    ///Gets the vector containing all the individuals of the population
    const population& get_pop() const noexcept {return m_pop;}

    ///Gets the vector containing all the individuals of the population
    population& get_pop() noexcept {return m_pop;}

    ///Returns ref to rng
    std::minstd_rand& get_rng() noexcept  {return m_rng;}

    ///Gets the number of ticks in the simulation
    int get_timestep() const noexcept {return m_sim_timesteps;}

    ///returns the value of m_death_is_active
    bool is_death_active() const noexcept{return m_death_is_active;}

    ///Resets the timesteps to 0
    void reset_timesteps() noexcept {m_sim_timesteps = 0;}

    ///Allows the first 50 cycles to select all types of inds, and then
    /// activates the select_only_spores_flag
    void select_for_spores() noexcept {m_select_for_spores = true;}

    ///Turns on the flag that signals that only spores are to be selected
    void select_only_spores() noexcept {m_pop.select_only_spores();}

    ///Sets the m_demo_sim to a new demographic_sim object
    void set_demo_sim( const demographic_sim& d_s) noexcept {m_demo_sim = d_s;}

    ///Sets the m_funders_success member to a new funders_success object
    void set_funders_success(const funders_success& f_s) noexcept {m_funder_success = f_s;}

    ///Sets the metaparameters of the simualtion given a metaparameter object
    void set_meta_param(const meta_param& m) noexcept {m_meta_param = m;}

    ///Ticks the counter of cycles up by oneù
    void tick_cycles() noexcept {++m_executed_cycles;}

    ///Ticks the number of timesteps counter up by one
    void tick_timesteps() noexcept {++m_sim_timesteps;}

private:

    ///Is death happening?
    bool m_death_is_active =  false;

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

    ///The timer that keeps track of how many timesteps
    ///have passed in a cycle
    int m_sim_timesteps = 0;

    ///Random number generator
    std::minstd_rand m_rng;

    ///Flag that signals that after 50 cycles of selection fo all inds to fund next generation
    /// switches to a regimen where only spores are selected
    bool m_select_for_spores = false;

};

///Adds a new vector of funders to the funders_success vector
void add_new_funders(simulation& s) noexcept;

///Sets the success of the funders fo the cycle
/// Based on the number of individuals
/// present at the end of the simulation
/// descended from that funder
void add_success_funders(simulation& s) noexcept;

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
void set_new_params(simulation& s, const env_param& e, const ind_param& i);

///Continues evolution of an already run evolution simulation with a given seed and change of freq
int continue_evo(int seed, int change_freq);

///Checks that the individuals in a population
/// have exactly the same grns of the funders of a particular generation
bool check_funder_equal_pop(const population &p, const funders_success &f, const int &generation);

///Creates the name for the before to last population file
/// where the before last funder object will be saved
/// This is then used to find the most succesful individuals in a sim
std::string create_before_last_pop_name(const simulation& s,const std::string& prefix = "");

///Creates a name for the file where the
/// run of  only the best individual
/// against random conditions is saved
std::string create_best_random_condition_name(const simulation& s, double amplitude);
std::string create_best_random_condition_name(double amplitude, int  change_freq, int seed);

///Creates a name for the file where the TEST run for random conditions is saved
std::string create_test_random_condition_name(const simulation& s, double amplitude);
std::string create_test_random_condition_name(double amplitude, int change_freq, int seed);

///Creates the name of a funders_success .csv file given
/// a simulation, or the simulation's seed and change_frequency
std::string create_funders_success_name(const simulation& s, std::string prefix = "", std::string suffix = "");
std::string create_funders_success_name(int seed, int change_freq, std::string prefix = "");

///Creates the name for the last population of an
/// evolutionary run so that it can be uploaded
/// for testing against random conditions
/// given a seed number and a change frequency
std::string create_last_pop_name(int seed, int change_freq, const std::string &prefix  = "");

///Creates the name for the last population of an
/// evolutionary run so that it can be uploaded
/// for testing against random conditions given a simulation
std::string create_last_pop_name(const simulation& s, const std::string& prefix = "");

///Creates the name for the file
/// where to save a vector of pairs
/// of randomly generated env_param and ind_param,
/// given the number of conditions,
/// the amplitude of the randomness
/// and the seed of the rng.
std::string create_name_vec_rand_cond(int n_of_conditions, double amplitude, int seed);

///Create prefix for simulation derived from using create_unif_extreme random conditions
std::string create_rand_extreme_prefix(int seq_index, int cond_per_seq, double amplitude);

///Creates a name for the file where
/// the deomgraphic of the evolving population
/// in random environment is saved
std::string create_sim_demo_name(const simulation& s, std::string prefix = "", std::string suffix = "");
std::string create_sim_demo_name(int seed, int change_freq, std::string prefix = "", std::string suffix = "");

///Creates a name for the file where the parameters
/// of the evolving population in random environment is saved
/// given a simulation
std::string create_sim_par_name(const simulation& s, std::string prefix = "", std::string suffix = "");

///Creates a name for the file where the parameters
/// of the evolving population in random environment is saved
/// given a seed and a frequency of change
std::string create_sim_par_name(int seed, int change_freq, std::string prefix = "");

///Creates a vector of random conditions given
/// parameter classes and how much wider
/// can be the oscillation of the new random conditions
std::vector<std::pair<env_param, ind_param>> create_rand_conditions_unif(const env_param& e,
                                                                         const ind_param& i,
                                                                         int n_rand_conditions,
                                                                         double amplitude,
                                                                         int seed);


///Creates a vector of random conditions given
/// parameter classes and how much wider
/// can be the oscillation of the new random conditions
std::vector<std::pair<env_param, ind_param>> create_rand_conditions_unif_extreme(const env_param& e,
                                                                                 const ind_param& i,
                                                                                 int n_rand_conditions,
                                                                                 double amplitude,
                                                                                 int seed = 123456);


///Creates a Matrix of random conditions given the number of columns and rows
/// (sequences and condition per sequences) and the amplitude,
/// each column of random condition is seeded with its index
std::vector<std::vector<std::pair<env_param, ind_param>>> create_rand_conditions_matrix(const env_param& ep,
                                                                                        const ind_param& ip,
                                                                                        int number_of_sequences,
                                                                                        int conditions_per_sequence,
                                                                                        double amplitude);


///Works like random matrix
/// but only samples from extremes of uniform distribution for paramter values
std::vector<std::vector<std::pair<env_param, ind_param>>> create_rand_conditions_matrix_extreme(const env_param& ep,
                                                                                                const ind_param& ip,
                                                                                                int number_of_sequences,
                                                                                                int conditions_per_sequence,
                                                                                                double amplitude);

///Creates a vector of random env param and ind param generated by changing
/// params of given env and ind param object by a given amplitude
std::vector<std::pair<env_param, ind_param>> create_vector_random_conditions(const env_param& e,
                                                                             const ind_param& i,
                                                                             double amplitude,
                                                                             int n_conditions = 50,
                                                                             int seed = 0);
///Selects a new population and places it in a new environment
void dispersal(simulation& s);

///Checks if death is active in simulation
bool death_is_active(const simulation& s);

///Returns a demographic cycle object storing data about a population
demographic_cycle demographics(const simulation &s, const env_param&e) noexcept;

///Runs a simulation for a given amount of cycles
/// stated by the simulation meta parameters
void exec(simulation& s) noexcept;

///Runs a simulation and changes the conditions
/// given in a random_conditions vector
/// The time of change of conditions is adjusted
/// so that all conditions in the vector
/// are experienced in the simulation
/// The number of condition has to be a dividend of the number of cycles
void exec_change(simulation& s, const std::vector<std::pair<env_param,ind_param>>& rand_conditions);

///Runs a cycle for a given amount of timesteps
/// stated by the simulation meta parameters
void exec_cycle(simulation& s) noexcept;

///Runs a cycle for a given amount of timesteps
/// OR until population reaches a limit
void exec_cycle(simulation& s) noexcept;

/// All the individuals feed from environment
void feeding(simulation& s);

///Individuals feed a proportion of total food in their cell not a fixed rate
void jordi_feeding(simulation& s);

///Find mean min and max of rand_cond vector
double mean_diff_coeff_rand_cond(const std::vector<std::pair<env_param,ind_param>>& rand_cond);
double find_min_diff_coeff_rand_cond(const std::vector<std::pair<env_param,ind_param>>& rand_cond);
double find_max_diff_coeff_rand_cond(const std::vector<std::pair<env_param,ind_param>>& rand_cond);

///Loads a simulation whose population is composed
/// exclusively of the best individual
/// in the before last generation of a
/// previously saved simulation's funders_success file
simulation load_best_ind_for_rand_cond(int seed, int change_freq);

///Loads random conditions given their filename
/// (usually automated with create_rand_cond_name).
std::vector<std::pair<env_param, ind_param>> load_random_conditions(const std::string& filename);

////Loads simulation (a sim_param and sim_demo and fuinders_succes object)
/// saved with a given seed and change freq parameters
/// Population is set to the last funders
simulation load_sim(int seed, int change_freq);

///Loads a simulation from another simulation, with
/// no funders_success or sim_demographic,
/// with same sim_param and with
/// pop initialized to BEFORE last fuinders in given sim
simulation load_sim_before_last_pop(int seed, int change_freq);

///Loads a simulation from another simulation, with
/// no funders_success or sim_demographic,
/// with same sim_param and with
/// pop initialized to last fuinders in given sim
simulation load_sim_last_pop(int seed, int change_freq, std::string prefix = "");

///Loads a simulation from the files produced from
/// an evolutionary run of the simulation
/// The population is not set
simulation load_sim_no_pop(int seed, int change_freq, std::string prefix = "");

/// Makes a copy of a simulation but erases
/// the funder_success and demographic sim vectors
simulation no_dem_and_fund_copy(const simulation& s);


///Stores ancestor_ID and GRN of funders of a cycle in funders_success
funders prepare_funders(const simulation& s);

///Sets the environment of a simulation to the environment that happened in a certain cycle
void reproduce_cycle_env(simulation&s, int cycle);

///Sets the individuals of the populaiton to be the
/// individuals of a specified generation of funders
/// And the enviromental conditions to be the conditions in that generation
void reproduce_cycle(simulation&s, int cycle);

///Resets the population of a simulation to the one at the start
/// of a specified generation
void recreate_generation(simulation& s, int cycle);

///Sets the environmental and individual conditions of a population
/// to the condition o the n_rand_cond element of the rand_cond vector
void reproduce_rand_cond(simulation&s, const std::vector<std::pair<env_param, ind_param>>& rand_cond, int n_rand_cond);

///Runs an evolution simulation
/// starting from the last population of a previous simulation
/// in a given random condition
int run_sim_evo_rand(double amplitude,
                     int change_frequency,
                     int num_of_sequences,
                     int conditions_per_seq,
                     int pop_max,
                     int seed,
                     int seq_index,
                     bool overwrite, bool death, bool eden);

///Runs an evolution simulation in a given random condition
funders_success run_evo_random_conditions(const simulation& rand_s,
                                          int number_of_sequences,
                                          int cond_per_seq,
                                          int seq_index,
                                          int pop_max,
                                          double amplitude,
                                          std::string prefix);

///Creates the reaction norm of the best individuall of a certain simulation (seed, change_freq)
/// over a range of food, energy and metabolite values (from 0 to max_*)
/// checking over all conditions witha given interval step
int run_reac_norm_best(int change_freq,
                       double max_food,
                       double max_energy,
                       double max_metabolite,
                       double step,
                       int seed,
                       bool overwrite);

///Runs the a population made of 100 clones of the best individual
/// of a given simulation(change_freq, seed) against n random condtions
int run_sim_best_rand(double amplitude,
                      int change_frequency,
                      int n_random_conditions,
                      int pop_max,
                      int seed,
                      bool overwrite);

///Runs the last population of funders of a given simulation(change_freq, seed)
/// against n random conditions
int run_sim_rand(double amplitude,
                 int change_frequency,
                 int n_random_conditions, int pop_max,
                 int seed,
                 bool overwrite);

///Createsa an evolution simulation with
/// the given parameters, seed and change_frequency
int run_sim_evo(const env_param& e,
                const ind_param& i,
                const meta_param& m,
                const pop_param& p,
                bool overwrite,
                bool death,
                bool select_for_spores);

///A standard simulation is created, with given parameters,
/// and run for boht evo and against random conditions
int run_standard(const env_param& e,
                 const ind_param &i,
                 const meta_param& m,
                 const pop_param& p,
                 double amplitude,
                 int change_frequency,
                 int n_random_conditions,
                 int pop_max,
                 int seed);

///Tests the first and the last population of rand_Evo_extreme runs
/// against
/// !!!ATTENTION!!! The number of random condition is hardcoded for now
demographic_sim run_test_extreme_rand_evo_beginning_end(int original_seed,
                                    int original_change,
                                    int cond_per_seq,
                                    int seq_index,
                                    double amplitude, bool death, bool eden);

///Runs a poplation from a simulation against a series of random conditions
demographic_sim test_against_random_conditions(const simulation &s,
                                               int n_number_rand_cond,
                                               int pop_max,
                                               double amplitude,
                                               std::string name);

///Resets a simulation to its initial conditions
void reset_sim(simulation& s) noexcept;

///Individuals read input from environment and determine their own phenotype
void response(simulation& s);



///Saves a vector of pairs of randomlly generated env_param and ind_param
void save_vector_of_rand_cond(const std::vector<std::pair<env_param, ind_param>>& rand_cond_v,
                              const std::string& filename);

///Save all necessary data at the end of the simulation
void save_data(const simulation& s, std::string prefix = "");

///Saves the funder that fund the before last cycle
funders save_before_last_pop(const simulation& s, std::string prefix);

///Saves the funder that fund the before last cycle
funders save_last_pop(const simulation& s, const std::string &prefix);

///saves the before last pop and the last pop,
/// a prefix can be added to the file names to
/// disinguish last pops coming from different kinds of simulations
std::vector<funders> save_two_last_pops(const simulation& s, const std::string &prefix = "");

///All individuals secrete metabolite into environment
void secretion_metabolite(simulation& s);

///Changes the demographic cycle object with a mroe recent one
void store_demographics(simulation &s) noexcept;

///Runs all necessary actions for a timestep but
/// check collision and resolves them only every n timesteps
int tick(simulation& s, int n_ticks = 0);

///Stores the demographic of the population in the simulation
/// at that point in time
demographic_sim update_demographics(const simulation& s) noexcept;


#endif // SIMULATION_H
