#ifndef META_PARAM_H
#define META_PARAM_H
#include <iostream>
#include <fstream>
#include <sstream>

class meta_param
{
public:
    meta_param(int n_cycles = 1,
               int cycle_duration = 200,
               int seed = 1,
               int change_frequency = 0,
               int pop_max = 1000000,
               int collision_check_interval = 1);

    ///Returns after how many generations
    ///  the parameters of  environment
    ///  and population will change
    int get_change_freq() const noexcept {return m_change_frequency;}

    ///Returns the number of timesteps after which collision are checked
    int get_collision_check_interval() const noexcept {return m_collision_check_interval; }

    ///Returns number of cycles for which the simulation will last
    int get_n_cycles() const noexcept {return m_n_cycles;}

    ///Returns number of ticks per cycle
    int get_cycle_duration() const noexcept {return m_cycle_duration;}

    ///Returns the max number of individual in a population
    int get_pop_max() const noexcept {return m_pop_max;}

    ///Returns the seed
    int get_seed() const noexcept {return  m_seed;}

private:

    ///Indicates how often the environment and population parameters will change
    int m_change_frequency;

    ///The number of timesteps after which collisions are checked
    int m_collision_check_interval;

    ///The number of timesteps executed in one colony cycle
    int m_cycle_duration;

    ///The number of cycles of funding and growing
    /// of individual colonies that will be simulated
    int m_n_cycles;

    ///The max number of individuals allowed in a population
    int m_pop_max;

    ///The seed with which the rng of the simulation will be seeded
    int m_seed;

};


//Shows metaparameters in terminal
std::ostream& operator<<(std::ostream& os, const meta_param& p);

//Instantiates a meta_param from terminal
std::ifstream& operator>> (std::ifstream& is, meta_param& p);

//Compares two instantiations of metaparameters to see if they are equal
bool operator==(const meta_param& lhs, const meta_param& rhs) noexcept;

//Loads metaparameters from a given file name
meta_param load_meta_parameters( const std::string& filename);

//Saves the parameters to a given file name
void save_meta_parameters( const meta_param& p, const std::string& filename);

void test_meta_param() noexcept;
#endif // META_PARAM_H
