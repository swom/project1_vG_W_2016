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
               double range_env_change = 0,
               double magnitude_env_change = 0,
               int seed = 1);

    ///Returns number of cycles for which the simulation will last
    int get_n_cycles() const noexcept {return m_n_cycles;}

    ///Returns number of ticks per cycle
    int get_cycle_duration() const noexcept {return m_cycle_duration;}

    ///Returns the range of environmental change
    double get_range_env_change() const noexcept {return m_range_env_change;}

    ///Returns the magnitude of environmental change
    double get_magnitude_env_change() const noexcept {return m_magnitude_env_change;}

    ///Returns the seed
    int get_seed() const noexcept {return  m_seed;}

private:

    ///The number of timesteps executed in one colony cycle
    int m_cycle_duration;

    ///The number of cycles of funding and growing
    /// of individual colonies that will be simulated
    int m_n_cycles;

    ///The minimum step of change from previous value for
    /// changing environmental paramters
    double m_magnitude_env_change;

    ///The maximum distance from the original
    /// initialization values to which env parameters
    /// can vary
    double m_range_env_change;

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
