#ifndef UTILITIES_H
#define UTILITIES_H
#include <iostream>
#include "random"
#include <vector>

///Compares two vectors of double with a certain tolerance to differences in precision
/// used mainly for operators that act on loaded objects
bool compare_with_tolerance(const std::vector<double>& lhs,const std::vector<double>& rhs);

///Checks from command line argument if simulation
/// takes parameters from command line, and which
void check_for_cmd_param(const std::vector<std::string>& args,
                         int& seed,
                         int& change_freq,
                         int& n_conditions,
                         int &replay_cycle,
                         double& amplitude,
                         bool &overwrite,
                         int &seed_rand_cond,
                         int &rand_cond_n);

///Creates a uniform distribution
std::uniform_real_distribution<double> create_unif_dist(double a, double b) noexcept;

///Creates a bernoulli distribution
std::bernoulli_distribution create_bernoulli_dist(double p) noexcept;

///Creates a normal distribution
std::normal_distribution<double> create_normal_dist(double m, double v);

///Checks if a file with a given name exists(not tested=
///Taken from:
///https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
bool exists (const std::string& name);

///Takes amplitude of oscillation of params from command line arguments
void take_amplitude_arg(const std::vector<std::string>& args, double &amplitude);

///Checks if the last argumetn of the command line is "--overwrite"
/// if it is then the simulations will be run even if they had already run and
/// will override the old saved files
void take_overwrite_arg(const std::vector<std::string>& args, bool& overwrite);

///Takes number of condition against which a population is tested
/// from command line arguments
void take_n_sequences_arg(const std::vector<std::string>& args, int& n_conditions);

///Takes change of frequency number from command line arguments
void take_change_freq_arg(const std::vector<std::string>& args, int &change_freq);

///Takes the index of the random conditions in the random condition vector
/// that needs to be replayed
void take_seq_index_arg(const std::vector<std::string>& args, int& rand_cond_n);

/// Takes the seed with which the random conditions vector was generated
/// and assigns it to seed_rand_cond
void take_cond_per_seq_arg(const std::vector<std::string>& args, int& cond_per_seq);

///Takes the number of the cycle that is gonna be visually replayed
void take_replay_cycle_arg(const std::vector<std::string>& args, int& replay_cycle);

///Takes seed number from command line arguments
void take_seed_arg(const std::vector<std::string>& args, int& seed);


#endif // UTILITIES_H
