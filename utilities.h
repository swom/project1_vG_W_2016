#ifndef UTILITIES_H
#define UTILITIES_H
#include <iostream>
#include "random"
#include <vector>


///Checks from command line argument if simulation
/// takes parameters from command line, and which
void check_for_cmd_param(const std::vector<std::string>& args,
                         int& seed,
                         int& change_freq,
                         int& n_conditions,
                         double& amplitude);

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

///Takes number of condition against which a population is tested
/// from command line arguments
void take_n_conditions_arg(const std::vector<std::string>& args, int& n_conditions);

///Takes change of frequency number from command line arguments
void take_change_freq_arg(const std::vector<std::string>& args, int &change_freq);

///Takes seed number from command line arguments
void take_seed_arg(const std::vector<std::string>& args, int& seed);


void test_utilities() noexcept;
#endif // UTILITIES_H
