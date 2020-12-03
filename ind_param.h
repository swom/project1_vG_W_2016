#ifndef IND_PARAM_H
#define IND_PARAM_H
#include <iostream>
#include <fstream>
#include <random>
#include <sstream>

class ind_param
{
public:
    ind_param(double radius  = 0.8,
              double treshold_energy = 10,
              double uptake_rate = 0.1,
              double uptake_rate_mean = 0.1,
              double uptake_rate_var = 0.01,
              double metabolic_rate = 0.01,
              double reproduction_prob = 0.5,
              double reproduction_prob_mean = 0.5,
              double reproduction_prob_var = 0.01,
              double spor_metabolic_rate = 0.5,
              double spor_metabolic_rate_mean = 0.5,
              double spor_metabolic_rate_var = 0.01,
              int transformation_time = 5,
              int transformation_time_mean = 5,
              int transformation_range = 1,
              double metab_secretion_rate = 1);

    ///gets the metabolic rate of an individual
    double get_metabolic_rate() const noexcept {return m_metabolic_rate;}

    ///Gets the rate of secretion of metabolite
    double get_metab_secr_rate() const noexcept{return m_metab_secr_rate;}

    ///gets the metabolic rate of a sporulating individual
    double get_spor_metabolic_rate() const noexcept {return m_spor_metabolic_rate;}

    ///Gets const ref to the range of all the possible sporulation metabolism values
    double get_spor_metabolic_rate_mean() const noexcept {return m_mean_spor_metabolic_rate;}

    ///Gets ref to the range of all the possible sporulation metabolism values
    double get_spor_metabolic_rate_var() const noexcept {return m_var_spor_metabolic_rate;}

    ///gets radius of individual
    double get_base_radius() const noexcept {return m_radius;}

    ///Gets the reproduction probability of an ind once sufficient energy is gathered
    double get_repr_prob() const noexcept {return m_repr_prob;}

    ///Gets const ref to the range of all the possible reproduction prob values
    double get_repr_prob_mean() const noexcept {return m_mean_repr_prob;}

    ///Gets const ref to the range of all the possible reproduction prob values
    double get_repr_prob_var() const noexcept {return m_var_repr_prob;}

    ///Gets the time it takes to transform into a spore
    int get_transformation_time() const noexcept {return m_transformation_time;}

    ///Gets the mean value of all possible transformation time values
    int get_transformation_time_mean() const noexcept {return m_mean_transformation_time;}

    ///Gets range of all the possible transformation/sporulation times
    int get_transformation_range() const noexcept
    {return m_transformation_range;}

    ///gets the food uptake_rate of an individual
    double get_uptake_rate() const noexcept {return m_uptake_rate;}

    ///Gets const ref to the range of all the possible uptake values
    double get_uptake_mean() const noexcept {return  m_mean_uptake_rate;}

    ///Gets const ref to the range of all the possible uptake values
    double get_uptake_var() const noexcept {return  m_var_uptake_rate;}

    ///gets the treshold energy at which an individual reproduces
    double get_treshold_energy() const noexcept {return m_treshold_energy;}

    ///Sets the radius, used only in m_relax object for managing hilbert_tree collisions
    void set_radius(float r) noexcept {m_radius = static_cast<double>(r);}

    ///Sets the reproduciton probability
    void set_repr_prob(double p) noexcept {m_repr_prob = p;}

    ///Sets the metabolic rate of sporulating individuals
    void set_spor_metabolic_rate(double m) noexcept {m_spor_metabolic_rate = m;}

    ///Sets the number of timesteps required to transform into a spore
    void set_transformation_time(int t) noexcept {m_transformation_time = t;}

    ///Sets the uptake rate of nutrient of an individual
    void set_uptake_rate(double u) noexcept {m_uptake_rate = u;}

private:

    ///The rate at which internal energy is depleted
    double m_metabolic_rate;

    ///The rate at which metabolite is secreted into the environment
    double m_metab_secr_rate;

    ///Radius of an individual, individuals are considered circular
    double m_radius;

    ///The current probability of reproduction once sufficient energy has been gathered
    double m_repr_prob;

    ///The mean of the normal distribution
    /// encompassing all possible values
    /// of reproduction probabilities
    double m_mean_repr_prob;

    /// The variance of the normal distribution
    /// encompassing all possible values
    /// of reproduction porbability
    double m_var_repr_prob;

    ///Metabolic rate for sporulating individuals
    double m_spor_metabolic_rate;

    ///The mean of the normal distribution
    /// encompassing all possible values
    /// of sporulation metabolic rate
    double m_mean_spor_metabolic_rate;

    /// The variance of the normal distribution
    /// encompassing all possible values
    /// of sporulation metabolic rate
    double m_var_spor_metabolic_rate;

    ///number of time steps the individual needs
    ///to go through to sporulate, if it is alive at
    ///m_transformation time + 1 it will become a spore
    int m_transformation_time;

    ///Mean of the possible sporulation times
    int m_mean_transformation_time;

    ///The possible range of sporulation time values
    /// mean - range, mean + range
    int m_transformation_range;

    ///The level of energy required to divide
    double m_treshold_energy;

    ///The current rate of food uptake, the conversion of food into energy is = 1
    double m_uptake_rate;

    ///The mean of the normal distribution encompassing all possible values
    /// of the rate of food uptake, the conversion of food into energy is = 1
    double m_mean_uptake_rate;

    /// The variance of the normal distribution encompassing all possible values
    /// of rate of food uptake, the conversion of food into energy is = 1
    double m_var_uptake_rate;

};

//Compares two instantiations of ind_param
bool operator==(const ind_param& lhs, const ind_param& rhs) noexcept;

//Prints the ind_param to terminal
std::ostream& operator<<(std::ostream& os, const ind_param& p);

//Initializes a instance p from a file stream
std::ifstream& operator>>(std::ifstream& is, ind_param& p);

///Returns a new ind_param that is a changed version of the
/// given ind_param
/// whith new values drawn from a normal distribution
/// with mean and variance as indicated by the m_mean* and m_var* members
ind_param change_ind_param_norm( ind_param i,  std::minstd_rand& rng);

///Returns a new ind_param that is a changed version of the
/// given ind_param
/// whith new values drawn from a uniform distribution
/// with mean as indicated by the m_mean* members
/// and range = m_mean* -/+ 3 * m_var* members
/// EXCEPT for m_transformation_time that uses m_mean -/+ m_range
/// without multiplication
ind_param change_ind_param_unif( ind_param i,  std::minstd_rand& rng);

//Initializes a instance p from a filename
ind_param load_ind_parameters( const std::string& filename);

///Returns a env_param whose variances are the variances
/// of the given env_param object multiplied
/// by a factor = amplitude
ind_param change_range_ind_param(const ind_param& i, double amplitude);

//Saves an instance of ind_param to a file name
void save_ind_parameters( const ind_param& p, const std::string& filename);

void test_ind_param() noexcept;

#endif // IND_PARAM_H
