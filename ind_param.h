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
              double uptake_range = 0.05,
              double metabolic_rate = 0.01,
              double reproduction_prob = 0.5,
              double reproduction_range = 0.25,
              double spor_metabolic_rate = 0.5,
              double spor_metab_range = 0.25,
              int transformation_time = 5,
              int transformation_range = 2,
              double metab_secretion_rate = 1,
              double range_on_change_ratio = 10);

    ///gets the metabolic rate of an individual
    double get_metabolic_rate() const noexcept {return m_metabolic_rate;}

    ///Gets the rate of secretion of metabolite
    double get_metab_secr_rate() const noexcept{return m_metab_secr_rate;}

    ///gets the metabolic rate of a sporulating individual
    double get_spor_metabolic_rate() const noexcept {return m_spor_metabolic_rate;}

    ///Gets const ref to the range of all the possible sporulation metabolism values
    const std::uniform_real_distribution<double>& get_spor_metabolic_range() const noexcept
    {return m_spor_metabolic_range;}

    ///Gets ref to the range of all the possible sporulation metabolism values
    std::uniform_real_distribution<double>& get_spor_metabolic_range() noexcept
    {return m_spor_metabolic_range;}

    ///gets radius of individual
    double get_radius() const noexcept {return m_radius;}

    ///Gets the ratio of the span of the range over the minimum change allowed
    /// of all the possible changeing values, except transformation time
    double get_range_on_change_ratio() const noexcept {return m_range_on_change_ratio;}

    ///Gets the reproduction probability of an ind once sufficient energy is gathered
    double get_repr_prob() const noexcept {return m_repr_prob;}

    ///Gets const ref to the range of all the possible reproduction prob values
    const std::uniform_real_distribution<double>& get_repr_range() const noexcept {return m_repr_prob_range;}

    ///Gets const ref to the range of all the possible reproduction prob values
    std::uniform_real_distribution<double>& get_repr_range() noexcept {return m_repr_prob_range;}

    ///Gets the time it takes to transform into a spore
    int get_transformation_time() const noexcept {return m_transformation_time;}

    ///Gets const ref to the range of all the possible transformation/sporulation times
    const std::uniform_int_distribution<int>& get_transformation_range() const noexcept
    {return m_transformation_range;}

    ///Gets const ref to the range of all the possible transformation/sporulation times
    std::uniform_int_distribution<int>& get_transformation_range() noexcept
    {return m_transformation_range;}

    ///gets the food uptake_rate of an individual
    double get_uptake_rate() const noexcept {return m_uptake_rate;}

    ///Gets const ref to the range of all the possible uptake values
    const std::uniform_real_distribution<double>& get_uptake_range() const noexcept {return  m_uptake_range;}

    ///Gets const ref to the range of all the possible uptake values
    std::uniform_real_distribution<double>& get_uptake_range() noexcept {return  m_uptake_range;}

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

    ///The ratio between total range span of a variable and the minimum change that
    /// that variable can undergo(e.g. = 10 means that the variable will change at least
    /// of a tenth of the delta between max and min of the range)
    /// P.S. this ratio is not applied to transformation time for now
    double m_range_on_change_ratio;

    ///Probability of reproduction once sufficient energy has been gathered
    double m_repr_prob;

    ///The range of possible reproduction probabilities
    std::uniform_real_distribution<double> m_repr_prob_range;

    ///Metabolic rate for sporulating individuals
    double m_spor_metabolic_rate;

    ///The range of possible sporulation metabolic rates
    std::uniform_real_distribution<double> m_spor_metabolic_range;

    ///number of time steps the individual needs
    ///to go through to sporulate, if it is alive at
    ///m_transformation time + 1 it will become a spore
    int m_transformation_time;

    ///The possible range of sporulation time values
    std::uniform_int_distribution<int> m_transformation_range;

    ///The level of energy required to divide
    double m_treshold_energy;

    ///The rate of food uptake, the conversion of food into energy is = 1
    double m_uptake_rate;

    ///The possible range of uptake values
    std::uniform_real_distribution<double> m_uptake_range;
};

//Compares two instantiations of ind_param
bool operator==(const ind_param& lhs, const ind_param& rhs) noexcept;

//Prints the ind_param to terminal
std::ostream& operator<<(std::ostream& os, const ind_param& p);

//Initializes a instance p from a file stream
std::ifstream& operator>>(std::ifstream& is, ind_param& p);

///Changes the individual parameters based on their ranges of change
ind_param change_ind_param(ind_param i, std::minstd_rand& rng);

///Changes the parameter of reproduction probability based on range and change
double change_ind_repr_prob(ind_param& i, std::minstd_rand& rng);

///Changes the parameter of sporulation metabolism based on range and change
double change_ind_spor_metab(ind_param& i, std::minstd_rand& rng);

///Changes the parameter of uptake rate based on range and change
double change_ind_uptake(ind_param& i, std::minstd_rand& rng);

//Initializes a instance p from a filename
ind_param load_ind_parameters( const std::string& filename);

//Saves an instance of ind_param to a file name
void save_ind_parameters( const ind_param& p, const std::string& filename);

void test_ind_param() noexcept;

#endif // IND_PARAM_H
