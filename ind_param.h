#ifndef IND_PARAM_H
#define IND_PARAM_H
#include <iostream>
#include <fstream>
#include <sstream>

class ind_param
{
public:
    ind_param(double radius  = 0.8,
              double treshold_energy = 10,
              double uptake_rate = 0.1,
              double metabolic_rate = 0.01,
              double spor_metabolic_rate = 0.5,
              int transformation_time = 5,
              double wiggle_room = 0.01,
              double metab_secretion_rate = 1);

    ///gets the metabolic rate of an individual
    double get_metabolic_rate() const noexcept {return m_metabolic_rate;}

    ///Gets the rate of secretion of metabolite
    double get_metab_secr_rate() const noexcept{return m_metab_secr_rate;}

    ///gets the metabolic rate of a sporulating individual
    double get_spor_metabolic_rate() const noexcept {return m_spor_metabolic_rate;}

    ///gets radius of individual
    double get_radius() const noexcept {return m_radius;}

    ///Gets the time it takes to transform into a spore
    int get_transformation_time() const noexcept {return m_transformation_time;}

    ///gets the food uptake_rate of an individual
    double get_uptake_rate() const noexcept {return m_uptake_rate;}

    ///gets the treshold energy at which an individual reproduces
    double get_treshold_energy() const noexcept {return m_treshold_energy;}

    ///Gets the wiglle room of an individual
    double get_wiggle_room() const noexcept {return m_wiggle_room;}

    ///Sets the radius, used only in m_relax object for managing hilbert_tree collisions
    void set_radius(float r) noexcept {m_radius = static_cast<double>(r);}

private:

    ///The rate at which internal energy is depleted
    double m_metabolic_rate;

    ///The rate at which metabolite is secreted into the environment
    double m_metab_secr_rate;

    ///Radius of an individual, individuals are considered circular
    double m_radius;

    ///Metabolic rate for sporulating individuals
    double m_spor_metabolic_rate;

    ///number of time steps the individual needs
    ///to go through to sporulate, if it is alive at
    ///m_transformation time + 1 it will become a spore
    int m_transformation_time;

    ///The level of energy required to divide
    double m_treshold_energy;

    ///The rate of food uptake, the conversion of food into energy is = 1
    double m_uptake_rate;

    ///The minimum amount of overlap necessary for a detection between two individuals
    /// to be detected, implemented to speed up collision management
    double m_wiggle_room;
};

//Compares two instantiations of ind_param
bool operator==(const ind_param& lhs, const ind_param& rhs) noexcept;

//Prints the ind_param to terminal
std::ostream& operator<<(std::ostream& os, const ind_param& p);

//Initializes a instance p from a file stream
std::ifstream& operator>>(std::ifstream& is, ind_param& p);

//Initializes a instance p from a filename
ind_param load_ind_parameters( const std::string& filename);

//Saves an instance of ind_param to a file name
void save_ind_parameters( const ind_param& p, const std::string& filename);

void test_ind_param() noexcept;

#endif // IND_PARAM_H
