#ifndef SIM_PARAMETERS_H
#define SIM_PARAMETERS_H
#include"env_changer.h"
#include"pop_param.h"
#include"meta_param.h"

#include <iostream>

class sim_param
{
public:

    sim_param(env_changer e, ind_param i, meta_param m, pop_param p);

    ///Gets const reference to population parameter
    const env_changer& get_env_changer() const noexcept {return m_env_changer;}

    ///Gets const reference to metaparameters
    const meta_param& get_meta_param() const noexcept {return  m_meta_param;}

    ///Gets const reference to population parameter
    const pop_param& get_pop_param() const noexcept {return m_pop_param;}

    ///Gets const reference to ind parameter
    const ind_param& get_ind_param() const noexcept {return m_ind_param;}


private:

    ///Parameters for the environment
    env_changer m_env_changer;

    ///Parameters pertaining individuals
    ind_param m_ind_param;

    ///Parameters concerning the simulation, such as duration, number of replicates etc
    meta_param m_meta_param;

    ///Parameters concerning the population of individuals
    pop_param m_pop_param;
};

//Shows sim_param in terminal
std::ostream& operator<<(std::ostream& os, const sim_param& p);

//Compares two instantiations of sim_param to see if they are equal
bool operator==(const sim_param& lhs, const sim_param& rhs) noexcept;

//Loads sim_param from a given file name
sim_param load_sim_parameters( const std::string& filename);

//Saves the sim parameters to a given file name
void save_sim_parameters( const sim_param& p, const std::string& filename);


#endif // SIM_PARAMETERS_H
