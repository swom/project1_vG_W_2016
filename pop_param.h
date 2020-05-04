#ifndef POP_PARAM_H
#define POP_PARAM_H
#include "ind_param.h"

class pop_param
{
public:
    pop_param(unsigned int start_pop_size = 1,
              unsigned int exp_new_pop_size = 1,
              double min_dist = 0.1,
              double mutation_prob = 0.01,
              double mutation_step = 0.1,
              double base_disp_prob = 0.01,
              double spore_advantage = 10.0,
              double reproduction_prob = 0.1,
              ind_param ind_parameters = ind_param());


    ///Returns the value of the variable m_base_fitness that indicates
    /// the basal fitness/dispersal probability of an individual
    double get_base_disp_prob() const noexcept {return m_base_disp_prob;}

    ///Returns the number of individuals that should be present
    ///  in the new funding population
    unsigned int get_exp_new_pop_size() const noexcept {return m_exp_new_pop_size;}

    ///Get ind parameters
    const ind_param& get_ind_param() const noexcept {return m_ind_param;}

    ///Get minimum distance between individuals at the start of the simulation
    double get_min_dist() const noexcept {return m_min_init_dist_btw_inds;}

    ///Gets the reference to the mutation probability distribution
    double get_mu_p() const noexcept {return m_mutation_prob;}

    ///Gets the reference to the mutation step distribution
    double get_mu_st() const noexcept {return m_mutation_step;}

    ///Gets the starting size of the population
    unsigned int get_pop_start_size() const noexcept {return m_start_pop_size;}

    ///Gets the probability of reproducing
    double get_repr_p() const noexcept {return m_repr_prob;}

    ///Returns the variable m_spore_advantage that indicates
    ///how many times a spore is more likely to get dispersed
    ///than the other phenotypes
    double get_spo_adv() const noexcept {return m_spore_advantage;}

private:

    ///The parameters of the individuals
    ind_param m_ind_param;

    ///The base dispersal probability of an individual(so when it is active or sporulating)
    double m_base_disp_prob;

    ///The expected size of a newly funded population
    unsigned int m_exp_new_pop_size;

    ///The minimum ditance between individuals when they are placed down at the start
    ///of a population cycle
    double m_min_init_dist_btw_inds;

    ///The probability of a mutation to happen
    double m_mutation_prob;

    ///The variance of the gaussian from which the size of the mutation is drawn
    double m_mutation_step;

    ///The probability that an individual will reproduce
    double m_repr_prob;

    ///The factor for which the dispersal probability is multiplied if the individual is a spore
    double m_spore_advantage;

    ///The initial starting size of the population
    unsigned int m_start_pop_size;
};

void test_pop_param() noexcept;

#endif // POP_PARAM_H
