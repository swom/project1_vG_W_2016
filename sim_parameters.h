#ifndef SIM_PARAMETERS_H
#define SIM_PARAMETERS_H


class sim_parameters
{
public:
    sim_parameters(unsigned int start_pop_size = 1, int exp_new_pop_size = 1, double min_dist = 0.1, int grid_side = 1,
                   double diff_coeff = 0.1, double init_food = 1.0, double mutation_prob = 0.01, double mutation_step = 0.1,
                   double base_disp_prob = 0.01, double spore_advantage = 10.0, double reproduction_prob = 0.1,
                   double metab_degradation_rate = 0.01);

    ///Returns the value of the variable m_base_fitness that indicates
    /// the basal fitness/dispersal probability of an individual
    double get_base_disp_prob() const noexcept {return m_base_disp_prob;}

    ///Returns the diffusion coefficient
    double get_diff_coeff() const noexcept {return  m_diff_coeff;}

    ///Returns the number of individuals that should be present in the new funding population
    int get_exp_new_pop_size() const noexcept {return m_exp_new_pop_size;}

    ///Returns the side of the grid, used to determine the total size of the grid(square)
    int get_grid_side() const noexcept {return m_grid_side;}

    ///Gets the initial food that will be provided in each grid_cell of the environment
    double get_init_food() const noexcept {return m_init_food;}

    ///Gets the degradation rate of the metabolite
    double get_metab_degr_rate() const noexcept {return m_metab_degradation_rate;}

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

    ///The base dispersal probability of an individual(so when it is active or sporulating)
    double m_base_disp_prob;

    /// The diffusion coefficient of substances in the grid
    double m_diff_coeff;

    ///The expected size of a newly funded population
    int m_exp_new_pop_size;

    ///The side of the grid
    int m_grid_side;

    ///The initial amount of food in each grid_cell at the start of a pop cycle
    double m_init_food;

    ///The minimum ditance between individuals when they are placed down at the start
    ///of a population cycle
    double m_min_init_dist_btw_inds;

    ///The rate at which metabolite degrades
    double m_metab_degradation_rate;

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

void test_sim_parameters() noexcept;

#endif // SIM_PARAMETERS_H
