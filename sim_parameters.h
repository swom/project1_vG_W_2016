#ifndef SIM_PARAMETERS_H
#define SIM_PARAMETERS_H


class sim_parameters
{
public:
    sim_parameters();

    ///Returns the number of individuals that should be present in the new funding population
    int get_exp_mew_pop_size() const noexcept {return m_exp_new_pop_size;}

    ///Returns the minimum distance between individuals
    ///when they are placed in space at the start of a
    ///population cycle
    double get_init_dist_btw_ind() const noexcept {return m_min_init_dist_btw_inds;}

    ///Returns the side of the grid, which is also used to determine
    /// the total size of the grid, since the grid is a square
    int get_grid_side() const noexcept {return m_grid_side;}

private:

    ///The expected size of a newly funded population
    int m_exp_new_pop_size;

    ///The minimum ditance between individuals when they are placed down at the start
    ///of a population cycle
    double m_min_init_dist_btw_inds;

    ///The side of the grid
    int m_grid_side;

    /// The diffusion coefficient of substances in the grid
    double m_diff_coeff;

    ///The rate at which metabolite degrades
    double m_metab_degradation_rate;

    ///The initial amount of food in each grid_cell at the start of a pop cycle
    double m_init_food;

    ///The probability of a mutation to happen
    double m_mutation_prob;

    ///The variance of the gaussian from which the size of the mutation is drawn
    double m_mutation_step;

    ///The base dispersal probability of an individual(so when it is active or sporulating)
    double m_base_disp_prob;

    ///The factor for which the dispersal probability is multiplied if the individual is a spore
    double m_spore_advantage;

    ///The probability that an individual will reproduce
    double m_repr_prob;
};

void test_sim_parameters() noexcept;

#endif // SIM_PARAMETERS_H
