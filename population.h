#ifndef POPULATION_H
#define POPULATION_H

#include "demographic_sim.h"
#include "funders_success.h"
#include "individual.h"
#include "pop_param.h"
#include "relaxation.hpp"
#include <utility>

class population
{
public:
    population(const pop_param &pop_parameters = pop_param(),
               const ind_param &ind_parameters = ind_param());


    ///Gets an inidividual at a certain index in the vector
    const individual& get_ind(int i) const
    {
        assert(i >= 0);
        assert(static_cast<unsigned int>(i)< m_pop.size());
        return m_pop[static_cast<unsigned int>(i)];
    }

    ///Gets an inidividual at a certain index in the vector and allows to change it
    individual& get_ind(int i)
    {
        assert(i >= 0);
        assert(static_cast<unsigned int>(i) < m_pop.size());
        return m_pop[static_cast<unsigned int>(i)];
    }

    ///Gets the vector containing all the individuals of the population
    const std::vector<individual>& get_new_v_ind() const noexcept {return m_new_pop_tmp_buffer;}

    ///Gets the vector containing all the individuals of the population
    std::vector<individual>& get_new_v_ind() noexcept {return m_new_pop_tmp_buffer;}

    ///Gets const ref to parameters
    const pop_param& get_param() const noexcept {return m_pop_param;}

    ///Gets const ref to parameters
    pop_param& get_param()  noexcept {return m_pop_param;}

    ///Gets the vector containing all the individuals of the population
    const std::vector<individual>& get_v_ind() const noexcept {return m_pop;}

    ///Gets the vector containing all the individuals of the population
    std::vector<individual>& get_v_ind() noexcept {return m_pop;}

    ///Gets the size of the population
    int get_pop_size() const noexcept {return static_cast<int>(m_pop.size());}

    ///Gets reference to Relaxation object that deals with collision management
    relaxation::Relaxation& get_relax() {return m_relax;};

    ///Gets the reference to m_rng
    auto& get_rng() noexcept {return  m_rng;}

    ///Sets the population individuals
    void set_pop_inds(const std::vector<individual>& v_inds) noexcept {m_pop = v_inds;}

private:
    ///The parameters of the population
    pop_param m_pop_param;

    ///The vector of individuals representing a population
    std::vector<individual> m_pop;

    ///The buffer vector on which the new population will be copied and that
    /// will be then swapped with the current population
    /// Used to reduce memory allocation time
    std::vector<individual> m_new_pop_tmp_buffer;

    ///The random number generator of pop(used for everything)
    std::minstd_rand m_rng;

    ///Hilbert tree collision manager
    relaxation::Relaxation m_relax;
};
///Checks if two populations are the same.
///  If they have the same parameters and the same individuals.
bool operator== (const population& lhs, const population& rhs);

///Checks if two populations are not the same.
/// If they don't have the same parameters or the same individuals.
bool operator!= (const population& lhs, const population& rhs);

///Active individuals lose energy due to their metabolism
void active_metabolism_pop(population &p);

///Checks if the entire population has been already drawn for funding the new population
bool all_ind_are_drawn(const population& s) noexcept;

///Appends a unique ID to the m_ancestor of each individual in the population
std::vector<individual> assign_ancestor_ID(const std::vector<individual>& p) noexcept;

///Counts the number of hexagonal layers necessary to place all individuals in hex pattern
int count_hex_layers(int pop_size)  noexcept ;

///Counts number of active individuals in a population
int count_active(const population& pop);

///Counts the number of spores in a population
int count_spores(const population& pop);

///Counts the number of sporulating individuals in a population
int count_sporulating(const population& pop);

///Calculate the distance between two cells given their positions
double calculate_distance(std::pair<double, double> lhs, std::pair<double, double> rhs) noexcept;

///Calculates angle between 3 positions(first position is the vertex of the angle)
///returns an angle in radiants betweem 0 and 2PI
double calc_angle_3_pos(std::pair<double, double> P1,
                        std::pair<double, double> P2,
                        std::pair<double, double> P3);

///Calculates the total displacement of each individual in the population if there are collisions
bool calc_tot_displ_pop(population& population);

///Creates a new vector of individual
///which is the same population with
//different params of the individuals in the population
std::vector<individual> change_inds(const population &p, const ind_param &new_ind_params);

///Removes dead inidviduals from population vector
void death(population& p) noexcept;

///Returns a demographic cycle object storing data about a population
demographic_cycle demographics(const population&p, const env_param&e) noexcept;

///The individuals in the vector are copied at the end of the pop vector
bool division(population& p) noexcept;

///Displaces the individuals after their total displacement
///has been calculated with calc_tot_disp_pop()
void displace_inds(std::vector<individual>& population) noexcept;

///Selects a new population of founders, swaps it with the old and positions them in hexagon patter
void fund_new_pop(population& p) noexcept;

///The new population becomes the actual population, the old population is cancelled
void place_new_pop(population& p) noexcept;

///finds the indexes of the individuals ready to divide
std::vector<int> get_dividing_individuals(const population& p) noexcept;

///Gets the excess energy from a referenced vector of index of individuals in population vector
std::vector<double> get_excess_energies(const population& p) noexcept;

///Gets an individual's energy
double get_ind_en(const population& p, int i);

//Pop
///Gets an individual's treshold energy
double get_ind_tr_en(const population& p, int i);

///Gets the position of an individual as a vector x,y
std::pair<double, double> get_ind_pos(int i);

///Gets distance in vector elements between two duaghters of same mother cell
std::vector<std::pair<int, int> > get_sisters_index_offset(const population& p) noexcept;

///Checks if any two cells are colliding, return an empty vector in this case
/// or a vecto of the indexes of the first two colliding cells
bool has_collision(population& p );

///Death without starvation risk
void jordi_death(population &p) noexcept;

/// Displaces colliding cells so that they do not collide anymore
int manage_static_collisions(population& p);

///All inidviduals lose energy due to metabolism
void metabolism_pop(population& p);

///Only spores lose energy due to metabolism
void spor_metabolism_pop(population &p);

///Calculates the modulus of the angles between all the individuals of a population
///and a given angle in radiants
std::vector<double> modulus_of_btw_ind_angles(population& p, double ang_rad);

///Checks if a mutation happens given the mutation probability of a simulation
bool mut_happens(population& p) noexcept;

///Gives back the mutation size based on initialization parameters
double mut_step() noexcept;

///Moves individuals that are perfectly on top of each other a little bit
///to allow correct displacement
void no_complete_overlap(population& p) noexcept;

///Check that in a certain range there is only one individual
bool only_one_focal(const std::vector<individual>::iterator first,
                    const std::vector<individual>::iterator last);

/// Places cells in an hexagonal grid
/// https://stackoverflow.com/questions/14280831/algorithm-to-generate-2d-magic-hexagon-lattice
void place_start_cells(population& p) noexcept;

///Returns the iterators to the range of individuals in a population
/// that could possibly be colliding with a focal individual along the x axis
std::pair<std::vector<individual>::iterator,std::vector<individual>::iterator>
possible_collisions_x(individual focal_ind, std::vector<individual>& population);

///Returns the iterators to the range of individuals that could
/// collide with a focal individual on the x_axis and as well on the y axis
/// !!! It changes the order of the elemnts in pop vector
std::pair<std::vector<individual>::iterator,std::vector<individual>::iterator>
possible_collisions_y(individual focal_ind,
                      std::vector<individual>::iterator first_x,
                      std::vector<individual>::iterator last_x,
                      std::vector<individual> &population);


///Returns a random angle between 0 and and 2PI (in radians)
double repr_angle(population& p) noexcept;

///Resets the drawn flag for all individuals in the new_population vector
void reset_drawn_fl_new_pop(population& p) noexcept;

///Resets the output nodes of the population so that the spores return active
void reset_output_nodes_pop(population &p) noexcept;

///Resets a population to its original parameters
void reset_pop(population& p) noexcept;

///Draws a 100 individual to fund the new population and puts them in m_new_pop
void select_new_pop(population& p);

///Places an individual of index i at position x,y
void set_ind_pos(individual& i, double x, double y);

///Places an individual of index i at position x,y
void set_ind_pos(individual& i, std::pair<double, double> pos);

///Sets and individual's energy
void set_ind_en(individual& i, double en);

///Sorts individuals in a given range of a vector by increasing x coordinate
void sort_inds_by_x_inc(std::vector<individual>::iterator start, std::vector<individual>::iterator end);

///Sets the population to be the funders of a certain generation
std::vector<individual> pop_from_funders(const funders_success &f_s,
                                         const demographic_sim& d_s,
                                         int generation);

///Normal death due to death rate
void senescence(population &p) noexcept;

///Kills individuals with energy = 0
void starvation(population &p) noexcept;

void test_population() noexcept;
#endif // POPULATION_H
